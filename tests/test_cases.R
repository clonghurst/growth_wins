library(tidyverse)
library(purrr)
# Below we write code to ensure the clustered win algorithm is correct

# We construct 6 treated tumors (T1–T6) and 6 control tumors (C1–C6)
# so that each pair (Ti, Ci) represents a distinct structural case
# for the clustered comparison logic in `hpc_clustered()`.
#
# Time grid: Day in {0, 5, 10, 15}
# Baseline: Baseline = 1 for all tumors
# Volumes: constant over time within each tumor
#   - "Low sAUC" tumor   -> Volume = 1 at all observed days
#   - "High sAUC" tumor  -> Volume = 2 at all observed days
#
# Survival / censoring coding:
#   - death = 1  -> tumor/animal dies at `last_day_obs`
#   - death = 0  -> tumor is censored at `last_day_obs`
#
# The 6 “golden” pairs (Ti, Ci) and expected behavior:
#
# 1) Case 1: Both die, different times  -> rule = "death_time"
#    - T1: death = 1, last_day_obs = 15
#    - C1: death = 1, last_day_obs = 10
#    - Interpretation:
#        Both tumors die; T1 dies later than C1.
#        => Treatment wins on survival time.
#        => Expected: rule == "death_time", score == +1
#
# 2) Case 2: Both die, same time  -> rule = "sAUC_equal_time"
#    - T2: death = 1, last_day_obs = 10, low sAUC  (Volume = 1)
#    - C2: death = 1, last_day_obs = 10, high sAUC (Volume = 2)
#    - Interpretation:
#        Both die on the same day; survival times are equal.
#        => Compare sAUC at the common death time (Day 10).
#        => T2 has lower sAUC, so treatment wins.
#        => Expected: rule == "sAUC_equal_time", score == +1
#
# 3) Case 3: Death vs censor, censor AFTER death  -> rule = "death_vs_censor_decided"
#    - T3: death = 0, last_day_obs = 15 (censored)
#    - C3: death = 1, last_day_obs = 10 (died earlier)
#    - Interpretation:
#        Control tumor dies at Day 10, treated tumor is still under observation
#        beyond that time (censored at Day 15).
#        => Survival alone decides the win (no sAUC comparison needed).
#        => Treatment wins on survival.
#        => Expected: rule == "death_vs_censor_decided", score == +1
#
# 4) Case 4: Death vs censor, censor EARLIER  -> rule = "sAUC_common_maxobs"
#    - T4: death = 1, last_day_obs = 10, low sAUC  (Volume = 1)
#    - C4: death = 0, last_day_obs = 5,  high sAUC (Volume = 2)
#    - Interpretation:
#        Control tumor is censored at Day 5, treated tumor dies at Day 10.
#        Censoring occurs on or before the death time, so survival times
#        are not strictly ordered in the “win” sense.
#        => Compare sAUC at the largest common observed day (Day 5).
#        => T4 has lower sAUC, so treatment wins.
#        => Expected: rule == "sAUC_common_maxobs", score == +1
#
# 5) Case 5: Both censored, same time  -> rule = "sAUC_equal_time"
#    - T5: death = 0, last_day_obs = 10, low sAUC  (Volume = 1)
#    - C5: death = 0, last_day_obs = 10, high sAUC (Volume = 2)
#    - Interpretation:
#        Both tumors are censored at Day 10; no deaths observed.
#        => Survival windows are equal.
#        => Compare sAUC at Day 10.
#        => T5 has lower sAUC, so treatment wins.
#        => Expected: rule == "sAUC_equal_time", score == +1
#
# 6) Case 6: Both censored, different times  -> rule = "sAUC_common_maxobs"
#    - T6: death = 0, last_day_obs = 15, high sAUC (Volume = 2)
#    - C6: death = 0, last_day_obs = 10, low sAUC  (Volume = 1)
#    - Interpretation:
#        Both tumors are censored, but with different follow-up lengths.
#        => Compare sAUC at the largest common observed day (Day 10).
#        => C6 has lower sAUC at Day 10, so control wins.
#        => Expected: rule == "sAUC_common_maxobs", score == -1
#
# These six pairs form a minimal 'golden dataset' for checking that
# the decision algorithm in `hpc_clustered()` correctly assigns
# the intended rule and win/loss direction for each canonical case.


# Helper: build one tumor's longitudinal record
make_tumor <- function(id, mouse, group,
                       last_day, death,
                       vol_value) {
  days <- seq(0, last_day, by = 5)
  data.frame(
    ID           = id,
    Mouse        = mouse,
    Group        = group,
    Day          = days,
    Volume       = vol_value,
    Baseline     = 1,          # scaled AUC
    death        = death,      # 1 = died, 0 = censored
    last_day_obs = last_day,   # used as "terminal time" in hpc_clustered
    stringsAsFactors = FALSE
  )
}


# Treated tumors (T1–T6)
trt_list <- list(
  # Case 1: both die, different times (T1 dies later)
  make_tumor("T1", "M1", "Trt", last_day = 15, death = 1, vol_value = 1),
  
  # Case 2: both die, same time, trt lower sAUC
  make_tumor("T2", "M2", "Trt", last_day = 10, death = 1, vol_value = 1),
  
  # Case 3: censored after control's death (trt wins on survival)
  make_tumor("T3", "M3", "Trt", last_day = 15, death = 0, vol_value = 1),
  
  # Case 4: dies later than censored control, but censor occurs earlier -> sAUC_common_maxobs
  make_tumor("T4", "M4", "Trt", last_day = 10, death = 1, vol_value = 1),
  
  # Case 5: both censored, same time, trt has lower sAUC
  make_tumor("T5", "M5", "Trt", last_day = 10, death = 0, vol_value = 1),
  
  # Case 6: both censored, different times, trt has higher sAUC
  make_tumor("T6", "M6", "Trt", last_day = 15, death = 0, vol_value = 2)
)

# Control tumors (C1–C6)
ctl_list <- list(
  # Case 1: dies earlier than T1
  make_tumor("C1", "M7", "Ctl", last_day = 10, death = 1, vol_value = 2),
  
  # Case 2: dies same time as T2, higher sAUC
  make_tumor("C2", "M8", "Ctl", last_day = 10, death = 1, vol_value = 2),
  
  # Case 3: dies at 10, trt tumor T3 censored at 15
  make_tumor("C3", "M9", "Ctl", last_day = 10, death = 1, vol_value = 1),
  
  # Case 4: censored earlier than T4's death (at 5)
  make_tumor("C4", "M10", "Ctl", last_day = 5, death = 0, vol_value = 2),
  
  # Case 5: both censored, same time, higher sAUC
  make_tumor("C5", "M11", "Ctl", last_day = 10, death = 0, vol_value = 2),
  
  # Case 6: censored earlier (10) with lower sAUC than T6
  make_tumor("C6", "M12", "Ctl", last_day = 10, death = 0, vol_value = 1)
)

test_dat <- bind_rows(trt_list, ctl_list)
res <- hpc_clustered(
  df             = test_dat,
  tumor_id_col   = "ID",
  cluster_id_col = "Mouse",
  group_col      = "Group",
  day_col        = "Day",
  vol_col        = "Volume",
  baseline_col   = "Baseline",
  death_col      = "death",
  last_col       = "last_day_obs",
  trt_label      = "Trt",
  ctl_label      = "Ctl",
  conf_level     = 0.95,
  correct        = "none"
)

# Look at the pairwise rules and scores
res$pairs %>%
  select(treat_id, control_id, rule, score) %>%
  arrange(treat_id, control_id)

golden <- res$pairs %>%
  filter(
    (treat_id == "T1" & control_id == "C1") |
      (treat_id == "T2" & control_id == "C2") |
      (treat_id == "T3" & control_id == "C3") |
      (treat_id == "T4" & control_id == "C4") |
      (treat_id == "T5" & control_id == "C5") |
      (treat_id == "T6" & control_id == "C6")
  ) %>%
  select(treat_id, control_id, rule, score)


