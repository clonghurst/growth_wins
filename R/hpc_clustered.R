# Author: Colin A. Longhurst
# Date: Nov 2025
# Website:  https://github.com/clonghurst/growth_wins (user guide w/examples)

# Hierarchical / Prioritized Estimator for Clustered Tumor Data
#
# hpc_clustered() implements the hierarchical and prioritized win-ratio
# estimator for clustered data (e.g., bilateral tumors within mouse) as
# described in Section 3 the main paper. It is the clustered analogue of
# hpc_independent(), allowing multiple tumors per animal and using
# cluster-robust inference at the mouse level.
#
# The function expects a long-format data frame with one row per
# tumor–day, containing:
#   - a unique tumor ID
#   - a cluster/animal ID (mouse)
#   - treatment group
#   - observation day
#   - tumor volume
#   - baseline volume
#   - death indicator (1 = death/event, 0 = censored/none)
#   - last observed day for that tumor
#
# For each tumor, the function:
#   * Identifies its baseline volume and last observed day (death day if
#     the tumor/animal dies, or last follow-up time otherwise).
#   * Restricts the longitudinal series to observed times up to the
#     last observed day (with a small tolerance).
#   * Computes a cumulative scaled area under the curve (sAUC) for
#     volume / baseline over time, storing the cumulative sAUC at each
#     observed day.

# Tumors are then compared between treatment and control arms using a
# pairwise decision rule that prioritizes survival time and then tumor
# volume:
#   * If both tumors die:
#       - The one with the later death time wins ("death_time").
#       - If death times are equal, compare sAUC at that time; lower sAUC
#         wins ("sAUC_equal_time").
#   * If one dies and the other is censored:
#       - If the censored tumor is strictly beyond the death
#         time, it wins on survival ("death_vs_censor_decided").
#       - Otherwise, compare sAUC at the largest common observed day
#         ("sAUC_common_maxobs").
#   * If both are censored:
#       - If follow-up times are equal, compare sAUC at
#         that common time ("sAUC_equal_time").
#       - If follow-up times differ, compare sAUC at the largest common
#         observed day ("sAUC_common_maxobs").
#   * Ties or missing sAUC at the comparison time are labeled with
#     tie/NA-specific rules (e.g. "tie_equal_time", "sAUC_common_maxobs_na").
#     This should not happen often.
#
# From these pairwise comparisons between treated and control tumors, the
# function:
#   * Builds matrices of treatment wins (S1), treatment losses (S2), and
#     ties (TIE).
#   * Counts total wins/losses/ties and decomposes wins/losses by rule
#     (death-based vs sAUC-based).
#   * Maps each tumor to its cluster/animal ID for subsequent
#     cluster-robust inference.
#
# The win ratio and its inference are obtained by calling aux_clustered_wr(),
# which uses the S1/S2 matrices and the cluster labels to:
#   * Estimate the win ratio Psi = kappa / (1 - kappa).
#   * Compute a cluster-robust standard error for log(WR).
#   * Construct confidence intervals at the requested conf_level.
#   * Optionally apply small-sample corrections (CR1 / CR1s) based on
#     the number of clusters and tumors.
#
# The returned object includes:
#   * win_ratio, its standard error, CI, and p-value (from aux_clustered_wr()).
#   * counts of wins/losses/ties and rule-specific win counts.
#   * the full pairwise comparison table (which rule decided each pair
#     and the corresponding score).
#   * the S1, S2, and TIE matrices and the ordering/cluster labels for
#     treated and control tumors.


hpc_clustered <- function(df,
                          tumor_id_col   = "ID",        # unique tumor id (e.g. "401:R")
                          cluster_id_col = "Mouse",     # animal id for clustering (e.g. "401")
                          group_col      = "Group",     # e.g. "JQ1", "Vehicle"
                          day_col        = "Day",       # numeric day
                          vol_col        = "Volume",    # tumor volume
                          baseline_col   = "Baseline",  # baseline volume at day 0
                          death_col      = "death",     # 1=death, 0=censored/none
                          last_col       = "last_day_obs", # last observed day for that tumor
                          trt_label      = "Treatment",
                          ctl_label      = "Control",
                          conf_level     = 0.95,
                          correct        = c("none", "CR1", "CR1s"),
                          tol            = 1e-12) {
  
  ## ---- 0) Argument validation ----------------------------------------------
  correct <- match.arg(correct)
  
  # df must be data.frame-like
  if (!is.data.frame(df)) {
    stop("Error: `df` must be a data.frame or tibble.")
  }
  
  if (missing(baseline_col)) {
    warning("Argument `baseline_col` was not specified; using default \"Baseline\". ",
            "It is recommended to explicitly pass `baseline_col` to avoid mistakes.")
  }
  if (missing(last_col)) {
    warning("Argument `last_col` was not specified; using default \"last_day_obs\". ",
            "It is recommended to explicitly pass `last_col` to avoid mistakes.")
  }
  
  # required columns present, report which are missing
  need <- c(tumor_id_col, cluster_id_col, group_col, day_col, vol_col,
            baseline_col, death_col, last_col)
  missing_cols <- setdiff(need, names(df))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Error: the following required column(s) are missing from `df`: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }
  
  # types / basic coding checks
  if (!is.numeric(df[[day_col]])) {
    stop(sprintf("Error: `%s` must be numeric (observation day/time).", day_col))
  }
  if (!is.numeric(df[[vol_col]])) {
    stop(sprintf("Error: `%s` must be numeric (tumor volume).", vol_col))
  }
  if (!is.numeric(df[[baseline_col]])) {
    stop(sprintf("Error: `%s` must be numeric (baseline volume).", baseline_col))
  }
  if (!is.numeric(df[[last_col]])) {
    stop(sprintf("Error: `%s` must be numeric (last observed day).", last_col))
  }
  if (!is.numeric(df[[death_col]])) {
    stop(sprintf("Error: `%s` must be numeric 0/1 (death indicator).", death_col))
  }
  bad_death <- setdiff(na.omit(unique(df[[death_col]])), c(0, 1))
  if (length(bad_death) > 0) {
    stop(sprintf("Error: `%s` must be coded 0/1. Found values: %s",
                 death_col, paste(bad_death, collapse = ", ")))
  }
  
  # conf_level and tol
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("Error: `conf_level` must be a single numeric in (0, 1).")
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0) {
    stop("Error: `tol` must be a single non-negative numeric (e.g., 1e-12).")
  }
  
  # package check
  if (!requireNamespace("MESS", quietly = TRUE)) {
    stop("Error: package 'MESS' is required. Install via install.packages('MESS').")
  }
  
  # group labels present
  present_groups <- unique(as.character(df[[group_col]]))
  need_groups <- c(trt_label, ctl_label)
  missing_groups <- setdiff(need_groups, present_groups)
  if (length(missing_groups) > 0) {
    stop(sprintf("Error: `%s` is missing required label(s): %s",
                 group_col, paste(missing_groups, collapse = ", ")))
  }
  
  # some checks around last-day values
  if (any(!is.finite(df[[last_col]]))) {
    stop(sprintf("Error: `%s` contains non-finite values.", last_col))
  }
  if (any(df[[last_col]] < 0, na.rm = TRUE)) {
    stop(sprintf("Error: `%s` contains negative values.", last_col))
  }
  
  ## rep data and helper funcs
  dn <- function(x) df[[x]]
  
  # Drop rows with missing observation day (cannot be ordered / used)
  df <- df[!is.na(dn(day_col)), , drop = FALSE]
  
  # Split into per-tumor longitudinal records
  split_by_tumor <- split(df, dn(tumor_id_col))
  
  ## Per-tumor cumulative sAUC vector (scaled by baseline) 
  # For each tumor: determine baseline, restrict to <= last_day, ensure Day 0, compute cumulative sAUC
  tumor_info <- lapply(split_by_tumor, function(d) {
    d <- d[order(d[[day_col]]), , drop = FALSE]
    
    # Baseline: prefer explicit column; fallback to volume nearest to Day 0
    base <- unique(d[[baseline_col]])
    base <- base[!is.na(base)]
    if (length(base) != 1L) {
      # fallback: closest-to-zero day’s volume
      base <- d[[vol_col]][which.min(abs(d[[day_col]] - 0))]
    }
    if (is.na(base) || !is.finite(base) || base <= 0) {
      stop("Error: Non-positive or missing baseline encountered for a tumor in `", baseline_col, "`.")
    }
    
    # Terminal observed day for this tumor (death day if dead; otherwise censoring/last follow-up)
    last <- as.numeric(unique(d[[last_col]])[1])
    if (!is.finite(last) || last < 0) {
      stop("Error: invalid `last_day_obs` detected (must be finite and non-negative).")
    }
    
    # Restrict to observed non-NA volumes up to (last + tol)
    d <- d[!is.na(d[[vol_col]]) & d[[day_col]] <= last + tol, , drop = FALSE]
    
    # Ensure a Day 0 row exists; add synthetic row if needed to anchor sAUC
    if (!any(abs(d[[day_col]] - 0) <= tol)) {
      d0 <- d[1, , drop = FALSE]
      d0[[day_col]] <- 0
      d0[[vol_col]] <- base
      d <- rbind(d0, d)
      d <- d[order(d[[day_col]]), , drop = FALSE]
    }
    
    # Deduplicate days: keep the last per day if duplicates
    d <- d[!duplicated(d[[day_col]], fromLast = TRUE), , drop = FALSE]
    
    # Build cumulative sAUC across observed times
    days <- as.numeric(d[[day_col]])
    svol <- as.numeric(d[[vol_col]]) / base
    cum_sauc <- vapply(seq_along(days), function(k) {
      MESS::auc(x = days[1:k], y = svol[1:k])
    }, numeric(1))
    names(cum_sauc) <- as.character(days)
    
    list(
      cluster  = as.character(unique(d[[cluster_id_col]])[1]),
      group    = as.character(unique(d[[group_col]])[1]),
      death    = as.integer(unique(d[[death_col]])[1]),
      last     = last,
      days     = days,
      cum_sauc = cum_sauc
    )
  })
  
  ## Determine arm membership 
  all_ids <- names(tumor_info)
  t_ids <- all_ids[vapply(tumor_info, function(x) identical(x$group, trt_label), logical(1))]
  c_ids <- all_ids[vapply(tumor_info, function(x) identical(x$group, ctl_label), logical(1))]
  if (length(t_ids) == 0L || length(c_ids) == 0L) {
    stop("Error: No treated or control tumors found. Check `group_col`, `trt_label`, `ctl_label`.")
  }
  
  ## Some utilities for pairwise decisions 
  # Largest common observed day up to a target 'upto' (default: min of terminal days)
  max_common_day <- function(tu, cu, upto = min(tu$last, cu$last)) {
    td <- tu$days[tu$days <= upto + tol]
    cd <- cu$days[cu$days <= upto + tol]
    common <- intersect(td, cd)
    if (length(common) == 0L) return(NA_real_)
    max(common)
  }
  # sAUC value at a specific observed day
  sAUC_at <- function(tu, d_star) {
    if (is.na(d_star)) return(NA_real_)
    tu$cum_sauc[[as.character(d_star)]]
  }
  
  ## Pairwise comparisons implementing the win strategy from the paper
  decide_pair <- function(tu, cu) {
    # Case A: both died -> later death wins; tie on time -> compare sAUC at that time
    if (tu$death == 1L && cu$death == 1L) {
      if (tu$last > cu$last + tol) return(list(score = +1L, rule = "death_time"))
      if (cu$last > tu$last + tol) return(list(score = -1L, rule = "death_time"))
      d  <- max_common_day(tu, cu, upto = tu$last)
      st <- sAUC_at(tu, d); sc <- sAUC_at(cu, d)
      if (is.na(st) || is.na(sc)) return(list(score = 0L, rule = "sAUC_equal_time_na"))
      if (st + tol < sc) return(list(score = +1L, rule = "sAUC_equal_time"))
      if (sc + tol < st) return(list(score = -1L, rule = "sAUC_equal_time"))
      return(list(score = 0L, rule = "tie_equal_time"))
    }
    
    # Case B: one died, one censored
    if (tu$death != cu$death) {
      dead <- if (tu$death == 1L) tu else cu
      cens <- if (tu$death == 0L) tu else cu
      # Censor strictly after death => censored wins on survival (decided on death vs censor)
      if (cens$last > dead$last - tol) { # this was the bug
        if (identical(cens, tu)) return(list(score = +1L, rule = "death_vs_censor_decided"))
        else                     return(list(score = -1L, rule = "death_vs_censor_decided"))
      }
      # Otherwise compare sAUC at the largest common observed day
      upto <- min(tu$last, cu$last)
      d  <- max_common_day(tu, cu, upto = upto)
      st <- sAUC_at(tu, d); sc <- sAUC_at(cu, d)
      if (is.na(st) || is.na(sc)) return(list(score = 0L, rule = "sAUC_common_maxobs_na"))
      if (st + tol < sc) return(list(score = +1L, rule = "sAUC_common_maxobs"))
      if (sc + tol < st) return(list(score = -1L, rule = "sAUC_common_maxobs"))
      return(list(score = 0L, rule = "tie_common_maxobs"))
    }
    
    # Case C: both censored (no deaths)
    if (tu$death == 0L && cu$death == 0L) {
      # Equal follow-up (within tol) -> compare sAUC at that equal time
      if (abs(tu$last - cu$last) <= tol) {
        d  <- max_common_day(tu, cu, upto = tu$last)
        st <- sAUC_at(tu, d); sc <- sAUC_at(cu, d)
        if (is.na(st) || is.na(sc)) return(list(score = 0L, rule = "sAUC_equal_time_na"))
        if (st + tol < sc) return(list(score = +1L, rule = "sAUC_equal_time"))
        if (sc + tol < st) return(list(score = -1L, rule = "sAUC_equal_time"))
        return(list(score = 0L, rule = "tie_equal_time"))
      } else {
        # Unequal follow-up -> compare at largest common observed day
        d  <- max_common_day(tu, cu)
        st <- sAUC_at(tu, d); sc <- sAUC_at(cu, d)
        if (is.na(st) || is.na(sc)) return(list(score = 0L, rule = "sAUC_common_maxobs_na"))
        if (st + tol < sc) return(list(score = +1L, rule = "sAUC_common_maxobs"))
        if (sc + tol < st) return(list(score = -1L, rule = "sAUC_common_maxobs"))
        return(list(score = 0L, rule = "tie_common_maxobs"))
      }
    }
    
    # Fallback (should be rare)
    list(score = 0L, rule = "indeterminate")
  }
  
  # see the test_cases.R 
  
  ##  Build S1/S2/TIE and the pairwise table for the aux variance function
  m <- length(t_ids); n <- length(c_ids)
  S1  <- matrix(0L, nrow = m, ncol = n, dimnames = list(t_ids, c_ids))  # treatment wins
  S2  <- matrix(0L, nrow = m, ncol = n, dimnames = list(t_ids, c_ids))  # treatment losses
  TIE <- matrix(0L, nrow = m, ncol = n, dimnames = list(t_ids, c_ids))  # ties
  
  rows <- vector("list", m * n); k <- 0L
  for (i in seq_len(m)) {
    tu <- tumor_info[[t_ids[i]]]
    for (j in seq_len(n)) {
      cu <- tumor_info[[c_ids[j]]]
      out <- decide_pair(tu, cu)
      if (out$score > 0L)      S1[i, j] <- 1L
      else if (out$score < 0L) S2[i, j] <- 1L
      else                     TIE[i, j] <- 1L
      k <- k + 1L
      rows[[k]] <- data.frame(
        treat_id        = t_ids[i],
        control_id      = c_ids[j],
        treat_cluster   = tu$cluster,
        control_cluster = cu$cluster,
        treat_death     = tu$death,
        control_death   = cu$death,
        treat_last      = tu$last,
        control_last    = cu$last,
        rule            = out$rule,
        score           = out$score,
        stringsAsFactors = FALSE
      )
    }
  }
  pairs <- do.call(rbind, rows)
  
  ##  Tally up the decisions
  n_wins   <- sum(S1)
  n_losses <- sum(S2)
  n_ties   <- sum(TIE)
  n_comp   <- m * n
  
  death_rules   <- c("death_time", "death_vs_censor_decided")
  is_sauc_rule  <- grepl("^sAUC_", pairs$rule)
  
  n_trt_wins_death <- sum(pairs$score == +1L & pairs$rule %in% death_rules)
  n_trt_wins_sAUC  <- sum(pairs$score == +1L & is_sauc_rule)
  n_ctl_wins_death <- sum(pairs$score == -1L & pairs$rule %in% death_rules)
  n_ctl_wins_sAUC  <- sum(pairs$score == -1L & is_sauc_rule)
  
  ## Note Cluster labels for sandwich variance estimator code
  t_clusters <- vapply(t_ids, function(id) tumor_info[[id]]$cluster, character(1))
  c_clusters <- vapply(c_ids, function(id) tumor_info[[id]]$cluster, character(1))
  
  ## Cluster-robust WR via aux_clustered_wr 
  est <- aux_clustered_wr(S1, S2,
                          t_clusters = t_clusters,
                          c_clusters = c_clusters,
                          conf_level = conf_level,
                          correct    = correct)
  
  ##  Things to return to the user
  c(est, list(
    counts = list(
      n_wins               = n_wins,
      n_losses             = n_losses,
      n_ties               = n_ties,
      n_comparisons        = n_comp,
      n_trt_wins_on_death  = n_trt_wins_death,
      n_trt_wins_on_sAUC_d = n_trt_wins_sAUC,
      n_ctl_wins_on_death  = n_ctl_wins_death,
      n_ctl_wins_on_sAUC_d = n_ctl_wins_sAUC
    ),
    pairs       = pairs,         # which rule decided each comparison
    S1 = S1, S2 = S2, TIE = TIE, # matrices used for WR estimation
    t_order     = t_ids, c_order = c_ids,
    t_clusters  = t_clusters, c_clusters = c_clusters
  ))
}
