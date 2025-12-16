# a function that calculates win/loss for multi group studies
# this is a dev function, as it needs to be wrapped in with the aux
# function for easy use for future users



pw_build_clustered_longitudinal <- function(
    df,
    id_col        = "Tumor_ID",      # unique tumor ID
    cluster_col   = "Cluster_ID",    # animal / mouse ID
    day_col       = "Day",
    vol_col       = "volume",
    baseline_col  = "Baseline",
    death_col     = "death",         # 1 = death, 0 = censored / no death
    last_col      = "last_day_obs",  # last observed day for that tumor
    subject_cols  = NULL,            # columns to keep for covariates (e.g. Drug1, Drug2, Group, LV)
    tol           = 1e-8
) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame / tibble.")
  }
  
  need <- c(id_col, cluster_col, day_col, vol_col,
            baseline_col, death_col, last_col)
  missing_cols <- setdiff(need, names(df))
  if (length(missing_cols)) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
  }
  
  if (!is.numeric(df[[day_col]]))
    stop(sprintf("`%s` must be numeric (observation day).", day_col))
  if (!is.numeric(df[[vol_col]]))
    stop(sprintf("`%s` must be numeric (tumor volume).", vol_col))
  if (!is.numeric(df[[baseline_col]]))
    stop(sprintf("`%s` must be numeric (baseline volume).", baseline_col))
  if (!is.numeric(df[[last_col]]))
    stop(sprintf("`%s` must be numeric (last observed day).", last_col))
  if (!is.numeric(df[[death_col]]))
    stop(sprintf("`%s` must be numeric 0/1 (death indicator).", death_col))
  
  bad_death <- setdiff(na.omit(unique(df[[death_col]])), c(0, 1))
  if (length(bad_death) > 0) {
    stop(sprintf("`%s` must be coded 0/1. Found: %s",
                 death_col, paste(bad_death, collapse = ", ")))
  }
  
  if (!requireNamespace("MESS", quietly = TRUE)) {
    stop("Package 'MESS' is required (for AUC). Install via install.packages('MESS').")
  }
  
  # Drop rows with missing day or volume
  df <- df[!is.na(df[[day_col]]) & !is.na(df[[vol_col]]), , drop = FALSE]
  
  # Per-tumor longitudinal info: days, cum sAUC, death, last, cluster
  split_by_tumor <- split(df, df[[id_col]])
  
  tumor_info <- lapply(split_by_tumor, function(d) {
    d <- d[order(d[[day_col]]), , drop = FALSE]
    
    # Baseline: prefer baseline_col; fallback to volume at time 0
    base <- unique(d[[baseline_col]])
    base <- base[!is.na(base)]
    if (length(base) != 1L) {
      base <- d[[vol_col]][which.min(abs(d[[day_col]] - 0))]
    }
    if (is.na(base) || !is.finite(base) || base <= 0) {
      stop("Non-positive or missing baseline encountered for a tumor.")
    }
    
    last <- unique(d[[last_col]])
    last <- last[!is.na(last)]
    if (length(last) != 1L) {
      stop("Each tumor must have a unique `last_day_obs`.")
    }
    last <- as.numeric(last)
    
    # Restrict to observed days up to last+tol
    d <- d[d[[day_col]] <= last + tol, , drop = FALSE]
    
    # Ensure day 0 exists; if not, add synthetic point at baseline
    if (!any(abs(d[[day_col]] - 0) <= tol)) {
      d0 <- d[1, , drop = FALSE]
      d0[[day_col]] <- 0
      d0[[vol_col]] <- base
      d <- rbind(d0, d)
      d <- d[order(d[[day_col]]), , drop = FALSE]
    }
    
    # Deduplicate days: keep last occurrence
    d <- d[!duplicated(d[[day_col]], fromLast = TRUE), , drop = FALSE]
    
    days <- as.numeric(d[[day_col]])
    svol <- as.numeric(d[[vol_col]]) / base
    
    cum_sauc <- vapply(seq_along(days), function(k) {
      MESS::auc(x = days[1:k], y = svol[1:k])
    }, numeric(1))
    names(cum_sauc) <- as.character(days)
    
    list(
      cluster  = as.character(unique(d[[cluster_col]])[1]),
      death    = as.integer(unique(d[[death_col]])[1]),
      last     = last,
      days     = days,
      cum_sauc = cum_sauc
    )
  })
  
  # Utility: largest common observed day up to `upto`
  max_common_day <- function(tu, cu, upto = min(tu$last, cu$last)) {
    td <- tu$days[tu$days <= upto + tol]
    cd <- cu$days[cu$days <= upto + tol]
    common <- intersect(td, cd)
    if (length(common) == 0L) return(NA_real_)
    max(common)
  }
  
  # sAUC at specific day (if observed)
  sAUC_at <- function(tu, d_star) {
    if (is.na(d_star)) return(NA_real_)
    tu$cum_sauc[[as.character(d_star)]]
  }
  
  # Pairwise decision rule: 1 = first wins, -1 = second wins, 0 = tie/NA
  decide_pair <- function(tu, cu) {
    # (A) both died
    if (tu$death == 1L && cu$death == 1L) {
      if (tu$last > cu$last + tol) return(list(score = +1L, rule = "death_time"))
      if (cu$last > tu$last + tol) return(list(score = -1L, rule = "death_time"))
      # equal death time -> compare sAUC at that time
      d  <- max_common_day(tu, cu, upto = tu$last)
      st <- sAUC_at(tu, d)
      sc <- sAUC_at(cu, d)
      if (is.na(st) || is.na(sc)) return(list(score = 0L, rule = "sAUC_equal_time_na"))
      if (st + tol < sc) return(list(score = +1L, rule = "sAUC_equal_time"))
      if (sc + tol < st) return(list(score = -1L, rule = "sAUC_equal_time"))
      return(list(score = 0L, rule = "tie_equal_time"))
    }
    
    # (B) one died, one censored
    if (tu$death != cu$death) {
      dead <- if (tu$death == 1L) tu else cu
      cens <- if (tu$death == 0L) tu else cu
      
      # censor strictly after death -> censored wins on survival
      if (cens$last > dead$last - tol) { # this was the bug!!!
        if (identical(cens, tu)) return(list(score = +1L, rule = "death_vs_censor_decided"))
        else                     return(list(score = -1L, rule = "death_vs_censor_decided"))
      }
      # otherwise compare sAUC at largest common observed day
      upto <- min(tu$last, cu$last)
      d  <- max_common_day(tu, cu, upto = upto)
      st <- sAUC_at(tu, d)
      sc <- sAUC_at(cu, d)
      if (is.na(st) || is.na(sc)) return(list(score = 0L, rule = "sAUC_common_maxobs_na"))
      if (st + tol < sc) return(list(score = +1L, rule = "sAUC_common_maxobs"))
      if (sc + tol < st) return(list(score = -1L, rule = "sAUC_common_maxobs"))
      return(list(score = 0L, rule = "tie_common_maxobs"))
    }
    
    # (C) both censored (no deaths)
    if (tu$death == 0L && cu$death == 0L) {
      # equal follow-up -> compare sAUC at that time
      if (abs(tu$last - cu$last) <= tol) {
        d  <- max_common_day(tu, cu, upto = tu$last)
        st <- sAUC_at(tu, d)
        sc <- sAUC_at(cu, d)
        if (is.na(st) || is.na(sc)) return(list(score = 0L, rule = "sAUC_equal_time_na"))
        if (st + tol < sc) return(list(score = +1L, rule = "sAUC_equal_time"))
        if (sc + tol < st) return(list(score = -1L, rule = "sAUC_equal_time"))
        return(list(score = 0L, rule = "tie_equal_time"))
      } else {
        # unequal follow-up -> compare at largest common observed day
        d  <- max_common_day(tu, cu)
        st <- sAUC_at(tu, d)
        sc <- sAUC_at(cu, d)
        if (is.na(st) || is.na(sc)) return(list(score = 0L, rule = "sAUC_common_maxobs_na"))
        if (st + tol < sc) return(list(score = +1L, rule = "sAUC_common_maxobs"))
        if (sc + tol < st) return(list(score = -1L, rule = "sAUC_common_maxobs"))
        return(list(score = 0L, rule = "tie_common_maxobs"))
      }
    }
    
    # should not reach here
    list(score = 0L, rule = "indeterminate")
  }
  
  # Build pairwise grid over all tumors (unordered pairs)
  all_ids  <- names(tumor_info)
  n_tumors <- length(all_ids)
  if (n_tumors < 2L) stop("Need at least 2 tumors to form pairwise comparisons.")
  
  n_pairs <- n_tumors * (n_tumors - 1L) / 2L
  rows <- vector("list", n_pairs)
  k <- 0L
  
  for (a in seq_len(n_tumors - 1L)) {
    for (b in seq.int(a + 1L, n_tumors)) {
      id_i <- all_ids[a]
      id_j <- all_ids[b]
      ti   <- tumor_info[[id_i]]
      tj   <- tumor_info[[id_j]]
      out  <- decide_pair(ti, tj)
      k <- k + 1L
      rows[[k]] <- data.frame(
        id_i      = id_i,
        id_j      = id_j,
        cluster_i = ti$cluster,
        cluster_j = tj$cluster,
        death_i   = ti$death,
        death_j   = tj$death,
        last_i    = ti$last,
        last_j    = tj$last,
        win_ij    = out$score,   # 1 = id_i wins, -1 = id_j wins, 0 = tie
        rule      = out$rule,
        stringsAsFactors = FALSE
      )
    }
  }
  pairs <- do.call(rbind, rows)
  
  # Subject-level data (one row per tumor) for covariate modeling
  keep_cols <- c(id_col, cluster_col, subject_cols)
  keep_cols <- intersect(keep_cols, names(df))
  subjects <- df[!duplicated(df[[id_col]]), keep_cols, drop = FALSE]
  rownames(subjects) <- NULL
  
  list(
    subjects    = subjects,
    pairs       = pairs,
    id_col      = id_col,
    cluster_col = cluster_col
  )
}

# set.seed(1)
# skeleton <- PW_clust_data_create(
#   number_cluster      = 10,
#   number_per_cluster  = 2,
#   longitudinal_schedule = 4
# )
# sim_dat <- volume_PW_clust_create(
#   skeleton_data  = skeleton,
#   re_cluster     = 0.2,
#   re_tumor       = 0.1,
#   residual_sd    = 300,
#   intercept      = 5.3,
#   linear_predictor = "+ 0.8*Day_scale +
#                       (-0.18*Drug1*Day_scale) +
#                       (-0.10*Drug2*Day_scale) +
#                       Cluster_RE + Tumor_RE",
#   lower_censor   = 20,
#   upper_censor   = 1500,
#   cluster_censor = TRUE
# )
# pw_obj <- pw_build_clustered_longitudinal(
#   df          = sim_dat,
#   id_col      = "Tumor_ID",
#   cluster_col = "Cluster_ID",
#   day_col     = "Day",
#   vol_col     = "volume",
#   baseline_col= "Baseline",
#   death_col   = "death",
#   last_col    = "last_day_obs",
#   subject_cols = c("Group", "Drug1", "Drug2")  # to use as covariates
# )
# 
# 
# head(sim_dat)

