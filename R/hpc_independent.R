# this file analyzes independent data according to the methods described in 
# Section 2.1 of the main paper
hpc_independent <- function(df,
                            group_col    = "Group",
                            mouse_col    = "Mouse",
                            day_col      = "Day",
                            vol_col      = "Volume",
                            baseline_col = "Baseline",
                            death_col    = "death",     # 1=death, 0=survive
                            sac_day_col  = "sac_day",   # day of death if death==1, else NA
                            trt_label    = "Treatment",
                            ctl_label    = "Control",
                            estimator_type = c("kappa","win ratio"),
                            p_val_type     = c("exact","asymptotic"),
                            comparison_details = FALSE,
                            study_end_day = NULL,       # if NULL, infer from survivors
                            conf_level    = 0.95,
                            ci_method     = c("asymptotic","percentile"),
                            correct       = c("none","CR1","CR1s"),
                            nBoot         = 2000,
                            seed          = NULL) {
  
  ## ---- 0) Basic argument validation ----------------------------------------
  # typed choices
  estimator_type <- match.arg(estimator_type)
  p_val_type     <- match.arg(p_val_type)
  ci_method      <- match.arg(ci_method)
  correct        <- match.arg(correct)
  
  # df must be data.frame-like
  if (!is.data.frame(df)) {
    stop("Error: `df` must be a data.frame or tibble.")
  }
  
  # required columns must exist (name-aware)
  req <- c(group_col, mouse_col, day_col, vol_col, baseline_col, death_col, sac_day_col)
  missing_cols <- setdiff(req, names(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Error: the following required column(s) are missing from `df`: %s",
                 paste(missing_cols, collapse = ", ")))
  }
  
  # column class sanity checks (lightweight)
  if (!is.numeric(df[[day_col]])) {
    stop(sprintf("Error: `%s` must be numeric (observation day/time).", day_col))
  }
  if (!is.numeric(df[[vol_col]])) {
    stop(sprintf("Error: `%s` must be numeric (tumor volume).", vol_col))
  }
  if (!is.numeric(df[[baseline_col]])) {
    stop(sprintf("Error: `%s` must be numeric (baseline volume).", baseline_col))
  }
  if (!is.numeric(df[[death_col]])) {
    stop(sprintf("Error: `%s` must be numeric 0/1 (death indicator).", death_col))
  }
  if (!is.null(study_end_day)) {
    if (!is.numeric(study_end_day) || length(study_end_day) != 1L || !is.finite(study_end_day) || study_end_day < 0) {
      stop("Error: `study_end_day` must be a single non-negative numeric value or NULL.")
    }
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1) {
    stop("Error: `conf_level` must be a single numeric in (0, 1).")
  }
  if (!is.logical(comparison_details) || length(comparison_details) != 1L) {
    stop("Error: `comparison_details` must be a single logical (TRUE/FALSE).")
  }
  if (!is.numeric(nBoot) || length(nBoot) != 1L || nBoot < 1 || !is.finite(nBoot)) {
    stop("Error: `nBoot` must be a single positive numeric (e.g., 2000).")
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed))) {
    stop("Error: `seed` must be a single numeric or NULL.")
  }
  
  # packages
  if (!requireNamespace("dplyr", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE) ||
      !requireNamespace("purrr", quietly = TRUE) ||
      !requireNamespace("MESS",  quietly = TRUE)) {
    stop("Error: Please install packages: dplyr, tidyr, purrr, MESS.")
  }
  dplyr <- asNamespace("dplyr"); tidyr <- asNamespace("tidyr"); purrr <- asNamespace("purrr")
  
  # check group labels present in data
  present_groups <- unique(as.character(df[[group_col]]))
  need_groups <- c(trt_label, ctl_label)
  missing_groups <- setdiff(need_groups, present_groups)
  if (length(missing_groups) > 0) {
    stop(sprintf("Error: `group_col` (%s) is missing required label(s): %s",
                 group_col, paste(missing_groups, collapse = ", ")))
  }
  
  # check death coding (warn if unexpected values)
  bad_death <- setdiff(na.omit(unique(df[[death_col]])), c(0, 1))
  if (length(bad_death) > 0) {
    stop(sprintf("Error: `%s` must be coded 0/1. Found values: %s",
                 death_col, paste(bad_death, collapse = ", ")))
  }
  
  ## ---- 1) Infer study end (T) if not provided ------------------------------
  # T is taken as max day among survivors (fallback to max observed day)
  if (is.null(study_end_day)) {
    T_candidates <- df |>
      dplyr::as_tibble() |>
      dplyr::group_by(.data[[mouse_col]], .data[[death_col]]) |>
      dplyr::summarise(max_day = max(.data[[day_col]], na.rm = TRUE), .groups = "drop") |>
      dplyr::filter(.data[[death_col]] == 0) |>
      dplyr::pull(max_day)
    if (length(T_candidates) == 0) T_candidates <- df[[day_col]]
    study_end_day <- max(T_candidates, na.rm = TRUE)
  }
  T_total <- as.numeric(study_end_day)
  
  ## ---- 2) Per-mouse sAUC and Q (scaled by Baseline) ------------------------
  # For each mouse: compute sAUC up to death day or T; then compute composite Q.
  per_mouse <- df |>
    dplyr::as_tibble() |>
    dplyr::mutate(
      !!group_col := as.character(.data[[group_col]]),
      !!mouse_col := as.character(.data[[mouse_col]])
    ) |>
    dplyr::group_by(.data[[mouse_col]]) |>
    dplyr::group_modify(~{
      d <- dplyr::arrange(.x, .data[[day_col]])
      
      grp   <- unique(d[[group_col]])[1]
      death <- as.integer(unique(d[[death_col]])[1])
      sac   <- unique(d[[sac_day_col]])[1]
      base  <- unique(d[[baseline_col]])[1]
      
      # Early baseline sanity: must be positive, finite
      if (is.na(base) || !is.finite(base) || base <= 0) {
        stop(sprintf("Error: Non-positive or missing baseline in `%s` for mouse '%s'.",
                     baseline_col, unique(d[[mouse_col]])[1]))
      }
      
      # terminal day = death day (if died) or T_total (if survived)
      last_day <- if (isTRUE(death == 1)) as.numeric(sac) else T_total
      
      d_use <- d |>
        dplyr::filter(!is.na(.data[[vol_col]]), .data[[day_col]] <= last_day) |>
        dplyr::distinct(.data[[day_col]], .keep_all = TRUE) |>
        dplyr::arrange(.data[[day_col]]) |>
        dplyr::mutate(svol = .data[[vol_col]] / base)
      
      days <- d_use[[day_col]]
      svol <- d_use[["svol"]]
      
      # sAUC via MESS::auc; define as 1 if only one observation (avoids 1/0 issues)
      sauc <- if (length(days) <= 1) 1 else MESS::auc(x = days, y = svol)
      sauc <- max(sauc, .Machine$double.eps)
      
      # Composite Q: prioritize survival, then 1/sAUC as tie-breaker
      Q <- if (death == 1L) as.numeric(sac) + 1 / sauc else (T_total + 1) + 1 / sauc
      
      dplyr::tibble(
        !!group_col := grp,
        !!mouse_col := d[[mouse_col]][1],
        death   = death,
        sac_day = ifelse(isTRUE(death == 1L), as.numeric(sac), NA_real_),
        last_day = last_day,
        sAUC    = sauc,
        Q       = Q
      )
    }) |>
    dplyr::ungroup()
  
  # Split arms; ensure both present after filtering
  Q_trt <- per_mouse |> dplyr::filter(.data[[group_col]] == trt_label)
  Q_ctl <- per_mouse |> dplyr::filter(.data[[group_col]] == ctl_label)
  if (nrow(Q_trt) == 0 || nrow(Q_ctl) == 0) {
    stop("Error: after processing, one of the arms has zero mice; check `trt_label`, `ctl_label`, or missing data.")
  }
  
  xQ <- Q_trt$Q; yQ <- Q_ctl$Q
  
  ## ---- 3) Wilcoxon p-value on Q -------------------------------------------
  wt <- suppressWarnings(stats::wilcox.test(
    xQ, yQ,
    exact = (p_val_type == "exact"),
    alternative = "two.sided",
    conf.int = FALSE
  ))
  p_value_wilcox_Q <- wt$p.value
  
  ## ---- 4) Pairwise S1/S2/TIE from Q and optional rule grid -----------------
  cmp  <- outer(xQ, yQ, `-`)
  S1   <- (cmp > 0) * 1L
  S2   <- (cmp < 0) * 1L
  TIE  <- (cmp == 0) * 1L
  m <- nrow(S1); n <- ncol(S1)
  
  pairs <- NULL
  details <- NULL
  if (comparison_details) {
    # Label each pair by how the decision would be made under the independent W1
    pairs <- tidyr::crossing(
      Q_trt |> dplyr::select(t_mouse = dplyr::all_of(mouse_col),
                             t_death = death, t_sac = sac_day,
                             t_last  = last_day, t_sAUC = sAUC, t_Q = Q),
      Q_ctl |> dplyr::select(c_mouse = dplyr::all_of(mouse_col),
                             c_death = death, c_sac = sac_day,
                             c_last  = last_day, c_sAUC = sAUC, c_Q = Q)
    ) |>
      dplyr::mutate(
        rule = dplyr::case_when(
          t_death != c_death ~ "death_priority",
          t_death == 1L & c_death == 1L & t_sac != c_sac ~ "death_time",
          TRUE ~ "sAUC_equal_time"
        ),
        score = dplyr::case_when(
          rule == "death_priority" & t_death == 0L & c_death == 1L ~ +1L,
          rule == "death_priority" & t_death == 1L & c_death == 0L ~ -1L,
          rule == "death_time"     & t_sac   >  c_sac              ~ +1L,
          rule == "death_time"     & c_sac   >  t_sac              ~ -1L,
          rule == "sAUC_equal_time" & t_sAUC + 1e-12 < c_sAUC      ~ +1L,
          rule == "sAUC_equal_time" & c_sAUC + 1e-12 < t_sAUC      ~ -1L,
          TRUE ~ 0L
        )
      )
    
    # Component counts (wins/losses by rule)
    details <- list(
      n_wins               = sum(S1),
      n_losses             = sum(S2),
      n_ties               = sum(TIE),
      n_comparisons        = m * n,
      n_trt_wins_on_death  = sum(pairs$score == +1L & pairs$rule %in% c("death_priority","death_time")),
      n_trt_wins_on_sAUC   = sum(pairs$score == +1L & pairs$rule == "sAUC_equal_time"),
      n_ctl_wins_on_death  = sum(pairs$score == -1L & pairs$rule %in% c("death_priority","death_time")),
      n_ctl_wins_on_sAUC   = sum(pairs$score == -1L & pairs$rule == "sAUC_equal_time")
    )
  }
  
  ## ---- 5) Point estimates from Q-comparisons -------------------------------
  kappa_hat <- sum(S1) / (m * n)
  wr_hat    <- kappa_hat / max(1 - kappa_hat, .Machine$double.eps)
  
  ## ---- 6) Confidence intervals --------------------------------------------
  ci_out <- NULL; se_out <- NULL; df_out <- NULL; corr_used <- correct
  t_clusters <- Q_trt[[mouse_col]]
  c_clusters <- Q_ctl[[mouse_col]]
  
  if (ci_method == "asymptotic") {
    if (estimator_type == "win ratio") {
      wr <- aux_clustered_wr(S1, S2, t_clusters, c_clusters,
                             conf_level = conf_level, correct = correct)
      ci_out  <- wr$ci
      se_out  <- wr$se_log_wr
      df_out  <- wr$df
    } else {
      ka <- aux_clustered_kappa(S1, t_clusters, c_clusters,
                                conf_level = conf_level, correct = correct)
      ci_out  <- ka$ci
      se_out  <- ka$se_kappa
      df_out  <- ka$df
    }
  } else {
    # Percentile bootstrap: resample mice within arm (cluster-level bootstrap)
    if (!is.null(seed)) set.seed(seed)
    B <- as.integer(nBoot)
    boot_stat <- numeric(B)
    for (b in seq_len(B)) {
      i_trt <- sample(seq_len(nrow(Q_trt)), nrow(Q_trt), replace = TRUE)
      i_ctl <- sample(seq_len(nrow(Q_ctl)), nrow(Q_ctl), replace = TRUE)
      xb <- Q_trt$Q[i_trt]; yb <- Q_ctl$Q[i_ctl]
      cmpb <- outer(xb, yb, `-`)
      S1b  <- (cmpb > 0) * 1L
      kapb <- sum(S1b) / (length(xb) * length(yb))
      boot_stat[b] <- if (estimator_type == "win ratio") {
        kapb / max(1 - kapb, .Machine$double.eps)
      } else kapb
    }
    alpha <- (1 - conf_level) / 2
    qs <- stats::quantile(boot_stat, probs = c(alpha, 1 - alpha), names = FALSE, type = 6)
    ci_out   <- setNames(qs, c("lower","upper"))
    corr_used <- "bootstrap_percentile"
    se_out   <- stats::sd(boot_stat)
  }
  
  estimate <- if (estimator_type == "win ratio") wr_hat else kappa_hat
  
  ## ---- 7) Return -----------------------------------------------------------
  out <- list(
    estimator_type   = estimator_type,
    estimate         = estimate,
    p_value_wilcox_Q = p_value_wilcox_Q,
    conf_int         = ci_out,
    conf_level       = conf_level,
    ci_method        = ci_method,
    correction       = corr_used,
    study_end_day    = T_total,
    tie_rate_Q       = sum(TIE) / (m * n),
    Q_table          = per_mouse
  )
  if (!is.null(se_out)) out$se <- se_out
  if (!is.null(df_out)) out$df <- df_out
  if (!is.null(details)) out$comparison_details <- details
  if (!is.null(pairs))   out$pairs <- pairs
  class(out) <- c("indep_W1_result","list")
  out
}
