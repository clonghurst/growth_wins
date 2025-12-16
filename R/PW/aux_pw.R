aux_pw <- function(
    subject_data,
    pair_data,
    covariates,
    id_col      = "Tumor_ID",
    id_i_col    = "id_i",
    id_j_col    = "id_j",
    cluster_col = NULL,
    control     = list(tol = 1e-8, maxit = 40),
    small_sample = FALSE,                    # independent case
    correction   = c("none", "CR1", "CR1s")  # clustered case
) {
  correction <- match.arg(correction)
  
  if (!is.data.frame(subject_data))
    stop("`subject_data` must be a data.frame with one row per tumor.")
  if (!is.data.frame(pair_data))
    stop("`pair_data` must be a data.frame produced by pw_build_clustered_longitudinal().")
  
  if (!all(c(id_col) %in% names(subject_data)))
    stop("`subject_data` must contain column `", id_col, "`.")
  if (!all(c(id_i_col, id_j_col, "win_ij") %in% names(pair_data)))
    stop("`pair_data` must contain `", id_i_col, "`, `", id_j_col, "`, and `win_ij`.")
  
  # ---- 1. design matrix Z (subject-level) -----------------------------------
  if (inherits(covariates, "formula")) {
    Z <- model.matrix(covariates, subject_data)
  } else {
    Z <- as.matrix(subject_data[, covariates, drop = FALSE])
    Z <- cbind("(Intercept)" = 1, Z)
  }
  # drop zero-variance columns (e.g. constant intercept)
  keep_col <- apply(Z, 2, var) > 0
  Z <- Z[, keep_col, drop = FALSE]
  
  p <- ncol(Z)
  n <- nrow(Z)
  if (p == 0L)
    stop("No varying covariates to estimate.")
  
  ids <- as.character(subject_data[[id_col]])
  
  # ---- 2. map pairs to row indices & build outcomes -------------------------
  i_ids <- as.character(pair_data[[id_i_col]])
  j_ids <- as.character(pair_data[[id_j_col]])
  win_ij <- as.integer(pair_data[["win_ij"]])   # 1, -1, or 0
  
  i_idx <- match(i_ids, ids)
  j_idx <- match(j_ids, ids)
  if (anyNA(i_idx) || anyNA(j_idx)) {
    stop("Some pairwise IDs not found in `subject_data`.")
  }
  
  # We treat "first wins" vs "first loses" as binary; ties dropped
  keep <- win_ij != 0L
  if (!any(keep)) stop("All pairwise comparisons are ties or missing.")
  
  i_idx  <- i_idx[keep]
  j_idx  <- j_idx[keep]
  win_ij <- win_ij[keep]
  
  # Binary response: 1 if first wins, 0 if first loses
  wins01 <- as.integer(win_ij == 1L)
  
  m_eff <- length(wins01)
  
  # Pairwise covariate differences
  dZ <- Z[i_idx, , drop = FALSE] - Z[j_idx, , drop = FALSE]   # m_eff × p
  
  # ---- 3. Newton–Raphson for β̂ ---------------------------------------------
  score_fn <- function(beta) {
    mu <- plogis(drop(dZ %*% beta))                      # P(first wins)
    colSums((wins01 - mu) * dZ)
  }
  hess_fn <- function(beta) {
    eta <- drop(dZ %*% beta)
    w   <- plogis(eta) * (1 - plogis(eta))
    t(dZ * w) %*% dZ
  }
  
  beta_hat <- rep(0, p)
  for (iter in seq_len(control$maxit)) {
    S <- score_fn(beta_hat)
    if (max(abs(S)) < control$tol) break
    H <- hess_fn(beta_hat)
    step <- tryCatch(solve(H, S), error = function(e) NA)
    if (anyNA(step))
      stop("Hessian singular – check collinearity or extreme separation.")
    beta_hat <- beta_hat + step
    if (max(abs(step)) < control$tol) break
  }
  converged <- iter < control$maxit
  if (!converged)
    warning("aux_pw: Newton–Raphson did not fully converge (", iter, " iterations).")
  
  # ---- 4. Influence functions (Mao & Wang style) ----------------------------
  eta_hat <- drop(dZ %*% beta_hat)
  mu_ij   <- plogis(eta_hat)
  w_den   <- mu_ij * (1 - mu_ij)
  
  Ahat <- - crossprod(dZ * w_den, dZ) / m_eff   # mean derivative
  invA <- solve(Ahat)
  
  pair_score <- ((wins01) - mu_ij) * dZ        # m_eff × p
  
  # subject-level kappa: accumulate pair contributions for BOTH i and j
  kappa <- matrix(0, n, p)
  for (k in seq_len(m_eff)) {
    kappa[i_idx[k], ] <- kappa[i_idx[k], ] + pair_score[k, ]
    kappa[j_idx[k], ] <- kappa[j_idx[k], ] + pair_score[k, ]
  }
  kappa <- kappa / (n - 1)                     # divide by (n-1)
  
  psi_mat <- -2 * kappa %*% t(invA)           # n × p
  
  # ---- 5. Variance: independent vs clustered -------------------------------
  G_clusters <- NA_integer_
  df         <- NA_integer_
  test_type  <- "z"
  
  if (is.null(cluster_col)) {
    # independent (original Mao–Wang variance)
    S  <- crossprod(psi_mat) / n
    infl <- if (isTRUE(small_sample)) n / (n - p - 1) else 1
    vcov <- S * infl / n
    
  } else {
    # cluster-robust at animal level
    if (!cluster_col %in% names(subject_data))
      stop("`subject_data` must contain `cluster_col` when cluster-robust variance is requested.")
    
    cl <- subject_data[[cluster_col]]
    if (length(cl) != n)
      stop("cluster_col must have length equal to nrow(subject_data).")
    
    cl <- droplevels(factor(cl))
    G  <- nlevels(cl)
    G_clusters <- G
    
    psi_clu <- rowsum(psi_mat, group = cl)      # G × p (sum over tumors in cluster)
    
    N_tot <- n
    Sigma <- crossprod(psi_clu) / (N_tot^2)     # (1/N^2) Σ_g Ψ_g Ψ_gᵀ
    
    if (correction == "CR1") {
      Sigma <- (G / (G - 1)) * Sigma
    } else if (correction == "CR1s") {
      Sigma <- (G / (G - 1)) * ((N_tot - 1) / (N_tot - p)) * Sigma
    }
    vcov <- Sigma
    
    if (correction != "none") {
      df        <- G - 1L
      test_type <- "t"
    }
  }
  
  se   <- sqrt(diag(vcov))
  stat <- beta_hat / se
  
  # ---- 6. p-values & summary table -----------------------------------------
  if (test_type == "z") {
    pvals <- 2 * pnorm(-abs(stat))
    tab <- cbind(
      Estimate    = beta_hat,
      `Std. Error` = se,
      `z value`   = stat,
      `Pr(>|z|)`  = pvals
    )
  } else {
    pvals <- 2 * (1 - pt(abs(stat), df = df))
    tab <- cbind(
      Estimate     = beta_hat,
      `Std. Error` = se,
      `t value`    = stat,
      `Pr(>|t|)`   = pvals
    )
  }
  rownames(tab) <- colnames(Z)
  
  out <- list(
    coef        = setNames(drop(beta_hat), colnames(Z)),
    vcov        = vcov,
    se          = se,
    table       = tab,
    iterations  = iter,
    converged   = converged,
    n_pairs     = m_eff,
    n_obs       = n,
    G_clusters  = G_clusters,
    correction  = if (is.null(cluster_col)) "none" else correction,
    df          = df,
    test_type   = test_type,
    clustered   = !is.null(cluster_col)
  )
  class(out) <- "ee_wins"
  out
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
# # 1) Build pairwise wins using the clustered longitudinal kernel
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
# subjects <- pw_obj$subjects   # one row per tumor
# pairs    <- pw_obj$pairs      # pairwise wins/losses grid
# 
# # 2) Fit proportional win-fraction model with cluster-robust variance
# fit_pw <- aux_pw(
#   subject_data = subjects,
#   pair_data    = pairs,
#   covariates   = ~ Drug1 * Drug2,#~ Group,   # or ~ Drug1 + Drug2, etc.
#   id_col       = "Tumor_ID",
#   cluster_col  = "Cluster_ID",
#   correction   = "CR1s"     # small-sample cluster correction
# )
# 
# fit_pw$table




