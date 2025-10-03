#Generic engine for clustered win ratio
aux_clustered_wr <- function(S1, S2,
                             t_clusters,  # length = nrow(S1), treated tumor -> mouse id
                             c_clusters,  # length = ncol(S1), control tumor -> mouse id
                             conf_level = 0.95,
                             correct = c("none", "CR1", "CR1s")) {
  correct <- match.arg(correct)
  
  # --- dimensions & basic estimates -----------------------------------------
  m <- nrow(S1); n <- ncol(S1)
  stopifnot(identical(dim(S1), dim(S2)))
  stopifnot(length(t_clusters) == m, length(c_clusters) == n)
  
  S1_i <- rowSums(S1);  S2_i <- rowSums(S2)
  S1_j <- colSums(S1);  S2_j <- colSums(S2)
  
  mn       <- m * n
  tau1_hat <- sum(S1) / mn
  tau2_hat <- sum(S2) / mn
  eps      <- .Machine$double.eps
  tau1_hat <- max(tau1_hat, eps)
  tau2_hat <- max(tau2_hat, eps)
  WR_hat   <- tau1_hat / tau2_hat
  
  # --- row-level influence contributions (already on 1/(mn) scale) ----------
  psi_trt_row <- cbind(
    (S1_i - n * tau1_hat) / mn,
    (S2_i - n * tau2_hat) / mn
  )
  psi_ctl_row <- cbind(
    (S1_j - m * tau1_hat) / mn,
    (S2_j - m * tau2_hat) / mn
  )
  
  # --- aggregate to cluster (mouse) level -----------------------------------
  psi_trt_clu <- rowsum(psi_trt_row, group = t_clusters)  # one row per treated mouse
  psi_ctl_clu <- rowsum(psi_ctl_row, group = c_clusters)  # one row per control mouse
  psi_all     <- rbind(psi_trt_clu, psi_ctl_clu)
  
  G     <- nrow(psi_all)        # independent clusters (mice)
  N_tot <- m + n                # total tumors (used in CR1s)
  K     <- 2                    # parameters (tau1, tau2) for the delta step
  
  # --- sandwich Sigma (on correct scale already) -----------------------------
  Sigma_CR <- crossprod(psi_all)
  
  # finite-sample options
  if (correct == "CR1") {
    Sigma_CR <- (G / (G - 1)) * Sigma_CR
  } else if (correct == "CR1s") {
    Sigma_CR <- (G / (G - 1)) * ((N_tot - 1) / (N_tot - K)) * Sigma_CR
  }
  
  # --- delta for log(WR) -----------------------------------------------------
  grad        <- c(1 / tau1_hat, -1 / tau2_hat)
  var_logWR   <- drop(t(grad) %*% Sigma_CR %*% grad)
  se_logWR    <- sqrt(var_logWR)
  logWR_hat   <- log(WR_hat)
  
  # --- CI and p-value --------------------------------------------------------
  if (correct == "none") {
    crit  <- qnorm(1 - (1 - conf_level) / 2)
    stat  <- logWR_hat / se_logWR
    p_val <- 2 * (1 - pnorm(abs(stat)))
    df    <- NA_integer_
  } else {
    df    <- G - 1L
    crit  <- qt(1 - (1 - conf_level) / 2, df = df)
    stat  <- logWR_hat / se_logWR
    p_val <- 2 * (1 - pt(abs(stat), df = df))
  }
  ci <- exp(logWR_hat + c(-1, 1) * crit * se_logWR)
  
  list(
    G_clusters   = G,
    n_tumors_trt = m,
    n_tumors_ctl = n,
    win_prob_trt = tau1_hat,
    win_prob_ctl = tau2_hat,
    win_ratio    = WR_hat,
    se_log_wr    = se_logWR,
    ci           = setNames(ci, c("lower", "upper")),
    p_value      = p_val,
    test_stat    = stat,
    correction   = correct,
    df           = df
  )
}

#aux_clustered_wr(S1 = S1, S2=S2,t_clusters = id_x,c_clusters = id_y,correct = 'none')
