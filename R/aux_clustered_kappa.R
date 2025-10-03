aux_clustered_kappa <- function(S1,
                                t_clusters,  # length = nrow(S1)
                                c_clusters,  # length = ncol(S1)
                                conf_level = 0.95,
                                correct = c("none", "CR1", "CR1s")) {
  # Asymptotic CI for kappa = P(T > C): same IF/sandwich with a 1D parameter
  correct <- match.arg(correct)
  m <- nrow(S1); n <- ncol(S1)
  mn <- m * n
  kappa_hat <- sum(S1) / mn
  eps <- .Machine$double.eps
  kappa_hat <- max(min(kappa_hat, 1 - eps), eps)
  
  S1_i <- rowSums(S1)
  S1_j <- colSums(S1)
  
  psi_trt_row <- (S1_i - n * kappa_hat) / mn
  psi_ctl_row <- (S1_j - m * kappa_hat) / mn
  
  psi_trt_clu <- rowsum(matrix(psi_trt_row, ncol = 1), group = t_clusters)
  psi_ctl_clu <- rowsum(matrix(psi_ctl_row, ncol = 1), group = c_clusters)
  psi_all     <- rbind(psi_trt_clu, psi_ctl_clu)  # G x 1
  
  G     <- nrow(psi_all)
  N_tot <- m + n
  K     <- 1
  
  Sigma_CR <- crossprod(psi_all)  # scalar
  if (correct == "CR1") {
    Sigma_CR <- (G / (G - 1)) * Sigma_CR
  } else if (correct == "CR1s") {
    Sigma_CR <- (G / (G - 1)) * ((N_tot - 1) / (N_tot - K)) * Sigma_CR
  }
  se_kappa <- sqrt(drop(Sigma_CR))
  
  if (correct == "none") {
    crit <- qnorm(1 - (1 - conf_level) / 2)
    df   <- NA_integer_
  } else {
    df   <- G - 1L
    crit <- qt(1 - (1 - conf_level) / 2, df = df)
  }
  ci <- kappa_hat + c(-1, 1) * crit * se_kappa
  ci <- pmin(pmax(ci, 0), 1)
  
  list(
    G_clusters   = G,
    n_tumors_trt = m,
    n_tumors_ctl = n,
    kappa        = kappa_hat,
    se_kappa     = se_kappa,
    ci           = setNames(ci, c("lower", "upper")),
    correction   = correct,
    df           = df
  )
}
