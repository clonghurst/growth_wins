library(tidyverse)
library(magrittr)
library(numDeriv)
pocl <- function(
    formula,
    data,
    ref = NULL,     # baseline factor level (e.g. "A")
    start = c(0,0,0),  # initial guess for (beta2, beta3, beta4)
    maxit = 10,
    tol   = 1e-7, # maybe make this lower
    trace = TRUE
) {
  hyper4_condstats_fast <- function(s, nA, nB, nC, nD, psiB, psiC, psiD) {
    #s = 40; nA = 50; nB = 10; nC =  5; nD = 35; psiB =  2.0; psiC =  0.5; psiD = 1.2
    
    # Vectorized and optimized version of hyper4_condstats
    
    # Handle edge case where s is larger than the sum of group sizes
    if (s > (nA + nB + nC + nD)) {
      return(list(mean = c(A = NaN, B = NaN, C = NaN, D = NaN),
                  var = matrix(NaN, nrow = 4, ncol = 4,
                               dimnames = list(c("A", "B", "C", "D"), c("A", "B", "C", "D")))))
    }
    
    # Pre-calculate log-combinations
    a_vals <- 0:min(nA, s)
    b_vals <- 0:min(nB, s)
    c_vals <- 0:min(nC, s)
    d_vals <- 0:min(nD, s)
    
    
    # Create all possible combinations
    combos <- data.table::CJ(a = a_vals, b = b_vals, c = c_vals, d = d_vals)
    # substantially faster than expand.grid()
    
    # Filter combinations based on the sum constraint
    combos <- combos[combos$a + combos$b + combos$c + combos$d == s, ]
    
    
    if (nrow(combos) == 0) {
      return(list(mean = c(A = NaN, B = NaN, C = NaN, D = NaN),
                  var = matrix(NaN, nrow = 4, ncol = 4,
                               dimnames = list(c("A", "B", "C", "D"), c("A", "B", "C", "D")))))
    }
    
    # Pre-calculate log-choose values for the valid ranges
    lchoose_a <- lchoose(nA, combos$a)
    lchoose_b <- lchoose(nB, combos$b)
    lchoose_c <- lchoose(nC, combos$c)
    lchoose_d <- lchoose(nD, combos$d)
    
    
    # Calculate log-numerators
    lnumer <- lchoose_a + lchoose_b + lchoose_c + lchoose_d +
      combos$b * log(psiB) + combos$c * log(psiC) + combos$d * log(psiD)
    
    # Log-sum-exp trick
    max_lnumer <- max(lnumer)
    sum_logw <- max_lnumer + log(sum(exp(lnumer - max_lnumer)))
    denom_log <- sum_logw
    
    # Calculate weights (probabilities)
    lw <- lnumer - denom_log
    w <- exp(lw)
    
    # Calculate means (vectorized)
    EA <- sum(combos$a * w)
    EB <- sum(combos$b * w)
    EC <- sum(combos$c * w)
    ED <- sum(combos$d * w)
    
    # Calculate second moments (vectorized)
    EA2 <- sum(combos$a^2 * w)
    EB2 <- sum(combos$b^2 * w)
    EC2 <- sum(combos$c^2 * w)
    ED2 <- sum(combos$d^2 * w)
    
    # Calculate cross-product moments (vectorized)
    EAB <- sum(combos$a * combos$b * w)
    EAC <- sum(combos$a * combos$c * w)
    EAD <- sum(combos$a * combos$d * w)
    EBC <- sum(combos$b * combos$c * w)
    EBD <- sum(combos$b * combos$d * w)
    ECD <- sum(combos$c * combos$d * w)
    
    # Var & Cov
    varA <- EA2 - EA^2
    varB <- EB2 - EB^2
    varC <- EC2 - EC^2
    varD <- ED2 - ED^2
    
    covAB <- EAB - EA * EB
    covAC <- EAC - EA * EC
    covAD <- EAD - EA * ED
    covBC <- EBC - EB * EC
    covBD <- EBD - EB * ED
    covCD <- ECD - EC * ED
    
    V <- matrix(c(
      varA, covAB, covAC, covAD,
      covAB, varB, covBC, covBD,
      covAC, covBC, varC, covCD,
      covAD, covBD, covCD, varD
    ), nrow = 4, byrow = TRUE,
    dimnames = list(c("A", "B", "C", "D"), c("A", "B", "C", "D")))
    
    list(
      mean = c(A = EA, B = EB, C = EC, D = ED),
      var = V
    )
  }
  # 1) Parse formula & data
  mf <- model.frame(formula, data)
  if (ncol(mf)!=2L) {
    stop("Need exactly one predictor + one ordinal response.")
  }
  resp_var <- names(mf)[1]
  grp_var  <- names(mf)[2]
  
  # ensure response is ordered factor
  if (!is.factor(mf[[resp_var]])) {
    mf[[resp_var]] <- factor(mf[[resp_var]], ordered=TRUE)
  }
  if (!is.ordered(mf[[resp_var]])) {
    mf[[resp_var]] <- ordered(mf[[resp_var]])
  }
  
  # ensure group is a factor w/ exactly 4 levels
  if (!is.factor(mf[[grp_var]])) {
    mf[[grp_var]] <- factor(mf[[grp_var]])
  }
  if (nlevels(mf[[grp_var]])!=4L) {
    stop("This function requires exactly 4 levels in the predictor factor.")
  }
  if (!is.null(ref)) {
    if (!(ref %in% levels(mf[[grp_var]]))) {
      stop("ref must be one of the factor levels.")
    }
    mf[[grp_var]] <- relevel(mf[[grp_var]], ref=ref)
  }
  g_levels <- levels(mf[[grp_var]])  # e.g. c("A","B","C","D")
  
  # 2) Aggregate counts
  agg <- mf %>%
    group_by(
      !!rlang::sym(resp_var),
      !!rlang::sym(grp_var)
    ) %>%
    summarise(count=n(), .groups="drop")
  
  wide <- agg %>%
    pivot_wider(
      names_from=!!rlang::sym(grp_var),
      values_from="count",
      values_fill=0
    ) %>%
    arrange(!!rlang::sym(resp_var))
  
  # rename => N_<level>
  for (lvl in g_levels) {
    old_nm <- lvl
    new_nm <- paste0("N_", lvl)
    wide <- rename(wide, !!new_nm := all_of(old_nm))
  }
  
  # 3) compute cumulative sums from top to bottom
  wide_rev <- wide %>% arrange(desc(!!rlang::sym(resp_var)))
  wide_rev <- wide_rev %>%
    mutate(
      cA = cumsum(.data[[paste0("N_", g_levels[1])]]),
      cB = cumsum(.data[[paste0("N_", g_levels[2])]]),
      cC = cumsum(.data[[paste0("N_", g_levels[3])]]),
      cD = cumsum(.data[[paste0("N_", g_levels[4])]])
    ) %>%
    arrange(!!rlang::sym(resp_var))
  
  k <- nrow(wide_rev)
  if (k<2) stop("Need at least 2 ordinal categories.")
  
  cdf <- wide_rev[1:(k-1), ] %>%
    mutate(
      S_A=cA,
      S_B=cB,
      S_C=cC,
      S_D=cD,
      S_tot=S_A+S_B+S_C+S_D
    )
  
  # group totals
  nA <- sum(wide[[paste0("N_", g_levels[1])]])
  nB <- sum(wide[[paste0("N_", g_levels[2])]])
  nC <- sum(wide[[paste0("N_", g_levels[3])]])
  nD <- sum(wide[[paste0("N_", g_levels[4])]])
  
  # We'll do 3 unknown parameters => (bB, bC, bD)
  beta_curr <- start
  
  # 4) define the unweighted score function:
  # U1= sum_j[S_B - mu_B], U2= sum_j[S_C - mu_C], U3= sum_j[S_D - mu_D].
  # We'll do a numeric derivative approach just like before.
  
  # helper => sum residuals in each dimension
  score_func <- function(bvec) {
    bB <- bvec[1]
    bC <- bvec[2]
    bD <- bvec[3]
    psiB <- exp(bB)
    psiC <- exp(bC)
    psiD <- exp(bD)
    
    sum_residB <- 0
    sum_residC <- 0
    sum_residD <- 0
    
    for (i in seq_len(nrow(cdf))) {
      s <- cdf$S_tot[i]
      sB<- cdf$S_B[i]
      sC<- cdf$S_C[i]
      sD<- cdf$S_D[i]
      
      st <- hyper4_condstats_fast(
        s, nA, nB, nC, nD,
        psiB, psiC, psiD
      )
      muB <- st$mean["B"]
      muC <- st$mean["C"]
      muD <- st$mean["D"]
      
      sum_residB <- sum_residB + (sB - muB)
      sum_residC <- sum_residC + (sC - muC)
      sum_residD <- sum_residD + (sD - muD)
    }
    
    c(sum_residB, sum_residC, sum_residD)
  }
  
  
  # We'll use numDeriv for the Jacobian
  
  for (it in seq_len(maxit)) {
    U <- score_func(beta_curr)
    J <- jacobian(score_func, beta_curr, method="simple")  # 3x3
    
    step <- tryCatch(solve(J, U), error=function(e) rep(0,3))
    beta_new <- beta_curr - step
    diff_val <- sqrt(sum((beta_new - beta_curr)^2))
    if (trace) {
      cat(sprintf(
        "Iter %d: bB=%.4f, bC=%.4f, bD=%.4f,  U=(%.2f,%.2f,%.2f), diff=%.2e\n",
        it, beta_new[1], beta_new[2], beta_new[3],
        U[1], U[2], U[3], diff_val
      ))
    }
    beta_curr <- beta_new
    if (diff_val<tol) break
  }
  
  # final => name them
  names(beta_curr) <- c(
    paste0("logOR(",g_levels[2]," vs ",g_levels[1],")"),
    paste0("logOR(",g_levels[3]," vs ",g_levels[1],")"),
    paste0("logOR(",g_levels[4]," vs ",g_levels[1],")")
  )
  list(
    estimate   = beta_curr,
    #std_error  = std_errors,  # Added
    iterations = it#,
    #varcov     = varcov_matrix 
  )
}
