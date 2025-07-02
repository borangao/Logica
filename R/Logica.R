#' @useDynLib Logica, .registration = TRUE
#' @import Rcpp
#' @import RcppArmadillo
#' @import susieR
#' @import dplyr
#' @import CompQuadForm
#' @import HDMT
NULL

# Reliance on susieR , CompQuadForm
#' Preprocess Data for Logica
#'
#' This function preprocesses data by removing SNPs based on diagnostic criteria and aligning SNPs across populations.
#'
#' @param z1 Data frame or matrix for population 1.
#' @param z2 Data frame or matrix for population 2.
#' @param R1 LD matrix for population 1.
#' @param R2 LD matrix for population 2.
#'
#' @return A list containing updated data and LD matrices.
#' @export
preprocess_data <- function(z1, z2, R1, R2) {
  # Function to iteratively remove SNPs based on diagnostic criteria
  remove_snps <- function(z, R) {
    rm_id <- 1
    while (length(rm_id) != 0) {
      diagnostic <- kriging_rss(z$Z, R)
      rm_id <- which(diagnostic$conditional_dist$logLR > 2)
      if (length(rm_id) != 0) {
        z <- z[-rm_id, ]
        R <- R[-rm_id, -rm_id]
      }
    }
    list(z = z, R = R)
  }
  
  # Process population 1
  result_pop1 <- remove_snps(z1, R1)
  z1 <- result_pop1$z
  R1 <- result_pop1$R
  
  # Process population 2
  result_pop2 <- remove_snps(z2, R2)
  z2 <- result_pop2$z
  R2 <- result_pop2$R
  
  # Find common SNPs
  common_SNP <- intersect(z1$SNP, z2$SNP)
  z1_update <- z1 %>% filter(SNP %in% common_SNP) %>% arrange(match(SNP, common_SNP))
  z2_update <- z2 %>% filter(SNP %in% common_SNP) %>% arrange(match(SNP, common_SNP))
  
  # Align covariance matrices
  R1_update <- R1[match(common_SNP, colnames(R1)), match(common_SNP, colnames(R1))]
  R2_update <- R2[match(common_SNP, colnames(R2)), match(common_SNP, colnames(R2))]
  
  # Check SNP alignment
  if (all(z1_update$SNP == z2_update$SNP) &&
      all(z1_update$SNP == colnames(R1_update)) &&
      all(z1_update$SNP == colnames(R2_update))) {
    cat("SNP pair up successful\n")
  } else {
    cat("SNP order mismatch\n")
  }
  
  # LD shrinkage estimation
  lambda1 <- estimate_s_rss(z1_update$Z, R1_update, r_tol = 1e-08, method = "null-mle")
  lambda2 <- estimate_s_rss(z2_update$Z, R2_update, r_tol = 1e-08, method = "null-mle")
  
  R1_update <- R1_update * (1 - lambda1) + diag(lambda1, nrow(R1_update))
  R2_update <- R2_update * (1 - lambda2) + diag(lambda2, nrow(R2_update))
  
  list(
    pop1 = list(z = z1_update, R = R1_update),
    pop2 = list(z = z2_update, R = R2_update)
  )
}

#' Run Logica Analysis
#'
#' This function runs the Logica pipeline using summary statistics and LD matrices.
#'
#' @param sumstat_1 Summary statistics for population 1.
#' @param sumstat_2 Summary statistics for population 2.
#' @param R1 LD matrix for population 1.
#' @param R2 LD matrix for population 2.
#' @param n1 Sample size for population 1.
#' @param n2 Sample size for population 2.
#' @param ... Additional parameters.
#'
#' @return A data frame with results.
#' @export
run_Logica <- function(sumstat_1, sumstat_2, R1, R2, n1, n2, z1_intercept = 1, z2_intercept = 1, fix_intercept = TRUE, initial_method = "Uni",screen = FALSE) {
  
  t0 <- Sys.time()

  chr<-unique(sumstat_1$CHROM)
  start<-min(sumstat_1$POS)
  end<-max(sumstat_1$POS)
  
  convert_size <- function(size_bytes) {
    if (size_bytes >= 1024^2) {
      size <- round(size_bytes / 1024^2, 2)
      unit <- "MB"
    } else if (size_bytes >= 1024) {
      size <- round(size_bytes / 1024, 2)
      unit <- "KB"
    } else {
      size <- size_bytes
      unit <- "bytes"
    }
    paste0(size, " ", unit)
  }
  
  start_readable <- convert_size(start)
  end_readable <- convert_size(end)
  
  
  
  z1<-sumstat_1$Z
  z2<-sumstat_2$Z
  num_SNP <- length(z1)
  # Estimate initial heritabilities
  if(initial_method=="MOM"){
    h2_1_initial = min(max(mean(z1^2-z1_intercept)/mean(apply(R1,2,function(x)sum(x^2))),10^-40),0.1)
    h2_2_initial = min(max(mean(z2^2-z2_intercept)/mean(apply(R2,2,function(x)sum(x^2))),10^-40),0.1)
  }else{
  
  h2_1_initial <- optim(10^-20, llk_h2_est, method = "Brent", sigma = z1_intercept, z = z1,
                        R = R1, n = n1, num_SNP = num_SNP,
                        lower = 10^-40, upper = 10^-1)$par
  
  h2_2_initial <- optim(10^-20, llk_h2_est, method = "Brent", sigma = z2_intercept, z = z2,
                        R = R2, n = n2, num_SNP = num_SNP,
                        lower = 10^-40, upper = 10^-1)$par
  
  }
  # Estimate null score and heritabilities
  est_score <- est_h2_null_score(
    h2_1_initial, h2_2_initial, z1_intercept, z2_intercept, 
    z_1 = z1, z_2 = z2,
    R_1 = R1, R_2 = R2,
    n_1 = n1, n_2 = n2, num_SNP = num_SNP, max_iter = 1000, fix_intercept = fix_intercept
  )
  
  # Compute Davies p-values for heritability tests
  davies_score_h2_1 <- davies(est_score$h2_1_score, est_score$h2_1_eigvals, acc = 1e-10)$Qq
  davies_score_h2_2 <- davies(est_score$h2_2_score, est_score$h2_2_eigvals, acc = 1e-10)$Qq
  davies_score_rho <- 2 * davies(abs(est_score$rho_score), sort(est_score$rho_eigvals), acc = 1e-10)$Qq
  
  t1 <- Sys.time()
  step1_mins <- as.numeric(difftime(t1, t0, units = "mins"))
  cat(sprintf("Step 1 Screening heritable region: %.3f minutes\n", step1_mins))

  if (screen) {
    return(data.frame(
      chr        = as.character(chr),
      start      = start_readable,
      end        = end_readable,
      Est_h2_1   = est_score$est[1],
      Est_h2_2   = est_score$est[2],
      P_h2_1     = davies_score_h2_1,
      P_h2_2     = davies_score_h2_2,
      P_rho_score= davies_score_rho,
       Time_step1  = step1_mins
    ))
  }

  # Initial rho estimate
  rho_initial <- if(min(est_score$est[1], est_score$est[2]) < 1e-6) {
    0
  } else {
    optim(0, llk_h2_rho, method = "Brent",
          h2_1 = est_score$est[1], h2_2 = est_score$est[2], sigma_1 = est_score$est[3], sigma_2 = est_score$est[4],
          n_1 = n1, n_2 = n2, num_SNP = num_SNP,
          z_1 = z1, z_2 = z2, R_1 = R1, R_2 = R2,
          lower = -1, upper = 1)$par
  }
  
  # Final estimation via PX-EM
  est_rho <- PX_EM_alt(
    est_score$est[1], est_score$est[2], rho_initial,
    est_score$est[3], est_score$est[4],
    z_1 = z1, z_2 = z2,
    R_1 = R1, R_2 = R2,
    n_1 = n1, n_2 = n2,
    num_SNP = num_SNP, n_iter = 100,
    fix_intercept = fix_intercept, fix_h2 = est_score$fix_h2
  )
  
  # Likelihood ratio test
  llk_test <- pchisq(2 * (est_rho$llk_alt - est_score$llk_null), df = 1, lower.tail = FALSE)
  
  t2 <- Sys.time()
  step2_mins <- as.numeric(difftime(t2, t1, units = "mins"))
  cat(sprintf("Step 2 Estimation of genetic correlation:    %.3f minutes\n", step2_mins))

  # 3. Total
  total_mins <- as.numeric(difftime(t2, t0, units = "mins"))
  cat(sprintf("Total (start → end):              %.3f minutes\n\n", total_mins))


  # Output results as named vector
  results <- data.frame(
    chr = as.character(chr),
    start = start_readable,
    end = end_readable,
    Est_h2_1 = est_rho$est[1],
    Est_h2_2 = est_rho$est[2],
    Est_cov = est_rho$est[3],
    Est_rho = est_rho$est[3]/sqrt(est_rho$est[1] * est_rho$est[2]),
    P_h2_1 = davies_score_h2_1,
    P_h2_2 = davies_score_h2_2,
    P_rho_score = davies_score_rho,
    P_rho_LLK = llk_test,
    Time_step1  = step1_mins,
    Time_step2  = step2_mins,
    Time_total  = total_mins
  )
  
  return(results)
}

#' Composite Null Hypothesis Testing
#'
#' This function runs the Logica pipeline using summary statistics and LD matrices.
#'
#' @param sumstat_1 Summary statistics for population 1.
#' @param sumstat_2 Summary statistics for population 2.
#' @param R1 LD matrix for population 1.
#' @param R2 LD matrix for population 2.
#' @param n1 Sample size for population 1.
#' @param n2 Sample size for population 2.
#' @param ... Additional parameters.
#'
#' @return A data frame with results.
#' @export
p_HDMT<-function(p1,p2){
  input_pvalues = as.matrix(cbind(p1,p2))
  nullprop<-null_estimation(input_pvalues)
  pmax_value<-pmax(input_pvalues[,1],input_pvalues[,2])
  p_null_cons<-nullprop$alpha10+nullprop$alpha01+nullprop$alpha00
  p_HDMT<-(nullprop$alpha10+nullprop$alpha01)/p_null_cons*pmax_value+nullprop$alpha00/p_null_cons*pmax_value^2
  return(p_HDMT)
}
 