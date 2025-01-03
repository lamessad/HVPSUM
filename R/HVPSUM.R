#' HVPSUM: Genetic Correlation
#'
#' This function applies the HVP model to disentangle horizontal pleiotropy effects from vertical pleiotropy using GWAS summary data. It corrects genetic correlation estimates derived from cross-trait LD-score regression (LDSC).

#' @param ldsc_env_path Character string. The path to the environment where LDSC is installed.
#' @param ldsc_exe_path Character string. The path to the LDSC executable.
#' @param ld_path Character string. The path to the reference LD blocks.
#' @param sumstats1 Character string. Path to the first summary statistics file (e.g., .sumstats.gz).
#' @param sumstats2 Character string. Path to the second summary statistics file (e.g., .sumstats.gz).
#' @param h21_data Data frame or vector. Precomputed heritability estimates for the first trait. 
#'   If `NULL`, it will be computed using LDSC if the LDSC paths are provided.
#' @param h22_data Data frame or vector. Precomputed heritability estimates for the second trait.
#'   If `NULL`, it will be computed using LDSC if the LDSC paths are provided.
#' @param gencov_data Data frame or vector. Precomputed genetic covariance estimates.
#'   If `NULL`, it will be computed using LDSC if the LDSC paths are provided.
#' @param tau Numeric value. The causal effect to adjust the genetic covariance and heritability.
#' @param intern Logical. If `TRUE`, the LDSC command is run internally.
#' 
#' @return A list containing the following:
#'   \item{means}{A named vector of means for heritability, genetic covariance, and correlation estimates.}
#'   \item{se}{A named vector of standard errors for the estimates.}
#'   \item{pvalues}{A named vector of p-values for the estimates.}
#' 
#' @details This r package performs the following steps:
#'   1. If `ldsc_env_path`, `ldsc_exe_path`, and `ld_path` are provided, it runs the LDSC analysis to generate 
#'      heritability and genetic covariance estimates based on summary statistics provided.
#'   2. If precomputed data is provided for `h21_data`, `h22_data`, and `gencov_data`, these are used directly.
#'   3. Calculates and returns the estimates, standard errors, and p-values for the corrected and uncorrected parameters.
#' 
#' @importFrom stats var pchisq
#' @importFrom utils read.table
#' @export
hvpsum <- function(ldsc_env_path = NULL, ldsc_exe_path = NULL, ld_path = NULL, sumstats1 = NULL, sumstats2 = NULL,
                   h21_data = NULL, h22_data = NULL, gencov_data = NULL,
                   tau, intern = TRUE) {
  
  # Helper Function: Run LDSC
  run_ldsc <- function(ldsc_env_path, ldsc_exe_path, ld_path, sumstats1, sumstats2) {
    
    ldsc_command <- paste(
      ldsc_env_path, ldsc_exe_path,  # Use the path to the environment and LDSC script
      "--rg", paste0(basename(sumstats1), ",", basename(sumstats2)),  # Paths to the summary statistics files
      "--ref-ld-chr", ld_path, # Path to reference LD blocks
      "--w-ld-chr", ld_path, # Path to weights LD blocks
      "--print-delete-vals", 
      "--out", paste0(gsub(".sumstats.gz", "", basename(sumstats1)), "_", gsub(".sumstats.gz", "", basename(sumstats2)))
    )
    system(ldsc_command,intern = intern )
    
    base_name <- (paste0(gsub(".sumstats.gz", "", sumstats1), "_", gsub(".sumstats.gz", "", sumstats2),sumstats1,"_",sumstats2))
    
    h21_path <- file.path(getwd(),paste0(base_name, ".hsq1.delete"))
    h22_path <- file.path(getwd(),paste0(base_name, ".hsq2.delete"))
    gencov_path <- file.path(getwd(),paste0(base_name, ".gencov.delete"))
    
    file_path <- return(list(h21_path = h21_path, h22_path = h22_path, gencov_path = gencov_path))
  }
  
  # Helper Function: Downstream Analysis
  downstream_analysis <- function(h21_data, h22_data, gencov_data, tau) {
    # Validate inputs
    if (is.null(h21_data) || is.null(h22_data) || is.null(gencov_data)) {
      stop("Precomputed data not provided. Please supply the necessary data.")
    }
    
    # Data preparation and calculations
    data <- data.frame(
      hsq.1 = as.numeric(h21_data$V1), 
      hsq.2 = as.numeric(h22_data$V1), 
      gencov = as.numeric(gencov_data$V1)
    )
    
    data$gencov_corrected <- data$gencov - tau * data$hsq.2
    data$hsq1c_corrected <- data$hsq.1 - tau^2 * data$hsq.2 - 2 * tau * data$gencov_corrected
    data$cor <- data$gencov / sqrt(data$hsq.1 * data$hsq.2)
    data$cor_corrected <- data$gencov_corrected / sqrt(data$hsq1c_corrected * data$hsq.2)
    
    # Mean, pseudovalues, and SE calculations
    blocks_est <- as.matrix(data)
    blocks_mean <- apply(blocks_est, 2, mean)
    
    
    pseudovalues <- matrix(0, nrow(blocks_est), ncol(blocks_est))
    for (j in 1:nrow(blocks_est)) {
      for (k in 1:ncol(blocks_est)) {
        pseudovalues[j, k] <- nrow(blocks_est) * blocks_mean[k] - (nrow(blocks_est) - 1) * blocks_est[j, k]
      }
    }
    
    
    se <- sqrt(apply(pseudovalues, 2, var) / nrow(blocks_est))
    names(se) <- names(blocks_mean)
    pvalues <- sapply(1:length(blocks_mean), function(i) pchisq((blocks_mean[i] / se[i])^2, 1, lower.tail = FALSE))
    
    return(list(means = blocks_mean, se = se, pvalues = pvalues))
  }
  
  # Main Logic
  if (!is.null(ldsc_env_path) && !is.null(ldsc_exe_path) && !is.null(ld_path)) {
    message("Running LDSC analysis...")
    Sys.setenv(PYTHONWARNINGS = "ignore::FutureWarning")
    file_paths <- run_ldsc(ldsc_env_path, ldsc_exe_path, ld_path, sumstats1, sumstats2)
    
    # Validate intermediate files
    if (!file.exists(file_paths$h21_path) || !file.exists(file_paths$h22_path) || !file.exists(file_paths$gencov_path)) {
      stop("Intermediate files from LDSC analysis not found.")
    }
    
    h21_data <- read.table(file_paths$h21_path)
    h22_data <- read.table(file_paths$h22_path)
    gencov_data <- read.table(file_paths$gencov_path)
  } else if (is.null(h21_data) || is.null(h22_data) || is.null(gencov_data)) {
    stop("Precomputed data or paths for LDSC analysis are missing.")
  }
  
  # Perform downstream analysis
  downstream_results <- downstream_analysis(h21_data, h22_data, gencov_data, tau)
  
  # Print results
  cat(sprintf("Heritability 1: Estimate = %.4f, SE = %.4f, P-value = %s\n", 
              downstream_results$means["hsq.1"], downstream_results$se["hsq.1"], 
              format.pval(downstream_results$pvalues["hsq.1"])))
  cat(sprintf("Heritability 1 Corrected: Estimate = %.4f, SE = %.4f, P-value = %s\n", 
              downstream_results$means["hsq1c_corrected"], downstream_results$se["hsq1c_corrected"], 
              format.pval(downstream_results$pvalues["hsq1c_corrected"])))
  cat(sprintf("Heritability 2: Estimate = %.4f, SE = %.4f, P-value = %s\n", 
              downstream_results$means["hsq.2"], downstream_results$se["hsq.2"], 
              format.pval(downstream_results$pvalues["hsq.2"])))
  cat(sprintf("Genetic Covariance: Estimate = %.4f, SE = %.4f, P-value = %s\n", 
              downstream_results$means["gencov"], downstream_results$se["gencov"], 
              format.pval(downstream_results$pvalues["gencov"])))
  cat(sprintf("Genetic Covariance Corrected: Estimate = %.4f, SE = %.4f, P-value = %s\n", 
              downstream_results$means["gencov_corrected"], downstream_results$se["gencov_corrected"], 
              format.pval(downstream_results$pvalues["gencov_corrected"])))
  cat(sprintf("Genetic Correlation: Estimate = %.4f, SE = %.4f, P-value = %s\n", 
              downstream_results$means["cor"], downstream_results$se["cor"], 
              format.pval(downstream_results$pvalues["cor"])))
  cat(sprintf("Genetic Correlation Corrected: Estimate = %.4f, SE = %.4f, P-value = %s\n", 
              downstream_results$means["cor_corrected"], downstream_results$se["cor_corrected"], 
              format.pval(downstream_results$pvalues["cor_corrected"])))
}
