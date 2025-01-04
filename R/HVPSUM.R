#' HVPSUM: Adjusting Genetic Correlation for Horizontal and Vertical Pleiotropy
#'
#' This function applies the HVP model to disentangle horizontal pleiotropy effects from vertical pleiotropy using GWAS summary statistics. It adjusts genetic correlation estimates obtained from cross-trait LD-score regression (LDSC).
#'
#' @param ldsc_env_path Character string. The path to the environment where LDSC is installed.
#' @param ldsc_exe_path Character string. The path to the LDSC executable file.
#' @param ld_path Character string. The path to the reference LD blocks.
#' @param sumstats_outcome Character string. The file path to the munged outcome summary statistics (e.g., a `.sumstats.gz` file).
#' @param sumstats_exposure Character string. The file path to the munged exposure summary statistics (e.g., a `.sumstats.gz` file).
#' @param h2_outcome Data frame or vector. Precomputed heritability block estimates from jackknife resampling for the outcome trait. 
#'   If `NULL`, these estimates will be computed using LDSC if the required LDSC paths are provided.
#' @param h2_exposure Data frame or vector. Precomputed heritability block estimates from jackknife resampling for the exposure trait. 
#'   If `NULL`, these estimates will be computed using LDSC if the required LDSC paths are provided.
#' @param gencov Data frame or vector. Precomputed genetic covariance block estimates from jackknife resampling. 
#'   If `NULL`, these estimates will be computed using LDSC if the required LDSC paths are provided.
#' @param tau Numeric. The causal effect estimate from a Mendelian Randomization (MR) method, used to adjust the genetic covariance and heritability estimates.
#' @param intern Logical. If `TRUE`, the LDSC command will be executed internally.
#' @param remove_temp_files Logical. Whether to delete temporary files after the analysis (default is TRUE).
#' 
#' @return A list containing the following:
#'   \item{means}{A named vector of means for heritability, genetic covariance, and correlation estimates.}
#'   \item{se}{A named vector of standard errors for the estimates.}
#'   \item{pvalues}{A named vector of p-values for the estimates.}
#' 
#' @importFrom stats var pchisq
#' @importFrom utils read.table
#' @export
hvpsum <- function(ldsc_env_path = NULL, ldsc_exe_path = NULL, ld_path = NULL, sumstats_outcome = NULL, sumstats_exposure = NULL,
                   h2_outcome = NULL, h2_exposure = NULL, gencov = NULL,
                   tau, intern = TRUE, remove_temp_files = TRUE ) {
  
  # Helper Function: Run LDSC
  run_ldsc <- function(ldsc_env_path, ldsc_exe_path, ld_path, sumstats_outcome, sumstats_exposure) {
    
    ldsc_command <- paste(
      ldsc_env_path, ldsc_exe_path,  # Use the path to the environment and LDSC script
      "--rg", paste0(basename(sumstats_outcome), ",", basename(sumstats_exposure)),  # Paths to the summary statistics files
      "--ref-ld-chr", ld_path, # Path to reference LD blocks
      "--w-ld-chr", ld_path, # Path to weights LD blocks
      "--print-delete-vals", 
      "--out", paste0(gsub(".sumstats.gz", "", basename(sumstats_outcome)), "_", gsub(".sumstats.gz", "", basename(sumstats_exposure)))
    )
    system(ldsc_command,intern = intern )
    
    base_name <- (paste0(gsub(".sumstats.gz", "", sumstats_outcome), "_", gsub(".sumstats.gz", "", sumstats_exposure),sumstats_outcome,"_",sumstats_exposure))
    
    h2_outcome_path <- file.path(getwd(),paste0(base_name, ".hsq1.delete"))
    h2_exposure_path <- file.path(getwd(),paste0(base_name, ".hsq2.delete"))
    gencov_path <- file.path(getwd(),paste0(base_name, ".gencov.delete"))
    
    return(list(h2_outcome_path = h2_outcome_path, h2_exposure_path = h2_exposure_path, gencov_path = gencov_path))
  }
  
  remove_files <- function(file_paths) {
    unlink(file_paths$h2_outcome_path)
    unlink(file_paths$h2_exposure_path)  
    unlink(file_paths$gencov_path)
    message("Temporary files deleted.")
  }
  downstream_analysis <- function(h2_outcome, h2_exposure, gencov, tau) {
    # Validate inputs
    if (is.null(h2_outcome) || is.null(h2_exposure) || is.null(gencov)) {
      stop("Precomputed data not provided. Please supply the necessary data.")
    }
    
    # Data preparation and calculations
    data <- data.frame(
      hsq.1 = as.numeric(h2_outcome$V1), 
      hsq.2 = as.numeric(h2_exposure$V1), 
      gencov = as.numeric(gencov$V1)
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
  
  file_paths <- NULL
  if (!is.null(ldsc_env_path) && !is.null(ldsc_exe_path) && !is.null(ld_path)) {
    message("Running LDSC analysis...")
    Sys.setenv(PYTHONWARNINGS = "ignore::FutureWarning")
    file_paths <- run_ldsc(ldsc_env_path, ldsc_exe_path, ld_path, sumstats_outcome, sumstats_exposure)
    
    # Validate intermediate files
    if (!file.exists(file_paths$h2_outcome_path) || !file.exists(file_paths$h2_exposure_path) || !file.exists(file_paths$gencov_path)) {
      stop("Intermediate files from LDSC analysis not found.")
    }
    
    h2_outcome <- read.table(file_paths$h2_outcome_path)
    h2_exposure <- read.table(file_paths$h2_exposure_path)
    gencov <- read.table(file_paths$gencov_path)
  } else if (is.null(h2_outcome) || is.null(h2_exposure) || is.null(gencov)) {
    stop("Precomputed data or paths for LDSC analysis are missing.")
  }
  
  # Perform downstream analysis
  downstream_results <- downstream_analysis(h2_outcome, h2_exposure, gencov, tau)
  file_paths<- file_paths
  cat(sprintf("%-32s: Estimate = %.4f, SE = %.4f, P-value = %s\n", 
              "h2 outcome", downstream_results$means["hsq.1"], 
              downstream_results$se["hsq.1"], 
              format.pval(downstream_results$pvalues["hsq.1"])))
  cat(sprintf("%-32s: Estimate = %.4f, SE = %.4f, P-value = %s\n", 
              "h2 outcome (corrected)", downstream_results$means["hsq1c_corrected"], 
              downstream_results$se["hsq1c_corrected"], 
              format.pval(downstream_results$pvalues["hsq1c_corrected"])))
  cat(sprintf("%-32s: Estimate = %.4f, SE = %.4f, P-value = %s\n", 
              "h2 exposure", downstream_results$means["hsq.2"], 
              downstream_results$se["hsq.2"], 
              format.pval(downstream_results$pvalues["hsq.2"])))
  cat(sprintf("%-32s: Estimate = %.4f, SE = %.4f, P-value = %s\n", 
              "Genetic Covariance", downstream_results$means["gencov"], 
              downstream_results$se["gencov"], 
              format.pval(downstream_results$pvalues["gencov"])))
  cat(sprintf("%-32s: Estimate = %.4f, SE = %.4f, P-value = %s\n", 
              "Genetic Covariance (corrected)", downstream_results$means["gencov_corrected"], 
              downstream_results$se["gencov_corrected"], 
              format.pval(downstream_results$pvalues["gencov_corrected"])))
  cat(sprintf("%-32s: Estimate = %.4f, SE = %.4f, P-value = %s\n", 
              "Genetic Correlation", downstream_results$means["cor"], 
              downstream_results$se["cor"], 
              format.pval(downstream_results$pvalues["cor"])))
  cat(sprintf("%-32s: Estimate = %.4f, SE = %.4f, P-value = %s\n", 
              "Genetic Correlation (corrected)", downstream_results$means["cor_corrected"], 
              downstream_results$se["cor_corrected"], 
              format.pval(downstream_results$pvalues["cor_corrected"])))

  if (remove_temp_files && !is.null(file_paths)) {
    remove_files(file_paths)
  }
}