#' HVPSUM: Genetic Correlation
#'
#' This function applies the HVP model to disentangle horizontal pleiotropy effects from vertical pleiotropy using GWAS summary data. It corrects genetic correlation estimates derived from cross-trait LD-score regression (LDSC).
#'
#' @param ldsc_env_path Path to the LDSC environment.
#' @param ldsc_exe_path Path to the LDSC executable script.
#' @param ld_path Path to the reference LD files.
#' @param sumstats1 Path to the first trait munged summary statistics file (outcome).
#' @param sumstats2 Path to the second trait munged summary statistics file(exposure).
#' @param h21block_path Path to the pre-computed jackknife blocks estimate  for the heritability of trait 1.
#' @param h22block_path Path to the pre-computed jackknife blocks estimate for the heritability of the trait 2.
#' @param gen_cov_path Path to the pre-computed jackknife blocks estimate for the genetic covariance between the traits.
#' @param tau Causal effect of the exposure on the outcome, as estimated by a robust Mendelian Randomization (MR) method.
#' @param intern Logical; if \code{TRUE}, captures the LDSC output internally when running the command. 
#'               Defaults to \code{TRUE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{Estimates}{Estimates for heritabilities, genetic covariance, and genetic correlation (corrected and uncorrected).}
#'   \item{SE}{Standard errors for each estimate.}
#'   \item{p-values}{P-values associated with each estimate.}
#' }
#'
#' @details 
#' The function operates in two modes:
#' 
#' \strong{1. LDSC Mode}: If \code{ldsc_env_path}, \code{ldsc_exe_path}, \code{ld_path}, \code{sumstats1}, and \code{sumstats2} are provided, 
#' the function runs LDSC analysis to compute jackknife blocks estimates for heritability and genetic covariance for the input summary statistics.
#' 
#' \strong{2. Pre-computed Mode}: If paths to jackknife blocks estimates (\code{h21block_path}, \code{h22block_path}, \code{gen_cov_path}) 
#' are provided, the function performs downstream analysis without running LDSC.
#' 
#' The downstream analysis computes the heritabilities, genetic covariance, and genetic correlation 
#' with and without corrections, along with their standard errors and p-values using pseudovalues.
#'
#' @examples
#' # Example using LDSC Mode
#' #hvpsum(ldsc_env_path = "path/to/env//ldsc/bin/python", ldsc_exe_path = "path/to/ldsc.py", 
#' #      ld_path = "path/to/ld", sumstats1 = "path/to/file1.sumstats.gz", 
#' #       sumstats2 = "path/to/file2.sumstats.gz", tau = 0.5, intern = TRUE)
#' 
#' # Example using Precomputed Mode
#' #hvpsum(h21block_path = "path/to/file1.hsq1.delete", h22block_path = "path/to/file2.hsq2.delete", 
#' #       gen_cov_path = "path/to/file3.gencov.delete", tau = 0.5, intern = FALSE)
#'
#' @importFrom stats var pchisq
#' @importFrom utils read.table
#' @export

hvpsum <- function(ldsc_env_path = NULL, ldsc_exe_path = NULL, ld_path = NULL,sumstats1=NULL, sumstats2=NULL,
h21block_path = NULL, h22block_path = NULL, gen_cov_path = NULL,
tau, intern = TRUE) {
  
 # Helper Function: Run LDSC
  run_ldsc <- function(ldsc_env_path, ldsc_exe_path, ld_path, sumstats1, sumstats2) {
    #   # Construct and execute LDSC command
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
    list(h21_path, h22_path, gencov_path)
    
    }

  downstream_analysis <- function(h21_path , h22_path, gencov_path, tau) {
   
    dat1<- read.table(h21_path)$V1
    dat2 <- read.table(h22_path)$V1
    dat3 <- read.table(gencov_path)$V1
    
   
    data <- data.frame("hsq.1"=dat1, "hsq.2"=dat2, "gencov"=dat3)
   
    data$gencov_corrected <- data$gencov - tau * data$hsq.2
    data$hsq1c_corrected <- data$hsq.1 - tau^2 * data$hsq.2 - 2 * tau * data$gencov_corrected
    data$cor <- data$gencov / sqrt(data$hsq.1 * data$hsq.2)
    data$cor_corrected <- data$gencov_corrected / sqrt(data$hsq1c_corrected * data$hsq.2)
    
   
    blocks_est <- as.matrix(data)
    blocks_mean <- apply(blocks_est, 2, mean)
    pseudo_values <- function(blocks_mean, blocks_est) {
      n_blocks <- nrow(blocks_est)
      p <- ncol(blocks_est)
      pseudovalues <- matrix(0, n_blocks, p)
      for (j in 1:n_blocks) {
        for (k in 1:p) {
          pseudovalues[j, k] <- n_blocks * blocks_mean[k] - (n_blocks - 1) * blocks_est[j, k]
        }
      }
      
      return(pseudovalues)
    }
    standard_error <- function(pseudovalues) {
      n_blocks <- nrow(pseudovalues)
      
      variance <- apply(pseudovalues, 2, var)
     
      standard_error <- sqrt(variance / n_blocks)
      return(standard_error)
    }
    pseudovalues <- pseudo_values(blocks_mean, blocks_est)
    se=standard_error(pseudovalues)
    names(se)=names(blocks_mean)
   
    pvalues <- sapply(seq_along(blocks_mean), function(i) {
      z_score <- blocks_mean[i] / se[i]
      pchisq((z_score)^2, 1,lower.tail=F)
    })
    
    
    list(
      means = blocks_mean,
      se = se,
      pvalues = pvalues
    )   
  }
  
  
  if (!is.null(ldsc_env_path) && !is.null(ldsc_exe_path) && !is.null(ld_path)) {
    
    message("Running LDSC analysis...")
    Sys.setenv(PYTHONWARNINGS = "ignore::FutureWarning")
    ldsc_results <- run_ldsc(ldsc_env_path, ldsc_exe_path, ld_path, sumstats1, sumstats2)
    
    h21block_path <- ldsc_results[[1]]
    h22block_path <- ldsc_results[[2]]
    gen_cov_path <- ldsc_results[[3]]
    
  } else if (!is.null(h21block_path) && !is.null(h22block_path) && !is.null(gen_cov_path)) {
    message("Using precomputed jacknife estimate for downstream analysis...")
  } else {
    stop("Please provide either Set 1 (LDSC inputs) or Set 2 (jacknife estimate paths).")
  }
  
  
  results <- downstream_analysis(h21block_path, h22block_path, gen_cov_path, tau)
  
 
  cat("Heritability 1: Estimate =", round(results$means["hsq.1"], 4),
      "SE =", round(results$se["hsq.1"], 4),
      "P-value =", format.pval(results$pvalues["hsq.1"]), "\n")
  cat("Heritability 1 Corrected: Estimate =", round(results$means["hsq1c_corrected"], 4),
      "SE =", round(results$se["hsq1c_corrected"], 4),
      "P-value =", format.pval(results$pvalues["hsq1c_corrected"]), "\n")
  cat("Heritability 2: Estimate =", round(results$means["hsq.2"], 4),
      "SE =", round(results$se["hsq.2"], 4),
      "P-value =", format.pval(results$pvalues["hsq.2"]), "\n")
  cat("Genetic Covariance: Estimate =", round(results$means["gencov"], 4),
      "SE =", round(results$se["gencov"], 4),
      "P-value =", format.pval(results$pvalues["gencov"]), "\n")
  cat("Genetic Covariance Corrected: Estimate =", round(results$means["gencov_corrected"], 4),
      "SE =", round(results$se["gencov_corrected"], 4),
      "P-value =", format.pval(results$pvalues["gencov_corrected"]), "\n")
  cat("Genetic Correlation: Estimate =", round(results$means["cor"], 4),
      "SE =", round(results$se["cor"], 4),
      "P-value =", format.pval(results$pvalues["cor"]), "\n")
  cat("Genetic Correlation: Estimate =", round(results$means["cor_corrected"], 4),
      "SE =", round(results$se["cor_corrected"], 4),
      "P-value =", format.pval(results$pvalues["cor_corrected"]), "\n")
  
   invisible(results)
}