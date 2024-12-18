#' Jackknife Delete Values from LDSC for heritability of trait 1
#'
#' This dataset contains the jackknife "delete values" derived from Linkage Disequilibrium Score Regression (LDSC).
#' These values represent the block-specific estimates heritability obtained when a particular block is excluded from the analysis (i.e., leave-one-block-out).
#'
#' In LDSC, the genome is divided into multiple blocks to account for correlations among nearby SNPs due to linkage disequilibrium (LD).
#' The jackknife procedure recalculates the parameter of interest by leaving out one block at a time. These block estimates are used to:
#' - Compute the jackknife mean, which is the average estimate across all blocks.
#' - Estimate the variance and standard error of the parameter by analyzing the variation among these block-specific estimates.
#'
#' @format A data frame with one column containing the block-specific parameter estimates obtained when each block is excluded:
#' \describe{
#'   \item{V1}{Numeric values of block-specific estimates.}
#' }
#' @usage data(dat1)
#' @keywords datasets
"dat1"