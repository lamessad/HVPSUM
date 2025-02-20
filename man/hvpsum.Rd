% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HVPSUM.R
\name{hvpsum}
\alias{hvpsum}
\title{HVPSUM: Adjusting Genetic Correlation for Horizontal and Vertical Pleiotropy}
\usage{
hvpsum(
  ldsc_env_path = NULL,
  ldsc_exe_path = NULL,
  ld_path = NULL,
  sumstats_outcome = NULL,
  sumstats_exposure = NULL,
  h2_outcome = NULL,
  h2_exposure = NULL,
  gencov = NULL,
  tau,
  intern = TRUE,
  remove_temp_files = TRUE
)
}
\arguments{
\item{ldsc_env_path}{Character string. The path to the environment where LDSC is installed.}

\item{ldsc_exe_path}{Character string. The path to the LDSC executable file.}

\item{ld_path}{Character string. The path to the reference LD blocks.}

\item{sumstats_outcome}{Character string. The file path to the munged outcome summary statistics (e.g., a \code{.sumstats.gz} file).}

\item{sumstats_exposure}{Character string. The file path to the munged exposure summary statistics (e.g., a \code{.sumstats.gz} file).}

\item{h2_outcome}{Data frame or vector. Precomputed heritability block estimates from jackknife resampling for the outcome trait.
If \code{NULL}, these estimates will be computed using LDSC if the required LDSC paths are provided.}

\item{h2_exposure}{Data frame or vector. Precomputed heritability block estimates from jackknife resampling for the exposure trait.
If \code{NULL}, these estimates will be computed using LDSC if the required LDSC paths are provided.}

\item{gencov}{Data frame or vector. Precomputed genetic covariance block estimates from jackknife resampling.
If \code{NULL}, these estimates will be computed using LDSC if the required LDSC paths are provided.}

\item{tau}{Numeric. The causal effect estimate from a Mendelian Randomization (MR) method, used to adjust the genetic covariance and heritability estimates.}

\item{intern}{Logical. If \code{TRUE}, the LDSC command will be executed internally.}

\item{remove_temp_files}{Logical. Whether to delete temporary files after the analysis (default is TRUE).}
}
\value{
A list containing the following:
\item{means}{A named vector of means for heritability, genetic covariance, and correlation estimates.}
\item{se}{A named vector of standard errors for the estimates.}
\item{pvalues}{A named vector of p-values for the estimates.}
}
\description{
This function applies the HVP model to disentangle horizontal pleiotropy effects from vertical pleiotropy using GWAS summary statistics. It adjusts genetic correlation estimates obtained from cross-trait LD-score regression (LDSC).
}
