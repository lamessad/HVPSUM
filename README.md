
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Overview of HVPSUM

<!-- badges: start -->

<!-- badges: end -->

**HVPSUM** is an R package designed to provide unbiased estimates of
heritability and genetic correlation by disentangling horizontal and
vertical pleiotropy using summary data. Unlike Linkage Disequilibrium
Score Regression (LDSC), the HVPSUM model adjusts for vertical
pleiotropy by incorporating causal effect estimate from Mendelian
Randomization (MR). This enables a clear distinction between:

- **Shared genetic influences** (horizontal pleiotropy), and
- **Genetic effects mediated through causal pathways** (vertical
  pleiotropy).

By addressing the limitations of conventional approaches, the HVPSUM
package offers a robust framework for analyzing complex genetic
relationships.

## Installation

You can install the development version of HVPSUM from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("lamessad/HVPSUM")
```

# Steps for the Horizontal and Vertical Pleiotropy (HVP) Model based on summary data

## 1. Munge the Summary Data (see LDSC for details)

## 2. Run MR

- Estimate the causal effect size (`tau`) of the exposure on the outcome
  using robust MR methods (e.g., MRLOVA).
- Use the **same summary dataset** to ensure consistency between the
  LDSC analysis and MR.

## 3. Apply the HVPSUM

This R package performs the following steps:

1.  **LDSC Mode**: If `ldsc_env_path`, `ldsc_exe_path`, `ld_path`, and
    munged `summary data` are provided, the package performs LDSC
    analysis to estimate heritability and genetic covariance for each
    block using jackknife resampling. The required inputs are:

    - `ldsc_env_path`: Path to the LDSC environment.
    - `ldsc_exe_path`: Path to the LDSC executable.
    - `ld_path`: Path to the reference LD data.
    - `sumstats_outcome`: Path to the munged summary statistics file for
      the outcome trait.
    - `sumstats_exposure`: Path to the munged summary statistics file
      for the exposure trait.

    When these paths are provided, the package will execute LDSC
    analysis and compute the necessary block estimates from the
    jackknife resampling.

    ``` r
    # Example using LDSC Mode
    #library(HVPSUM)
    #hvpsum(ldsc_env_path = "/data/alh-admlda/anaconda3/envs/ldsc/bin/python", 
    #        ldsc_exe_path = "/data/alh-admlda/ldsc/ldsc.py",                   
    #        ld_path = "/data/alh-admlda/ldsc/eur_w_ld_chr/", 
    #        sumstats_outcome = "/data/alh-admlda/DM2.sumstats.gz", 
    #        sumstats_exposure = "/data/alh-admlda/mets.sumstats.gz", 
    #        tau = 0.21)
    ```

2.  **Pre-computed Mode**: If pre-computed data are provided for
    `h2_outcome`, `h2_exposure`, and `gencov`, these are used directly
    to bypass the LDSC analysis. The required inputs are:

    - `h2_outcome`: Pre-computed heritability estimates for each block
      for the outcome.
    - `h2_exposure`: Pre-computed heritability estimates for each block
      for the exposure.
    - `gencov`: Pre-computed genetic covariance estimates for each
      block.

    The outputs include estimates, standard errors, and p-values for
    corrected and uncorrected parameters.

    ``` r
    #Example using Pre-computed Mode

    library(HVPSUM)
    data("dat1")
    data("dat2")
    data("dat3")
    head(dat1)
    #>           V1
    #> 1 0.07696392
    #> 2 0.07642106
    #> 3 0.07556504
    #> 4 0.07559329
    #> 5 0.07599179
    #> 6 0.07575593
    head(dat2)
    #>           V1
    #> 1 0.04595247
    #> 2 0.04651566
    #> 3 0.04587754
    #> 4 0.04633127
    #> 5 0.04651247
    #> 6 0.04634499
    head(dat3)
    #>           V1
    #> 1 0.04297810
    #> 2 0.04307912
    #> 3 0.04293479
    #> 4 0.04297755
    #> 5 0.04302241
    #> 6 0.04294928
    tau <- 0.21
    hvpsum(h2_outcome = dat1, h2_exposure = dat2, gencov = dat3, tau =tau)   
    #> h2 outcome                      : Estimate = 0.0753, SE = 0.0086, P-value = < 2.22e-16
    #> h2 outcome (corrected)          : Estimate = 0.0594, SE = 0.0076, P-value = 5.9497e-15
    #> h2 exposure                     : Estimate = 0.0460, SE = 0.0065, P-value = 1.1615e-12
    #> Genetic Covariance              : Estimate = 0.0428, SE = 0.0056, P-value = 3.203e-14
    #> Genetic Covariance (corrected)  : Estimate = 0.0331, SE = 0.0052, P-value = 2.6609e-10
    #> Genetic Correlation             : Estimate = 0.7269, SE = 0.0796, P-value = < 2.22e-16
    #> Genetic Correlation (corrected) : Estimate = 0.6339, SE = 0.1025, P-value = 6.2006e-10
    ```

## Contact

Please contact Lamessa Amente (<lamessa.amente@mymail.unisa.edu.au>) or
Hong Lee (<hong.lee@unisa.edu.au>) if you have any queries.
