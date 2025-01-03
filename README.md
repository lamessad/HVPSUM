
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
    munged `sumstat` are provided, it runs the LDSC analysis to generate
    heritability and genetic covariance estimates. The required inputs
    are:

    - `ldsc_env_path`: Path to the LDSC environment.
    - `ldsc_exe_path`: Path to the LDSC executable.
    - `ld_path`: Path to the LD data.
    - `sumstats1`: Path to the munged summary data of outcome.
    - `sumstats2`: Path to the munged summary data of exposure.

    When these paths are provided, the package will execute LDSC
    analysis and compute the necessary delete values.

    ``` r
    # Example using LDSC Mode
    # library(HVPSUM)
    # hvpsum(ldsc_env_path = "/data/alh-admlda/anaconda3/envs/ldsc/bin/python", 
    #        ldsc_exe_path = "./ldsc.py", 
    #        ld_path = "/data/alh-admlda/ldsc/eur_w_ld_chr/", 
    #        sumstats1 = "DM2.sumstats.gz", 
    #        sumstats2 = "mets.sumstats.gz", 
    #        tau = 0.21)
    ```

2.  **Pre-computed Mode**: If pre-computed data are provided for
    `h21_data`, `h22_data`, and `gencov_data`, these are used directly
    to bypass the LDSC analysis. The required inputs are:

    - `h21_data`: Pre-computed heritability block estimates for trait 1
      (outcome).
    - `h22_data`: Pre-computed heritability block estimates for trait 2
      (exposure).
    - `gencov_data`: Pre-computed genetic covariance block estimates.

    The outputs include estimates, standard errors, and p-values for
    corrected and uncorrected parameters.

    ``` r
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
    hvpsum(h21_data = dat1, h22_data = dat2, gencov_data = dat3, tau = tau)
    #> Heritability 1: Estimate = 0.0753, SE = 0.0000, P-value = < 2.22e-16
    #> Heritability 1 Corrected: Estimate = 0.0594, SE = 0.0000, P-value = < 2.22e-16
    #> Heritability 2: Estimate = 0.0460, SE = 0.0000, P-value = < 2.22e-16
    #> Genetic Covariance: Estimate = 0.0428, SE = 0.0000, P-value = < 2.22e-16
    #> Genetic Covariance Corrected: Estimate = 0.0331, SE = 0.0000, P-value = < 2.22e-16
    #> Genetic Correlation: Estimate = 0.7269, SE = 0.0000, P-value = < 2.22e-16
    #> Genetic Correlation Corrected: Estimate = 0.6339, SE = 0.0000, P-value = < 2.22e-16
    ```

## Contact

Please contact Lamessa Amente (<lamessa.amente@mymail.unisa.edu.au>) or
Hong Lee (<hong.lee@unisa.edu.au>) if you have any queries.
