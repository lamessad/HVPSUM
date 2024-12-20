
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

By addressing the limitations of conventional approaches, the HVP
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

## 1. Munge the Summary Data data (see LDSC for the detail)

## 2. Run MR

- Estimate the causal effect size (`tau`) of the exposure on the outcome
  using robust MR methods (e.g., MRLOVA).
- Use the **same summary dataset** to ensure consistency between the
  LDSC analysis and MR.

## 3. Apply the HVPSUM

- This package operates in two distinct modes, and the required
  **inputs** depend on the selected mode:

- 1.  `**LDSC Mode**`: In this mode, the package performs LDSC analysis
      to compute delete values for the input summary statistics. The
      following inputs must be provided:

  - `ldsc_env_path`: Path to the LDSC environment.
  - `ldsc_exe_path`: Path to the LDSC executable.
  - `ld_path`: Path to the LD data.
  - `sumstats1`: Path to the `munge` summary data of outcome.
  - `sumstats2`: Path to the `munge` summary data of outcome.

When these paths are provided, the package will execute LDSC analysis
and compute the necessary delete values.

- 2.  `**Pre-computed Mode**`: In this mode, the package performs the
      analysis without running LDSC in this rpackage. Instead, it uses
      pre-computed delete values from LDSC computed else where. The
      following inputs must be provided:

  - `h21block_path`: Path to the pre-computed `hertability` blocks
    estemates for trait 1(outcome).
  - `h22block_path`: Path to the pre-computed `hertability` blocks
    estemates for trait 2(exposure).
  - `gen_cov_path`: Path to the pre-computed `genetic covariance` for
    blocks.

If these paths are provided, the package will bypass the LDSC analysis
and directly perform the analysis using the provided data. - Outputs
includes both LDSC and corrected for heritability, genetic covaraince,
and genetic correlation.

``` r
library(HVPSUM)
data("dat1")#
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
#results <- hvpsum(h21block_path = "path/to/file1.hsq1.delete", h22block_path = "path/to/file2.hsq2.delete", gen_cov_path = "path/to/file3.gencov.delete",tau=tau)
#results
#Heritability 1: Estimate = 0.0991 SE = 0.0309 P-value = 0.001333
#Heritability 1 Corrected: Estimate = 0.1005 SE = 0.0313 P-value = 0.0013128
#Heritability 2: Estimate = 0.0753 SE = 0.0086 P-value = < 2.22e-16
#Genetic Covariance: Estimate = -0.0045 SE = 0.0049 P-value = 0.35232
#Genetic Covariance Corrected: Estimate = 0.0113 SE = 0.005 P-value = 0.023646
#Genetic Correlation: Estimate = -0.0523 SE = 0.0588 P-value = 0.37434
#Genetic Correlation: Estimate = 0.1299 SE = 0.0553 P-value = 0.01876
```

## Contact

Please contact Lamessa Amente (<lamessa.amente@mymail.unisa.edu.au>) or
Hong Lee (<hong.lee@unisa.edu.au>) if you have any queries.
