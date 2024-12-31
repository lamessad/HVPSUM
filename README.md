
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

- This r package performs the following steps:
  - 1.  If `ldsc_env_path`, `ldsc_exe_path`, and `ld_path` are provided,
        it runs the LDSC analysis to generate heritability and genetic
        covariance estimates.
  - 2.  If precomputed data is provided for `h21_data`, `h22_data`, and
        `gencov_data`, these are used directly.
  - 3.  It calculates corrected heritability, genetic covariance, and
        correlation, along with their associated p-values.
  - 4.  It outputs the estimates, standard errors, and p-values for the
        corrected and uncorrected quantities.

- 1.  `LDSC Mode`: In this mode, the package performs LDSC analysis to
      compute delete values for the input summary statistics. The
      following inputs must be provided:

  - `ldsc_env_path`: Path to the LDSC environment.
  - `ldsc_exe_path`: Path to the LDSC executable.
  - `ld_path`: Path to the LD data.
  - `sumstats1`: Path to the `munge` summary data of outcome.
  - `sumstats2`: Path to the `munge` summary data of outcome.

When these paths are provided, the package will execute LDSC analysis
and compute the necessary delete values.

``` r
@examples
# Example using LDSC Mode
#hvpsum(ldsc_env_path = "path/to/env//ldsc/bin/python", ldsc_exe_path = "path/to/ldsc.py", 
#      ld_path = "path/to/ld", sumstats1 = "path/to/file1.sumstats.gz", 
#       sumstats2 = "path/to/file2.sumstats.gz", tau = 0.5, intern = TRUE)
```

- 2.  `Pre-computed Mode`: In this mode, the package performs the
      analysis without running LDSC in this rpackage. Instead, it uses
      pre-computed delete values from LDSC computed else where. The
      following inputs must be provided:

  - `h21_data`: pre-computed `hertability` blocks estemates for trait
    1(outcome).
  - `h22_data`: pre-computed `hertability` blocks estemates for trait
    2(exposure).
  - `gencov_data`: pre-computed `genetic covariance` estimates for
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
tau <- as.numeric(0.21)
results <- hvpsum(h21_data = dat1, h22_data = dat2, gencov_data = dat3,tau=tau)
#> Heritability 1: Estimate = 0.0753, SE = NA, P-value = < 2.22e-16
#> Heritability 1 Corrected: Estimate = 0.0594, SE = NA, P-value = 5.9497e-15
#> Heritability 2: Estimate = 0.0460, SE = NA, P-value = 1.1615e-12
#> Genetic Covariance: Estimate = 0.0428, SE = NA, P-value = 3.203e-14
#> Genetic Covariance Corrected: Estimate = 0.0331, SE = NA, P-value = 2.6609e-10
#> Genetic Correlation: Estimate = 0.7269, SE = NA, P-value = < 2.22e-16
#> Genetic Correlation Corrected: Estimate = 0.6339, SE = NA, P-value = 6.2006e-10
results
#> NULL
```

## Contact

Please contact Lamessa Amente (<lamessa.amente@mymail.unisa.edu.au>) or
Hong Lee (<hong.lee@unisa.edu.au>) if you have any queries.
