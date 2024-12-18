
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
relationships. \## Installation

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

- 2.  `Pre-computed Mode`: In this mode, the package performs the
      analysis without running LDSC in this rpackge. Instead, it uses
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
#library(HVPSUM)
#data(c("dat1", "dat2", "dat3"))
#head(dat1)
#head(dat2)
#head(dat3)
#tau <- 0.21
#results <- hvpsum(h21block_path = "dat1", h22block_path = "dat2", gen_cov_path = "dat3",tau=tau)
#results
```

## Contact

Please contact Lamessa Amente (<lamessa.amente@mymail.unisa.edu.au>) or
Hong Lee (<hong.lee@unisa.edu.au>) if you have any queries.
