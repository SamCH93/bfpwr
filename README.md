# Closed-Form Power and Sample Size Calculations for Bayes Factors

This repository contains 

1. `./package` The R package **bfpwr** to perform power and sample size
   calculations for Bayes factor analysis

2. `./paper` Code and data to reproduce result from the paper: *Pawel, S., Held.
   L. (2025). Closed-Form Power and Sample Size Calculations for Bayes Factors.
   The American Statistician. <https://doi.org/10.1080/00031305.2025.2467919>*

To cite our work, use the following BibTeX reference

```BibTeX
@article{PawelHeld2025,
  year = {2025},
  author = {Samuel Pawel and Leonhard Held},
  title = {Closed-Form Power and Sample Size Calculations for {Bayes} Factors},
  journal = {The American Statistician},
  doi = {10.1080/00031305.2025.2467919}
}
```

## Reproducing the paper with Docker

Make sure to have Docker and Make installed, then run `make docker` from the
root directory of this git repository. This will install all necessary
dependencies. RStudio Server can then be opened from a browser
(<http://localhost:8787>), and the R scripts in `./paper`, e.g., `bfssd.R` which
contains all code for the results from the paper, can be rerun (make sure to set
the working directory to `./paper` when running R interactively).
