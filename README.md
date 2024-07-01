# Closed-Form Power and Sample Size Calculations for Bayes Factors

This repository contains 

1. `./package` The R package **bfpwr** to perform power and sample size calculations for
   planned Bayes factor analysis

2. `./paper` Code and data to reproduce result from the preprint: *Pawel, S.,
   Held. L. (2024). Closed-Form Power and Sample Size Calculations for Bayes
   Factors. <https://doi.org/10.48550/arxiv.2406.19940>*

To cite our work, use the following BibTeX reference

```BibTeX
@article{Pawel2024,
  year = {2024},
  author = {Samuel Pawel and Leonhard Held},
  title = {Closed-Form Power and Sample Size Calculations for {Bayes} Factors},
  doi = {10.48550/arxiv.2406.19940},
  note = {Preprint}
}
```

## Reproducing the paper with Docker

Make sure to have Docker and Make installed, then run `make docker` from the
root directory of this git repository. This will install all necessary
dependencies. RStudio Server can then be opened from a browser
(<http://localhost:8787>), and the R scripts in `./paper` can be rerun (make
sure to set the working directory to `./paper` when running R interactively).
