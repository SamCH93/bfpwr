# Closed-Form Power and Sample Size Calculations for Bayes Factors

This directory contains code and data related to the preprint

   Pawel, S., Held, L. (2024). Closed-Form Power and Sample Size Calculations
   for Bayes Factors.


## Reproducibility

The results can be reproduced by installing the necessary R packages

``` r
## CRAN packages
pkgs <- c("BayesRep", "lamW", "xtable", "remotes", "knitr")
install.packages(pkgs)

## GitHub package
remotes::install_github("https://github.com/SamCH93/bfpwr", subdir = "package")
```

and then rerunning the code in `paper/bfssd.R`. To recompile the manuscript make
sure to have LaTeX installed (tested only with TeX Live 2022/dev/Debian) and
then run

``` sh
make
```

which should produce `paper/bfssd.pdf`. The R and R package versions that were
used when the paper was successfully compiled before submission are visible in
the following output

``` r
sessionInfo()

#> R version 4.4.0 (2024-04-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 22.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
#> LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=de_CH.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=de_CH.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=de_CH.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=de_CH.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: Europe/Zurich
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] BayesRep_0.42.2 lamW_2.2.3      xtable_1.8-4    bfpwr_0.1      
#> [5] knitr_1.46     
#> 
#> loaded via a namespace (and not attached):
#> [1] compiler_4.4.0     tools_4.4.0        Rcpp_1.0.12        xfun_0.43         
#> [5] RcppParallel_5.1.7

cat(paste(Sys.time(), Sys.timezone(), "\n"))

#> 2024-05-27 14:16:12.421251 Europe/Zurich 
```
