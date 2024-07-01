# Closed-Form Power and Sample Size Calculations for Bayes Factors

This directory contains code and data related to the preprint

   Pawel, S., Held, L. (2024). Closed-Form Power and Sample Size Calculations
   for Bayes Factors. <https://doi.org/10.48550/arxiv.2406.19940>


## Reproducibility

The results can be reproduced by installing the necessary R packages

``` r
## CRAN packages
pkgs <- c("BayesRep", "lamW", "xtable", "remotes", "knitr")
install.packages(pkgs)

## GitHub packages
remotes::install_github("SamCH93/bfpwr", subdir = "package")
remotes::install_github("nicebread/BFDA", subdir = "package")
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
#> [1] parallel  stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] BFDA_0.5.0        ggplot2_3.5.1     dplyr_1.1.4       doParallel_1.0.17
#>  [5] iterators_1.0.14  foreach_1.5.2     BayesRep_0.42.2   lamW_2.2.4       
#>  [9] xtable_1.8-4      bfpwr_0.1         knitr_1.47       
#> 
#> loaded via a namespace (and not attached):
#>  [1] utf8_1.2.4         generics_0.1.3     stringi_1.8.4      lattice_0.22-5    
#>  [5] digest_0.6.35      magrittr_2.0.3     grid_4.4.0         fastmap_1.2.0     
#>  [9] sfsmisc_1.1-18     plyr_1.8.9         Matrix_1.6-5       hypergeo_1.2-13   
#> [13] deSolve_1.40       gridExtra_2.3      promises_1.3.0     mgcv_1.9-1        
#> [17] doRNG_1.8.6        fansi_1.0.6        scales_1.3.0       truncnorm_1.0-9   
#> [21] abtest_1.0.1       codetools_0.2-19   cli_3.6.2          shiny_1.8.1.1     
#> [25] rlang_1.1.4        contfrac_1.1-12    munsell_0.5.1      splines_4.4.0     
#> [29] withr_3.0.0        plotrix_3.8-4      tools_4.4.0        reshape2_1.4.4    
#> [33] colorspace_2.1-0   httpuv_1.6.15      rngtools_1.5.2     mime_0.12         
#> [37] vctrs_0.6.5        R6_2.5.1           qgam_1.3.4         lifecycle_1.0.4   
#> [41] stringr_1.5.1      MASS_7.3-60.2      pkgconfig_2.0.3    later_1.3.2       
#> [45] RcppParallel_5.1.7 pillar_1.9.0       gtable_0.3.5       glue_1.7.0        
#> [49] Rcpp_1.0.12        TeachingDemos_2.13 xfun_0.44          tibble_3.2.1      
#> [53] tidyselect_1.2.1   htmltools_0.5.8.1  nlme_3.1-165       elliptic_1.4-0    
#> [57] compiler_4.4.0 

cat(paste(Sys.time(), Sys.timezone(), "\n"))
#> 2024-06-28 15:03:59.985729 Europe/Zurich
```
