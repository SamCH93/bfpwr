## set R version (https://hub.docker.com/r/rocker/verse/tags)
FROM rocker/verse:4.4.1

## set up directories
WORKDIR home/rstudio
RUN mkdir /home/rstudio/paper
COPY package /home/rstudio/package

## install R packages from CRAN the last day of the specified R version
RUN install2.r --error --skipinstalled --ncpus -1 \
    BayesRep lamW xtable remotes knitr ggplot2 dplyr && \
    R -e 'remotes::install_github("nicebread/BFDA@1fd2bf7ae2d06ceee97deb4e573dcd90a7e0e2a7", subdir = "package")' && \
    R CMD INSTALL  package/out/bfpwr_0.1.3.tar.gz
