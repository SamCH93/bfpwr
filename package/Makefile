## ----------------------------------------------------------------------------
## Document/Build/Check/Install R package (depends on devtools and roxygen2)
## Author: Samuel Pawel
## adapted from Manuela Ott, Sebastian Meyer, Florian Gerber
## ----------------------------------------------------------------------------

PACKAGE = bfpwr
VERSION = 0.1.6
TAR = $(PACKAGE)_$(VERSION).tar.gz

all: build

## update DESCRIPTION file
description:
	sed -i -r -- 's/^Version:.*/Version: '$(VERSION)'/g' DESCRIPTION ;      
	sed -i -r -- 's/^Date:.*/Date: '`date +'%F'`'/g' DESCRIPTION ;
	R -e 'roxygen2::roxygenize()'

## generate documentation with roxygen2
document: description
	R -e 'devtools::document()'

## build package manual
manual: document
	R -e 'devtools::build_manual(path = "out/")'

## build package
$(TAR): manual
	R -e 'devtools::build(path = "out/")'

## build package
build: $(TAR)

## install package
install: $(TAR)
	R -e 'devtools::install_local(path = "out/$(TAR)")'

## perform normal check
check: $(TAR)
	R -e 'devtools::check_built(path = "out/$(TAR)", cran = FALSE)'

## perform cran check
cran: $(TAR)
	R -e 'devtools::check_built(path = "out/$(TAR)", cran = TRUE, remote = TRUE)'

## perform only tests
test:
	R -e 'tinytest::build_install_test(pkgdir = ".")'

## build vignettes
vignettes:
	R -e 'devtools::build_vignettes(pkg = ".")'

.PHONY: all document manual build install check cran description test vignettes
