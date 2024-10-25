all: pdf code supplement clean

# names of .Rnw .bib. and .bst files
FILE=bfssd
BIB=bibliography
BIBSTYLE=apalikedoiurl

## generate .tex file from .Rnw file with knitr
tex: $(FILE).Rnw $(BIB).bib $(BIBSTYLE).bst
	Rscript -e "knitr::knit('$(FILE).Rnw')" --vanilla

## generate .pdf file from .tex file
pdf: tex
	pdflatex $(FILE)
	bibtex $(FILE)
	pdflatex $(FILE)
	pdflatex $(FILE)

## generate .pdf file from .tex (without first running the tex recipe)
pdf2:
	pdflatex $(FILE)
	bibtex $(FILE)
	pdflatex $(FILE)
	pdflatex $(FILE)

supplement: supplement.Rnw appendix.Rnw tex
	Rscript -e "knitr::knit2pdf('supplement.Rnw')" --vanilla

## extract R code
code: $(FILE).Rnw
	Rscript -e "knitr::purl('$(FILE).Rnw')" --vanilla

## use texdiff to show difference to previous version
texdiff: tex old/$(FILE)V1.tex old/ld.cfg
	-latexdiff -c old/ld.cfg old/$(FILE)V1.tex $(FILE).tex > $(FILE)-changes.tex
	-pdflatex --interaction nonstopmode $(FILE)-changes
	-bibtex $(FILE)-changes
	-pdflatex --interaction nonstopmode $(FILE)-changes
	-pdflatex --interaction nonstopmode $(FILE)-changes

## clean
clean:  
	-rm $(FILE).aux  $(FILE).blg  $(FILE).log  $(FILE).tex  $(FILE).bbl \
	$(FILE).out $(FILE).brf $(FILE).synctex.gz \
	supplement.aux  supplement.blg  supplement.log  supplement.tex  supplement.bbl \
	supplement.out supplement.brf supplement.synctex.gz \
	$(FILE)-changes.aux  $(FILE)-changes.blg  $(FILE)-changes.log \
	$(FILE)-changes.tex  $(FILE)-changes.bbl $(FILE)-changes.out \
	$(FILE)-changes.brf $(FILE)-changes.synctex.gz
	-rm -r ./figure/
