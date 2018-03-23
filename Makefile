objects := $(wildcard R/*.R) DESCRIPTION
version := $(shell grep "Version" DESCRIPTION | sed "s/Version: //")
pkg := $(shell grep "Package" DESCRIPTION | sed "s/Package: //")
tar := $(pkg)_$(version).tar.gz
checkLog := $(pkg).Rcheck/00check.log
citation := inst/CITATION
yr := $(shell date +"%Y")
dt := $(shell date +"%Y-%m-%d")

rmd := $(wildcard vignettes/*.Rmd)
vignettes := $(patsubst %.Rmd,%.html,$(rmd))
cprt := COPYRIGHT


.PHONY: check
check: $(checkLog)

.PHONY: build
build: $(tar)

.PHONY: install
install: $(tar)
	R CMD INSTALL $(tar)

.PHONY: preview
preview: $(vignettes)

.PHONY: pkgdown
pkgdown:
	@mkdir -p docs/inst/bib/; cp inst/bib/$(pkg).bib docs/inst/bib/;
	@echo "added/updated the package bib file"
	Rscript -e "library(methods); pkgdown::build_site();"

$(tar): $(objects)
	@if [ "$$(uname)" == "Darwin" ];\
	then echo "remeber to update date and version number";\
	else make -s updateMeta;\
	fi;\
	Rscript -e "library(methods); devtools::document();";
	R CMD build .

$(checkLog): $(tar)
	R CMD check --as-cran $(tar)

vignettes/%.html: vignettes/%.Rmd
	Rscript -e "library(methods); rmarkdown::render('$?')"


## update copyright year in HEADER, R script and date in DESCRIPTION
.PHONY: updateMeta
updateMeta: $(objects) $(cprt)
	@echo "Updating date, version, and copyright year"
	@sed -i "s/Copyright (C) 2015-*[0-9]*/Copyright (C) 2015-$(yr)/" $(cprt)
	@for Rfile in R/*.R; do \
	if ! grep -q 'Copyright (C)' $$Rfile;\
	then cat $(cprt) $$Rfile > tmp;\
	mv tmp $$Rfile;\
	fi;\
	sed -i "s/Copyright (C) 2015-*[0-9]*/Copyright (C) 2015-$(yr)/" $$Rfile;\
	done;
	@sed -i "s/Date: [0-9]\{4\}-[0-9]\{1,2\}-[0-9]\{1,2\}/Date: $(dt)/" DESCRIPTION
	@sed -i "s/version [0-9]\.[0-9]\.[0-9]\(\.[0-9][0-9]*\)*/version $(version)/" $(citation)
	@sed -i "s/20[0-9]\{2\}/$(yr)/" $(citation)

## make tags
.PHONY: TAGS
TAGS:
	Rscript -e "utils::rtags(path = 'R', ofile = 'TAGS')"
	gtags

.PHONY: clean
clean:
	rm -rf *~ */*~ *.Rhistroy *.tar.gz src/*.so src/*.o *.Rcheck/ .\#*
