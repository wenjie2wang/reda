pkg = reda
header = HEADER

Rpkg: Rd build  
	make check 
	make INSTALL

Rd: 
	Rscript -e "library(methods); devtools::document();" 

build:  
	R CMD build ../$(pkg)

check: $(pkg)_*.tar.gz
	R CMD check --as-cran $(pkg)_*.tar.gz

INSTALL: $(pkg)_*.tar.gz
	R CMD INSTALL --build $(pkg)_*.tar.gz

## update copyright year in HEADER, R script and date in DESCRIPTION
updateHeader:
	yr=$$(date +"%Y");\
	sed -i "s/Copyright (C) [0-9]\{4\}/Copyright (C) $$yr/" $(header);\
# add HEADER file if there is no header
	for Rfile in R/*.R; do \
	if ! grep -e 'Copyright (C)' $$Rfile ;\
	then cat $(header) $$Rfile > tmp ;\
	mv tmp $$Rfile;\
	fi;\
	yr=$$(date +"%Y");\
	sed -i "s/Copyright (C) [0-9]*/Copyright (C) $$yr/" $$Rfile;\
	done;\
	dt=$$(date +"%Y-%m-%d");\
	sed -i "s/Date: [0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}/Date: $$dt/" DESCRIPTION;

clean: 
	rm -rf *~ */*~ */*.Rd *.Rhistroy NAMESPACE *.tar.gz *.Rcheck/ .\#* 
