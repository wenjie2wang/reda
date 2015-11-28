pkg = reda

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

## update copyright year in each R source and date in DESCRIPTION
updateYear: 
	yr=$$(date +"%Y");\
	for f in R/*.R; do sed -i "s/Copyright (C) [0-9]\{4\}/Copyright (C) $$yr/" $$f; done;\
	dt=$$(date +"%Y-%m-%d");\
	sed -i "s/Date: [0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}/Date: $$dt/" DESCRIPTION;

clean: 
	rm -rf *~ */*~ */*.Rd *.Rhistroy NAMESPACE *.tar.gz *.Rcheck/ .\#* 
