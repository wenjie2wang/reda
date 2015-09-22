Rpkg: Rd build  
	make check 
	make INSTALL

Rd: 
	Rscript -e "library(methods); devtools::document();" 

build:  
	R CMD build ../reda

check: reda_*.tar.gz
	R CMD check --as-cran reda_*.tar.gz

INSTALL: reda_*.tar.gz
	R CMD INSTALL --build reda_*.tar.gz

clean: 
	rm -rf *~ */*~ */*.Rd *.Rhistroy NAMESPACE *.tar.gz *.Rcheck/
