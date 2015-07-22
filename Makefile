all: Rpkg 

Rpkg: 
	make build
	make check
	make INSTALL

build:  ../heart
	R CMD build -v; Rscript -e "Rd2roxygen::rab(install=TRUE)" ../heart

check: heart_*.tar.gz
	R CMD check --as-cran heart_*.tar.gz

INSTALL: heart_*.tar.gz
	R CMD INSTALL --build heart_*.tar.gz

clean: 
	rm -rf *~ .*~ */*.Rd NAMESPACE *.tar.gz *.Rcheck/
