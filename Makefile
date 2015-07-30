# prep:
#	rm -fr ../heart; \
#	mkdir ../heart; \
#	cp -pr data ../heart; \
#	cp -p DESCRIPTION ../heart; \
#	cp -pr R ../heart; \
#	cp -pr vignettes ../heart; \
#	cd ../heart; \
#	Rscript -e "library(methods); devtools::document();"; \

Rpkg: Rd build  
	make check 
	make INSTALL

Rd: 
	Rscript -e "library(methods); devtools::document();" 

build:  
	R CMD build ../heart

check: heart_*.tar.gz
	R CMD check --as-cran heart_*.tar.gz

INSTALL: heart_*.tar.gz
	R CMD INSTALL --build heart_*.tar.gz

clean: 
	rm -rf *~ .*~ */*.Rd NAMESPACE *.tar.gz *.Rcheck/
