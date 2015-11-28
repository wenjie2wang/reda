# R Package reda 

## To-do

* ~~Understand the origin code, add some necessary comments and 
	document every function with Roxygen2.~~

* ~~Rewrite plot function by using ggplot2 if necessary.~~

* ~~Add S4 class functions for print and summary.~~

* ~~test spline rate function by simulation.~~

* Polish the documentation for every function in a better organization.

* Test by using testthat package.  Write vignettes.

* Tests based on datasets from other R packages for 
	recurrent event analysis, such as 
	[survrec](http://cran.r-project.org/web/packages/survrec/index.html),
	[etc](http://cran.r-project.org/web/views/Survival.html).

* Setup repository on Github for issue or bug report.

* Add more handy functions for recurrent event analysis.

* Implement M-spline to replace B-spline bases
  for rate function (I-spline for MCF).

* Write up manuscript for JSS.

## Notes

* 'heart' is renamed by 'rateReg' since the model is
based on counts and rate function.
Therefore, in 'Survr', 'r' represents 'rate'.

* Later on, function fitting model based on gap times can be named as
'gapsReg' (or 'gapReg'). Response formula function can be named as 'Survg',
where 'g' means 'gaps' and which would have different data checking procedure
with 'Survr'.

* All the reference in the package manual follows
  [APA style](http://www.apastyle.org/learn/faqs/index.aspx).

