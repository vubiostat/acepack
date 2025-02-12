acepack
=======

`acepack` is an [R](https://www.r-project.org) package that provides two
nonparametric methods for multiple regression transform selection.

The first, Alternating Conditional Expectations (ACE), 
is an algorithm to find the fixed point of maximal
correlation, i.e. it finds a set of transformed response variables that
maximizes R^2 using smoothing functions [see Breiman, L., and J.H. Friedman.
1985. "Estimating Optimal Transformations for Multiple Regression and
Correlation". Journal of the American Statistical Association. 80:580-598. 
<doi:10.1080/01621459.1985.10478157>].

Also included is the Additivity Variance Stabilization (AVAS) method which works
better than ACE when correlation is low [see Tibshirani, R.. 1986. "Estimating 
Transformations for Regression via Additivity and Variance Stabilization".
Journal of the American Statistical Association. 83:394-405. 
<doi:10.1080/01621459.1988.10478610>].

A good introduction to these two methods is in chapter 16 of
Frank Harrell's "Regression Modeling Strategies" in the Springer Series in Statistics.

History
===============

This package is based on public domain S and FORTRAN code for AVAS by 
Tibshirani, and on FORTRAN code for ACE from Statlib, written by Spector
and Friedman.

The FORTRAN code has been edited to use double precision, for
compatibility with R, and the R code and documentation for ace() have been
added by Thomas Lumley, based on that for avas().

Shawn Garbett has refactored with the assistance of ChatGPT to F90 and 
cleaned up the R interface to current standards.

