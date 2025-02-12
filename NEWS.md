# *News*

## 1.6.1

* Date 02-12-2025
* [Fixed bug in weight specification when x was a vector](https:///github.com/vubiostat/acepack/issues/15)
* When ierr isn't 0 calls a function that defaults to warning.

## 1.6.0

* Date: 02-05-2025
* [Added independence test using ACE as provided by Holzmann & Klar.](https://github.com/vubiostat/acepack/issues/11)
* [Added formula interface to ACE.](https://github.com/vubiostat/acepack/issues/6)
* [Added S3 plot/summary/print to ace and AVAS output.](https://github.com/vubiostat/acepack/issues/7)
* [Added Checking of ACE error codes.](https://github.com/vubiostat/acepack/issues/13)

## 1.5.2

* Date: 01-27-2025
* [Fixing parallel make build issue and type checking on scratch memory.](https://github.com/vubiostat/acepack/issues/9)

## 1.5.1

* Date: 01-21-2025
* [Fixing errors detected by Ripley in future version of R.](https://github.com/vubiostat/acepack/issues/8)

## 1.5.0

* Date: 01-05-2025
* Complete refactor into F90 assisted by ChatGPT.
* Merged all globals to a single parameter set to control routines and provided calls to set them in the FORTRAN.

## 1.4.2

* Date: 03-29-2018
* Updated email address of maintainer.

## 1.4.1

* Date: 10-27-2016
* Converted all routines to "implicit none" by request for Prof Brian Ripley. Added reference 

## 1.4.0

* Date: 10-13-2016
* Jonathan Baron (Penn) has transferred maintainer status to Shawn Garbett (Vanderbilt). A shoutout and thanks to Dr. Baron's work for the last several years. 

## 1.3.x

Older notes. Unknown versions.

* 10-6-2016 Shawn Garbett cleaned up fortran6 warnings to maintain CRAN status. Stated precision of exported variables now matches internal computations. An extremely unlikely but possible division by zero path is prevented.
* 11-23-2014 Removed non-standard files from top directory.
* 4-20-2013 Brian Ripley fixed a Fortran error.
* 4-5-2012 Added namespace and removed stray print in avas.
* 7-4-2010 Fixed options circ, cat, and mon, in both ace and avas, so that they now can apply to the dependent variable, as specified previously in both the help page and the fortran code. Colin McCullogh did most of the work. Frank Harrell also reported this bug.
* Fixed the checks on the options so that they apply to the correct dimension. Previously circ and mon were not working as described. Thanks to Frank Harrell.
* Expanded the help pages to make them clearer and provide more examples. --Jon Baron
