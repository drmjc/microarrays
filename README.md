microarrays
===========
An R package containing a collection of microarray and genomic helper code. 

Description
===========
An R package containing a collection of microarray and genomic helper code. Lots of code for raw microarray data and algorithm result file I/O, conversion, plots, and summarisation. 

Most of this was developed from 2007-2010, so may be a little out of date. caveat emptor.

Installation
============
Note there are quite a few dependencies, listed in the DESCRIPTION file, which will probably need installing. I think this is the complete list. If you get an error, then try installing that package, either through install.packages, biocLite, or install_github.

    # install.packages("devtools")
    library(devtools)
    install_github("excelIO", "drmjc")
    install_github("mjcbase", "drmjc")
    install_github("mjcgraphics", "drmjc")
    install_github("mjcstats", "drmjc")
    install_github("mjcaffy", "drmjc")
    install_github("microarrays", "drmjc")

Usage
=====
Extensive package documentations is available via:

	library(microarrays)
	?microarrays

Exporting Illumina data to GEO
==============================
Take a look at ?LumiBatch2GEOarchive
