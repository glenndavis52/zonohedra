

### Version 0.5-0 [Sep 1, 2025]

* added new function `as.mesh3d.zonohedron()`, which is derived from `rgl::as.mesh3d()`
* moved package `rgl` from Suggests to Imports
* in the User Guide, replaced an animated GIF with a WEBM video
* in Suggests, dropped package `gifski` and added `av` and `base64enc`
* in function `plot2trans()`, fixed a bug when type='e'


### Version 0.4-0 [Feb 1, 2025]

* for function `boundarypgramdata()`, improved correctness and slightly changed the return value
* for function `inside()`, fixed error on the man page
* added new item - `collapsetosimple` - to matroid list
* fixed a compilation warning


### Version 0.3-0 [July 7, 2024]

* improved the wording of three vignettes
* in `inside2trans()`, made calculation of linking number faster
* in `matroid.matrix()`, rank-deficient matrix was not always detected
* added global option `zonohedra.stoponerror`
* fixed `rchk` compilation issues
* fixed `noRemap` compilation issues


### Version 0.2-2 [May 30, 2023]

* fixed some *Additional issues* for variants clang-UBSAN, gcc-UBSAN, and M1mac


### Version 0.2-1 [May 27, 2023]

* initial version on CRAN
