\encoding{UTF-8}

\name{genlist}

\alias{genlist}
\alias{classics.genlist}
\alias{colorimetry.genlist}

\docType{data}

\title{
zonohedra generators useful for testing and plotting
}

\description{
\tabular{ll}{
\code{classics.genlist} \tab 13 classic zonohedra generators\cr
\code{colorimetry.genlist} \tab 4 sets of Color Matching Functions (each set is a 3xN matrix)\cr
}}

\format{
Each is an S3 class \bold{genlist} object organized as a list of 3xN matrices (N varies).
The list must have names, preferably short names or abbreviations.
Each matrix can have optional attributes 
\code{"shortname"} and \code{"fullname"}
which are useful when printing with \code{\link{print.genlist}()}.
}

\note{
Making these S3 class \bold{genlist} makes it possible to easily
print a short summary using \code{\link{print.genlist}()}.

For colorimetry.genlist[[2]] a few remarks are in order.
These generators come from the xyz CIE color matching functions of 1931,
from 360 to 830 nm with 1 nm step.
From 699 to 830 nm, the angles between the generators only differ by a few microradians,
and it apparent that the designers tapered all 3 functions identically in that nm range.
For an illustration of this in the chromaticity domain, see \cite{Burns}, Figure 10.
When the zonohedron is constructed from these 132 generators, with the default options,
all these generators a 'collapsed' into a single one.
In the original matroid these 132 \emph{points} form a \emph{multiple group},
and in the simplified matroid they are collapsed to a single \emph{point}, labeled with 699.
}

\references{
ASTM E 308 - 01.
Standard Practice for Computing the Colors of Objects by Using the CIE System.
Table 1 

Scott A Burns.
\bold{The location of optimal object colors with more than two transitions}.
Color Research & Application.
Vol. 46.
No. 6.
pp 1180-1193.
2021.

Günther Wyszecki and W.S. Stiles.
\bold{Color Science : Concepts and Methods, Quantitative Data and Formulae}.
Second Edition.
Wiley-Interscience. 1982.
Table I(3.3.1).
pp. 723-735.
}

\source{
David Eppstein.
\bold{Zonohedra and Zonotopes}.
\cr
\url{https://www.ics.uci.edu/~eppstein/junkyard/ukraine/ukraine.html}

Colour & Vision Research Laboratory.
University College London.
\url{http://www.cvrl.org}
}

\seealso{
\code{\link{print.genlist}()}
}

\examples{
# get the names of 3 sets of color matching functions
names(colorimetry.genlist)
# [1] "xyz1931.5nm" "xyz1931.1nm" "lms2000.1nm"


# print zonohedra metrics associated with 3 sets of color matching functions
colorimetry.genlist
#                            fullname generators vertices  edges facets      area     volume pointed
# xyz1931.5nm  xyz at 5nm step                81     5100  10146   5048  1582.722   4070.345    TRUE
# xyz1931.1nm  xyz at 1nm step               471   112910 225720 112812 39586.707 509434.149    TRUE
# lms2000.1nm  lms at 1nm step               441   146642 292860 146220 22736.652 181369.085    TRUE
# ciexyzjv.5nm xyz at 5nm step (1978)         90     8012  16020   8010  1553.535   3951.899    TRUE
}

\keyword{datasets}
