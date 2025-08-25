\name{zonoseg-invert}

\alias{invert.zonoseg}

\title{
invert points in a zonoseg
}

\description{
For points in a zonoseg, find points in the unit cube that map to those points.
}

\usage{
\method{invert}{zonoseg}( x, z, tol=0, ... )
}


\arguments{
\item{x}{
a \code{zonoseg} object as returned by the constructor \code{zonoseg()}
}

\item{z}{a numeric M-vector}

\item{tol}{points that are within \code{tol} of a boundary point are
taken to be that point}

\item{...}{not used}
}


\details{
For a point in the interior of the zonoseg,
there are infinitely many points in the cube that map to it.
This function tries to find one with the fewest number of non-zero
components.
}


\value{
\code{invert.zonoseg()} returns a \code{data.frame} with M rows and these columns:

\item{z}{the given point}
\item{pcube}{a point in the unit cube that maps to \code{z}. For the 2 boundary points, \code{pcube} is a vertex. If \code{z} is outside the zonoseg, \code{pcube} is all \code{NA}, and a warning message is printed}

If the names of \code{z} are unique,
they are copied to the row names of the output.
The column names are copied from the ground set of the associated matroid.
\cr\cr
In case of error, the function returns \code{NULL}.
}


\seealso{
\code{\link{zonoseg}()}
}

\examples{
zono1 = zonoseg( c(1,-2,3,0,-3,-4) )
zono1
# generators:        6 -- 3 negative, 2 positive, and 1 loops.
#
# segment:                   [-9,4]
#      value pcube.1 pcube.2 pcube.3 pcube.4 pcube.5 pcube.6
# zmin    -9       0       1       0       0       1       1
# zmax     4       1       0       1       0       0       0
#
# 2-transition subsegment:   [-8,3]
#             value source.1 source.2 source.3 source.4 source.5 source.6
# tmin-2trans    -8        1        1        0        0        1        1
# tmax-2trans     3        0        0        1        1        0        0


z = c( 0, -3*pi, pi, 2*pi, getsegment(zono1) )

invert( zono1, z )
#           z   pcube.1   pcube.2   pcube.3   pcube.4   pcube.5   pcube.6
# 1  0.000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
# 2 -9.424778        NA        NA        NA        NA        NA        NA
# 3  3.141593 1.0000000 0.0000000 0.7138642 0.0000000 0.0000000 0.0000000
# 4  6.283185        NA        NA        NA        NA        NA        NA
# 5 -9.000000 0.0000000 1.0000000 0.0000000 0.0000000 1.0000000 1.0000000
# 6  4.000000 1.0000000 0.0000000 1.0000000 0.0000000 0.0000000 0.0000000
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
