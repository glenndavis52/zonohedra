\name{boundarypgramdata}

\alias{boundarypgramdata}

\title{
compute data about specific parallelograms in the boundary of a zonohedron
}

\description{
The boundary of a zonohedron is the union of parallelograms,
where some of them may be facets, and some may be tiles
in the standard tiling of more complex facets.
The edges of each parallelogram are given by a pair of distinct generators.
If a zonohedron has \eqn{n} of these generators,
then there are \eqn{n(n{-}1)/2} such pairs.
However, if the two generators are multiples of each other, or 0,
the parallelogram is degenerate and does not count.
For each pair of generators,
there are 2 parallelograms which are antipodal to each other.
The total number of parallelograms is \eqn{n(n{-}1)}.
This function computes data about one parallelogram from this antipodal pair.

}

\usage{
boundarypgramdata( x, gndpair, cube=FALSE )
}


\arguments{
\item{x}{
a \bold{zonohedron} object as returned by the constructor \code{zonohedron()}
}

\item{gndpair}{
an Mx2 integer matrix.
Each row of \code{gndpair} must contain a pair of points
in the ground set of the matroid of the zonohedron \code{x}.
If the 1st point is less that the 2nd point, then the standard parallelogram is used.
If the 2nd point is less than the 1st point, then the antipodal parallelogram is used.
If one of the points is not in the ground set of the matroid,
it is a silent error and all the returned data is set to \code{NA}.
If the two vector generators are multiples of each other, or 0,
then the parallelogram is degenerate and all the returned data is set to \code{NA}.
\cr
\code{gndpair} can also be a numeric vector that can be converted
to such a matrix, by row.

}

\item{cube}{
if \code{TRUE}, then a point of the cube that maps to the center
of the given parallelogram is returned, see \bold{Value}.
}
}



\value{
\code{boundarypgramdata()} returns a \code{data.frame}
with M rows and these columns:

\item{gndpair}{the given \code{gndpair}}

\item{hyperplaneidx}{
the index of the hyperplane in the simplified matroid
of \code{x} that contains \code{gndpair}
}

\item{center}{the center of the standard or antipodal parallelogram.
The centers of the standard and antipodal parallelograms add to white.
}

\item{transitions}{
the number of transitions in \code{pcube} - a point in the \eqn{n}-cube
that maps to \code{center},
where \code{n} is the number of generators of \code{x}.
This is a positive even integer.
}

And if \code{cube} is \code{TRUE}, then this column is added:

\item{pcube}{
a point in the \eqn{n}-cube that maps to \code{center}.
The sum of the spectra of the standard and the antipodal spectra is identically 1.
}

If a row of \code{gndpair} has an invalid pair,
the other columns are filled with \code{NA}s.

In case of global error, the function returns \code{NULL}.
}



\seealso{
\code{\link{zonohedron}()}
}

\examples{
zono =  zonohedron( colorimetry.genlist[[2]] )
boundarypgramdata( zono, c(570,608,  608,570, 400,450,  650,700,  650,720, 700,720,  650,900) )
##   gndpair.1 gndpair.2 hyperplaneidx     center.x     center.y     center.z transitions
## 1       570       608         49284  34.01432310  23.49690880   0.03214207           8
## 2       608       570         49284  72.85114639  83.36000830 106.86010920           8
## 3       400       450         12831   9.89612333   0.57529647  49.17990701           2
## 4       650       700             1   4.58023729   1.69316773   0.00000000           2
## 5       650       720             1   4.70484309   1.73816516   0.00000000           2
## 6       700       720            NA           NA           NA           NA          NA
## 7       650       900            NA           NA           NA           NA          NA

# In rows 1 and 2, the ground pairs are swapped, so the hyperlane index remains the same
# but the parallelograms are antipodal; the sum of their centers is the white point.
# Row 3 is a parallelogram facet, which is the usual situation.
# In rows 4 and 5, since generators for ground points 700 and 720 are multiples,
# the hyperplane index is the same. Both parallelograms are in a tiling of a non-trivial facet.
# In row 6, since the generators are multiples, the parallelogram is degenerate.
# In row 7, the point 900 is not in the ground set, so the parallelogram is undefined.
}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
