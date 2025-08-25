\name{invertboundary}

\alias{invertboundary}
\alias{invertboundary.zonohedron}

\title{
invert points on the boundary of a zonohedron
}

\description{
A \emph{zonohedron} \eqn{Z} is the image of a linear map
\eqn{[0,1]^n \to Z \subset \bold{R}^3}, from the n-cube to 3D space.
For a point on the boundary of the zonohedron, this function
computes a point in the unit cube that maps to it.
All coordinates of the point in the cube are 0 or 1, 
except for two of them.
The point is not necessarily unique.
}

\usage{
\method{invertboundary}{zonohedron}( x, point, tol=5.e-14 )
}


\arguments{
\item{x}{
a \bold{zonohedron} object as returned by the constructor \code{zonohedron()}
}

\item{point}{
Mx3 matrix with points on the boundary of \code{x} in the rows.
Such a matrix typically is returned by
\code{\link{raytrace.zonohedron}()}
or
\code{\link{section.zonohedron}()}.
\code{point} can also be a numeric vector that can be converted
to such a matrix, by row.
}

\item{tol}{
points that are not within \code{tol} of the boundary are
skipped, see \bold{Details}
}

}


\details{
Given the boundary point, the function determines the facet
that contains it.
The \code{pcube} coordinates of the \emph{base vertex} of this facet are
all 0 or 1, and fairly easy to determine.
If the facet is a parallelogram, the other two coordinates are
fairly easy to determine too.
If the facet is a zonogon with \code{K} generators, with \code{K>2},
then the unknown \code{K} coordinates are calculated
using \code{\link{invert.zonogon}()}.
\cr
Because of floating point behaviour, coordinates can be slightly
negative or slightly more than 1.
After the calculation, they are clamped to [0,1].
}


\value{
\code{invertboundary.zonohedron()} returns a \code{data.frame}
with M rows and these columns:

\item{point}{the given boundary point}

\item{distance}{
signed distance to the boundary of \code{x};
for successful inversion its absolute value is \eqn{\le} \code{tol}
}

\item{facetidx}{index of the facet pair that contains the point}

\item{sign}{sign of the facet pair; either +1 or -1}

\item{pcube}{
a point in the unit n-cube that maps to the given boundary point;
all coordinates of \code{pcube} are 0 or 1,
except for 2 of them.
}

\item{transitions}{
the number of transitions in \code{pcube} - a non-negative even integer
}

If a point \code{point} cannot be inverted, e.g. because
\code{distance} is too large, the other columns are all \code{NA}.  

If the row names of \code{point} are unique,
they are copied to the row names of the output.
The column names of \code{pcube} are copied from the ground set of the
associated matroid.
\cr\cr
In case of global error, the function returns \code{NULL}.
}

\seealso{
\code{\link{zonohedron}()},
\code{\link{section.zonohedron}()},
\code{\link{raytrace.zonohedron}()},
\code{\link{invert.zonogon}()}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
