\name{support}

\alias{support}
\alias{support.zonotope}

\title{
the support function for a zonotope
}

\description{
Compute the classical support function for a zonotope.
It also computes a point on the boundary where the linear functional
is maximized, and the dimension of the face where the supporting
hyperplane intersects the zonotope.
}

\usage{
\method{support}{zonotope}( x, direction, tol=5.e-15 )
}


\arguments{
\item{x}{
a zonotope object -  a \bold{zonohedron}, a \bold{zonogon},
or a \bold{zonoseg}
}

\item{direction}{
an NxM matrix with N directions in the rows.
If \code{x} is a \bold{zonohedron}, M must be 3.
If \code{x} is a \bold{zonogon}, M must be 2.
If \code{x} is a \bold{zonoseg}, M must be 1.
\code{direction} can also be a vector that can be converted to
such a matrix, by row.
The direction is normal to the supporting hyperplane
}

\item{tol}{
the tolerance for determining whether the supporting hyperplane
intersects a face with positive dimension.
This does not affect the value of the support function.
For a \bold{zonoseg}, \code{tol} is ignored.
}

}



\value{
The function returns a \code{data.frame} with N rows and these columns:

\item{direction}{the given direction}
\item{value}{the value of the support function of \code{x}, in the given direction}
\item{argmax}{a point on the boundary of \code{x} where the functional max is taken.  This point is the center of the face where the supporting hyperplane intersects the zonotope.}
\item{dimension}{of the face where the supporting hyperplane intersects the zonotope. 0 means a vertex, 1 means an edge, and 2 means a 2-face.}

If \code{direction} is 0, the other columns are \code{NA}.
If the rownames of \code{direction} are unique,
they are copied to the row names of the output.
\cr
In case of error, the function returns \code{NULL}.
}


\references{
\bold{Wikipedia - Support function}\cr
\url{https://en.wikipedia.org/wiki/Support_function}
}

\seealso{
\code{\link{zonoseg}()}
}



