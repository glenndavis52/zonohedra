\name{lintransform}

\alias{lintransform}
\alias{lintransform.matroid}
\alias{lintransform.zonohedron}
\alias{lintransform.zonogon}


\title{
linear transformations of zonotopes, and vector matroids
}

\description{
These functions perform straightforward linear transformations
on the generators of a zonotope,
and the column vectors of a vector matroid
}

\usage{
\method{lintransform}{zonohedron}( x, W )
\method{lintransform}{zonogon}( x, W )
\method{lintransform}{matroid}( x, W )
}


\arguments{
\item{x}{
\code{x} can be a vector matroid object, as returned from the constructor \code{\link{matroid}()} that takes a matrix as input.
}

\item{W}{
An invertible matrix that matches the rank of \code{x}.
This invertibility is verified.
\code{W} can also be a scalar; it is then replaced by that scalar
multiplied by the identity matrix of the appropriate rank.
}

}


\value{
If \code{x} is a zonohedron (or zonogon),
\code{lintransform(x)} returns the zonohedron (or zonogon)
whose generators are the generators of \code{x} with the matrix
\code{W} applied on the left side.
This function is optimized - it is \emph{not} necessary
to transform the generators and start all over again.
\cr\cr
If \code{x} is a vector matroid, \code{lintransform(x)} returns the matroid
whose generators are the generators of \code{x} with the matrix
\code{W} applied on the left side.
If \code{x} is a matroid, but *not* a vector matroid,
it returns the original matroid and prints a warning message.
\cr\cr
In case of error, e.g. invalid \code{W},
the function prints an error message and returns \code{NULL}.
}



\references{
Matroid - Wikipedia.
\url{https://en.wikipedia.org/w/index.php?title=Matroid&oldid=1086234057}

}


\seealso{
\code{\link{rank}()}, 
\code{\link{matroid}()}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
