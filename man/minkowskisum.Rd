\name{minkowskisum}

\alias{minkowskisum}
\alias{minkowskisum.zonotope}
\alias{\%+\%}
\alias{\%+\%.zonotope}

\title{
Minkowski sum of Two zonotopes
}

\description{
A zonotope can be viewed as a Minkowski sum of line segments,
with one endpoint at 0.
Therefore, the Minkowski sum of two zonotopes (in the same dimension)
is also a zonotope.
}

\usage{
\S3method{minkowskisum}{zonotope}( zono1, zono2, e0=0, e1=1.e-6, e2=1.e-10, ground=NULL, ... )

\S3method{\%+\%}{zonotope}(zono1,zono2)
}


\arguments{
\item{zono1}{
a zonotope object - a \bold{zonohedron}, a \bold{zonogon}, or a \bold{zonoseg}
}
\item{zono2}{
a zonotope object with the same dimension as \code{zono1}
}
\item{e0}{see \code{\link{zonohedron}()}}

\item{e1}{see \code{\link{zonohedron}()}}

\item{e2}{see \code{\link{zonohedron}()}}

\item{ground}{the ground set of the returned zonotope.
If \code{ground} is \code{NULL}, it is set to the ground set of \code{zono1}
followed by the ground set of \code{zono2} translated sufficiently
to not intersect that of \code{zono1}.
}

\item{...}{not used}
}

\value{
\code{minkowskisum()} returns a zonotope of the same dimension
as \code{zono1} and \code{zono2}.\cr
\code{\%+\%} is a more convenient binary operator that calls
\code{minkowskisum()}, but without the flexibility of
the extra arguments.
\cr
In case of error, the function returns \code{NULL}.
}

\details{
After verifying that \code{zono1} and \code{zono2} are the same dimension,
it takes the 2 matrices, \code{cbind}s them,
and passes the new matrix to the appropriate constructor,
along with the other arguments.
There are no special optimizations.
}



\references{
\bold{Zonohedron - Wikipedia}.\cr
\url{https://en.wikipedia.org/wiki/Zonohedron}.
}


\seealso{
\code{\link{zonohedron}()},
\code{\link{zonogon}()},
\code{\link{zonoseg}()}
}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
