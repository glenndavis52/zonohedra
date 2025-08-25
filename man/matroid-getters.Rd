\name{matroid-getters}

\alias{gethyperplane}
\alias{getground}
\alias{getmatrix}
\alias{getloop}
\alias{getmultiple}

\alias{gethyperplane.matroid}
\alias{getground.matroid}
\alias{getmatrix.matroid}
\alias{getloop.matroid}
\alias{getmultiple.matroid}


\title{
matroid get functions
}

\description{
get some important members of a matrix
}

\usage{
\method{getground}{matroid}( x )

\method{gethyperplane}{matroid}( x )

\method{getmatrix}{matroid}( x )

\method{getloop}{matroid}( x )

\method{getmultiple}{matroid}( x )
}


\arguments{
\item{x}{
a matroid object, as returned from the constructor
\code{\link{matroid}()}
}
}

\value{
\code{getground()} returns an vector of positive integers in strictly
increasing order = the ground set of the matroid \code{x}.
\cr\cr
\code{gethyperplane()} returns a list of vectors of positive integers 
= the hyperplanes of the matroid.
\cr
If \code{x} is the simplification of an "original matroid",
the \code{"lmdata"} attribute of the returned list is set to the
\emph{loop} and \emph{multiple group} data of the "original hyperplanes".
These hyperplanes can be recovered using \code{\link{unsimplify}()}.
\cr
If \code{x} was constructed from a matrix,
these hyperplanes are sorted in decreasing order by length.
The non-trivial hyperplanes come first, followed by the trivial hyperplanes.
A hyperplane is \emph{trivial} iff it is independent in the matroid.
For a matroid of rank 3, a hyperplane is trivial iff it has 2 points.
\cr\cr
\code{getmatrix()} returns the matrix passed to the
\code{\link{matroid.matrix}()} constructor,
or \code{NULL} if the list constructor was used.
The column names are labeled with the ground set.
\cr\cr
\code{getloop()} returns an integer vector with the loops of \code{x}.
If \code{x} is simple, it is the empty vector.
\cr\cr
\code{getmultiple()} returns a list of integer vectors - the
multiple groups of \code{x}.
If \code{x} is simple, it is the empty list.
}


\references{
Matroid - Wikipedia.
\url{https://en.wikipedia.org/w/index.php?title=Matroid&oldid=1086234057}
}


\seealso{
\code{\link{rank}()}, 
\code{\link{simplify}()},
\code{\link{unsimplify}()},
\code{\link{getsimplified}()},
}

\examples{
# construct a classic matroid with 7 points, but assign an unusual ground set
mat = matroid( classics.genlist[['TRD']], ground=11:17 )

getmatrix( mat )

##       11 12 13 14       15       16       17
##  [1,]  1  1  1  1 1.732051 0.000000 0.000000
##  [2,]  1  1 -1 -1 0.000000 1.732051 0.000000
##  [3,]  1 -1  1 -1 0.000000 0.000000 1.732051
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
