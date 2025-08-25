\name{zonotope-getters}


\alias{getmatroid}
\alias{getcenter}

\alias{getmatrix.zonotope}
\alias{getmatroid.zonotope}
\alias{getcenter.zonotope}


\title{
zonotope get functions
}

\description{
get some important members of a zonotope
}

\usage{
\method{getmatrix}{zonotope}( x )

\method{getmatroid}{zonotope}( x )

\method{getcenter}{zonotope}( x )
}


\arguments{
\item{x}{
a zonotope object -  a \bold{zonohedron}, a \bold{zonogon}, or a \bold{zonoseg}
}
}

\value{
\code{getmatrix()} returns the matrix originally used to construct
the zonotope \code{x}.
\cr\cr
\code{getmatroid()} returns the matroid (possibly nonsimple)
constructed from the matrix
\cr\cr
\code{getcenter()} returns the center of the zonotope;
which is also the center of radial symmetry.
If \code{x} is an object-color solid, the center corresponds
to the 50\% graypoint.  For the whitepoint multiply by 2.
}



\seealso{
\code{\link{zonohedron}()},
\code{\link{zonogon}()},
\code{\link{zonoseg}()}
}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
