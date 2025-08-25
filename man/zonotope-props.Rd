\name{zonotope-props}

\alias{is_salient}
\alias{is_pointed}

\alias{is_salient.zonotope}
\alias{is_pointed.zonotope}

\title{zonotope properties}

\description{
Get some important boolean properties of a zonotope.
\cr\cr
\emph{pointed} means that 0 is a vertex of the zonotope.
\emph{salient} means that 0 is in the boundary of the zonotope.
So pointed implies salient, but not the reverse.
\cr\cr
A zonotope has an associated \emph{convex cone} - allow the
coefficients of the generators to be any non-negative numbers.
For convex cones, \emph{pointed} means that the cone
is in an open linear halfspace (except for 0).
And \emph{salient} means that the cone is in a
closed linear halfspace (the cone may contain a line).
\cr\cr
In terms of generators (of both zonotopes and convex cones),
\emph{pointed} means that the generators
are in an open linear halfspace (except for 0 generators).
And \emph{salient} means that the generators
are in a closed linear halfspace.
}

\usage{
\method{is_pointed}{zonotope}( x )

\method{is_salient}{zonotope}( x )
}


\arguments{
\item{x}{
a zonotope object -  a \bold{zonohedron}, a \bold{zonogon}, or a \bold{zonoseg}
}
}

\value{
\code{TRUE} or \code{FALSE}
}

\details{
For a \bold{zonohedron}, if 0 is in the interior of an edge or a facet,
then the zonohdron is salient but not pointed.
\cr
For a \bold{zonogon}, if 0 is in the interior of an edge,
then the zonogon is salient but not pointed.
\cr
For a \bold{zonoseg}, both \emph{pointed} and \emph{salient} are equivalent
to 0 being a boundary point.
And this is equivalent to all the non-zero generators having
the same sign  (all negative or all positive).
}



\references{
\bold{Zonohedron - Wikipedia}.\cr
\url{https://en.wikipedia.org/wiki/Zonohedron}
}


\seealso{
\code{\link{zonohedron}()},
\code{\link{zonogon}()},
\code{\link{zonoseg}()}
}

\examples{
zono1 = zonoseg( c(1,-2,3,0,-3,-4) )

is_pointed( zono1 )
# [1] FALSE

is_salient( zono1 )
# [1] FALSE
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
