\name{inside}

\alias{inside}
\alias{inside.zonotope}

\title{
test points for being inside a zonotope
}

\description{
Test points for being inside a zonotope.
The boundary points are considered to be inside.
}

\usage{
\method{inside}{zonotope}( x, p )
}


\arguments{
\item{x}{
a \bold{zonotope} object - a \bold{zonohedron}, \bold{zonogon},
or \bold{zonoseg}
}

\item{p}{
an NxM numeric matrix, where M is the dimension of the zonotope.
The points to be tested are in the rows.
\code{p} can also be a numeric vector that can be converted
to such a matrix, by row.
}
}


\details{
The given zonotope is viewed as the intersection of \emph{slabs};
there is a slab for each hyperplane in the simplified matroid.
For each slab a signed distance to boundary of the slab
is computed.
For points outside the slab the distance is positive,
for points on the boundary, the distance is 0,
and for points in the interior of the slab the distance is negative.
The distance to the zonotope is computed as the
maximum over all these slab distances,
and the \emph{critical hyperplane} index is recorded.
A point is inside iff the zonotope distance \eqn{\le} 0.
}


\value{
\code{inside.zonotope()} returns a \code{data.frame} with N rows and these columns:

\item{p}{the given point}
\item{distance}{the signed distance from the point to the zonotope.
When \code{distance} < 0, the point is in the interior, and \code{distance} is the
true Euclidean distance.
For boundary points, \code{distance} is 0.
When \code{distance} > 0, the point is outside, and \code{distance} may be larger
than the true Euclidean distance, so it is really a pseudo-distance.}

\item{inside}{whether the point is inside the zonotope. \code{inside = distance<=0}.}  

\item{idxhyper}{the index of the \emph{critical hyperplane} in the
simplified matroid. This is the index of the slab where the maximum slab distance was taken.
For a \bold{zonoseg} there is only 1 hyperplane (the empty set)
so this is always 1.
}

If the row names of \code{p} are unique,
they are copied to the row names of the output.
\cr
In case of error, the function returns \code{NULL}.
}

\seealso{
\code{\link{inside2trans}()}
}


\examples{
zono1 = zonoseg( c(1,-2,3,0,-3,-4) )

getsegment(zono1)
# [1] -9 4

p = c( 0, -3*pi, pi, 2*pi, getsegment(zono1) )

inside( zono1, p )
#           p inside   distance idxhyper
# 1  0.000000   TRUE -4.0000000        1
# 2 -9.424778  FALSE  0.4247780        1
# 3  3.141593   TRUE -0.8584073        1
# 4  6.283185  FALSE  2.2831853        1
# 5 -9.000000   TRUE  0.0000000        1
# 6  4.000000   TRUE  0.0000000        1
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
