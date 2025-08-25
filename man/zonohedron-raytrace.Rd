\name{zonohedron-raytrace}

\alias{raytrace}
\alias{raytrace.zonohedron}

\title{
compute the intersection of a ray, based in the interior of
a zonohedron, and the boundary of that zonohedron
}

\description{
The \emph{open ray} with basepoint \eqn{\bold{b}}
and non-zero direction \eqn{\bold{d}} is the set of the form
\eqn{\bold{b} + t \bold{d}} where \eqn{t > 0}.
\cr\cr
This function computes the intersection of an open ray
and the boundary of a zonohedron \eqn{Z}.
The basepoint is normally required to be in the interior
of \eqn{Z}, but an exception is made if the basepoint is 0,
and on the boundary of \eqn{Z},
and the direction points into the interior of \eqn{Z}.
In these two cases the intersection of the open ray
and the boundary of \eqn{Z} is unique.
In the second case, the basepoint is also allowed to be
the sum of all the generators - the so-called \emph{white point}
of \eqn{Z}.
}

\usage{
\method{raytrace}{zonohedron}( x, base, direction, invert=FALSE, plot=FALSE, ... )
}

\arguments{
\item{x}{
a \bold{zonohedron} object as returned by the constructor \code{zonohedron()}
}

\item{base}{
a numeric 3-vector - the basepoint of all the rays.
\code{base} must either be in the interior of \code{x},
or 0 or the \emph{white point} and on the boundary of \code{x}.
}

\item{direction}{
a numeric Mx3 matrix with M non-zero directions in the rows.
The basepoint and these directions define M rays.
\cr
\code{direction} can also be a numeric vector that can be converted
to such a matrix, by row.
}

\item{invert}{
if \code{TRUE}, then compute a point in the unit cube that maps to the point on the boundary
of \code{x} and add it as a column in the returned \code{data.frame}
}

\item{plot}{
if \code{TRUE}, the computed rays, up to the boundary,
are \emph{added} to an existing
plot of the zonohedron \code{x}, see \code{\link{plot.zonohedron}()}.
The segments are drawn in the color red.
If there is no open plot, a warning is issued.
}

\item{...}{not used}
}


\details{
If \code{plot} is \code{TRUE}, the rays are drawn
with \code{rgl::segments3d()}.
}


\value{
\code{raytrace.zonohedron()} returns a \code{data.frame} with M rows and these columns:

\item{base}{
the given basepoint - this is the same in every row
}

\item{direction}{
the given direction
}

\item{facetidx}{
the index of the facet (a zonogon) where ray exits the zonohedron
}

\item{sign}{
of the facet, either +1 or -1
}

\item{tmax}{
ray parameter of the intersection with the exit facet, always positive
}
\item{point}{
the point on the boundary; the intersection of the ray and the facet
}
\item{timetrace}{
the computation time, in seconds
}

And if \code{invert} is \code{TRUE}, then these columns are added:
\item{distance}{signed distance to the boundary of \code{x}}
\item{pcube}{a point in the unit cube that maps to \code{boundary}}
\item{transitions}{the number of transitions in \code{pcube} - a non-negative even integer}

If \code{base} and \code{direction} in a row cannot be
processed, the rest of the row is \code{NA}.

If the row names of \code{direction} are unique,
they are copied to the row names of the output.
\cr\cr
In case of error, the function returns \code{NULL}.
}

\note{
The package \pkg{rgl} is required for 3D plotting.
}


\seealso{
\code{\link{zonohedron}()},
\code{\link{plot.zonohedron}()},
\code{\link{section.zonohedron}()},
\code{\link{invertboundary}()},
\code{\link{raytrace.zonogon}()}
}

\examples{
#   make a regular prism, a regular 20-gon extruded 1 unit along z-axis
rp10 = regularprism( 10 )

#   make 7 random directions
set.seed(0)
dir = matrix(rnorm(7*3),7,3)

#  use basepoint in the interior of the zonohedron
raytrace( rp10, c(0.5,0.5,0.5), dir )
#   base.1 base.2 base.3  direction.1  direction.2  direction.3 facetidx sign      tmax   ...
# 1    0.5    0.5    0.5  1.262954285 -0.294720447 -0.299215118        1    1 1.6710386   ...
# 2    0.5    0.5    0.5 -0.326233361 -0.005767173 -0.411510833        1    1 1.2150348   ...
# 3    0.5    0.5    0.5  1.329799263  2.404653389  0.252223448        6   -1 0.8724774   ...
# 4    0.5    0.5    0.5  1.272429321  0.763593461 -0.891921127        1    1 0.5605877   ...
# 5    0.5    0.5    0.5  0.414641434 -0.799009249  0.435683299        1   -1 1.1476226   ...
# 6    0.5    0.5    0.5 -1.539950042 -1.147657009 -1.237538422        1    1 0.4040279   ...
# 7    0.5    0.5    0.5 -0.928567035 -0.289461574 -0.224267885        1    1 2.2294766   ...

#  use basepoint 0 on the boundary of the zonohedron
#  note that only 2 directions point into the interior
raytrace( rp10, c(0,0,0), dir )
#   base.1 base.2 base.3  direction.1  direction.2  direction.3 facetidx sign     tmax   ...
# 1      0      0      0  1.262954285 -0.294720447 -0.299215118       NA   NA       NA   ...
# 2      0      0      0 -0.326233361 -0.005767173 -0.411510833       NA   NA       NA   ...
# 3      0      0      0  1.329799263  2.404653389  0.252223448        6   -1 1.128580   ...
# 4      0      0      0  1.272429321  0.763593461 -0.891921127       NA   NA       NA   ...
# 5      0      0      0  0.414641434 -0.799009249  0.435683299        1   -1 2.295245   ...
# 6      0      0      0 -1.539950042 -1.147657009 -1.237538422       NA   NA       NA   ...
# 7      0      0      0 -0.928567035 -0.289461574 -0.224267885       NA   NA       NA   ...
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

