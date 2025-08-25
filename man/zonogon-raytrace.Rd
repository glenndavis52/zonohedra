\name{zonogon-raytrace}

\alias{raytrace.zonogon}

\title{
compute the intersection of a ray, based in the interior of
a zonogon, and the boundary of that zonogon
}

\description{
The \emph{open ray} with basepoint \eqn{\bold{b}}
and non-zero direction \eqn{\bold{d}} is the set of the form
\eqn{\bold{b} + t \bold{d}} where \eqn{t > 0}.
\cr\cr
This function computes the intersection of an open ray
and the boundary of a zonogon \eqn{Z}.
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
\method{raytrace}{zonogon}( x, base, direction, plot=FALSE, ... )
}

\arguments{
\item{x}{
a \bold{zonogon} object as returned by the constructor \code{zonogon()}
}

\item{base}{
a numeric 2-vector - the basepoint of all the rays.
\code{base} must either be in the interior of \code{x},
or 0 or the \emph{white point} and on the boundary of \code{x}.
}

\item{direction}{
a numeric Mx2 matrix with M non-zero directions in the rows.
The basepoint and these directions define M rays.
\cr
\code{direction} can also be a numeric vector that can be converted
to such a matrix, by row.
}

\item{plot}{
if \code{TRUE}, the computed rays, up to the boundary,
are \emph{added} to an existing
plot of the zonogon \code{x}, see \code{\link{plot.zonogon}()}.
The segments are drawn in the color red.
If there is no open plot, a warning is issued.
}

\item{...}{not used}
}


\value{
\code{raytrace.zonogon()} returns a \code{data.frame} with M rows and these columns:

\item{base}{
the given basepoint - this is the same in every row
}

\item{direction}{
the given direction
}

\item{facetidx}{
the index of the facet (an edge) where ray exits the zonogon
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

If \code{base} and \code{direction} in a row cannot be
processed, the rest of the row is \code{NA}.

If the row names of \code{direction} are unique,
they are copied to the row names of the output.
\cr\cr
In case of error, the function returns \code{NULL}.
}


\seealso{
\code{\link{zonogon}()},
\code{\link{plot.zonogon}()},
\code{\link{section.zonohedron}()}
}

\examples{
#   make a zonogon with 5 generators
pz20 = polarzonogon( 20, 5 )

#   make 4 random directions
set.seed(0)
dir = matrix(rnorm(4*2),4,2)

#  use basepoint in the interior of the zonogon
raytrace( pz20, c(0.5,0.5), dir )
#   base.1 base.2 direction.1 direction.2 facetidx sign      tmax boundary.1 boundary.2    timetrace
# 1    0.5    0.5   1.2629543   0.4146414        4   -1 2.0503073 3.08944438 1.35014236 7.680000e-05
# 2    0.5    0.5  -0.3262334  -1.5399500        1   -1 0.3246859 0.39407664 0.00000000 4.649995e-05
# 3    0.5    0.5   1.3297993  -0.9285670        2   -1 0.4868719 1.14744192 0.04790678 4.310103e-05
# 4    0.5    0.5   1.2724293  -0.2947204        2   -1 0.9354693 1.69031851 0.22429808 4.149997e-05


#  use basepoint at 0 - on the boundary of the zonogon
raytrace( pz20, c(0,0), dir )
#   base.1 base.2 direction.1 direction.2 facetidx sign     tmax boundary.1 boundary.2 timetrace
# 1      0      0   1.2629543   0.4146414        4   -1 2.192481  2.7690037  0.9090936 0.0001216
# 2      0      0  -0.3262334  -1.5399500       NA   NA       NA         NA         NA        NA
# 3      0      0   1.3297993  -0.9285670       NA   NA       NA         NA         NA        NA
# 4      0      0   1.2724293  -0.2947204       NA   NA       NA         NA         NA        NA
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
