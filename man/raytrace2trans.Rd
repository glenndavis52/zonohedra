\name{raytrace2trans}

\alias{raytrace2trans}

\title{
compute the intersection of a ray and the 2-transition surface
associated with a zonohedron
}

\description{
The \emph{open ray} with basepoint \eqn{\bold{b}}
and non-zero direction \eqn{\bold{d}} is the set of the form
\eqn{\bold{b} + t \bold{d}} where \eqn{t > 0}.
\cr\cr
This function computes the intersection of an open ray
and the 2-transition surface associated with a zonohedron.
The linking number of the surface and \eqn{\bold{b}} must be \eqn{\pm 1}.
This is verified at the beginning, and if not true, then it is an error.
The linking number condition implies that an intersection exists for every
ray based at \eqn{\bold{b}}.
Note also that the condition implies that \eqn{\bold{b}} is not on the surface.
For discussion of uniqueness, see \bold{Details}.
For the definition of \emph{linking number} see
\href{../doc/transitions.html}{The 2-Transition Subcomplex and the 2-Transition Surface}.

The 2-transition surface is a union of parallelograms.
The surface is symmetric about the center of the zonhedron,
so each parallelogram has an antipodal parallelogram.
Each parallelogram is specified by an ordered pair of
distinct generators from the
\emph{simplified} matroid associated with the zonohedron.
Thus, if there are \eqn{N} generators, there are \eqn{N(N-1)} parallelograms.
Swapping the generators of a parallelogram changes it
to the antipodal parallelogram.

The 2-transition surface has two \emph{poles} -
the point 0 and the sum of all the generators.
It is OK for the ray to pass through one of these poles.
}

\usage{
raytrace2trans( x, base, direction, invert=FALSE, plot=FALSE, tol=1.e-12, ...  )
}

\arguments{
\item{x}{
a \bold{zonohedron} object as returned by the constructor \code{zonohedron()}
}

\item{base}{
a numeric 3-vector - the basepoint of all the rays.
The surface must be strictly starshaped at \code{base},
and this is verified.
}

\item{direction}{
a numeric Mx3 matrix with M non-zero directions in the rows.
The basepoint and these directions define M rays.
\cr
\code{direction} can also be a numeric vector that can be converted
to such a matrix, by row.
}

\item{invert}{
if \code{TRUE}, then compute a point in the unit cube that maps to the point on the 2-transition surface associated with \code{x},
and add it as a column in the returned \code{data.frame}
}

\item{plot}{
if \code{TRUE}, the computed rays, up to the boundary,
are \emph{added} to an existing
plot of the zonohedron \code{x}, see \code{\link{plot.zonohedron}()}.
The segments are drawn in the color red.
If there is no open 3D plot, a warning is issued.
}

\item{tol}{
the tolerance for being strictly starshaped,
and for intersection with a \emph{pole}.
}

\item{...}{not used}
}


\details{
The function is designed for the situation when the intersection
of the ray and the surface exists and is unique.
This is guaranteed for all ray directions \eqn{\bold{d}}
when the surface is strictly starshaped at \eqn{\bold{b}}.
This condition is checked at the beginning of the function,
and if false then a warning is issued that the intersection point
may not be unique.
For the definition of \emph{strictly starshaped} see
\href{../doc/transitions.html}{The 2-Transition Subcomplex and the 2-Transition Surface}.

For finding a parallelogram of intersection,
a brute-force search is used; all parallelograms are searched until
the first one that intersects the ray is found.
To speed things up, the 3D problem is reduced to 2D,
and the search is programmed in plain C.

If \code{plot} is \code{TRUE}, the rays are drawn in red
using \code{rgl::segments3d()} and \code{rgl::points3d()}.
}


\value{
\code{raytrace2trans()} returns a \code{data.frame} with M rows and these columns:

\item{base}{
the given basepoint - this is the same in every row
}

\item{direction}{
the given direction
}

\item{gndpair}{
the 2 generators of the parallelogram that the ray intersects, 
taken from the ground set of the simplified matroid.
If the ray passes through a pole, both of these are \code{NA}.
}

\item{alpha}{
the 2 coordinates of the intersection point within the parallelogram
}

\item{tmax}{
ray parameter of the intersection with the parallelogram, always positive
}

\item{point}{
the point on the surface; the intersection of the ray and the parallelogram
}

\item{iters}{
the number of parallelograms searched, until the desired one was found.
If the ray intersects a pole, this is 0.
}

\item{timetrace}{
the computation time for the given ray, in seconds.
This does not include the initial preprocessing time.
}

And if \code{invert} is \code{TRUE}, then this column is added:
\item{pcube}{a point in the unit cube that maps to \code{point}.  This point in the cube always has 2 transitions.}


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
\code{\link{section2trans}()},
\code{\link{raytrace.zonohedron}()}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
