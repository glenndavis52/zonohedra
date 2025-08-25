\name{zonohedron-section}

\alias{section}
\alias{section.zonohedron}

\title{
compute the intersection of a plane and the boundary of a zonohedron
}

\description{
Generically, a plane intersects the boundary of a
zonohedron in a convex polygon.
Computing that polygon is the chief goal of this function.
\cr
For a supporting plane, the intersection is a face of the zonohedron,
but in this function only one point of intersection is computed and returned.
}

\usage{
\method{section}{zonohedron}( x, normal, beta, tol=1.e-10, plot=FALSE, ... )
}


\arguments{
\item{x}{
a \bold{zonohedron} object as returned by the constructor \code{zonohedron()}
}

\item{normal}{
a non-zero numeric 3-vector - the normal of all the planes
}

\item{beta}{
a numeric M-vector of line-constants.
The equation of the k'th plane k is: \code{<x,normal> = beta[k]}.
}

\item{tol}{
a small positive number, used as the tolerance for the plane
being considered a supporting plane
}

\item{plot}{
if \code{TRUE}, the polygons formed by the the intersection
of the planes and the boundary of the zonohedron are
\emph{added} to an existing 3D plot of the zonohedron \code{x},
see \code{\link{plot.zonohedron}()}.
The polygons are drawn in red.
}

\item{...}{not used}
}


\details{
Given a plane, the function finds all the facets of the zonohedron
that intersect the plane.
For each such facet it computes a single point of intersection
on the boundary of the facet.
For the parallelograms, the computation is done in a C function;
and for zonogon facets with 3 or more generators,
the computation is done in \code{\link{section.zonogon}()}.
Orientation is handled carefully so that no point appears twice.
The facets are not processed in order around the boundary,
so these points are in no particular order.
They are put in polygon order by sorting them by angle
around a suitable "diameter" of the zonohedron.
}


\value{
\code{section.zonohedron()} returns a list of length M
(=\code{length(beta)}),
and the i'th item in the list is a data frame with these columns:

\item{point}{
a Px3 matrix with the P points of the i'th polygon in the rows.
If the plane does not intersect the zonohedron, then P=0
and the matrix has 0 rows.
If the plane is a supporting plane, the polygon is degenerate
and P=1 and the matrix has 1 row.
The row names of \code{section} are the indexes of the facets
that contain the vertices of the polygon; see \bold{Details}.
}

\item{hyperidx}{
index of a hyperplane that contains the given point
}

\item{sign}{
The sign specifying which of the 2 facets (selected or antipodal)
contains the given point. The value is +1 or -1.
}


The names of the list are readable strings that contain
\code{normal} and \code{beta[i]}.

In case of error, the function returns \code{NULL}.
}

\note{
The package \pkg{rgl} is required for 3D plotting.
}


\seealso{
\code{\link{zonohedron}()},
\code{\link{plot.zonohedron}()},
\code{\link{section.zonogon}()}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
