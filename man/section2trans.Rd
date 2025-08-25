\name{section2trans}

\alias{section2trans}

\title{
compute the intersection of a plane and the
2-transition surface associated with a zonohedron
}

\description{
In general, the 2-transition surface may be highly non-convex,
possibly with self-intersections.
The intersection of a plane and the 2-transition surface is
a union of polygons,
possibly with self-intersections and intersecting each other.
This function computes one of those polygons.
If there are other polygons, it issues a warning
and does not try to compute them.
}

\usage{
section2trans( x, normal, beta, invert=FALSE, plot=FALSE, tol=1.e-12, ... )
}


\arguments{
\item{x}{
a \bold{zonohedron} object as returned by the constructor \code{zonohedron()}
}

\item{normal}{
a non-zero numeric 3-vector - the normal of all the planes
}

\item{beta}{
a numeric M-vector of plane constants.  
The equation of the k'th plane k is: \code{<x,normal> = beta[k]}.
}

\item{invert}{
if \code{TRUE}, then compute a point in the unit cube that maps to the point on the 2-transition surface,
and add it as a column in the returned \code{data.frame}
}

\item{plot}{
if \code{TRUE}, the polygons formed by the the intersection 
of the planes and the 2-transition surface.
\emph{added} to an existing 3D plot of the zonohedron \code{x},
see \code{\link{plot.zonohedron}()}.
The polygons are drawn in red.
}

\item{tol}{
a small positive number, used as the tolerance for the plane
intersecting the interior of each parallelogram, see \bold{Details}.
}

\item{...}{not used}
}


\details{
The function is designed for the situation when the
intersection of a plane and the surface is a single polygon.

Given a plane, the function finds all the parallelograms of the surface
whose interiors intersect the plane.
Each intersection is a line segment.
For each parallelogram it associates one of the endpoints of the segment.
The parallelograms are put in polygon order by picking an
arbitrary one as the starting point,
and then "marching" from one to the next using
the canonical parallelogram adjacency relation.
After returning to the starting point, if there are other parallelograms
remaining, it means that there are other polygons
in the section and a warning is issued.
}

\value{
\code{section2trans()} returns a list of length M
(=\code{length(beta)}),
and the i'th item in the list is a data frame with these columns:

\item{point}{
a Px3 matrix with the P points of the i'th polygon in the rows.
If the plane does not intersect the 2-transition surface, then P=0
and the matrix has 0 rows.
The row names of \code{point} are the indexes of the facets
that contain the vertices of the polygon; see \bold{Details}.
}

\item{gndpair}{
the 2 indexes from the ground set that generates the parallelogram
containing \code{point}.
See \bold{Details} for a description of the "marching parallelogram"
procedure.
}

And if \code{invert} is \code{TRUE}, then this column is added:
\item{pcube}{a point in the unit cube that maps to \code{point}.  This point in the cube always has 2 transitions.}

The names of the returned list are readable strings that contain
\code{normal} and \code{beta[i]}.

In case of error, the function returns \code{NULL}.
}



\note{
The package \pkg{rgl} is required for 3D plotting.
}



\seealso{
\code{\link{zonohedron}()},
\code{\link{plot.zonohedron}()},
\code{\link{section.zonohedron}()},
\code{\link{raytrace2trans}()}
}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
