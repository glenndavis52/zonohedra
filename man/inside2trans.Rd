\name{inside2trans}

\alias{inside2trans}


\title{
test points against a 2-transition surface
}

\description{
This function tests points for being inside the 2-transition surface
associated with a zonohedron.
}

\usage{
inside2trans( x, p )
}


\arguments{
\item{x}{
a \bold{zonohedron} object
}

\item{p}{
an Nx3 numeric matrix. The points to be tested are in the rows.
\code{p} can also be a numeric vector that can be converted
to such a matrix, by row.
}
}


\details{
If the surface has no self-intersections,
the the definition of whether a point \code{p} is
"inside" is fairly straightforward:
it is where the linking number of
\code{p} and the surface is non-zero.
In fact, if it is non-zero then it must be +1 or -1.
The \emph{linking number} is analogous the \emph{winding number}
in 2D, for more discussion see \bold{Note}.

Unfortunately, there is currently no test for whether the
surface \emph{has} self-intersections,
For a bad surface with self-intersections,
the linking number might be any integer.
Since there is no such test,
we simply use the same non-zero linking number rule always.

The computed \code{linkingnumber} is returned so that the user
can apply the non-zero rule, or the even-odd rule,
as appropriate for their situation.
These 2 rules are analogous to the two winding number rules
used for polygons in computer graphics,
see \bold{Point in polygon}.

The case where a point is \emph{on} the surface
(i.e. the \code{distance} to the surface is 0) is problematic.
The \code{linkingnumber} is then undefined,
and we currently set \code{inside} to be undefined as well.
Thus \code{inside} should be interpreted as \emph{strictly inside}.
However, in some situations, the user may want to consider
\code{inside} to be \code{TRUE} in this problematic case.
Or the user may want to consider points that are within
a very small epsilon of the surface,
where roundoff might have occurred, to have \code{inside=FALSE}
or \code{inside=NA}.
So the both the computed \code{linkingnumber} and \code{distance} are
returned so the user can use them
to make their own definition of what "inside" means.
}


\value{
\code{inside2trans()} returns a \code{data.frame}
with N rows and these columns:

\item{p}{the given point}

\item{distance}{the distance from the point to the surface.
This is the true Euclidean distance,
and not a "pseudo-distance" as in the case of \code{\link{inside}()}.
If the point is on the surface, the \code{distance} should be 0
up to numerical precision.
}

\item{linkingnumber}{the linking number of the point and the surface.
If the point is \emph{on} the surface (\code{distance==0}),
the (mathematical) linking number is undefined,
and the computed \code{linkingnumber} is \code{NA} (integer).
}

\item{inside}{whether the point is inside the surface; a logical.
This is currently set to \code{linkingnumber != 0}.
If the \code{linkingnumber} is \code{NA} (integer),
then \code{inside} is \code{NA} (logical).
}

\item{timecalc}{the time to do the calculations, in seconds}

If the row names of \code{p} are unique,
they are copied to the row names of the output.
\cr
In case of error, the function returns \code{NULL}.
}

\note{
The standard definition of the
\emph{linking number} of a point and a surface uses intersections with rays,
see the vignette
\href{../doc/transitions.html}{The 2-Transition Subcomplex and the 2-Transition Surface}
for the precise definition.
This is fine in theory, but in practice does not handle well
the case when the ray intersects the boundary of a parallelogram.
So this function uses an integral formula for the degree of a \emph{linking map}
that reduces to summing the signed area of a lot of spherical triangles, see
\bold{Spivak} p. 75 and
\bold{Guillemin and Pollack} p. 188.
}

\references{
\bold{Guillemin, Victor and Alan Pollack}.
\emph{Differential Topology}.
Prentice-Hall.
1974.

\bold{Point in polygon --- Wikipedia, The Free Encyclopedia}.
\url{https://en.wikipedia.org/w/index.php?title=Point_in_polygon&oldid=1139808558}.
2023.

\bold{Spivak, Michael}.
\emph{A Comprehensive Introduction to Differential Geometry}.
Volume 1.
3rd edition.
Publish or Perish.
1999.
}

\seealso{
\code{\link{inside}()}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
