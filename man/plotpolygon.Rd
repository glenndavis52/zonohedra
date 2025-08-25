\name{plotpolygon}

\alias{plotpolygon}

\title{
plot the \emph{generator polygon} associated with a pointed zonohedron
}

\description{
A zonohedron is \emph{pointed} iff there is a vector \bold{n}
so the inner product of all the zonohedron generators with \bold{n} is positive.
In other terminology, it is \emph{pointed} iff there is an open halfspace
that contains all the generators.
\cr
When \bold{n} exists, a neighborhood of 0 can be cut by a plane
orthogonal to \bold{n} and the intersection is a polygon.
Since \bold{n} is not unique, the polygon is only unique up
to a 2D projective transformation.
}

\usage{
plotpolygon( x, normal=NULL, points=TRUE, labels=TRUE )
}

\arguments{
\item{x}{
a \bold{zonohedron} object as returned by the constructor \code{zonohedron()}.
It must be pointed.
}

\item{normal}{
the vector \bold{n} to use - a non-zero numeric vector of length 3.
If it is given, the validity is checked and if invalid it is an error.
\cr
If it is \code{NULL}, a few canonical normals are first tested for validity.
If they are invalid, then a valid one is computed.
}

\item{points}{
If \code{TRUE} then draw the vertices of the polygon
}

\item{labels}{
If \code{TRUE} then draw labels,
taken from the ground set of the simplified matroid, near the vertices
}

}


\details{
The selected normal vector \bold{n} is added to the title of the plot.
}


\value{
The function returns \code{TRUE}; or \code{FALSE} in case of error.
}


\seealso{
\code{\link{zonohedron}()},
\code{\link{plot.zonohedron}()}
}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
