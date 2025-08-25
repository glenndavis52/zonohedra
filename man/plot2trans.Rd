\name{plot2trans}

\alias{plot2trans}

\title{
plot the 2-transition surface associated with a zonohedron
}

\description{
The 2-transition surface has the topology of a sphere and is contained in the zonohedron.
All the facets are parallelograms.
The surface is centrally symmetric, with the same center as the zonohedron.
The surface may have self-intersections.
}

\usage{
plot2trans(  x, type='ef', ecol='black', econc=FALSE,
                 fcol='yellow', falpha=0.5, level=NULL,
                 normals=FALSE, both=TRUE, bgcol="gray40", add=FALSE, ... )
}

\arguments{
\item{x}{
a \bold{zonohedron} object as returned by the constructor \code{zonohedron()}
}

\item{type}{
a character string with what parts to draw.
If \code{type} contains \code{'e'}, then draw the edges.
If \code{type} contains \code{'f'}, then draw filled facets.
If \code{type} contains \code{'p'}, then draw points
at the centers of the facets.
}


\item{ecol}{
The color to use when drawing the edges.
}

\item{econc}{
If \code{TRUE} then draw the concave edges in red,
and with extra thickness
}

\item{fcol}{
The color to use when drawing the facets.
}

\item{falpha}{
The opacity to use when drawing the facets.
}

\item{level}{
An integer vector which is a subvector of \code{0:(M-2)},
where \code{M} is the number of simplified generators.
Only the facets and edges at the specified levels are drawn.
When \code{level=NULL} then \emph{all} facets and edges are drawn.
This argument does not affect the drawing of points.
}

\item{normals}{
If \code{TRUE} then draw the unit facet normals.
}

\item{both}{
if \code{FALSE} then draw only one half of the centrally symmetric surface.
Otherwise draw both halves (the default).
}

\item{bgcol}{
the background color
}

\item{add}{
If \code{TRUE} then add to the current 3D plot.
If there is no current 3D plot, it is an error.
}

\item{...}{not used}
}


\details{
Facets and regular edges are drawn with \code{rgl::quads3d()}.
Concave edges are drawn with \code{rgl::segments3d()}.
Points are drawn with \code{rgl::points3d()}.
}


\value{
The function returns \code{TRUE}; or \code{FALSE} in case of error.
}


\note{
The package \pkg{rgl} is required for 3D plots.
A large black point is drawn at 0, 
a large white point at the "white point",
and a 50\% gray point at the center.
A line from the black point to the white point is also drawn.
}


\seealso{
\code{\link{zonohedron}()},
\code{\link{plothighertrans}()},
\code{\link{plot.zonohedron}()}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
