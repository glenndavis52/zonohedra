\name{zonohedron-plot}

\alias{plot.zonohedron}

\title{
plot a zonohedron
}

\description{
Plot a \bold{zonohedron} object in 3D, with many options.
}

\usage{
\method{plot}{zonohedron}( x, type='e', pcol=NULL, ecol=NULL, ewd=3, etcol=NA,
        fcol=NULL, falpha=1, normals=FALSE, bgcol="gray40", both=TRUE, ... )
}


\arguments{
\item{x}{
a \bold{zonohedron} object as returned by the constructor \code{zonohedron()}
}

\item{type}{
a string of letter with what parts to draw.
If \code{type} contains an \code{'p'}, then draw a point at the center
of each facet.
If \code{type} contains an \code{'e'}, then draw the edges.
If \code{type} contains an \code{'f'}, then draw filled facets.
}

\item{pcol}{
The color to use when drawing points.
It can be a vector of 2 colors, and then when \code{both} is \code{TRUE},
the first color is used for one half,
and the second color is used for the antipodal half.
When \code{pcol} is \code{NULL},
it is set to \code{c('black','red')}.
}

\item{ecol}{
A vector of colors to use when drawing the edges.
Let N be the number of \emph{simplified} generators of the zonohedron.
Each edge is parallel to exactly one of the generators,
so this divides the edges into N \emph{zones}, or \emph{belts}.
\code{ecol} can be a vector of N colors, one for each zone.
If \code{ecol} is shorter than N, it is extended to length N using the last color.
If \code{ecol} is longer than N, the extra colors are ignored.
If \code{ecol} is \code{NULL}, it is set to \code{rainbow(N)}.
}

\item{ewd}{
width of the edges, in pixels
}

\item{etcol}{
color of the \emph{tiling edges}, for the standard tiling
of the facets by parallelograms.
This only applies to facets that are \emph{not} parallelograms.
The default \code{etol=NA} means do not draw these edges.
}

\item{fcol}{
A vector of colors to use when drawing the facets.
The 1st color is used for parallelograms,
the next color for hexagons, etc.
For facets with more edges than colors available, the last color is used.
If \code{fcol} is \code{NULL}, it is set to
\code{c( 'blue', 'red', 'yellow', 'green', 'orange', 'purple' )}.
}

\item{falpha}{
opacity of the facets
}


\item{normals}{
if \code{TRUE} then draw an outward-pointing unit normal from each facet
}

\item{bgcol}{
the background color
}


\item{both}{
if \code{FALSE} then draw only one half of the centrally symmetric boundary.
Otherwise draw both halves.
This affects points, edges, and facets.
}


\item{...}{not used}
}


\details{
Points are drawn with \code{rgl::points3d()}.
Edges are drawn with \code{rgl::segments3d()}.
Edges of the tiles are drawn with \code{rgl::quads3d()}.
Facets are drawn with \code{rgl::quads3d()};
facets with more than 4 edges are split into trapezoids.
Facet normals are drawn with \code{rgl::arrow3d()}.
}


\value{
The function returns \code{TRUE}; or \code{FALSE} in case of error.
}


\note{
The package \pkg{rgl} is required for 3D plots.
A large black point is drawn at 0,
a 50\% gray point at the center,
and a large white point at the "white point" (which is 2*center).

A line from the black point to the white point is also drawn.
}


\seealso{
\code{\link{zonohedron}()},
\code{\link{spherize.zonotope}()}
}




% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
