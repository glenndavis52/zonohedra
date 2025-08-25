\name{plothighertrans}

\alias{plothighertrans}

\title{
plot abundant and deficient parallelograms
}

\description{
The 2-transition surface associated with a zonohedron
is a topological sphere and is contained in the zonohedron.
The surface is centrally symmetric, with the same center as the zonohedron.
The surface may have self-intersections.
For this function, the surface is required to be strictly starshaped at
the center.
For the definition of \emph{strictly starshaped} see the vignette
\href{../doc/transitions.html}{The 2-Transition Subcomplex and the 2-Transition Surface}.


The 2-transition surface is a union of parallelograms.
Each parallelogram has a unit normal that defines a linear functional.
\cr
If a 2-transition parallelogram is in the \emph{interior} of the zonohedron
then the functional is not maximized on the parallelogram,
and there is a corresponding similar parallelogram
on the boundary of the zonohedron where the functional *is* maximized.
The first parallelogram (in the surface) is called \emph{deficient} because
the functional is not maximized,
and the second parallelogram (in the boundary)
is called \emph{abundant} because the number of corresponding
transitions across this parallelogram is more than 2.
\cr
If the 2-transition parallelogram is on the boundary,
then it is called \emph{coincident}.
The coincident parallelograms are ignored and
not drawn in this function.

Because of this 1-1 correspondence between deficient parallelograms (in the 2-transition surface)
and the abundant parallelograms (in the boundary of the zonohedron),
the area of these two surfaces are the same.
}

\usage{
plothighertrans( x, abalpha=1, defcol='green', defalpha=0, ecol=NA,
                  connections=FALSE, bgcol="gray40", both=TRUE, ...  )
}


\arguments{
\item{x}{
a \bold{zonohedron} object as returned by the constructor \code{zonohedron()}
}

\item{abalpha}{
The opacity to use when drawing the abundant parallelograms.
If \code{abalpha=0} then they are not drawn.
}

\item{defcol}{
The color to use when drawing the deficient parallelograms
}

\item{defalpha}{
The opacity to use when drawing the deficient parallelograms.
If \code{defalpha=0} (the default), then they are not drawn.
}

\item{ecol}{
The color to use when drawing the edges.
If \code{ecol=NA} (the default), then they are not drawn.
}

\item{connections}{
If \code{TRUE} then draw segments between centers of the
deficient parallelograms in the 2-transition surface,
and centers of the the corresponding abundant parallelograms
in the zonohedron boundary
}

\item{bgcol}{
the background color
}

\item{both}{
if \code{FALSE} then draw only one half of the centrally symmetric boundary.
Otherwise draw both halves.
This affects edges and parallelograms.
}

\item{...}{not used}
}

\section{WARNING}{
This function currently only works when the 2-transition surface
is starshaped at the center.
This excludes many of the classic zonohedra.
}


\details{
Connections are drawn with \code{rgl::segments3d()} and 
\code{rgl::points3d()} .
parallelograms and edges are drawn with \code{rgl::quads3d()}.
Both parallelograms and edges are drawn unlit (\code{lit=FALSE}).
The parallelograms are colored by the number of transitions
using the color codes in \bold{Burns}, up to 10 transitions.

A large black point is drawn at 0, 
a large white point at the "white point",
and a 50\% gray point at the center.
A line from the black point to the white point is also drawn.
}

\value{
The function returns \code{TRUE} when successful; 
or \code{FALSE} in case of error.
}


\note{
The package \pkg{rgl} is required for 3D plots.
}

\references{
Scott A Burns.
\bold{The location of optimal object colors with more than two transitions}.
Color Research & Application.
Vol. 46.
No. 6.
pp 1180-1193.
2021.
}

\seealso{
\code{\link{zonohedron}()},
\code{\link{plot.zonohedron}()},
\code{\link{plot2trans}()}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
