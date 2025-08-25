\name{zonogon-plot}

\alias{plot.zonogon}

\title{
plot a zonogon
}

\description{
Plot a \bold{zonogon} object, with many options.
}

\usage{
\method{plot}{zonogon}( x, orientation=TRUE, normals=FALSE, elabels=FALSE,
                            tiling=FALSE,  tlabels=FALSE,
                            trans2=FALSE, trans2type='both', ... )
}


\arguments{
\item{x}{
a \bold{zonogon} object as returned by the constructor \code{zonogon()}
}

\item{orientation}{
if \code{TRUE} then draw the edges with orientation arrows.
Otherwise just draw unoriented line segments.
}

\item{normals}{
if \code{TRUE} then draw an outward-pointing unit normal on each edge
}

\item{elabels}{
if \code{TRUE} then label each edge with its generator
}

\item{tiling}{
if \code{TRUE} then draw the standard tiling of the zonogon by parallelograms
}

\item{tlabels}{
if \code{TRUE} then label each parallelogram in the tiling with
its generators.
If \code{tiling} is \code{FALSE} then this is ignored.
}


\item{trans2}{
if \code{TRUE} then draw the image of the 2-transition subcomplex
of the unit cube \eqn{[0,1]^n}, in the color blue.
\cr
\code{trans2} can also be an integer 2-vector defining a range
of levels of the subcomplex, where the \emph{level} of a vertex
of the \eqn{n}-cube is the number of 1s.
Both integers should be between \eqn{0} and \eqn{n}.
}

\item{trans2type}{
which part of the 2-transition subcomplex to draw.
It can be \code{'BP'} for \emph{bandpass} (aka Type 1),
\code{'BS'} for \emph{bandstop} (aka Type 2),
or \code{'both'} for both.
}

\item{...}{not used}
}


\details{
A white dot is plotted at the center of the zonogon.
A suitable is title is added above the plot.
If the zonogon was returned from \code{\link{spherize.zonotope}()}
the string \code{"[spherized]"} is added to the title.
}


\value{
The function returns \code{TRUE}; or \code{FALSE} in case of error.
}



\seealso{
\code{\link{zonogon}()},
\code{\link{spherize.zonotope}()}
}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
