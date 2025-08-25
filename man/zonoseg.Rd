\name{zonoseg}

\alias{getsegment}
\alias{getsegment2trans}

\alias{zonoseg}
\alias{getsegment.zonoseg}
\alias{getsegment2trans.zonoseg}
\alias{print.zonoseg}

\title{
zonoseg construction, querying, and printing
}

\description{
Construct a zonoseg from a numeric matrix with one row.
\cr\cr
A \emph{zonoseg} ("zonotope" + "segment") is my own personal term
for a 1-dimensional zonotope.
I could not find an alternative term.
It is a linear image of the unit cube \eqn{[0,1]^n} in the real numbers,
and a compact segment of reals.
The order of the generators has no effect on the zonoseg.
\cr\cr
The image of the \emph{2-transition subcomplex} of \eqn{[0,1]^n}
is a compact subsegment of the zonoseg.
The order of the generators affects this subsegment in a major way.
}

\usage{
zonoseg( mat, e0=0, ground=NULL )

\method{getsegment}{zonoseg}( x )

\method{getsegment2trans}{zonoseg}( x )

\method{print}{zonoseg}( x, ... )
}


\arguments{
\item{mat}{
a numeric matrix with 1 row whose entries determine the zonoseg.
One or more entries must be non-zero.
It is OK to have both positive and negative entries.
\cr
\code{mat} can also be a numeric vector which is then
converted to a matrix with 1 row.
}


\item{e0}{
threshold for an entry of \code{mat} to be considered 0.
Since the default is \code{e0=0},
by default an entry must be exactly 0 to become a loop in the
associated matroid.
}

\item{ground}{
The \emph{ground set} of the associated matroid of rank 1 -
an integer vector in strictly increasing order, or \code{NULL}.
\cr
When \code{ground} is \code{NULL}, it is set to \code{1:ncol(mat)}.
If \code{ground} is not \code{NULL}, \code{length(ground)} must be equal to \code{ncol(mat)}.
The point \code{ground[i]} corresponds to the \emph{i'th} column of \code{mat}.
}

\item{x}{
a \code{zonoseg} object as returned by \code{zonoseg()}
}

\item{...}{not used}
}


\details{
A \bold{zonoseg} object is a list with only 3 items:
the associated matroid, the endpoints of the segment,
and endpoints of the \emph{2-transition subsegment}.
\cr\cr
\code{print.zonoseg()} prints some information about the generators,
and the endpoints of the segment plus the 2 vertices of the unit cube
that map to these endpoints.
It prints similar data for the 2-transition subsegment.
Finally, it prints data on the associated matroid.
}


\value{
\code{zonoseg()} returns a list with S3 class \code{'zonoseg'}.
In case of error, e.g. invalid \code{mat},
the function prints an error message and returns \code{NULL}.
\cr\cr
\code{getsegment()} and \code{getsegment2trans()}
return numeric 2-vectors - the min and max endpoints of the corresponding
segments.
\cr\cr
\code{print.zonoseg()} returns \code{TRUE} or \code{FALSE}.
}



\note{
The \emph{ground set} of positive integers should not be too sparse;
otherwise performance may suffer.
}

\references{
\bold{Matroid - Wikipedia}.\cr
\url{https://en.wikipedia.org/w/index.php?title=Matroid&oldid=1086234057}
}


\seealso{
\code{\link{rank}()}
}

\examples{
zono1 = zonoseg( c(1,-2,3,0,-3,-4) )
zono1

# generators:        6 -- 3 negative, 2 positive, and 1 loops.
#
# segment:                   [-9,4]
#      value pcube.1 pcube.2 pcube.3 pcube.4 pcube.5 pcube.6
# zmin    -9       0       1       0       0       1       1
# zmax     4       1       0       1       0       0       0
#
# 2-transition subsegment:   [-8,3]
#             value source.1 source.2 source.3 source.4 source.5 source.6
# tmin-2trans    -8        1        1        0        0        1        1
# tmax-2trans     3        0        0        1        1        0        0
#
# matroid:
# ground set:           6 points   {1 2 3 4 5 6}
# hyperplanes:          1     {4}
# rank:                 1
# loops:                1   {4}
# multiple groups:      1     {1 2 3 5 6}
# uniform:              FALSE
# paving:               TRUE
# simple:               FALSE
# This matroid is constructed from a 1x6 real matrix.
#      1  2 3 4  5  6
# [1,] 1 -2 3 0 -3 -4
#
# The summary of the simplified matroid is:
#     ground set:           1 points   {1}
#                       Point 1 corresponds to the multiple group {1 2 3 5 6} in the original ...
#      hyperplanes:          1     {}
#      rank:                 1
#      loops:                0   {}
#      multiple groups:      0   {}
#      uniform:              TRUE
#      paving:               TRUE
#      simple:               TRUE
#      This matroid is constructed from a 1x1 real matrix.
#           1+...+6
#      [1,]     -13


## so the 2-transition subsegment is a proper subset of the zonoseg
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
