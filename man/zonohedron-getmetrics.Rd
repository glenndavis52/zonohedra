% File validate.Rd

\name{getmetrics.zonohedron}

\alias{getmetrics}
\alias{getmetrics.zonohedron}
\alias{print.zonohedron}
\alias{summary.zonohedron}

\title{Get Important Metrics about a Zonohedron, and Print and Summarize Them}

\description{
Get some important zonohedron metrics; for most some computation is needed.
\cr\cr
The \code{print()} function prints nicely formatted
facts about a zonohedron, including its matroid.
\cr\cr
The \code{summary()} function prints a single-line summary,
formatted as a row in a data frame.
}

\usage{
\S3method{getmetrics}{zonohedron}( x )

\S3method{print}{zonohedron}( x, trans2=FALSE, matroid=TRUE, ... )

\S3method{summary}{zonohedron}( object, ... )
}

\arguments{
\item{x}{a \code{zonohedron} object as returned by the constructor \code{\link{zonohedron}()}}

\item{trans2}{if \code{TRUE} then print extra metrics on the 2-transition surface associated with \code{x}}

\item{matroid}{if \code{TRUE} then print extra information about the matroid associated with \code{x}}

\item{object}{a \code{zonohedron} object as returned by the constructor \code{\link{zonohedron}()}}

\item{...}{for \code{print()} further arguments are ignored;
for \code{summary()} the further arguments can be *more* zonohedron objects, which are summarized by adding more rows to the same data frame; see \bold{Examples}.
}
}


\value{
\code{getmetrics.zonohedron()} returns a list with these items:

\item{vertices}{the number of vertices}
\item{edges}{the number of edges}
\item{facets}{the number of facets (2D faces); all of them are zonogons}
\item{area}{the sum of the areas of all the facets}
\item{volume}{as a polytope}

All of these are always positive.
\cr\cr
\code{print.zonohedron()} returns \code{TRUE} or \code{FALSE}.
\cr\cr
\code{summary.zonohedron()} returns a data frame, see \bold{Examples}.
}    

\seealso{
\code{\link{genlist}},
\code{\link{zonohedron}()},
}

\examples{
zono = zonohedron( classics.genlist[['BD']] )
zono
# zonohedron:
# fullname:                          Bilinski dodecahedron 
# generators (original):             4 
# generators with multiples:         0 
# generators (simplified):           4 
# number of facets:                  12  [6 antipodal facet-pairs]
# facets that contain 0:             4    { 1 3 4 6 }
# number of edges:                   24 
# center:                            0.809017 2.118034 1.309017 
# pointed:                           TRUE 
# salient:                           TRUE 
# area:                              38.83282 
# volume:                            16.94427 
# 
# matroid:
# ground set:           4 points   {1 2 3 4}
# hyperplanes:          6     {1 2}  {1 3}  {1 4}  {2 3}  {2 4}  {3 4}
# rank:                 3
# loops:                0   {}
# multiple groups:      0   {}
# uniform:              TRUE
# paving:               TRUE
# simple:               TRUE
# This matroid is constructed from a 3x4 real matrix.
#             1        2        3         4
# [1,] 1.000000 1.618034 0.000000 -1.000000
# [2,] 1.618034 0.000000 1.000000  1.618034
# [3,] 0.000000 1.000000 1.618034  0.000000


summary( zono )
#                fullname generators vertices edges facets     area   volume
# 1 Bilinski dodecahedron          4       14    24     12 38.83282 16.94427

zono4 = zonohedron( classics.genlist[['RI']] )
zono7 = zonohedron( classics.genlist[['TO']] )
summary( zono, zono4, zono7 )
#                fullname generators vertices edges facets     area   volume
# 1 Bilinski dodecahedron          4       14    24     12 38.83282 16.94427
# 2   rhombic icosahedron          5       22    40     20 64.72136 42.36068
# 3  truncated octahedron          6       24    36     14 53.56922 32.00000
}

