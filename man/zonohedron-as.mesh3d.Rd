\name{zonohedron-as.mesh3d}

\alias{as.mesh3d}
\alias{as.mesh3d.zonohedron}

\title{
convert the boundary of a zonohedron to an \code{rgl::mesh3d} object
}

\description{
\code{mesh3d} is an S3-class object defined in package \pkg{rgl}.
Polygons in the 3D mesh can be either triangles or quadrilaterals, or both.
In the case of a zonohedron, only quadrilaterals are used.
The function \code{as.mesh3d.zonohedron()} is derived from \code{rgl::as.mesh3d()}.

After a \code{mesh3d} object is constructed,
it can be plotted by \code{rgl::dot3d()}, \code{rgl::wire3d()} or \code{rgl::shade3d()}.
\code{rgl::addNormals()} can be used to add normal vectors to the vertices of the mesh.
An object can be exported to a text file using \code{rgl::writeOBJ()}, but the vertices are not shared.
There are many more methods in \pkg{rgl}.
}

\usage{
\method{as.mesh3d}{zonohedron}( x, fcolor=NULL, falpha=1, codes=FALSE, ... )
}


\arguments{
\item{x}{
a \bold{zonohedron} object as returned by the constructor \code{zonohedron()}
}

\item{fcolor}{
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

\item{codes}{
if \code{TRUE} then an \code{NxV} raw matrix of vertex codes is added to the returned list.
All values in the matrix are 0 or 1.
\code{N} is the number of simplified generators of \code{x}, and \code{V} is the number of vertices.
These codes in the columns of the matrix
correspond to the vertices of the cube \eqn{[0,1]^N} that map to the vertices of \code{x}.
}

\item{...}{not used}
}


\details{
On the boundary of the zonohedron, a parallelogram facet with 2 generators and 4 edges
is turned into a single quadrilateral in the mesh.
A hexagon facet with 3 generators and 6 edges is divided into 2 quadrilaterals.
A zonogon facet with \eqn{m} generators is divided into \eqn{m{-}1} quadrilaterals.

Let \eqn{V} be the #(vertices), \eqn{E} be the #(edges), and \eqn{Q} be the #(quadrilaterals).
Since an edge belongs to exactly 2 quadrilaterals, it easily follows that \eqn{E = 4Q/2 = 2Q}.
The boundary of the zonohedron is homeomorphic to the 2-sphere \eqn{S^2} whose Euler characteristic is 2.
Therefore \eqn{V - E + Q = 2 \implies V = Q + 2}.

The function \code{rgl::as.mesh3d.default()} has a \code{tolerance} argument
used to determine whether two or more vertices should be merged.
We avoid that tolerance by using the fact, as pointed out by \cite{Heckbert}, that each
vertex of a zonohedron has a unique "binary code", whose length is the number of simplified generators.
In this function the binary code storage is implemented in C++ using
\code{std::map<  std::vector<unsigned char>, int >}.
}


\value{
The function returns an object of class \code{c( "mesh3d", "shape3d" )}.
It is a list with items \code{vb}, \code{ib}, \code{material}, and  \code{meshColor}.
\code{vb} is a \code{4xV} double-precision matrix, where \code{V} is the number of vertices;
the 4th row is all 1s.
\code{ib} is an \code{4xQ} integer matrix with indexes into \code{vb},
where \code{Q} is the number of quadrilaterals.

\code{material} is a list with items \code{color} and \code{alpha}.
\code{color} is a character vector of length \code{Q}, and \code{alpha} is set to \code{falpha}.
\code{meshColor} is set to \code{"faces"}.

If \code{codes=TRUE}, then a \code{NxV} raw matrix of vertex codes is added to the returned list.
See above for details on this matrix.
WARNING: Use \code{codes=TRUE} with caution, since the matrix of codes can be large.
}


\note{
If the package \pkg{rgl} is not attached, this function can be called as \code{rgl::as.mesh3d(x)}.
If it \emph{is} attached, the namespace prefix \code{rgl::} is not needed.
}

\references{
Paul Heckbert.
\bold{An Efficient Algorithm for Generating Zonohedra}.
3-D Technical Memo 11.
24 February 1985.
Computer Graphics Lab.
New York Institute of Technology
}

\seealso{
\code{\link{zonohedron}()},
\code{\link[rgl:as.mesh3d]{rgl::as.mesh3d}()},
\code{\link[rgl:shade3d]{rgl::shade3d}()},
\code{\link[rgl:writeOBJ]{rgl::writeOBJ}()}
}

\examples{
#  view a zonohedron using rgl::shade3d()
zono12 = zonohedron( classics.genlist[[12]] )
mesh   = rgl::as.mesh3d( zono12 )
rgl::shade3d( mesh, back='culled' )


cube = zonohedron( classics.genlist[[ 'C' ]] )
mesh = rgl::as.mesh3d( cube )
class( mesh ) = "list"   # we bypass print.mesh3d() and see the contents of the mesh
mesh

##  $vb
##       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
##  [1,]    0    0    1    1    1    1    0    0
##  [2,]    0    1    1    0    1    0    0    1
##  [3,]    0    0    0    0    1    1    1    1
##  [4,]    1    1    1    1    1    1    1    1
##  
##  $ib
##       [,1] [,2] [,3] [,4] [,5] [,6]
##  [1,]    1    4    1    8    5    6
##  [2,]    2    6    7    7    3    4
##  [3,]    3    7    8    6    2    3
##  [4,]    4    1    2    5    8    5
##  
##  $material
##  $material$color
##  [1] "blue" "blue" "blue" "blue" "blue" "blue"
##  
##  $material$alpha
##  [1] 1
##  
##  $meshColor
##  [1] "faces"
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
