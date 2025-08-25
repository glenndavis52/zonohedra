\name{zonogon-section}

\alias{section.zonogon}

\title{
compute the intersection of a line and the boundary of a zonogon
}

\description{
Generically, a line intersects the boundary of a zonogon in 2 points.
Computing those 2 points is the chief goal of this function.
\cr
For a supporting line, the intersection is a face of the zonogon,
but in this function only one point of intersection is computed and returned.
}

\usage{
\method{section}{zonogon}( x, normal, beta, tol=1.e-10, plot=FALSE, ... )
}


\arguments{
\item{x}{
a \bold{zonogon} object as returned by the constructor \code{zonogon()}
}

\item{normal}{
a non-zero numeric 2-vector - the normal of all the lines
}

\item{beta}{
a numeric M-vector of line-constants.
The equation of the k'th line k is: \code{<x,normal> = beta[k]}.
}

\item{tol}{
a small positive number, used as the tolerance for the line
being considered a supporting line
}

\item{plot}{
if \code{TRUE},  the line segments formed by the the intersection
of the lines and the zonogon are \emph{added} to an existing
plot of the zonogon \code{x}, see \code{\link{plot.zonogon}()}.
The segments are drawn in dashed linestyle and the color red.
\code{boundary1} and \code{boundary2} are plotted as points.
If there is no open plot, a warning is issued.
}

\item{...}{not used}
}


\value{
\code{section.zonogon()} returns a \code{data.frame} with M rows and these columns:

\item{normal}{the given normal vector - this is the same in every row}

\item{beta}{
the given line constant
}

\item{boundary1}{
the 1st intersection point - a 2-vector
}
\item{boundary2}{
the 2nd intersection point - a 2-vector
}

Regarding orientation, if \code{normal} is considered "north" then
\code{boundary1} is on the "west" and \code{boundary2} is on the "east".

If a line is a supporting line of the zonogon,
then \code{boundary1} is some point in the boundary face (vertex or edge),
and \code{boundary2} is \code{NA}.
If a line does not intersect the zonogon,
both \code{boundary1} and \code{boundary2} are \code{NA}.


If the names of \code{beta} are unique,
they are copied to the row names of the output.
\cr\cr
In case of error, the function returns \code{NULL}.
}




\seealso{
\code{\link{zonogon}()},
\code{\link{plot.zonogon}()},
\code{\link{section.zonohedron}()}
}

\examples{
#   make a zonogon with 5 generators
pz20 = polarzonogon( 20, 5 )

section( pz20, normal=c(1,1), beta=-1:5 )
#   normal.1 normal.2 beta   boundary1.1   boundary1.2 boundary2.1 boundary2.2
# 1        1        1   -1            NA            NA          NA          NA
# 2        1        1    0 -2.220446e-16  0.000000e+00          NA          NA
# 3        1        1    1  2.452373e-01  7.547627e-01   1.0000000   0.0000000
# 4        1        1    2  6.203838e-01  1.379616e+00   1.7547627   0.2452373
# 5        1        1    3  1.095537e+00  1.904463e+00   2.3796162   0.6203838
# 6        1        1    4  1.674729e+00  2.325271e+00   2.9044629   1.0955371
# 7        1        1    5  2.420068e+00  2.579932e+00   3.3252706   1.6747294
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
