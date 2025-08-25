\name{zonogon-invert}

\alias{invert}
\alias{invert.zonogon}

\title{
invert points in a zonogon
}

\description{
A \emph{zonogon} \eqn{Z} is the image of a linear map
\eqn{[0,1]^n \to Z \subset \bold{R}^2}, from the n-cube to the plane.
For a point in a zonogon, this function
finds a point in the unit cube that maps to it.
There are infinitely many such points in general (unless \eqn{n=2}),
but this function picks a specific point using the \emph{standard tiling},
see \bold{Details}.
}

\usage{
\method{invert}{zonogon}( x, z, tol=0, plot=FALSE, ... )
}


\arguments{
\item{x}{
a \bold{zonogon} object as returned by the constructor \code{zonogon()}
}

\item{z}{
a numeric Mx2 matrix, with M points in the rows.
\code{z} can also be a numeric vector that can be converted to such
a matrix, by row.
}

\item{tol}{
points that are within \code{tol} of the boundary are
still processed, see \bold{Details}
}

\item{plot}{
if \code{TRUE} then the points in \code{z} are \emph{added} to an existing
plot of the zonogon \code{x}, using a red X symbol,
see \code{\link{plot.zonogon}()}.
If there is no plot open, a warning is issued.
}

\item{...}{not used}
}


\details{
The given points are first tested for being inside the zonogon,
using \code{\link{inside}()} and the given \code{tol}.
If any are outside, a warning is issued.
When the corresponding point \code{pcube} is computed,
it is clamped to the unit cube,
so the inversion error may be as large as \code{tol}.
\cr\cr
Inversion is not unique in general.
For this function, the \emph{standard tiling} of the zonogon by parallelograms
is computed; it is an example of a \emph{zonotopal tiling}.
It is a \emph{regular zonotopal tiling} because it arises from the
projection of a zonohedron onto the plane, see \cite{Ziegler}.
The function \code{\link{plot.zonogon}()} has an option to plot this tiling.
Given the point \code{z}, the function determines a parallelogram
that contains the point.
The \code{pcube} coordinates of the \emph{base} of this parallelogram are
all 0 or 1, and the coordinates of \code{z}
\emph{within} the parallelogram are in [0,1].
Thus, all coordinates of \code{pcube} are 0 or 1,
except possibly for 2 of them.

}


\value{
\code{invert.zonogon()} returns a \code{data.frame} with M rows and these columns:

\item{z}{the given point}
\item{pcube}{
a point in the unit cube that maps to \code{z}.
Every \code{pcube} has all coordinates 0 or 1,
except possibly for the 2 given by \code{hyper}, see \bold{Details}.
}
\item{hyper}{
the 2 indexes of the generators of the parallelogram
that contains \code{z}, in the simplified matroid.
These 2 coordinates in \code{pcube} are not 0 or 1 in general.
}
\item{hyperidx}{the index of the parallelogram that contains \code{z}}

If a point \code{z} cannot be inverted, the other columns are
all \code{NA}, and a warning message is printed.

If the row names of \code{z} are unique,
they are copied to the row names of the output.
The column names of \code{pcube} are copied from the ground set of the
associated matroid.
\cr\cr
In case of error, the function returns \code{NULL}.
}



\references{
Ziegler, G.M.
\bold{Lectures on Polytopes}.
Graduate Texts in Mathematics.
Springer New York.
2007.
}

\seealso{
\code{\link{zonogon}()},
\code{\link{inside}()},
\code{\link{plot.zonogon}()}
}

\examples{
#   make a zonogon with 5 generators
pz20 = polarzonogon( 20, 5 )

#   make 7 random points in the zonogon
set.seed(0)
pcube = matrix( runif(5*7), 5, 7 )
z = t( getmatrix(pz20) \%*\% pcube )

#  invert these 7 points back to the cube
invert( pz20, z )
#         z.1       z.2    pcube.1    pcube.2    pcube.3    pcube.4    pcube.5
# 1 2.0676319 1.6279807 0.00000000 0.70030526 1.00000000 1.00000000 0.01553241
# 2 2.4031738 1.9658035 0.00000000 0.96572153 1.00000000 1.00000000 0.28450140
# 3 0.9230336 1.0885446 0.00000000 0.00000000 0.39548689 1.00000000 0.04948838
# 4 2.5242122 1.7395069 0.16540765 1.00000000 1.00000000 1.00000000 0.03542132
# 5 2.2598725 1.0601592 0.38111324 1.00000000 1.00000000 0.20192029 0.00000000
# 6 1.1387813 1.2636700 0.00000000 0.00000000 0.65250505 1.00000000 0.07478012
# 7 1.6315341 1.0777737 0.00000000 0.64210923 1.00000000 0.36039509 0.00000000

#   hyper.1 hyper.2 hyperidx
# 1       2       5        7
# 2       2       5        7
# 3       3       5        9
# 4       1       5        4
# 5       1       4        3
# 6       3       5        9
# 7       2       4        6
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
