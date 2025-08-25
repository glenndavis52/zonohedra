\name{zonohedron}

\alias{zonohedron}
\alias{polarzonohedron}
\alias{regularprism}


\title{
zonohedron construction
}

\description{
Construct a zonohedron from a numeric matrix with 3 rows.
Also construct some special zonohedra useful for testing.
}

\usage{
zonohedron( mat, e0=0, e1=1.e-6, e2=1.e-10, ground=NULL )

polarzonohedron( n, m=n, height=pi, ground=NULL )

regularprism( n, m=n, axis=c(0,0,1), ground=NULL )

}


\arguments{
\item{mat}{
a numeric 3xM matrix, where 3 \eqn{\le} M.
The matrix must have rank 3 (verified).
The M columns are the generators of the zonohedron.
}

\item{e0}{
threshold for a column of \code{mat} to be considered 0,
in the \eqn{L^\infty} norm.
Since the default is \code{e0=0},
by default a column must be exactly 0 to be considered 0.
}

\item{e1}{
threshold, in a pseudo-angular sense,
for non-zero column vectors to be multiples of each other,
and thus members of a group of multiple (aka parallel) points in the associated matroid.
It OK for a column to be a negative multiple of another.
}

\item{e2}{
threshold, in a pseudo-angular sense, for the planes
spanned by pairs of column vectors to be considered coincident,
and thus the columns to be in the same hyperplane of the associated matroid.
}

\item{ground}{
The \emph{ground set} of the associated matroid of rank 3 -
an integer vector in strictly increasing order, or \code{NULL}.
\cr
When \code{ground} is \code{NULL}, it is set to \code{1:ncol(mat)}.
If \code{ground} is not \code{NULL}, \code{length(ground)} must be equal to \code{ncol(mat)}.
The point \code{ground[i]} corresponds to the \emph{i'th} column of \code{mat}.
}

\item{n}{
an integer \eqn{\ge} 3.
The generators are computed as \code{n} equally spaced points on
a circle.
See \bold{Details} for more on this circle.
}

\item{m}{an integer with 2 \eqn{\le} \code{m} \eqn{\le} \code{n}.
When  \code{m }< \code{n}, only the first \code{m}
points are used as generators of the zonohedron.
}

\item{height}{the z value at the apex of the zonohedron,
which is the sum of all the generators.
The z value of all the generators is set to make this happen.
When \code{height=pi}, as \eqn{n \to \infty} the zonohedron
converges to the interior of the surface of revolution
of the curve \eqn{x = sin(z)} for \eqn{z \in [0,\pi]},
see \bold{Chilton & Coxeter}.
}

\item{axis}{
the axis of the regular prism. It must be a 3-vector with z value non-zero.}
}


\details{
In \code{zonohedron()}, the contruction of the \emph{zones} (or \emph{belts}) is optimized by following
the procedure in \cite{Heckbert}.
The key step is sorting face normals that all lie on a great circle
of the unit sphere.
\cr\cr
For \code{polarzonohedron()} the circle is centered at
(0,0,\code{height/n}) and parallel to the xy-plane.
The radius is \code{height/n}.
\cr\cr
For \code{regularprism()} the circle is the unit circle in the xy-plane.
The 3-vector \code{axis} is added as column \code{m+1} of the matrix.
The returned zonohedron is the Minkowski sum of a zonogon and
the line segment defined by \code{axis}.
If \code{m} < \code{n}, the zonogon may not be regular.
\cr\cr
Both of these functions are useful for testing.
They load the matrix \code{mat} and pass it to \code{zonohedron()}.
}

\value{
\code{zonohedron()} and \code{polarzonohedron()} return a list with S3 class \code{'zonohedron'}.
In case of error, e.g. invalid \code{mat},
the functions print an error message and returns \code{NULL}.
}

\note{
The \emph{ground set} of positive integers should not be too sparse;
otherwise performance may suffer.
}

\references{
B. L. Chilton and H. S. M. Coxeter.
\bold{Polar Zonohedra}.
The American Mathematical Monthly.
Vol 70. No. 9.
pp. 946-951.
1963.

Paul Heckbert.
\bold{An Efficient Algorithm for Generating Zonohedra}.
3-D Technical Memo 11.
24 February 1985.
Computer Graphics Lab.
New York Institute of Technology
}


\seealso{
\code{\link{zonohedron}()},
\code{\link{zonoseg}()},
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
