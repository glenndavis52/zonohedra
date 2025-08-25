\name{zonogon}

\alias{zonogon}
\alias{polarzonogon}


\title{
zonogon construction
}

\description{
Construct a zonogon from a numeric matrix with 2 rows.
}

\usage{
zonogon( mat, e0=0, e1=1.e-6, ground=NULL )

polarzonogon( n, m=n, ground=NULL )

}



\arguments{
\item{mat}{
a numeric 2xM matrix, where 2 \eqn{\le} M.
The matrix must have rank 2 (verified).
The M columns are the generators of the zonogon.
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

\item{ground}{
The \emph{ground set} of the associated matroid of rank 2 -
an integer vector in strictly increasing order, or \code{NULL}.
\cr
When \code{ground} is \code{NULL}, it is set to \code{1:ncol(mat)}.
If \code{ground} is not \code{NULL}, \code{length(ground)} must be equal to \code{ncol(mat)}.
The point \code{ground[i]} corresponds to the \emph{i'th} column of \code{mat}.
}

\item{n}{
an integer \eqn{\ge} 3.
The generators are computed as \code{n} equally spaced points on
the unit circle, starting at (1,0).
}

\item{m}{an integer with 2 \eqn{\le} \code{m} \eqn{\le} \code{n}.
When  \code{m} < \code{n}, only the first \code{m}
points are used as generators of the zonogon.
}
}


\details{
\code{polarzonogon()} is useful for testing.
The term \emph{polar zonogon} is my own, and based on
the \emph{polar zonohedron} in \emph{Chilton & Coxeter}.
It it loads the matrix \code{mat} and passes it to \code{zonogon()}.
When \code{m=n} the zonogon is a regular 2n-gon.
When \code{m<n} the zonogon is a has 2\code{m} vertices,
but is not necessarily regular.
The generators correspond to the n'th-roots of unity.
}


\value{
\code{zonogon()} and \code{polarzonogon()} return a list with S3 class \code{'zonogon'}.
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
}


\seealso{
\code{\link{zonohedron}()},
\code{\link{zonoseg}()},
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
