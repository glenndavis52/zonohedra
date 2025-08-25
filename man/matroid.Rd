\name{matroid}

\alias{matroid}
\alias{matroid.matrix}
\alias{matroid.list}



\title{
matroid construction
}

\description{
Construct a matroid from a matrix, or from explicit list of hyperplanes
}

\usage{
\method{matroid}{matrix}( x, e0=0, e1=1.e-6, e2=1.e-10, ground=NULL, ... )

\method{matroid}{list}( x, ground=NULL, ... )
}


\arguments{
\item{x}{
\code{x} can be a numeric matrix with 3, 2, or 1 rows whose columns determine the matroid.
The matrix must be either square or "wide", i.e. more columns than rows.
The matrix must be \emph{full-rank}, i.e. the rank must be equal to the number of rows, which is then the rank of the constructed matroid.
Such a matroid is often called a \emph{column matroid} or \emph{vector matroid}.
\cr\cr
\code{x} can also be a list of vectors of positive integers,
which are thought of as sets, and are the hyperplanes of the matroid.
The hyperplanes are checked that they satisfy the matroid hyperplane axioms.
The rank of the constructed matroid is determined automatically,
and must be 3, 2, or 1.
}

\item{ground}{
The \emph{ground set} of the matroid -
a vector of positive integers in strictly increasing order.
\cr\cr
When \code{x} is a matrix,
\code{length(ground)} must be equal to \code{ncol(x)}.
The point \code{ground[i]} corresponds to the \emph{i'th} column of \code{x}.
If \code{ground} is \code{NULL},
the column names of \code{x} are converted to such a vector if possible.
If this is not possible, \code{ground} is set to \code{1:ncol(x)}.
\cr\cr
When \code{x} is a list, every set in the list must be a subset of \code{ground}.
If \code{ground} is \code{NULL}, it is set to the union of all
the sets in \code{x}.
For technical reasons, when the rank is 1, \code{ground} is required
and cannot be \code{NULL}, see \bold{Details}.
}

\item{e0}{
threshold, in the \eqn{L^{\infty}} norm, for a column vector of
\code{x} to be considered 0,
and thus that the corresponding point in the matroid is a loop.
% This tolerance is always used.
Since the default is \code{e0=0}, 
by default a column vector must be exactly 0 to become a loop.
}

\item{e1}{
threshold, in a pseudo-angular sense, for column vectors to be multiples
of each other,
and thus members of a group of multiple (aka parallel) points in the matroid.
This tolerance is only used when the rank is 2 or 3.
}

\item{e2}{
threshold, in a pseudo-angular sense, for the planes
spanned by pairs of column vectors to be considered coincident,
and thus the columns to be in the same hyperplane of the matroid.
This tolerance is used when the rank is 3.
}

\item{...}{not used}
}


\details{
It was mentioned above that the tolerances \code{e1} and \code{e2} are
\emph{pseudo-angular}.
Specifically, vectors are normalized to the \eqn{L^2} unit sphere and the
distance between them is computed in the \eqn{L^{\infty}} norm.
\cr\cr
Matroids are well-known to have many cryptomorphic definitions,
e.g. independent sets, bases, circuits, rank function, closure operator,
flats, and hyperplanes.  See \bold{Matroid - Wikipedia}.
In this package, matroids can only be constructed from hyperplanes,
but there are functions
\code{\link{rank}()} and \code{\link{is_independent}()}
that can be used \emph{after} construction.
\cr\cr
Checking that the hyperplanes satisfy the matroid hyperplane axioms
is made easier by the fact that all simple matroids of rank 3 or less
are \emph{paving matroids}, see \bold{Paving Matroid - Wikipedia}
\cr\cr
Rank 1 matroids are extremely simple - the loops form the
single hyperplane (possibly empty), and the non-loops form
a multiple group.
If \code{ground=NULL} the non-loops are unknown, so this is why
\code{ground} is required when the rank is 1.
}


\value{
\code{matroid()} returns an object with S3 class \bold{'matroid'}.
\cr
In case of error, e.g. invalid \code{x} or computed hyperplanes,
the function prints an error message and returns \code{NULL}.
}



\note{
The \emph{ground set} of positive integers should not be too sparse;
otherwise performance may suffer.
\cr\cr
When \code{x} is a matrix with 3 rows,
it may happen that the computed hyperplanes do not satisfy the axioms for a matroid.
In that case, the user will be prompted to try reducing tolerance \code{e2}.
Getting the expected hyperplanes may require some \emph{a priori} knowledge of the expected hyperplanes.
For best results, the matrix should be given with maximum precision.
}

\references{
\bold{Matroid - Wikipedia}.\cr
\url{https://en.wikipedia.org/w/index.php?title=Matroid&oldid=1086234057}


\bold{Paving Matroid - Wikipedia}.\cr
\url{https://en.wikipedia.org/w/index.php?title=Paving_matroid&oldid=1021966244}
}


\seealso{
\code{\link{rank}()}, 
\code{\link{simplify}()},
\code{\link{getsimplified}()}
}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
