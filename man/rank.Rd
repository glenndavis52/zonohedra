\name{rank}

\alias{rank}
\alias{is_independent}

\title{
Rank and Independence
}

\description{
calculate the rank of any subset of a matroid, or determine whether any subset is independent
}

\usage{
rank( x, subs )

is_independent( x, subs )
}


\arguments{
\item{x}{
a matroid object, as returned from the constructor \code{\link{matroid}()}
}
\item{subs}{
a list of integer vectors, representing subsets of the ground set of \code{x}.
\code{subs} can also be an integer vector, which is put into a
list of length 1.
}
}

\value{
\code{rank(x,subs)} returns an integer vector with the same length
as the list \code{subs}.
The i'th value is the rank of the i'th set in \code{subs}.
If a set is not a subset of the ground set of \code{x}, the value is \code{NA}, and a warning message is printed.
\cr\cr
\code{is_independent(x,subs)} returns a logical vector with the same length
as the list \code{subs}.
The i'th value is the independence of i'th set in \code{x}.
It is equal to \code{TRUE} iff the rank of the subset is equal to
the cardinality of the subset.
\cr\cr
For both functions the names are copied from input to output.
}


\references{
Matroid - Wikipedia.
\url{https://en.wikipedia.org/w/index.php?title=Matroid&oldid=1086234057}

}


\seealso{
\code{\link{matroid}()}
}

\examples{
# make a matroid with rank 3
mat = matroid( classics.genlist[['RT']] )


# the ground set itself should have rank 3
rank( mat, getground(mat) )
##  [1] 3


# single points should have rank 1  (there are no loops)
rank( mat, as.list(getground(mat)) )
##  [1] 1 1 1 1 1 1

# all hyperplanes should have rank 2
rank( mat, gethyperplane(mat) )
##  [1] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2


# a point not in the ground set should have rank NA
# and the emtpy set should have rank 0
rank( mat, list(100L,integer(0)) )
##  1 of 1 subsets are not a subset of ground.
##  [1] NA 0
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.


