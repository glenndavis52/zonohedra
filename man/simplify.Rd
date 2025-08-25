\name{simplify}

\alias{simplify}
\alias{simplify.list}
\alias{unsimplify}
\alias{unsimplify.list}
\alias{getsimplified}
\alias{getsimplified.matroid}

\title{
simplify and unsimplify
}

\description{
A \emph{simple matroid} has no loops and no multiple groups.
Simplification is the process of removing all loops,
and every point except one from each multiple group.
The result is a simple matroid.
\cr
The functions below simplify a matroid, or an explicit list of hyperplanes.
\cr
The hyperplanes can be \emph{unsimplified} if the original loops
and multiple groups are known.
}

\usage{
\method{getsimplified}{matroid}( x, ... )

\method{simplify}{list}( x, ground=NULL, ... )

\method{unsimplify}{list}( x, loop=NULL, multiple=NULL, ground=NULL, ... )
}


\arguments{
\item{x}{
\code{x} can be a matroid object, as returned from the constructor \code{\link{matroid}()}
\cr
\code{x} can also be a list of vectors of positive integers,
which are thought of as sets.
All must be subsets of \code{ground}.
They do not have to satisfy the matroid hyperplane axioms.
For a definition of \emph{loop} and \emph{multiple group} in this case,
see \bold{Details}.
}

\item{ground}{
The \emph{ground set} of the sets in \code{x}
(both \code{simplify()} and \code{unsimplify()}).
It must be a vector of positive integers in strictly increasing order
(not verified).
If \code{ground} is \code{NULL}, it is set to the union of
the sets in \code{x}.
}

\item{loop}{
a vector of positive integers, the loops, to add to the list \code{x};
\code{loop} must be disjoint from \code{ground} (verified).
If \code{loop} is \code{NULL}, the function looks for \code{loop}
in the attribute data \code{attr(x,'lmdata')}.
If there is no such attribute, \code{loop} is set to the empty set.
}

\item{multiple}{
a list of vectors of positive integers, the multiple groups,
to add to the list \code{x};
these groups must be pairwise disjoint and disjoint from \code{loop} (not verified).
Each group must intersect the ground set in exactly one point (verified).
If \code{multiple} is \code{NULL}, the function looks for \code{multiple}
in the attribute data \code{attr(x,'lmdata')}.
If there is no such attribute, \code{multiple} is set to the empty list.
}

\item{...}{not used}
}


\details{
First consider the case when \code{x} is a list of vectors of positive integers.
Each vector represents a subset of the ground set.
They are not required to satisfy the hyperplane axioms,
but by abuse of language we will call them hyperplanes in this paragraph.
A \emph{loop} is a point (an integer in the ground set)
that is in every hyperplane.
Imagine now that all loops have been removed.
Say that two points \eqn{p} and \eqn{q} are \emph{multiples} iff
for every hyperplane \eqn{H}, \eqn{p \in H} iff \eqn{q \in H}.
This is an equivalence relation, and the \emph{multiple groups} are the
equivalence classes with more than one point.
For computation it is convenient to think of a boolean \emph{incidence matrix}.
There is a column for each point in the ground set,
and a row for each hyperplane.
An entry is \code{TRUE} iff the point is in the hyperplane.
A \emph{loop} is then a column of all \code{TRUE}s.
A \emph{multiple group} is a maximal set of duplicate columns.
This is basically how \code{simplify()} is implemented,
except with optimizations that avoid computing the very large incidence matrix.
\cr\cr
Now consider the case when \code{x} is a matroid object.
When \code{x} was constructed, the simplification of \code{x} was
computed (with help from the *previous* \code{simplify()}) and stored
as a member of \code{x} (unless \code{x} was already simple).
So in this case \code{getsimplified(x)} does not do any real work
and only takes microseconds.
\cr\cr
These functions are accelerated with C/C++.
}


\value{
If \code{x} is a matroid, \code{getsimplified(x)} returns \code{x}
when \code{x} is simple, and a member of \code{x} when \code{x} is not simple.
It does not do any real work.
\cr\cr
If \code{x} is a list, \code{simplify(x)} returns a list
of the same length, but with all loops removed,
and every point except one from each multiple group removed.
The integer that remains is the smallest one in the group.
The order of the sets is preserved.
It also sets the \code{'lmdata'} attribute of the returned list
to a list of 2 objects -
the loop and multiple group data found in \code{x}.
\cr\cr
If \code{x} is a list, \code{unsimplify(x)} returns a list
of the same length, but with the loops and multiples added back.
The order of the sets is preserved.
\cr\cr
In case of error, e.g. invalid \code{x} etc.,
the function prints an error message and returns \code{NULL}.
}



\references{
Matroid - Wikipedia.
\url{https://en.wikipedia.org/w/index.php?title=Matroid&oldid=1086234057}
}


\seealso{
\code{\link{rank}()},
\code{\link{matroid}()}
}

\examples{
# an example using simplify.list() and unsimplify.list()
# get the matrix for CIE XYZ at 5 nm step size
mat3x81 = colorimetry.genlist[[1]]

# create the matroid
mat5 = matroid( mat3x81 )

#  test for simplicity
is_simple(mat5)
##  [1] FALSE

# get the list of hyperplanes, and simplify
hyper = gethyperplane( mat5 )
hypersimple = simplify( hyper )

# print the loop and multiple data found
attr(hypersimple,'lmdata')

# unsimplify and compare to the originals
# the list attr(hypersimple,'lmdata') is 'secretly' used in unsimplify()
identical( unsimplify(hypersimple), hyper )
##  [1] TRUE
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
