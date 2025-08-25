\name{transitionsdf}

\alias{transitionsdf}

\title{
summarize the number of transitions and associated data,
over all parallelograms in the boundary of a zonohedron
}

\description{
The 2-transition surface is a union of parallelograms.
For this function, the surface is required to be strictly starshaped at
the center.
For the definition of \emph{strictly starshaped} see
\href{../doc/transitions.html}{The 2-Transition Subcomplex and the 2-Transition Surface}.

Each parallelogram has a unit normal that defines a linear functional.
\cr
If the 2-transition parallelogram is in the \emph{interior} of the zonohedron
then the functional is not maximized on the parallelogram,
and there is a corresponding similar parallelogram
on the boundary of the zonohedron where the functional *is* maximized.
The first parallelogram (in the surface) is called \emph{deficient} because
the functional is not maximized,
and the second parallelogram (in the boundary)
is called \emph{abundant} because the number of corresponding
transitions across this parallelogram is more than 2.
The difference between the functional values is called the \emph{deficit}.
\cr
If the 2-transition parallelogram is on the boundary,
then it is called \emph{coincident}.
It is also called \emph{non-deficient} and the deficit is 0.
}

\usage{
transitionsdf( x, trans2=TRUE )
}

\arguments{
\item{x}{
a \bold{zonohedron} object as returned by the constructor \code{zonohedron()}.
The 2-transition surface must be strictly starshaped.
}

\item{trans2}{
if \code{TRUE}, then include metrics on the non-deficient (coincident)
parallelograms, with 2 transitions.
This is always the first row of the returned data frame.
\cr
if \code{FALSE}, then data on the non-deficient parallelograms
is not included, and the returned data frame only has data
on the deficient parallelograms, with more than 2 transitions.
}

}



\value{
\code{transitionsdf()} returns a \code{data.frame} with a row
for each number of transitions found,
plus a final row with totals on appropriate columns.
The columns are:

\item{transitions}{
the number of transitions, a positive even integer, in increasing order.
}

\item{parallelograms}{
the number of parallelograms with the given number of transitions
}

\item{area}{
the min and max of the area of the
parallelograms with the given number of transitions
}

\item{area.sum}{
the total area of the parallelograms with the given number of transitions
}

\item{deficit}{
the min and max of the deficit of the
parallelograms with the given number of transitions.
When there are 2 transitions the deficit should be exactly 0,
but is usually slightly non-0 due to truncation.
When there are more than 2 transitions the deficit is positive.
}

\item{example}{
the 2 generators (from the ground set of the simplified matroid)
of the parallelogram with the maximum area
}

In case of error, the function returns \code{NULL}.
}

\note{
Because of the 1-1 correspondence between similar parallelograms,
the surface areas of the 2-transition surface
and the boundary of the zonohedron are equal.
}


\seealso{
\code{\link{zonohedron}()},
\code{\link{plot.zonohedron}()}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
