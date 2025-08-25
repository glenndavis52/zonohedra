\name{grpDuplicated}
\alias{grpDuplicated}
\alias{grpDuplicated.default}
\alias{grpDuplicated.matrix}
\title{
Grouping by duplicated elements
}
\description{
\code{grpDuplicated()} is a generic function that takes an indexed set
of "elements", and outputs an integer vector with the same length.
The "elements" can be components of a vector,
or the row vectors or column vectors of a matrix.
In the output vector, a component is 0 if and only if the corresponding
element is unique.
When the element is unique, it forms a \emph{singleton group}.
Output components have equal positive integer values
if and only if the corresponding elements are identical to each other.
These elements form a \emph{non-singleton group},
and the positive integer is called the \emph{group number}.


The number of singleton groups is equal to #(zeros),
which is equal to the #(elements) - #(duplicated elements).\cr
The number of non-singleton groups is equal to max(output vector).\cr 
The number of all groups is equal to #(zeros) + max(output vector).
}

\usage{
\method{grpDuplicated}{default}( x, \dots )
	
\method{grpDuplicated}{matrix}( x, MARGIN=1, \dots )
}
\arguments{
  \item{x}{a vector or matrix of atomic mode \code{"numeric"}, \code{"integer"}, \code{"logical"}, \code{"complex"}, \code{"character"} or \code{"raw"}.  
  }
  
  \item{MARGIN}{an integer scalar, the matrix margin to be held fixed, as in \code{\link{apply}}. \code{MARGIN=1} means that it looks for duplicated rows, and \code{MARGIN=2} means that it looks for duplicated columns.  Other values are invalid.}  

  \item{\dots}{arguments for particular methods.}
}

\details{
The implementation is based on \code{std::unordered_map} in C++11,
which uses a hash-table.
}

\value{
The return value is an integer vector with all elements ranging from 0 to \code{K}, where \code{K} is the number of non-singleton groups.\cr
For vector \code{x} the elements are the vector components,
and the output is the same length as the input.\cr
For a matrix \code{x} with \code{MARGIN=1}, the elements are the rows of
the matrix and the output has length \code{nrow(x)}.\cr
For a matrix \code{x} with \code{MARGIN=2}, the elements are the columns of
the matrix and the output has length \code{ncol(x)}.\cr
The \code{'ngroups'} attribute of the returned vector is set to an integer 3-vector.
The 1st component is the total number of groups,
the 2nd component is the number of singleton groups,
and the 3rd component is the number of non-singleton groups \code{K}.
}

\author{
Long Qu and Glenn Davis
}

\note{
The templated C++ function that does the real work is taken from the package \bold{uniqueAtomMat}
by Long Qu,
but the returned vector is slightly modified by Glenn Davis.\cr
}

\source{
\url{https://github.com/cran/uniqueAtomMat/}\cr
The package \bold{uniqueAtomMat} was removed from CRAN by its author Long Qu.
}


\examples{
set.seed(0)

#   test a numeric vector
x = rnorm(7)
y = rnorm(5)
grpDuplicated( c(x,y,rev(x)) )
##  [1] 7 6 5 4 3 2 1 0 0 0 0 0 1 2 3 4 5 6 7
##  attr(,"ngroups")
##  [1] 12  5  7

# test a numeric matrix, both rows and columns
A = matrix( rnorm(3*7), 3, 7 )
B = matrix( rnorm(3*5), 3, 5 )

#   the columns of cbind(A,B,A) have the duplicates one would expect
grpDuplicated( cbind(A,B,A), MARGIN=2 )
##  [1] 1 2 3 4 5 6 7 0 0 0 0 0 1 2 3 4 5 6 7
##  attr(,"ngroups")
##  [1] 12  5  7

# but the rows of cbind(A,B,A) are unique
grpDuplicated( cbind(A,B,A), MARGIN=1 )
##  [1] 0 0 0
##  attr(,"ngroups")
##  [1] 3 3 0
}


