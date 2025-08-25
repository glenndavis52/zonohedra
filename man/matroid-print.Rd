% File validate.Rd

\name{print}

\alias{print}
\alias{print.matroid}

\title{Print Basic Facts about a Matroid}

\description{
The function prints a nicely formatted summary of a matroid,
including the ground set, the rank, loops, multiple groups,
and some boolean properties.
It prints the number of hyperplanes, broken down by their size.
If it is a vector matroid, and its matrix is not too large, it prints
that matrix.
If the matroid is not simple, it also prints the simplified matroid.
}

\usage{
\S3method{print}{matroid}( x, ... )
}

\arguments{
\item{x}{a \code{matroid} object as returned by the constructor \code{\link{matroid}()}}

\item{...}{further arguments ignored, but required by the generic \code{print()}}
}


\value{
The function returns \code{TRUE} or \code{FALSE}.
}    

\seealso{
\code{\link{matroid}()}
}

