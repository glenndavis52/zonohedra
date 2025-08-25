\name{symmetrize}

\alias{symmetrize}
\alias{symmetrize.zonotope}

\title{
symmetrize a zonotope
}

\description{
The input is a zonotope whose matroid is simple.
The function adds new generators that creates a new zonotope
that is a translate of the original, and has center of symmetry at 0.
}

\usage{
\method{symmetrize}{zonotope}( x, e0=0, e1=1.e-6, e2=1.e-10, ... )
}


\arguments{
\item{x}{
a zonotope object -  a \bold{zonohedron}, a \bold{zonogon},
or a \bold{zonoseg}.
The matroid of this zonotope must be simple.
}

\item{e0}{see \code{\link{zonohedron}()}}

\item{e1}{see \code{\link{zonohedron}()}}

\item{e2}{see \code{\link{zonohedron}()}}

\item{...}{not used}
}



\details{
Each generator \code{g} (a column of the matrix)
is replace by 2 generators: \code{g/2} and \code{-g/2}.
The new set of generators correponds to a \emph{star} at 0,
from Sec 2-8 of \bold{Coxeter}.
\cr
The new ground points are obtained by translating the original
ground points by the their maximum.
}


\value{
The function returns a zonotope that is a translate of the original,
and has center of symmetry at 0.
\cr
In case of error, the function returns \code{NULL}.
}



\references{
Coxeter, H.S.M.
\bold{Regular Polytopes}.
Dover Publications.
1973.
}

\seealso{
\code{\link{zonohedron}()},
\code{\link{zonogon}()},
\code{\link{zonoseg}()}
}


