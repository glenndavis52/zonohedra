\name{spherize}

\alias{spherize}
\alias{spherize.zonotope}

\title{
spherize a zonotope
}

\description{
The input is a zonotope with a best-fit ellipsoid
(or ellipse for a \bold{zonogon})
with axes that may have very different lengths.
The function computes a \emph{spherizing matrix} \code{W}, and then
transforms the zonotope so its boundary is close to a sphere.
}

\usage{
\method{spherize}{zonotope}( x, method='ZCA', ... )
}


\arguments{
\item{x}{
a zonotope object -  a \bold{zonohedron}, a \bold{zonogon},
or a \bold{zonoseg}.
}

\item{method}{
for computing the matrix \code{W}, either \code{'ZCA'} or \code{'PCA-COR'}.
Matching is partial and case-insensitive.
}

\item{...}{not used}
}



\details{
The 2 methods are taken from \bold{Kessy, et. al.}.
}


\value{
After computing the matrix \code{W},
the function return \code{lintransform(x,W)}.
The \code{"sphering"} attribute is set to \code{W}.
\cr
If \code{x} is a 1D \bold{zonoseg}, sphering is not really possible,
so the function prints a warning message and returns \code{x}.
In case of error, the function returns \code{NULL}.
}



\references{
Agnan Kessy, Alex Lewin, Korbinian Strimmer.
\bold{Optimal whitening and decorrelation}.\cr
\url{https://arxiv.org/abs/1512.00809}
v4 2016.
}

\seealso{
\code{\link{lintransform}()}
}

