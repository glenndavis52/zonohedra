% File zonohedra-package.Rd
\docType{package}
\name{zonohedra-package}
\alias{zonohedra-package}

\title{zonhedra package}



\description{
This package deals with \emph{zonohedra}, which are zonotopes of dimension 3.
It also handles \emph{zonogons} (2D zonotopes)
and \emph{zonosegs} (1D zonotopes).

The term \emph{zonoseg} ("zonotope" + "segment") is my own personal term;
I could not find an alternative term.
It is a linear image of the unit cube \eqn{[0,1]^n} in the real numbers,
and a compact segment of reals.
}


\section{S3 classes}{
\tabular{cc}{
\code{Z}   \tab \code{class(Z)}  \cr
\bold{zonohedron} \tab \code{"zonohedron" "zonotope" "list"}  \cr
\bold{zonogon}    \tab \code{"zonogon" "zonotope" "list"}  \cr
\bold{zonoseg}    \tab \code{"zonoseg" "zonotope" "list"}  \cr
}

For example, the \code{section()} returns very diffferent things
for a zonohedron and a zonogon, and so
\code{section.zonohedron()} and \code{section.zonogon()}
are coded and documented separately.
A section for a zonoseg does not make sense,
so \code{section.zonoseg()} is undefined.
}


\section{Terminology}{
For a convex polytope, a \emph{supporting hyperplane} is a hyperplane
that intersect the polytope's boundary but \emph{not} its interior.
\cr
A zonotope is a convex polytope.
A zonohedron has supporting planes,
and a zonogon has supporting lines.
\cr\cr
In the package \pkg{zonohedra},
a \emph{zonotope} mean a zonotope of dimension 3, 2, or 1.
\cr\cr
A \emph{face} of a zonotope is the intersection of the boundary
of the zonotope with some supporting hyperplane.
A \emph{d-face} is a face of dimension \emph{d}.
So a \emph{0-face} is a \emph{vertex},
and a \emph{1-face} is an \emph{edge}.
\cr\cr
A \emph{facet} of a zonotope is a face whose dimension is
1 less than the dimension of the zonotope.
A facet is a maximal proper face.
\cr\cr
A zonohedron has 0-faces (vertices), 1-faces (edges),
and 2-faces (facets).
\cr\cr
A zonogon has 0-faces (vertices) and 1-faces (edges).
Since the dimension of an edge is 1 less than the
dimension of the zonogon, an edge of a zonogon is also a facet
of a zonogon.
}



\keyword{package}

