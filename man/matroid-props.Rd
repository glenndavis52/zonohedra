\name{matroid-props}

\alias{is_simple}
\alias{is_uniform}
\alias{is_paving}
\alias{is_simple.matroid}
\alias{is_uniform.matroid}
\alias{is_paving.matroid}


\title{
matroid properties
}

\description{
get some important boolean properties of a matrix,
see \bold{Matroid - Wikipedia} for the definitions.
}

\usage{
\method{is_simple}{matroid}( x )

\method{is_uniform}{matroid}( x )

\method{is_paving}{matroid}( x )
}


\arguments{
\item{x}{
a matroid object, as returned from the constructor \code{\link{matroid}()}
}
}

\value{
\code{is_simple()} returns a logical.
A matroid is \emph{simple} iff it has no loops and no multiple groups.
\cr\cr
\code{is_uniform()} returns a logical.
A matroid is \emph{uniform} iff all the hyperplanes have the same size,
which is the rank-1.
\cr\cr
\code{is_paving()} returns a logical.
For the definition of \emph{paving} see \bold{Paving Matroid - Wikipedia}.
This property is important because the hyperplane axioms
are fairly easy to check.
}




\references{
Matroid - Wikipedia.
\url{https://en.wikipedia.org/w/index.php?title=Matroid&oldid=1086234057}

Paving Matroid - Wikipedia.
\url{https://en.wikipedia.org/w/index.php?title=Paving_matroid&oldid=1021966244}
}


\seealso{
\code{\link{matroid}()}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
