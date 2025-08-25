\name{zonogon-getmetrics}

\alias{getmetrics.zonogon}
\alias{print.zonogon}

\title{
get important metrics about a zonogon,
and print basic facts about a zonogon.
}

\description{
Get some important zonogon metrics; some computation is used.
Also print the data, and more.
}

\usage{
\method{getmetrics}{zonogon}( x )

\method{print}{zonogon}( x, ... )
}


\arguments{
\item{x}{
a \bold{zonogon} object
}
\item{...}{not used}
}

\value{
\code{getmetrics()} returns a list with these items:

\item{vertices}{the number of vertices}
\item{perimeter}{the sum of the lengths of the all edges}
\item{area}{as a polygon}

All of these are always positive.

\code{print.zonogon()} returns \code{TRUE} or \code{FALSE}.
}

\details{
\code{print.zonogon()} prints some basic information about the zonogon,
and the associated matroid.
}

\seealso{
\code{\link{zonogon}()},
\code{\link{getmetrics.zonohedron}()}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
