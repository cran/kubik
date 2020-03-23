\name{15_solve_method}
\alias{solve.chs}
\title{Solve Method}
\description{Solve method, for chs objects.}
\usage{
\method{solve}{chs}(a, b, \dots, to.list=FALSE)
}
\arguments{
\item{a}{A chs object.}
\item{b}{A vector of y values.}
\item{to.list}{If true, return list.\cr Ignored, if b has length two or greater.}
\item{\dots}{Other arguments, for the roots.chs.eval function.}
}
\value{
If to.list is false (the default), and length (b) is one:\cr
It returns a single numeric vector.

Otherwise:\cr
It returns a list of numeric vectors, one for each value in b.
}
\details{
This function is a wrapper for roots.chs.eval.\cr
(It calls roots.chs.eval with modified cy values).

It computes x values, where the spline's y value is equal to the values in b.
}
\references{Please refer to the help page for \link{chs} for background information and references.}
\seealso{\link{chs}, \link{roots.chs}}
\examples{
#control points
cx <- 1:4
cy <- c (-4, -1, 1, 4)

#cubic hermite spline
#(with function object, and default slopes)
f <- chs (cx, cy)

#solve for x, given y=2
solve (f, 2)
}
