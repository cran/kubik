\name{rd_5_solve_method}
\alias{solve.chs}
\title{Solve Method}
\description{Solve method, for chs objects.}
\usage{
\method{solve}{chs}(a, b, to.list, \dots)
}
\details{
Refer to chs.roots(), for more information.

The chs solve method computes x values, where the spline is equal to
the y values in b.

In principle, it just calls chs.roots.eval(), with modified cy values.
}
\value{
If length(b) is one, and to.list is false (the default):\cr
It returns a single numeric vector.

Otherwise:\cr
It returns a list of numeric vectors, one for each value in b.
}
\arguments{
\item{a}{A chs object.}
\item{b}{A vector of y values.}
\item{to.list}{Refer to the value section.}
\item{\dots}{Other arguments, for chs.roots.eval().}
}
\examples{
#control points
cx = 1:4
cy = c (-4, -1, 1, 4)

#cubic Hermite spline
#(with function object, and default slopes)
f = chs (cx, cy)

#solve for x, for y=2
solve (f, 2)
}
\seealso{\link{chs.tangents}, \link{chs.slopes}, \link{chs.roots}, \link{chs.roots.eval}}
