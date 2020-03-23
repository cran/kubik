\name{21_direct_chs_evaluation}
\alias{chs.eval}
\alias{chs.derivative.eval}
\alias{chs.integral.eval}
\title{Direct CHS Evaluation}
\description{Functions for evaluating cubic Hermite splines, and their derivatives and indefinite integrals, directly.}
\usage{
chs.eval (cx, cy, cb, x, \dots,
    outside = c (NA, NA) )

chs.derivative.eval (cx, cy, cb, x, \dots,
    outside = c (NA, NA) )

chs.integral.eval (cx, cy, cb, x, \dots,
    outside = c (NA, NA), constant=0)
}
\arguments{
\item{cx}{Vector of unique ascending control point x values.\cr (Preferably, equally-spaced).}
\item{cy}{Vector of control point y values.}
\item{cb}{Vector of control slopes.}
\item{x}{Vector of x values, where the spline is evaluated at.}
\item{outside}{A vector of length two, giving the value of the spline outside the control points.}
\item{constant}{Constant term.}
\item{\dots}{.}
}
\details{
Refer to the help page for \link{chs}, for more information.\cr
(The functions described in that help page are similar to these functions).

These functions (with a .eval suffix) evaluate cubic Hermite splines, and their derivatives and indefinite integrals, without using function objects, and with minimal error checking.\cr
Alternatively, you can use function objects, which are likely to be more convenient, in most cases.

\bold{chs.eval}\cr
Evaluate cubic Hermite splines.

\bold{chs.deritivative.eval}\cr
Evaluate (exact) derivatives of cubic Hermite splines.

\bold{chs.integral.eval}\cr
Evaluate indefinite integrals of cubic Hermite splines.
}

\references{Please refer to the help page for \link{chs} for background information and references.}
\seealso{\link{chs}}
\examples{
#control points
cx <- 1:4
cy <- c (-4, -1, 1, 4)

#control slopes
cb <- chs.slopes (cx, cy)

#evaluate
#(without function object)
chs.eval (cx, cy, cb, 3.5)
}
