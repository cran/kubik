\name{rd_7_spline_function_evaluation}
\alias{chs.eval}
\alias{chs.derivative.eval}
\alias{chs.integral.eval}
\title{Spline Function Evaluation}
\description{Functions for evaluating cubic Hermite splines (without function objects).}
\usage{
chs.eval (nc, cx, cy, cb, x, outside = c (NA, NA) )
chs.derivative.eval (nc, cx, cy, cb, x, outside = c (NA, NA) )
chs.integral.eval (nc, cx, cy, cb, x, outside = c (NA, NA), constant=0)
}
\details{
Refer the first help page, on chs.tangents() and chs.slopes(), for background information.

These functions (with a .eval suffix) evaluate cubic Hermite splines (without function objects), and their derivatives and integrals, with little error checking. Alternatively, you can use the other spline functions (with no suffix, and with function objects), which may be more convenient, in some cases.

\bold{chs.eval()}\cr
Evaluate cubic Hermite splines.

\bold{chs.deritivative.eval()}\cr
Evaluate (exact) derivatives of cubic Hermite splines.

\bold{chs.integral.eval()}\cr
Evaluate indefinite integrals of cubic Hermite splines.

Each of these functions requires two or more control points, and control slopes.

Note that exact derivatives aren't necessarily smooth.
}
\arguments{
\item{nc}{Number of control points.}
\item{cx}{Vector of unique ascending control point x values (preferably, equally spaced).}
\item{cy}{Vector of control point y values.}
\item{cb}{Vector of control slopes.}
\item{x}{Vector of x values, where the spline is evaluated at.}
\item{outside}{A vector of length two, giving the value of the spline outside the control points.}
\item{constant}{Constant term.}
}
\examples{
#number of control points
nc = 4

#control points
cx = 1:4
cy = c (-4, -1, 1, 4)

#control slopes
cb = chs.slopes (cx, cy)

#evaluate
#(without function object)
chs.eval (nc, cx, cy, cb, 3.5)
}
\seealso{\link{chs.tangents}, \link{chs.slopes}, \link{chs}, \link{chs.derivative}, \link{chs.integral}, \link{approx.chs.derivative}}
