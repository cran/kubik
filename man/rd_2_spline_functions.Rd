\name{rd_2_spline_functions}
\alias{chs}
\alias{chs.derivative}
\alias{chs.integral}
\alias{approx.chs.derivative}
\title{Spline Functions}
\description{Functions for constructing cubic Hermite splines (with functions objects).}
\usage{
chs (cx, cy, cb, correction=TRUE, outside = c (NA, NA) )
chs.derivative (cx, cy, cb, correction=TRUE, outside = c (NA, NA) )
chs.integral (cx, cy, cb, correction=TRUE, outside = c (NA, NA), constant=0)

approx.chs.derivative (cx, cy, cb, correction=TRUE, outside = c (NA, NA), nth=1)
}
\details{
Refer the first help page, on chs.tangents() and chs.slopes(), for background information.

All of these functions are constructors for (S3-based) function objects representing cubic Hermite splines, and their derivatives and integrals. The resulting function objects can be printed, plotted and evaluated (for x). Alternatively, you can use the other spline functions (with a .eval suffix, and without function objects), which may be more efficient, in some cases.

\bold{chs()}\cr
Construct a chs object, for cubic Hermite splines.

\bold{chs.derivative()}\cr
Construct a chs.derivative object, for (exact) first derivatives of cubic Hermite splines.

\bold{chs.integral()}\cr
Construct a chs.integral object, for indefinite integrals of cubic Hermite splines.

\bold{approx.chs.derivative()}\cr
Construct a chs object, for (new) cubic Hermite splines representing smooth approximations of first, second or higher derivatives of (initial) cubic Hermite splines.

Each of these objects requires two or more control points, and optionally their control slopes. The control points and control slopes become attributes of the resulting function objects, with the exception of chs objects returned by approx.chs.derivative(), which have new cy and cb values.

Note that exact derivatives aren't necessarily smooth.

Also note that if the control slopes are omitted, then they're computed using the same method as chs.slopes().
}
\value{
All the resulting function objects are functions of the form:

function (x) = \{...\}

Where x is a numeric vector.

And they have class and other attributes.
}
\arguments{
\item{cx}{Vector of unique ascending control point x values (preferably, equally spaced).}
\item{cy}{Vector of control point y values.}
\item{cb}{Optional vector of control slopes.}
\item{correction}{If the control slopes (cb) aren't supplied and correction is true, then a correction method is used when computing them.}
\item{outside}{A vector of length two, giving the value of the spline outside the control points.}
\item{constant}{Constant term.}
\item{nth}{Which derivative to approximate.}
}
\examples{
#control points
cx = 1:4
cy = c (-4, -1, 1, 4)

#cubic Hermite spline
#(with function object, and default slopes)
f = chs (cx, cy)

#evaluate
#(with function object)
f (3.5)
}
\seealso{\link{chs.tangents}, \link{chs.slopes}, \link{chs.eval}, \link{chs.derivative.eval}, \link{chs.integral.eval}}
