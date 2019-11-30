\name{01_spline_tangents}
\alias{chs.tangents}
\alias{chs.slopes}
\title{Spline Tangents}
\description{
Functions for computing (possibly) suitable tangents or slopes.

Note that the details section contains background information, relevant to this package.
}
\usage{
chs.tangents (cx, cy, correction=TRUE)
chs.slopes (cx, cy, correction=TRUE)
}
\details{
Here, a cubic Hermite spline is defined by a sequence of:\cr
(nc) control points and (nc) control slopes.

Each control point has an x (cx) value, y (cy) value and a related slope (cb) value.

This results in (nc - 1) spline segments.\cr
Each spline segment has a pair of surrounding control points, and a pair of surrounding control tangents, however, in general, only the slopes of the tangents are required.

chs.tangents() computes a (possibly) suitable matrix of tangents, and chs.slopes() computes a (possibly) suitable vector of slopes.

If there's only one spline segment:\cr
The initial surrounding control slopes are set equal to the secant's slope.

And if there's two or more spline segments:\cr
Each initial slope, except for the outermost slopes, is computed by fitting a quadratic to three consecutive control points, and then computing the derivative of the quadratic. The derivative is evaluated at the second (i.e. mid) control point. The outermost slopes are set equal to the secants' slopes.

By default, a correction method is used, that prevents optima from appearing within each spline segment, except where either of the following applies:

(1) The spline segment has a constant value (so, the sign of the secant's slope is zero) and the surrounding control slopes have opposite signs.\cr
(2) The secant's slope has the opposite sign to either of the surrounding control slopes.

Where the secant refers to the secant line between the surrounding control points.

If the cy values are ascending or descending, it also results in slopes, that constrain the spline to be monotonically increasing or monotonically decreasing.

Note that it's possible to use chs.slopes() to compute the slopes, and then modify them as desired, which may be particularly useful where some slopes are known, or for the first and last slopes.}
\arguments{
\item{cx}{Vector of unique ascending control point x values (preferably, equally spaced).}
\item{cy}{Vector of control point y values.}
\item{correction}{If true, use a correction method.}
}
\examples{
#control points
cx = 1:4
cy = c (-4, -1, 1, 4)

#control tangents
#(intercepts and slopes)
chs.tangents (cx, cy)

#control slopes
chs.slopes (cx, cy)
}
\references{
Help pages from the stats package, for splinefun(), splinefunH() and spline().

Wikipedia pages "Cubic Hermite spline" and "Monotone cubic interpolation".

Fritsch, F.N. & Carlson, R.E. (1980). Monotone Piecewise Cubic Interpolation. SIAM Journal on Numerical Analysis, 17 (2). doi:10.1137/0717021.
}
