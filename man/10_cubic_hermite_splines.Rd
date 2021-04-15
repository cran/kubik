\name{10_cubic_hermite_splines}
\alias{chs}
\alias{chs.derivative}
\alias{chs.integral}
\alias{approx.chs.derivative}
\alias{incr.chs}
\alias{decr.chs}
\title{Cubic Hermite Splines}
\description{Functions for constructing cubic Hermite splines (as S4-based function objects), and their derivatives and indefinite integrals.}
\usage{
chs (cx = 1:length (cy), cy, cb, \dots,
    constraints = chs.constraints (, \dots),
    transform=FALSE, outside = c (NA, NA), init.method)

chs.derivative (cx = 1:length (cy), cy, cb, \dots,
    constraints = chs.constraints (, \dots),
    transform=FALSE, outside = c (NA, NA), init.method)

chs.integral (cx = 1:length (cy), cy, cb, \dots,
    constraints = chs.constraints (, \dots),
    transform=FALSE, outside = c (NA, NA), init.method,
    constant=0)

approx.chs.derivative (cx = 1:length (cy), cy, cb, \dots,
    constraints = chs.constraints (, \dots),
    transform=FALSE, outside = c (NA, NA), init.method,
    apply.constraints.to=0, nth=1, trim=TRUE)

#call chs, setting the constraints arg
incr.chs (\dots)
decr.chs (\dots)
}
\arguments{
\item{cx}{
	Vector of unique ascending control point x values.\cr
	(Preferably, equally-spaced).}
\item{cy}{
	Vector of control point y values.}
\item{cb}{
	Vector of control slopes.\cr
	Optional, unless transform is true.}
\item{constraints}{
	A CHS.Constraints object.\cr
	Can be NULL or NA, for no constraints.\cr
	By default (if transform is false), ignored if cb supplied.\cr 
	Except for approximate derivatives).\cr\cr
	Refer to the help file for \link{chs.constraints} for more information.}
\item{transform}{Logical value, if true, transform cb (the user-supplied slopes), subject to constraints.}
\item{outside}{A vector of length two, giving the value of the spline outside the control points.}
\item{init.method}{
	Initialization method for slopes.\cr
	Possible values are "SL" and "SQ".\cr
	Ignored, if cb supplied.\cr
	(Except for approximate derivatives).\cr\cr
	Refer to the help file for \link{chs.slopes} for more information.}
\item{constant}{Constant term.}
\item{apply.constraints.to}{Integer, between 0 and nth, which derivative to apply constraints to.}
\item{nth}{Nonnegative integer, which derivative to approximate.}
\item{trim}{
	Logical value, if true, remove nth control points, from both the beginning and end of the spline.\cr
	For example, if nth=2, two control points are removed from each end.\cr\cr
	This should be true (the default), if non-default constraints are used.}
\item{\dots}{
	Alternative (possibly simpler) way of specifying constraints.\cr
	Ignored, if the constraints argument is set.\cr\cr
	This (specifying constraints in dots) should not be used inside packages, because the use of dots may change.}
}
\details{
These functions are constructors for function objects representing cubic Hermite splines, and their derivatives and indefinite integrals. 

The resulting function objects can be printed, plotted and evaluated (for x).\cr
Alternatively, you can use the spline eval functions (with a .eval suffix), without using function objects, which may be more efficient, in some cases.

\bold{chs}\cr
Construct a CHS object, for a cubic Hermite spline.

\bold{chs.derivative}\cr
Construct a CHSD object, for an (exact) first derivative of a cubic Hermite spline.

\bold{chs.integral}\cr
Construct a CHSI object, for an (indefinite) integral of a cubic Hermite spline.

\bold{approx.chs.derivative}\cr
Construct an ACHSD object (which is also a CHS object), for a (new) cubic Hermite spline representing a smooth approximation of the first, second or higher derivatives of an (initial) cubic Hermite spline.

Each of these objects requires two or more control points, and optionally their control slopes.

Note that currently, some functions may fail to work, for very small or very large input values.

Expanding on the constraints and transform arguments, there are three valid cases:\cr
(1) The slopes are omitted, in which case, default slopes are used subject to any constraints.\cr
(2) Slopes are provided and transform=FALSE (the default), in which case, user-supplied slopes are used.\cr
(3) Slopes are provided and transform=TRUE, in which case, the user-supplied slopes are transformed subject to any constraints, and the transformed slopes are used.

The approximate derivative is computed by constructing an initial cubic Hermite spline, then for each derivative, the cy values are set equal to the slopes, and new slopes are computed. This is unlikely to give a good result unless a large number of equally-spaced control points are used.

Also note that:\cr
(1) Exact derivatives aren't necessarily smooth.\cr
(2) In the chs.derivative and chs.integral objects, constraints and slopes apply to the original function, not the derivatives or integrals of that function.\cr
(3) In contrast, the outside values apply to the resulting derivatives or integrals.\cr
(4) If it's not possible to compute slopes that satisfy constraints, an error is generated.\cr
(5) In approximate derivatives, user-supplied constraints only apply to one iteration (determined by apply.constraints.to), with default constraints used for other iterations.\cr
(6) Approximate derivatives may give poor approximations in the tail regions, and applying non-default constraints to approximate derivatives may be problematic if tail regions are included.\cr
(7) Currently, constraints don't apply to outside values, however, this may be changed in the future.\cr
(8) If the control slopes are omitted, then they're computed using the same method as chs.slopes.

Refer to the help pages for \link{chs.constraints} and \link{chs.slopes} for more information on constraints and slopes.

Note that the help page on constraints, also discusses the relationship between the spline's shape, and its slopes. And the help page on slopes also discusses the initialization method.
}
\value{
Self-referencing S4-based function objects.

Refer to \link{Runtime Function Objects}.
}
\references{
Help pages from the stats package, for the splinefun, splinefunH and spline functions.

Wikipedia pages "Cubic Hermite spline" and "Monotone cubic interpolation".

Fritsch, F.N. & Carlson, R.E. (1980). Monotone Piecewise Cubic Interpolation. SIAM Journal on Numerical Analysis, 17 (2). doi:10.1137/0717021.
}
\seealso{\link{chs.constraints}, \link{apply.chs.constraints}, \link{chs.slopes}, \link{chs.eval}}
\examples{
############
#example (1)
############
#control points
cx <- 1:4
cy <- c (-4, -1, 1, 4)

#cubic hermite spline
#(with function object, and default slopes)
f <- chs (cx, cy)

#plot
#(with control points)
plot (f, control.points=TRUE)

#evaluate
f (3.5)

#add point to plot
points (3.5, f (3.5), pch=16, cex=2, col="blue")

############
#example (2)
############
#control points
#(sine wave)
cx <- seq (-2 * pi, 2 * pi, length.out=200)
cy <- sin (cx)

#cubic hermite spline, and approximate 4th derivative
f0 <- chs (cx, cy)
f4 <- approx.chs.derivative (cx, cy, nth=4)

#plot
#(approx 4th derivative in blue)
plot (f0)
lines (f4, col="blue")

#evaluate
#(results should be close to one)
f0 (pi / 2)
f4 (pi / 2)
}
