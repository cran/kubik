\name{16_slopes_and_tangents}
\alias{chs.slopes}
\alias{chs.tangents}
\alias{apply.chs.constraints}
\title{Slopes and Tangents}
\description{Functions for computing (possibly) suitable slopes and tangents.}
\usage{
chs.slopes (cx, cy, \dots, constraints = chs.constraints (, \dots),
    init.method)
chs.tangents (cx, cy, \dots, constraints = chs.constraints (, \dots),
    init.method)

apply.chs.constraints (cx, cy, cb, \dots,
    constraints = chs.constraints (, \dots) )
}
\arguments{
\item{cx}{Vector of unique ascending control point x values.\cr (Preferably, equally-spaced).}
\item{cy}{Vector of control point y values.}
\item{cb}{Vector of control slopes.}
\item{constraints}{A chs.constraints object.\cr Can be NULL or NA, for no constraints.\cr\cr Refer to the help file for \link{chs.constraints} for more information.}
\item{init.method}{Initialization method for slopes.\cr Possible values are "SL" and "SQ".}
\item{\dots}{Alternative (possibly simpler) way of specifying constraints.\cr Ignored if a chs.constraints object supplied.\cr\cr This (specifying constraints in dots) should not be used inside packages, because the use of dots may change.}
}
\value{
chs.slopes returns a (possibly) suitable vector of slopes.

chs.tangents returns a (possibly) suitable two-column matrix, giving the intercepts and slopes.

appl.chs.constraints returns a transformed vector of slopes.
}
\details{
These functions work by computing a initial vector of slopes, and then trying to transform them to satisfy constraints.

Currently, there are two methods that may be used to initialize the slopes:\cr
SL, Simple secant-based method with near-linear tails.\cr
SQ, Modified secant based-method with quadratic tails.

If the initialization method isn't specified by the user, it defaults to "SQ", if the control points are equally-spaced (or near-equally-spaced), and it defaults to "SL", otherwise.

If there's only one spline segment:\cr
(Regardless of of what initialization method is used).\cr
The initial slopes are set equal to the secant's slope.

Where the secant refers to the secant line between the surrounding control points.

And if there's two more spline segments:\cr
The second slope is computed from the (double-interval) secant line between the first and third control points, the third slope is computed from the (double-interval) secant line between the second and fourth control points, and so on, up to the second to last slope, which is computed from the (double-interval) secant line between the third-to-last and last control points.

For the SL method, the the outermost slopes are set equal to the adjacent secants' slopes, resulting in near-linear tails. For the SQ method, the first slope is computed by fitting a quadratic polynomial to the first three control points and evaluating the derivative at the first point, and the last slope is computed by fitting a quadratic polynomial to the last three control points and evaluating the derivative at the last point.

The SQ method should only be used for equally-spaced control points.

Note that while secants may sound linear, computing slopes from secants gives similar results to computing all the slopes via quadratic polynomials.\cr
(And precedence is given to the SQ method, because the level of resulting curvature, is relatively consistent over the entire spline, including the tails).

The apply.chs.constraints function tries to transform a vector of slopes to satisfy a set of constraints.

If a transformation isn't possible, an error is generated.\cr
(This applies to all three functions).

Note that it's possible to use chs.slopes to compute slopes, and then modify them as desired.
}
\references{Please refer to the help page for \link{chs} for background information and references.}
\seealso{\link{chs}, \link{chs.constraints}}
\examples{
#control points
cx <- 1:4
cy <- c (-4, -1, 1, 4)

#control slopes
chs.slopes (cx, cy)

#control tangents
#(intercepts and slopes)
chs.tangents (cx, cy)
}
