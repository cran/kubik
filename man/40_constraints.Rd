\name{40_constraints}
\alias{CHS.Constraints-class}
\alias{chs.constraints}
\title{Constraints}
\description{Constraints for cubic Hermite splines.}
\usage{
chs.constraints (correction=TRUE, \dots,
    increasing=FALSE, decreasing=FALSE,
    flim = c (fmin, fmax),
        fmin = -Inf, fmax = +Inf)
}
\arguments{
\item{correction}{Logical value, if true, apply the main correction method.}
\item{increasing}{Logical value, if true, a non-decreasing function.}
\item{decreasing}{Logical value, if true, a non-increasing function.}
\item{flim}{A ascending numeric length-two vector giving minimum and maximum allowed function values.}
\item{fmin, fmax}{
	Numeric scalars.\cr
	Alternative way, of specifying flim.}
\item{\dots}{
	Should be empty.\cr
	An error is produced, if not.}
}

\details{
A CHS.Constraints object defines constraints for cubic Hermite splines.

Such constraints, only apply to computation of control slopes.\cr
Not the control points.

So, the control point positions need to be suitable.\cr
(i.e. You can't apply an increasing constraint to decreasing control points).

The main correction method is derived from the Fritsch-Carlson algorithm, and prevents local optima from appearing within each spline segment, except where either of the following applies:\cr
(1) The sign of the secant's slope is zero, and the surrounding control slopes have opposite signs.\cr
(2) The secant's slope has the opposite sign to either of the surrounding control slopes.

Where the secant refers to the secant line between the surrounding control points.

In general, if the control points are equally-spaced (in x), are ascending or descending (in y), and the default slopes are used, then the main correction method should be sufficient to construct a monotonically increasing or decreasing function. Even so, I recommend you apply increasing or decreasing constraints, if you need an increasing or decreasing function.

Currently, both monotonicity and bound constraints, will set some slopes to zero, if the initial slopes cause the spline to break the constraints. This is not an optimal transformation, and hopefully it will be improved in the future.

Note that (re: bounded splines) this package allows control points to be exactly equal to flim, and I can not guarantee that very small floating point errors won't occur, when evaluating such splines.

If very small floating point errors are problematic, then users can either:\cr
(1) Ensure there's a relatively small margin between the control points and the flim values.\cr
(2) Apply a correction method after evaluating splines. 

Note that signs can be zero, hence not having the same sign, doesn't necessarily mean having opposite signs.

Also note that, given a single spline segment:\cr
(1) If the surrounding control slopes have the same sign, a pair of local optima may or may not be present.\cr
(2) If the surrounding control slopes have opposite signs, a single local optimum should be present.\cr
(3) If the surrounding control tangent lines are on the same side of the secant line, a single local inflection point may or may not be present.\cr
(4) If the surrounding control tangent lines are on opposite sides of the secant line, a single local inflection point should be present.
}
\value{
A CHS.Constraints object.

This is an S4 object, designed for use within the \link{chs} function.
}
\references{Please refer to the help page for \link{chs} for background information and references.}
\seealso{\link{chs}, \link{apply.chs.constraints}, \link{chs.slopes}}
\examples{
#constraints
#f (x) bounded within [0, 1]
constraints <- chs.constraints (flim = c (0, 1) )

#control points
#(control points close to bound-constraints)
cx <- 1:6
cy <- c (0.001, 0.999, 0.999, 0.001, 0.001, 0.999)

#default and bounded splines
f.def <- chs (cx, cy)
f.bnd <- chs (cx, cy, constraints=constraints)

#plot
#(bounded spline in blue)
plot (f.def)
abline (h = c (0, 1), lty=2, col="grey")
lines (f.bnd, col="blue")
}
