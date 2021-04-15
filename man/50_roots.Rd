\name{50_roots}
\alias{root.chs}
\alias{argmin.chs}
\alias{argmax.chs}
\alias{argflex.chs}
\alias{roots.chs}
\alias{argmins.chs}
\alias{argmaxs.chs}
\alias{argflexs.chs}
\alias{chs.roots.derivative}
\title{Roots}
\description{Compute roots of cubic Hermite splines, using function objects.}
\usage{
root.chs (sf, \dots, include.implied.roots=TRUE, warning=TRUE)
argmin.chs (sf, \dots, include.implied.roots=TRUE, warning=TRUE)
argmax.chs (sf, \dots, include.implied.roots=TRUE, warning=TRUE)
argflex.chs (sf, \dots, include.implied.roots=TRUE, warning=TRUE)

roots.chs (sf, \dots, include.implied.roots=TRUE, warning=TRUE)
argmins.chs (sf, \dots, include.implied.roots=TRUE, warning=TRUE)
argmaxs.chs (sf, \dots, include.implied.roots=TRUE, warning=TRUE)
argflexs.chs (sf, \dots, include.implied.roots=TRUE, warning=TRUE)

chs.roots.derivative (sf, \dots, include.implied.roots=TRUE, warning=TRUE)
}
\arguments{
\item{sf}{A CHS object.}
\item{include.implied.roots}{If true, include implied roots.}
\item{warning}{If true, give a warning if there's relevant constant/level (or linear) sections.}
\item{\dots}{Ignored.}
}
\details{
These functions compute roots of cubic Hermite splines, using function objects.\cr
Alternatively, you can you can use the spline root eval functions (with a .eval suffix), without using function objects, which may be more efficient, in some cases.

By default, the functions include implied roots (defined below) and remove spurious inflection points (also, defined below).

Here, an implied root refers to a root that coincides with a constant/level section.\cr
(Or coincides with a linear section, for inflection points).

Here, a constant/level section refers to one or more consecutive spline segments, with a constant value.

Here, a spurious inflection point refers to an inflection point within a quadratic-suggestive section.

And here, a quadratic-suggestive section refers to one or more consecutive spline segments, where all the control tangents within that section, are on the same side of the secant lines between consecutive control points.

If implied roots are included, then their value is equal to the midpoint of their combined interval.

Notes:\cr
(1) Optima and inflection points exclude the endpoints.\cr
(2) Currently, these functions may fail to find roots, if they're extremely close to the control points. Except, the functions should find optima in well-defined splines (with no more than one optimum per spline segment), regardless of where the optima are located.\cr
(3) The chs.roots.derivative function is likely to be deprecated in the future, and replaced by a better function.\cr
}
\value{
root.chs returns a single real root, if a single real root exists, otherwise, it returns an error.\cr
roots.chs returns all the real roots, which will be a zero-length numeric vector, if there are no real roots.

argmin.chs returns a single global argmin, if a single global argmin exists, otherwise, it returns an error.\cr
argmins.chs returns all argmins, which will be a zero-length numeric vector, if there are no minima.

argmax.chs and argmaxs.chs are the same as argmin.chs and argmins.chs, except they're for argmax(s) rather than argmin(s).\cr
argflex.chs() and argflexs.chs are similar to the functions above, except that they're for inflection points.

Unlike the other functions here, chs.roots.derivative returns a two-column matrix.\cr
The first column gives the roots of the derivative.\cr
And the second column gives a numeric value equal to the sign of the second derivative, so:\cr
1 (for minima), -1 (for maxima) and 0 (for stationary inflection points).
}
\references{Please refer to the help page for \link{chs} for background information and references.}
\seealso{\link{chs}, \link{solve.CHS}, \link{roots.chs.eval}}
\examples{
#control points
cx <- 1:4
cy <- c (1, -1, -1, 1)

#cubic hermite spline
#(with function object, and default slopes)
f <- chs (cx, cy)

#roots
roots.chs (f)
argmins.chs (f)
}
