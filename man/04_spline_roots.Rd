\name{04_spline_roots}
\alias{chs.root}
\alias{chs.argmin}
\alias{chs.argmax}
\alias{chs.argflex}
\alias{chs.roots}
\alias{chs.argmins}
\alias{chs.argmaxs}
\alias{chs.argflexs}
\alias{chs.roots.derivative}
\title{Spline Roots}
\description{Compute roots of spline functions (with function objects).}
\usage{
chs.root (f, include.implied.roots=TRUE, warning=TRUE)
chs.argmin (f, include.implied.roots=TRUE, warning=TRUE)
chs.argmax (f, include.implied.roots=TRUE, warning=TRUE)
chs.argflex (f, include.implied.roots=TRUE, warning=TRUE,
    all.inflection.points=FALSE)

chs.roots (f, include.implied.roots=TRUE, warning=TRUE)
chs.argmins (f, include.implied.roots=TRUE, warning=TRUE)
chs.argmaxs (f, include.implied.roots=TRUE, warning=TRUE)
chs.argflexs (f, include.implied.roots=TRUE, warning=TRUE,
    all.inflection.points=FALSE)

chs.roots.derivative (f, include.implied.roots=TRUE, warning=TRUE,
    all.inflection.points=FALSE)
}
\details{
Refer the first help page, on chs.tangents() and chs.slopes(), for background information.

These functions compute roots of cubic Hermite splines (with function objects). Alternatively, you can you can use the other spline root functions (with a .eval suffix, and without function objects), which may be more efficient, in some cases.

chs.root() returns a single real root, if a single real root exists, otherwise, it returns an error. chs.roots() returns all the real roots, which will be a zero-length numeric vector, if there are no real roots.

chs.argmin() returns a single global argmin, if a single global argmin exists, otherwise, it returns an error. chs.argmins() returns all argmins, which will be a zero-length numeric vector, if there are no minima.

chs.argmax() and chs.argmaxs() are the same as chs.argmin() and chs.argmins(), except they're for argmax(s) rather than argmin(s). chs.argflex() and chs.argflexs() are similar to the functions above, except that they're for inflection points.

Unlike the other functions here, chs.roots.derivative() returns a two-column matrix.

The first column gives the roots of the derivative, and the second column gives a numeric value equivalent to the sign of the second derivative, so 1 (for minima), -1 (for maxima) and 0 (for stationary inflection points).

Here, an implied root refers to a root that coincides with a level section (or a linear section, for inflection points), and a level section refers to one or more adjacent spline segments, with a constant value. If include.implied.roots is true (the default), implied roots are included, otherwise, they're not.

If implied roots are included, then their value is equal to the midpoint of their interval.

Here, a spurious inflection point refers to an inflection point within a quadratic-suggestive section, and a quadratic-suggestive section refers to one or more adjacent spline segments, where all the control tangents within that section, are on the same side of the secant lines. If all.inflection.points is true, spurious inflection points are included, otherwise (the default), they're not.

Note that the chs.roots.derivative() function contain this option, however, I haven't tested if it's possible for an inflection point to be both stationary and spurious.

Also note that these functions may fail to find optima and inflection points, if they're extremely close to the control points. However, optima and non-spurious inflection points, shouldn't be extremely close to the control points, unless the control slopes are extremely large relative to the distances between the control points.

Other notes:

(1) If a level section (or relevant linear section, as above) contains two or more spline segments, then the midpoint of the combined interval is used.\cr
(2) In principle, linear sections are quadratic-suggestive sections, however, this is of little relevance.\cr
(3) Quadratic-suggestive sections may contain spline segments, where the secant line is equal to one of the surrounding control tangents.\cr
(4) (Possible) optima and inflection points at the outermost control points are ignored.
}
\arguments{
\item{f}{A chs object.}
\item{include.implied.roots}{If true, include implied roots.}
\item{warning}{If true, give a warning if there's relevant level (or linear) sections.}
\item{all.inflection.points}{If true, include spurious inflection points.}
}
\examples{
#control points
cx = 1:4
cy = c (1, -1, -1, 1)

#cubic Hermite spline
#(with function object, and default slopes)
f = chs (cx, cy)

#roots
#(with function object)
chs.roots (f)
chs.argmins (f)
}
\seealso{\link{chs.tangents}, \link{chs.slopes}, \link{chs.root.eval}, \link{chs.argmin.eval}, \link{chs.argmax.eval}, \link{chs.argflex.eval}, \link{chs.roots.eval}, \link{chs.argmins.eval}, \link{chs.argmaxs.eval}, \link{chs.argflexs.eval}, \link{chs.roots.derivative.eval}}
