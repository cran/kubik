\name{07_spline_root_evaluation}
\alias{chs.root.eval}
\alias{chs.argmin.eval}
\alias{chs.argmax.eval}
\alias{chs.argflex.eval}
\alias{chs.roots.eval}
\alias{chs.argmins.eval}
\alias{chs.argmaxs.eval}
\alias{chs.argflexs.eval}
\alias{chs.roots.derivative.eval}
\title{Spline Roots}
\description{Compute roots of spline functions (without function objects).}
\usage{
chs.root.eval (nc, cx, cy, cb, include.implied.roots=TRUE, warning=TRUE)
chs.argmin.eval (nc, cx, cy, cb, include.implied.roots=TRUE, warning=TRUE)
chs.argmax.eval (nc, cx, cy, cb, include.implied.roots=TRUE, warning=TRUE)
chs.argflex.eval (nc, cx, cy, cb, include.implied.roots=TRUE, warning=TRUE,
    all.inflection.points=FALSE)

chs.roots.eval (nc, cx, cy, cb, include.implied.roots=TRUE, warning=TRUE)
chs.argmins.eval (nc, cx, cy, cb, include.implied.roots=TRUE, warning=TRUE)
chs.argmaxs.eval (nc, cx, cy, cb, include.implied.roots=TRUE, warning=TRUE)
chs.argflexs.eval (nc, cx, cy, cb, include.implied.roots=TRUE, warning=TRUE,
    all.inflection.points=FALSE)

chs.roots.derivative.eval (nc, cx, cy, cb,
    include.implied.roots=TRUE, warning=TRUE,
    all.inflection.points=FALSE)
}
\details{
Refer to chs.root(), chs.argmin() and chs.argmax(), etc, for more information.

These functions (with a .eval suffix) compute roots of cubic Hermite splines (without function objects), with little error checking. Alternatively, you can you can use the other spline root functions (with no suffix, and with function objects), which may be more convenient, in some cases.
}
\arguments{
\item{nc}{Number of control points.}
\item{cx}{Vector of unique ascending control point x values (preferably, equally spaced).}
\item{cy}{Vector of control point y values.}
\item{cb}{Vector of control slopes.}
\item{include.implied.roots}{If true, include implied roots.}
\item{warning}{If true, give a warning if there's relevant level (or linear) sections.}
\item{all.inflection.points}{If true, include spurious inflection points.}
}
\examples{
#number of control points
nc = 4

#control points
cx = 1:4
cy = c (1, -1, -1, 1)

#control slopes
cb = chs.slopes (cx, cy)

#roots
#(without function object)
chs.roots.eval (nc, cx, cy, cb)
chs.argmins.eval (nc, cx, cy, cb)
}
\seealso{\link{chs.tangents}, \link{chs.slopes}, \link{chs.root}, \link{chs.argmin}, \link{chs.argmax}, \link{chs.argflex}, \link{chs.roots}, \link{chs.argmins}, \link{chs.argmaxs}, \link{chs.argflexs}, \link{chs.roots.derivative}}
