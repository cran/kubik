\name{22_direct_chs_root_evaluation}
\alias{root.chs.eval}
\alias{argmin.chs.eval}
\alias{argmax.chs.eval}
\alias{argflex.chs.eval}
\alias{roots.chs.eval}
\alias{argmins.chs.eval}
\alias{argmaxs.chs.eval}
\alias{argflexs.chs.eval}
\alias{chs.roots.derivative.eval}
\title{Direct CHS Root Evaluation}
\description{Compute roots of cubic Hermite splines, directly.}
\usage{
root.chs.eval (cx, cy, cb, \dots,
    include.implied.roots=TRUE, warning=TRUE)
argmin.chs.eval (cx, cy, cb, \dots,
    include.implied.roots=TRUE, warning=TRUE)
argmax.chs.eval (cx, cy, cb, \dots,
    include.implied.roots=TRUE, warning=TRUE)
argflex.chs.eval (cx, cy, cb, \dots,
    include.implied.roots=TRUE, warning=TRUE)

roots.chs.eval (cx, cy, cb, \dots,
    include.implied.roots=TRUE, warning=TRUE)
argmins.chs.eval (cx, cy, cb, \dots,
    include.implied.roots=TRUE, warning=TRUE)
argmaxs.chs.eval (cx, cy, cb, \dots,
    include.implied.roots=TRUE, warning=TRUE)
argflexs.chs.eval (cx, cy, cb, \dots,
    include.implied.roots=TRUE, warning=TRUE)

chs.roots.derivative.eval (cx, cy, cb, \dots,
    include.implied.roots=TRUE, warning=TRUE)
}
\arguments{
\item{cx}{Vector of unique ascending control point x values.\cr (Preferably, equally-spaced).}
\item{cy}{Vector of control point y values.}
\item{cb}{Vector of control slopes.}
\item{include.implied.roots}{If true, include implied roots.}
\item{warning}{If true, give a warning if there's relevant constant/level (or linear) sections.}
\item{\dots}{Ignored.}
}
\details{
Refer to the help page for \link{roots.chs}, for more information.\cr
(The functions described in that help page are similar to these functions).

These functions (with a .eval suffix) compute roots of cubic Hermite splines, without using function objects, and with minimal error checking.\cr
Alternatively, you can use function objects, which are likely to be more convenient, in most cases.
}
\references{Please refer to the help page for \link{chs} for background information and references.}
\seealso{\link{chs}, \link{roots.chs}}
\examples{
#control points
cx <- 1:4
cy <- c (1, -1, -1, 1)

#control slopes
cb <- chs.slopes (cx, cy)

#roots
#(without function object)
roots.chs.eval (cx, cy, cb)
argmins.chs.eval (cx, cy, cb)
}
