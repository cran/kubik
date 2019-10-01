\name{rd_6_spline_quadraticization}
\alias{chs.bilinearize}
\title{Spline Quadraticization}
\description{
A function to make some spline sections more quadratic-like, and possibly minimize the number of spurious inflection points.

Please use caution.\cr
This function has had less design work and testing than other functions in this package.

Please email me, if you would like to discuss this subject further.
}
\usage{
chs.bilinearize (cx, cy, cb, correction=TRUE, at,
    distance = factor * mean (diff (cx) ),
    factor=0.125)
}
\details{
Refer the first help page, on chs.tangents() and chs.slopes(), for background information.\cr
Refer to the spline roots page, including chs.argflexs(), for more information, about linear and quadratic-suggestive sections.

chs.bilinearize() replaces the control points and control slopes, specified by the at parameter, with a pair of closely spaced control points, such that the two original surrounding spline segments are replaced by three new spline segments, linear-something-linear.

Also, this can be used to approximate angular splines, by using a small distance parameter.

The at parameter, which can be a vector, needs to have unique non-adjacent values, not equal to 1 or nc.

The distance parameter determines how far from the original control point, each new control point is, on the x axis. So, the distance between each resulting pair of control points, will be double this value. It needs to be positive, and small enough, such that the resulting cx values are still unique and ascending.

Note that if the control slopes are omitted, then they're computed using the same method as chs.slopes().
}
\value{
A named list with cx, cy and cb vectors.
}
\arguments{
\item{cx}{Vector of unique ascending control point x values (preferably, equally spaced).}
\item{cy}{Vector of control point y values.}
\item{cb}{Optional vector of control slopes.}
\item{correction}{If the control slopes (cb) aren't supplied and correction is true, then a correction method is used when computing them.}
\item{at}{Vector of which control points to replace.}
\item{distance}{Scalar (or vector with the same length as at), giving the distances from the original control points.}
\item{factor}{.}
}
\examples{
require (zeallot)

cx = 1:5
cy = c (0, 1, 2, 2.1, 2.2)
c (cx2, cy2, cb2) \%<-\% chs.bilinearize (cx, cy, at=3)
f1 = chs (cx, cy)
f2 = chs (cx2, cy2, cb2)
plot (f1, TRUE)
lines (f2, TRUE, col="blue")
}
\seealso{\link{chs.tangents}, \link{chs.slopes}, \link{chs}, \link{chs.argflexs} }
