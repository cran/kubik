\name{03_print_and_plot_methods}
\alias{print.k.spline}
\alias{plot.k.spline}
\alias{lines.k.spline}
\alias{points.k.spline}
\title{Print and Plot Methods}
\description{Print and plot methods, for chs and related objects.}
\usage{
\method{print}{k.spline}(x, \dots)

\method{plot}{k.spline}(x, with.points=FALSE, \dots, pch=16)
\method{lines}{k.spline}(x, with.points=FALSE, \dots, pch=16)
\method{points}{k.spline}(x, \dots, pch=16)
}
\arguments{
\item{x}{A k.spline object, which is the superclass of chs, chs.derivative and chs.integral objects.}
\item{with.points}{If true, plot the control points (for chs objects), or their corresponding points (for chs.derivative and chs.integral objects).}
\item{\dots}{Other arguments.}
\item{pch}{.}
}
\examples{
#control points
cx = 1:4
cy = c (-4, -1, 1, 4)

#cubic Hermite spline
#(with function object)
f = chs (cx, cy)

#print the object
#(with auto-printing)
f

#plot the object
plot (f, TRUE)
}
\seealso{\link{chs.tangents}, \link{chs.slopes}, \link{chs}, \link{chs.derivative}, \link{chs.integral}, \link{approx.chs.derivative}}

