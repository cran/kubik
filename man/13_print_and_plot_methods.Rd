\name{13_print_and_plot_methods}
\alias{print.kspline}
\alias{plot.kspline}
\alias{lines.kspline}
\alias{points.kspline}
\title{Print and Plot Methods}
\description{Print and plot methods, for chs-related objects.}
\usage{
\method{print}{kspline}(x, \dots)

\method{plot}{kspline}(x, \dots, control.points=FALSE)
\method{lines}{kspline}(x, \dots, control.points=FALSE)
\method{points}{kspline}(x, \dots)
}
\arguments{
\item{x}{A kspline object, which is the superclass of chs, chs.derivative and chs.integral objects.}
\item{control.points}{If true, plot the control points (for chs objects), or their corresponding points (for chs.derivative and chs.integral objects).}
\item{\dots}{Other arguments.}
}
\details{
The print method calls intoo::object.summary, if the intoo package is on the search path.

Otherwise, it calls the default print method.
}
\references{Please refer to the help page for \link{chs} for background information and references.}
\seealso{\link{chs}}
\examples{
#control points
cx <- 1:4
cy <- c (-4, -1, 1, 4)

#cubic hermite spline
#(with function object, and default slopes)
f <- chs (cx, cy)

#plot the object
#(with control points)
plot (f, control.points=TRUE)
}
