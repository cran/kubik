\name{31_plot_methods}
\alias{kub.plotf.KSpline}
\alias{kub.linesf.KSpline}
\alias{kub.pointsf.KSpline}
\title{Plot Methods}
\description{Plot methods, for chs-related objects.}
\usage{
\method{kub.plotf}{KSpline}(sf, \dots, control.points=FALSE)
\method{kub.linesf}{KSpline}(sf, \dots, control.points=FALSE)
\method{kub.pointsf}{KSpline}(sf, \dots)
}
\arguments{
\item{sf}{A KSpline object, which is the superclass of CHS, CHSD and CHSI objects.}
\item{control.points}{If true, plot the control points (for CHS objects), or their corresponding points (for CHSD and CHSI).}
\item{\dots}{Other arguments.}
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
