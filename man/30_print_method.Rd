\name{30_print_method}
\alias{kub.printf.KSpline}
\title{Print and Plot Methods}
\description{Print method, for chs-related objects.}
\usage{
\method{kub.printf}{KSpline}(sf, \dots)
}
\arguments{
\item{sf}{A KSpline object, which is the superclass of CHS, CHSD and CHSI objects.}
\item{\dots}{Ignored.}
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
f
}
