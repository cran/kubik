\name{00_bridging_functions}
\alias{print.KSpline}
\alias{plot.KSpline}
\alias{lines.KSpline}
\alias{points.KSpline}
\alias{kub.printf}
\alias{kub.plotf}
\alias{kub.linesf}
\alias{kub.pointsf}
\title{Bridging Functions}
\description{Bridging functions for printing/plotting function objects.}
\usage{
\method{print}{KSpline}(x, \dots)
\method{plot}{KSpline}(x, \dots)
\method{lines}{KSpline}(x, \dots)
\method{points}{KSpline}(x, \dots)

kub.printf (\dots)
kub.plotf (\dots)
kub.linesf (\dots)
kub.pointsf (\dots)
}
\arguments{
\item{x}{The first argument.} 
\item{\dots}{Arguments for specific methods.}
}
\details{
In general, do not call kub.printf, etc.\cr
Call print and plot instead.\cr
(These call kub.printf and kub.plotf).
}
