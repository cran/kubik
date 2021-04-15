\name{20_runtime_function_objects}
\alias{Runtime Function Objects}
\alias{KSpline-class}
\alias{CHS-class}
\alias{CHSD-class}
\alias{CHSI-class}
\alias{ACHSD-class}
\alias{kub.function.object}
\title{Runtime Function Objects}
\description{
Hard-coded functions, representing (runtime) S4-based function objects.

DO NOT CALL THESE FUNCTIONS.\cr
CALL A CONSTRUCTOR, WHICH SHOULD RETURN A FUNCTION OBJECT.}
\usage{
kub.function.object (x)
}
\arguments{
\item{x}{
	A numeric vector (or array), giving the evaluation points.\cr
	This can include nonfinite values.}
}
\value{
All function objects return a numeric vector (or numeric array).

The dimensions of the input, determine the dimensions of the output.
}
