#kubik: Cubic Hermite Splines and Related Foot Finding Methods
#Copyright (C), Abby Spurdle, 2019 to 2021

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.info.CHS   = c ("Cubic Hermite Spline", "f (x)")
.info.CHSD  = c ("Derivative of CHS", "f' (x)")
.info.CHSI  = c ("Indefinite Integral of CHS ", "F (x)")
.info.ACHSD = c ("Approx nth Derivative of CHS", "f* (x)")

setClass ("KSpline", contains="function",
	slots = c (
		.class.info="character",
		nc="integer",
		cx="numeric",
		cy="numeric",
		cb="numeric",
		outside="numeric"
	) )

setClass ("CHS.Constraints",
	slots = c (
		.bounded="logical",
		correction="logical",
		increasing="logical",
		decreasing="logical",
		flim="numeric"
	) )

setClass ("CHS", contains="KSpline")
setClass ("CHSD", contains="KSpline")
setClass ("CHSI", contains="KSpline", slots = c (constant="numeric") )
setClass ("ACHSD", contains="CHS")

setMethod ("show", "KSpline", function (object) print (object) )

kub.printf = function (...) UseMethod ("kub.printf")
kub.plotf = function (...) UseMethod ("kub.plotf")
kub.linesf = function (...) UseMethod ("kub.linesf")
kub.pointsf = function (...) UseMethod ("kub.pointsf")

print.KSpline = function (x, ...) kub.printf (x, ...)
plot.KSpline = function (x, ...) kub.plotf (x, ...)
lines.KSpline = function (x, ...) kub.linesf (x, ...)
points.KSpline = function (x, ...) kub.pointsf (x, ...)

.THAT = function ()
{	sf = sys.function (-1)
	attributes (sf)
}

kub.printf.KSpline = function (sf, ...)
{	cat ("<S4-Based Function Object>\n")
	cat (sf@.class.info [1], "\n", sep="")
	cat ("    ", sf@.class.info [2], "\n", sep="")
}

#function objects
#(for rd files)
kub.function.object = function (x) 0

.arg.error = function (...)
{	expr = format ( (sys.call (-1) ) )
	n = length (list (...) )
	if (n > 0)
	{	warning ("unsupported args, check arg names")

		cat ("call with unsupported args:\n")
		print (expr)
		cat ("check for incorrect argument names\n")
		cat ("check for unnamed non-leading arguments\n")
		cat ("e.g. in chs (x, y, slopes, constraints), constraints\n")
	}
}
