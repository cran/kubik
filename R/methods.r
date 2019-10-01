#kubik: Cubic Hermite Splines
#Copyright (C), Abby Spurdle, 2019

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

print.k.spline = function (x, ...)
	object.summary (x, ...)

plot.k.spline = function (x, with.points=FALSE, ..., pch=16)
{	f = x

	. = attributes (f)
	x = seq (.$cx [1], .$cx [.$nc], length.out=200)
	y = f (x)
	plot (x, y, type="l", ...)
	if (with.points)
		points (.$cx, f (.$cx), pch=pch)
}

lines.k.spline = function (x, with.points=FALSE, ..., pch=16)
{	f = x

	. = attributes (f)
	x = seq (.$cx [1], .$cx [.$nc], length.out=200)
	lines (x, f (x), ...)
	if (with.points)
		points (.$cx, f (.$cx), pch=pch, ...)
}

points.k.spline = function (x, ..., pch=16)
{	f = x

	. = attributes (f)
	points (.$cx, f (.$cx), pch=pch, ...)
}
