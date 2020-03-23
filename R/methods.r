#kubik: Cubic Hermite Splines and Related Optimization Methods
#Copyright (C), Abby Spurdle, 2020

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

print.kspline = function (x, ...)
{	if ("package:intoo" %in% search () )
	{	str = "object.summary (x, ...)"
		eval (str2lang (str) )
	}
	else
		print.default (x, ...)
}

plot.kspline = function (x, ..., control.points=FALSE)
{	f = x

	. = attributes (f)
	x = seq (.$cx [1], .$cx [.$nc], length.out=200)
	y = f (x)
	plot (x, y, type="l", ...)
	if (control.points)
		points (f, ...)
}

lines.kspline = function (x, ..., control.points=FALSE)
{	f = x

	. = attributes (f)
	x = seq (.$cx [1], .$cx [.$nc], length.out=200)
	lines (x, f (x), ...)
	if (control.points)
		points (f, ...)
}

.points.kspline = function (f, ..., pch=16)
{	. = attributes (f)
	points (.$cx, f (.$cx), pch=pch, ...)
}

points.kspline = function (x, ...)
	.points.kspline (x, ...)
