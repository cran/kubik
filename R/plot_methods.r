#kubik: Cubic Hermite Splines and Related Foot Finding Methods
#Copyright (C), Abby Spurdle, 2019 to 2021

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

kub.plotf.KSpline = function (sf, ..., control.points=FALSE)
{	. = attributes (sf)
	x = seq (.$cx [1], .$cx [.$nc], length.out=200)
	y = sf (x)
	plot (x, y, type="l", ...)
	if (control.points)
		.kpoints (sf, ...)
}

kub.linesf.KSpline = function (sf, ..., control.points=FALSE)
{	. = attributes (sf)
	x = seq (.$cx [1], .$cx [.$nc], length.out=200)
	lines (x, sf (x), ...)
	if (control.points)
		.kpoints (sf, ...)
}

kub.pointsf.KSpline = function (sf, ...)
	.kpoints (sf, ...)

.kpoints = function (sf, ..., pch=16)
{	. = attributes (sf)
	points (.$cx, sf (.$cx), pch=pch, ...)
}
