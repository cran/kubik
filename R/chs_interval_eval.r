#kubik: Cubic Hermite Splines and Related Foot Finding Methods
#Copyright (C), Abby Spurdle, 2019 to 2021

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.interval.eval = function (cx1, cx2, cy1, cy2, cb1, cb2, x)
{	dx = cx2 - cx1
	cb1 = dx * cb1
	cb2 = dx * cb2
	p = .params (cy1, cy2, cb1, cb2)
	x = (x - cx1) / dx

	sum (p * x ^ (0:3) )
}

.interval.derivative.eval = function (cx1, cx2, cy1, cy2, cb1, cb2, x)
{	dx = cx2 - cx1
	cb1 = dx * cb1
	cb2 = dx * cb2
	p = .params.derivative (cy1, cy2, cb1, cb2)
	x = (x - cx1) / dx

	sum (p * x ^ (0:2) ) / dx
}

.interval.integral.a2b = function (cx1, cx2, cy1, cy2, cb1, cb2)
{	dx = cx2 - cx1
	cb1 = dx * cb1
	cb2 = dx * cb2
	p = .params.integral (cy1, cy2, cb1, cb2)

	dx * sum (p)
}

.interval.integral.a2x = function (cx1, cx2, cy1, cy2, cb1, cb2, x)
{	dx = cx2 - cx1
	cb1 = dx * cb1
	cb2 = dx * cb2
	p = .params.integral (cy1, cy2, cb1, cb2)
	x = (x - cx1) / dx

	dx * sum (p * x ^ (1:4) )
}
