#kubik: Cubic Hermite Splines and Related Optimization Methods
#Copyright (C), Abby Spurdle, 2020

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.params = function (cy1, cy2, cb1, cb2)
{	P0 = cy1
	P1 = cb1
	P2 = 3 * cy2 - 3 * cy1 - 2 * cb1 - cb2
	P3 = 2 * cy1 - 2 * cy2 + cb1 + cb2

	c (P0, P1, P2, P3)
}

.params.derivative = function (cy1, cy2, cb1, cb2)
{	p = .params (cy1, cy2, cb1, cb2)[2:4]
	(1:3) * p
}

.params.integral = function (cy1, cy2, cb1, cb2)
{	p = .params (cy1, cy2, cb1, cb2)
	1 / (1:4) * p
}

.params.2nd = function (cy1, cy2, cb1, cb2)
{	p = .params (cy1, cy2, cb1, cb2)[3:4]
	c (2, 6) * p
}

.params.3rd = function (cy1, cy2, cb1, cb2)
{	p = .params (cy1, cy2, cb1, cb2)[4]
	6 * p
}
