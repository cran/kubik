#kubik: Cubic Hermite Splines and Related Foot Finding Methods
#Copyright (C), Abby Spurdle, 2019 to 2021

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.chs.eval = function (is.derivative, cx, cy, cb, x, outside)
{	nc = length (cx)
	nx = length (x)
	y = numeric (nx)
	for (i in seq_len (nx) )
	{	if (is.na (x [i]) )
			y [i] = NA
		else if (x [i] < cx [1])
			y [i] = outside [1]
		else if (x [i] > cx [nc])
			y [i] = outside [2]
		else
		{	nI = sum (cx <= x [i])
			if (is.derivative)
			{	if (cx [nI] == x [i])
					y [i] = cb [nI]
				else
					y [i] = .interval.derivative.eval (cx [nI], cx [nI + 1], cy [nI], cy [nI + 1], cb [nI], cb [nI + 1], x [i])
			}
			else
			{	if (cx [nI] == x [i])
					y [i] = cy [nI]
				else
					y [i] = .interval.eval (cx [nI], cx [nI + 1], cy [nI], cy [nI + 1], cb [nI], cb [nI + 1], x [i])
			}
		}
	}
	if (is.array (x) )
		dim (y) = dim (x)
	y
}

chs.eval = function (cx, cy, cb, x, ..., outside = c (NA, NA) )
	.chs.eval (FALSE, cx, cy, cb, x, outside)

chs.derivative.eval = function (cx, cy, cb, x, ..., outside = c (NA, NA) )
	.chs.eval (TRUE, cx, cy, cb, x, outside)

chs.integral.eval = function (cx, cy, cb, x, ..., outside = c (NA, NA), constant=0)
{	nc = length (cx)
	nx = length (x)
	if (nx == 0)
		y = numeric ()
	else
	{	areas = numeric (nc - 1)
		for (i in 1:(nc - 1) )
			areas [i] = .interval.integral.a2b (cx [i], cx [i + 1], cy [i], cy [i + 1], cb [i], cb [i + 1])
		
		y = rep (0, nx)
		for (i in 1:nx)
		{	if (is.na (x [i]) )
				y [i] = NA
			else if (x [i] < cx [1])
				y [i] = outside [1]
			else if (x [i] > cx [nc])
				y [i] = outside [2]
			else
			{	nI = sum (cx <= x [i])
				y [i] = constant
				if (nI > 1)
					y [i] = y [i] + sum (areas [1:(nI - 1)])
				if (cx [nI] != x [i])
					y [i] = y [i] + .interval.integral.a2x (cx [nI], cx [nI + 1], cy [nI], cy [nI + 1], cb [nI], cb [nI + 1], x [i])
			}
		}
	}
	if (is.array (x) )
		dim (y) = dim (x)
	y
}
