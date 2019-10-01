#kubik: Cubic Hermite Splines
#Copyright (C), Abby Spurdle, 2019

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

#level section -> -1
.interval.roots = function (cx1, cx2, cy1, cy2, cb1, cb2, include.lower, include.upper)
{	if (cy1 == 0 && cy2 == 0 && cb1 == 0 && cb2 == 0)
		list (-1, numeric () )
	else
	{	dx = cx2 - cx1
		cb1 = dx * cb1
		cb2 = dx * cb2
		p = .params (cy1, cy2, cb1, cb2)

		v = .cubic.roots (p, cy1, cy2, include.lower, include.upper)
		if (v [[1]] > 0)
			v [[2]] = .bt (cx1, cx2, dx, v [[2]])
		v
	}
}

#level section -> -1
.interval.roots.derivative = function (cx1, cx2, cy1, cy2, cb1, cb2)
{	if (cy1 == 0 && cy2 == 0 && cb1 == 0 && cb2 == 0)
		list (-1, numeric (), numeric () )
	else
	{	dx = cx2 - cx1
		cb1 = dx * cb1
	    	cb2 = dx * cb2
		p = .params (cy1, cy2, cb1, cb2)
		p1 = .params.derivative (cy1, cy2, cb1, cb2)
		p2 = .params.2nd (cy1, cy2, cb1, cb2)

		v = .quadratic.roots (p1, FALSE, FALSE)
		if (v [[1]] > 0)
		{	x = .bt (cx1, cx2, dx, v [[2]])
			y2 = p2 [1] + p2 [2] * v [[2]]
			if (v [[1]] == 1)
			{	y = sum (p * v [[2]] ^ (0:3) ) 
				if (y >= cy1 && y <= cy2)
					list (1, x, 0)
				else
					list (1, x, sign (y2) )
	
			}
			else
				list (2, x, sign (y2) )
		}
		else
			list (v [[1]], v [[2]], numeric () )
	}
}

#level (and *linear*) sections -> -1
.interval.argflexs = function (cx1, cx2, cy1, cy2, cb1, cb2)
{	dx = cx2 - cx1
	sec = (cy2 - cy1) / dx
	if (cy1 == 0 && cy2 == 0 && cb1 == 0 && cb2 == 0)
		list (-1, numeric (), numeric () )
	else if (sec == cb1 && sec == cb2)
		list (-1, numeric () )
	else
	{   	cb1 = dx * cb1
	    	cb2 = dx * cb2
		p2 = .params.2nd (cy1, cy2, cb1, cb2)

		v = .linear.root (p2, FALSE, FALSE)
		if (v [[1]] > 0)
			v [[2]] = .bt (cx1, cx2, dx, v [[2]])
		v
	}
}

.bt = function (cx1, cx2, dx, x)
{	y = cx1 + dx * x
	y [x == 0] = cx1
	y [x == 1] = cx2
	y
}
