#kubik: Cubic Hermite Splines
#Copyright (C), Abby Spurdle, 2019

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.roots = function (x, include.lower, include.upper)
{	x = .within (x, include.lower, include.upper)
	list (length (x), x)
}

.sympoly = function (p)
{	text = paste ("function (x)", p [1], "+", p [2], "* x +", p [3], "* x ^ 2 +", p [4], "* x  ^ 3")
	eval (parse (text=text) )
}

.linear.root = function (p, include.lower, include.upper)
{	if (p [1] == 0 && p [2] == 0)
		list (-1, numeric () )
	else if (p [2] == 0)
		list (0, numeric () )
	else
		.roots (.linear.root.ext (p), include.lower, include.upper)
}

.quadratic.roots = function (p, include.lower, include.upper)
{	if (p [1] == 0 && p [2] == 0 && p [3] == 0)
		list (-1, numeric () )
	else if (p [2] == 0 && p [3] == 0)
		list (0, numeric () )
	else
	{	if (p [3] == 0)
			x = .linear.root.ext (p)
		else
			x = .quadratic.roots.ext (p)
		.roots (x, include.lower, include.upper)
	}
}

.cubic.roots = function (p, cy1, cy2, include.lower, include.upper)
{	if (p [1] == 0 && p [2] == 0 && p [3] == 0 && p [4] == 0)
		list (-1, numeric ())
	else if (p [2] == 0 && p [3] == 0 && p [4] == 0)
		list (0, numeric () )
	else
	{	if (p [3] == 0 && p [4] == 0)
			x = .linear.root.ext (p)
		else if (p [4] == 0)
			x = .quadratic.roots.ext (p)
		else
			x = .cubic.roots.ext (p, cy1, cy2)
		.roots (x, include.lower, include.upper)
	}
}

.linear.root.ext = function (p)
	-p [1] / p [2]

.quadratic.roots.ext = function (p)
{	delta = p [2] ^ 2 - 4 * p [1] * p [3]
	if (delta == 0)
		p [2] / (-2 * p [3])
	else if (delta > 0)
	{	x1 = (-p [2] + sqrt (delta) ) / (2 * p [3])
		x2 = (-p [2] - sqrt (delta) ) / (2 * p [3])
		sort (c (x1, x2) )
	}
	else
		numeric ()
}

.cubic.roots.ext = function (p, cy1, cy2)
{	f = .sympoly (p)
	p1 = (1:3) * p [2:4]
	r1 = .within (.quadratic.roots.ext (p1), FALSE, FALSE)
	nr1 = length (r1)
	y = c (cy1, f (r1), cy2)
	if (nr1 == 0)
		.uniroot.2 (f, c (0, 1), y, TRUE)
	else if (nr1 == 1)
	{	x1 = .uniroot.2 (f, c (0, r1), y [1:2], FALSE)
		x2 = .uniroot.2 (f, c (r1, 1), y [2:3], TRUE)
		c (x1, x2)
	}
	else
	{	x1 = .uniroot.2 (f, c (0, r1 [1]), y [1:2], FALSE)
		x2 = .uniroot.2 (f, r1, y [2:3], FALSE)
		x3 = .uniroot.2 (f, c (r1 [2], 1), y [3:4], TRUE)
		c (x1, x2, x3)
	}
}

.within = function (x, include.lower, include.upper)
{	if (length (x) == 0)
		numeric ()
	else
	{	x [x > -0.000001 & x < 0.000001] = 0
		x [x > 0.999999 & x < 1.000001] = 1
		if (include.lower)
		{	if (include.upper)
				x [x >= 0 & x <= 1]
			else
				x [x >= 0 & x < 1]
		}
		else
		{	if (include.upper)
				x [x > 0 & x <= 1]
			else
				x [x > 0 & x < 1]
		}
	}
}

.uniroot.2 = function (f, x, y, include.upper)
{	if (y [1] == 0)
		x [1]
	else if (y [2] == 0)
	{	if (include.upper)
			x [2]
		else
			numeric ()
	}
	else
	{	if (sum (sign (y) ) == 0)
			uniroot (f, x)$root
		else
			numeric ()
	}
}

.quad.params = function (cx, cy)
	solve (cbind (1, cx, cx ^ 2), cy)

.quad.params.2 = function (cx, cy, cxi, cbi)
{	r1 = cbind (1, cx, cx ^ 2)
	r2 = c (0, 1, 2 * cxi)
	solve (rbind (r1, r2), c (cy, cbi) )
}

.quads.params = function (nc, cx, cy)
{	p = matrix (0, nc, 3)
	for (i in 2:(nc - 1) )
	{	I = (i - 1):(i + 1)
		p [i,] = .quad.params (cx [I], cy [I])
	}
	p [1,] = p [2,]
	p [nc,] = p [nc - 1,]
	p
}
