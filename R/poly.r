#kubik: Cubic Hermite Splines and Related Foot Finding Methods
#Copyright (C), Abby Spurdle, 2019 to 2021

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

.roots.2 = function (x)
{	x [x > -0.000001 & x < 0] = 0
	x [x > 1 & x < 1.000001] = 1
	x = x [x >= 0 & x <= 1]
	if (length (x) != 1)
		stop ("root finding algorithm failed")
	list (1, x)
}

.sympoly = function (p)
{	text = paste ("function (x)", p [1], "+", p [2], "* x +", p [3], "* x ^ 2 +", p [4], "* x  ^ 3")
	eval (str2lang (text) )
}

.linear.root = function (p, include.lower, include.upper)
{	ap = abs (p)
	tol = 1e-10

	if (ap [1] < tol && ap [2] < tol)
		list (-1, numeric () )
	else if (ap [2] < tol)
		list (0, numeric () )
	else
		.roots (.linear.root.ext (p), include.lower, include.upper)
}

.quadratic.roots = function (p, include.lower, include.upper, root.expected=FALSE)
{	ap = abs (p)
	tol = 1e-10

	if (ap [1] < tol && ap [2] < tol && ap [3] < tol)
		list (-1, numeric () )
	else if (ap [2] < tol && ap [3] < tol)
		list (0, numeric () )
	else
	{	if (ap [3] < tol)
			x = .linear.root.ext (p)
		else
			x = .quadratic.roots.ext (p)
		if (root.expected)
			.roots.2 (x)
		else
			.roots (x, include.lower, include.upper)
	}
}

.cubic.roots = function (p, cy1, cy2, include.lower, include.upper)
{	ap = abs (p)
	tol = 1e-10

	if (ap [1] < tol && ap [2] < tol && ap [3] < tol && ap [4] < tol)
		list (-1, numeric ())
	else if (ap [2] < tol && ap [3] < tol && ap [4] < tol)
		list (0, numeric () )
	else
	{	if (ap [3] < tol && ap [4] < tol)
			x = .linear.root.ext (p)
		else if (ap [4] < tol)
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
	r1 = .quadratic.roots.ext (p1)
	r1 = r1 [r1 > 0 & r1 < 1]
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

.quad.params = function (nc, cx, cy)
{	if (nc == 3)
	{	p = solve (cbind (1, cx, cx ^ 2), cy)
		rbind (p, p)
	}
	else
	{	Ia = 1:3
		Ib = (nc - 2):nc
		p1 = solve (cbind (1, cx [Ia], cx [Ia] ^ 2), cy [Ia])
		p2 = solve (cbind (1, cx [Ib], cx [Ib] ^ 2), cy [Ib])
		rbind (p1, p2)
	}
}
