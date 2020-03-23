#kubik: Cubic Hermite Splines and Related Optimization Methods
#Copyright (C), Abby Spurdle, 2020

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.likely.2b.equally.spaced = function (cx)
{	dx = diff (cx)
	diff (range (dx) ) < min (dx) / 200
}

chs.constraints = function (correction=TRUE, ..., bounds, increasing=FALSE, decreasing=FALSE)
{	if (missing (bounds) )
	{	bounded = c (FALSE, FALSE)
		bounds = c (NA, NA)
	}
	else
	{	bounded = is.finite (bounds)
		if (diff (bounds) <= 0)
			stop ("diff (bounds) <= 0")
	}
	if (increasing && decreasing)
		stop ("increasing and decreasing constraints")
	v = list (correction=correction, bounded=bounded, bounds=bounds, increasing=increasing, decreasing=decreasing)
	.EXTEND (v, "chs.constraints")
}

chs.slopes = function (cx, cy, ..., constraints = chs.constraints (, ...), init.method)
{	nc = .test.cps (cx, cy)
	.chs.slopes (nc, cx, cy, constraints, init.method)
}

.chs.slopes = function (nc, cx, cy, constraints, init.method)
{	sec = diff (cy) / diff (cx)
	cb = .slopes.ext (nc, cx, cy, sec, init.method)
	.apply.ext (nc, cx, cy, sec, cb, constraints)
}

chs.tangents = function (cx, cy, ..., constraints = chs.constraints (, ...), init.method)
{	cb = chs.slopes (cx, cy, constraints, init.method)
	cbind (intercept = cy - cb * cx, slope=cb)
}

apply.chs.constraints = function (cx, cy, cb, ..., constraints = chs.constraints (, ...) )
{	nc = .test.cps (cx, cy)
	.apply.chs.constraints (nc, cx, cy, cb, constraints)
}

.apply.chs.constraints = function (nc, cx, cy, cb, constraints)
{	sec = diff (cy) / diff (cx)
	.apply.ext (nc, cx, cy, sec, cb, constraints)
}

.slopes.ext = function (nc, cx, cy, sec, init.method)
{	if (nc == 2)
		cb = c (sec, sec)
	else
	{	if (missing (init.method) )
		{	if (.likely.2b.equally.spaced (cx) )
				init.method = "SQ"
			else
				init.method = "SL"
		}
		cb = numeric (nc)
		for (i in 2:(nc - 1) )
			cb [i] = (cy [i + 1] - cy [i - 1]) / (cx [i + 1] - cx [i - 1])
		if (init.method == "SL")
		{	cb [1] = sec [1]
			cb [nc] = sec [nc - 1]
		}
		else if (init.method == "SQ")
		{	ps = .quad.params (nc, cx, cy)
			cb [1] = ps [1, 2] + 2 * ps [1, 3] * cx [1]
			cb [nc] = ps [2, 2] + 2 * ps [2, 3] * cx [nc]
		}
		else
			stop ('init.method needs to be "SL" or "SQ"')
	}
	cb
}

.apply.ext = function (nc, cx, cy, sec, cb, constraints)
{	if (is.null (constraints) || (is.vector (constraints) && length (constraints) == 1 && is.na (constraints) ) )
		NULL
	else if (inherits (constraints, "chs.constraints") )
	{	if (constraints$bounded [1] && any (cy < constraints$bounds [1]) )
			stop ("some cy < bounds [1]")
		if (constraints$bounded [2] && any (cy > constraints$bounds [2]) )
			stop ("some cy > bounds [2]")
		if (constraints$increasing && any (diff (cy) < 0) )
			stop ("increasing constraint needs non-decreasing cy")
		if (constraints$decreasing && any (diff (cy) > 0) )
			stop ("decreasing constraint needs non-increasing cy")

		ns = nc - 1
		mono = (constraints$increasing || constraints$decreasing)

		if (constraints$correction)
			cb = .apply.ext.2 (ns, sec, cb)
		if (any (constraints$bounded) || mono)
		{	I = rep (FALSE, ns)
			for (i in 1:ns)
			{	K = c (i, i + 1)
				x = .chs.roots.derivative.eval (cx [K], cy [K], cb [K], TRUE, FALSE, FALSE)
				x = x [x [,2] == 1 | x [,2] == -1, 1]
				if (length (x) > 0)
				{	if (mono)
						I [i] = TRUE
					else
					{	y = chs.eval (cx [K], cy [K], cb [K], x)
						I [i] = (
							(constraints$bounded [1] && any (y < constraints$bounds [1]) ) ||
							(constraints$bounded [2] && any (y > constraints$bounds [2]) ) )
					}
				}
			}
			I = which (I)
			cb [I] = cb [I + 1] = 0
		}
	}
	else
		stop ("constraints needs to be NULL, NA or chs.constraints object")
	cb
}

.apply.ext.2 = function (ns, sec, cb)
{	sign.sec = sign (sec)
	sign.cb = sign (cb)
	for (i in 1:(ns) )
	{	if (sec [i] == 0)
		{	if (sign.cb [i] + sign.cb [i + 1] != 0)
			{	cb [i] = 0;
				cb [i + 1] = 0;
			}
		}
		else if (sign.sec [i] + sign.cb [i] != 0 && sign.sec [i] + sign.cb [i + 1] != 0)
		{	if (cb [i] / sec [i] > 2.999999)
				cb [i] = 2.999999 * sec [i]
			if (cb [i + 1] / sec [i] > 2.999999)
				cb [i + 1] = 2.999999 * sec [i]
		}
	}
	cb
}
