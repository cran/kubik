#kubik: Cubic Hermite Splines and Related Optimization Methods
#Copyright (C), Abby Spurdle, 2020

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.test.cps = function (cx, cy, cb)
{	nx = length (cx)
	ny = length (cy)
	if (nx < 2)
		stop ("needs 2 or more control points")
	if (! (all (is.finite (cx) ) && all (is.finite (cy) ) ) )
		stop ("cx, cy and cb need to be finite")
	else
	{	if (nx != length (unique (cx) ) )
			stop ("needs unique cx values")
		if (any (diff (cx) < 0) )
			stop ("needs ascending cx values")
	}
	if (missing (cb) )
		nb = nx
	else
	{	nb = length (cb)
		if (! (all (is.finite (cb) ) ) )
			stop ("cx, cy and cb need to be finite")
	}
	if (nx != ny || nx != nb)
		stop ("length of cx, cy and cb need to be equal")
	nx
}

chs = function (cx, cy, cb, ..., constraints = chs.constraints (, ...), transform=FALSE, outside = c (NA, NA), init.method)
{	nc = .test.cps (cx, cy, cb)
	f = function (x)
	{	. = .THAT ()
		with (.,
			chs.eval (cx, cy, cb, x, outside=outside) )
	}
	cb = .chs.slopes.trans (nc, cx, cy, cb, constraints, transform, init.method)
	.EXTEND (f, c ("chs", "kspline"),
		nc=nc, cx=cx, cy=cy, cb=cb, outside=outside)
}

chs.derivative = function (cx, cy, cb, ..., constraints = chs.constraints (, ...), transform=FALSE, outside = c (NA, NA), init.method)
{	nc = .test.cps (cx, cy, cb)
	f = function (x)
	{	. = .THAT ()
		with (.,
			chs.derivative.eval (cx, cy, cb, x, outside=outside) )
	}
	cb = .chs.slopes.trans (nc, cx, cy, cb, constraints, transform, init.method)
	.EXTEND (f, c ("chs.derivative", "kspline"),
		nc=nc, cx=cx, cy=cy, cb=cb, outside=outside)
}

chs.integral = function (cx, cy, cb, ..., constraints = chs.constraints (, ...), transform=FALSE, outside = c (NA, NA), init.method, constant=0)
{	nc = .test.cps (cx, cy, cb)
	f = function (x)
	{	. = .THAT ()
		with (.,
			chs.integral.eval (cx, cy, cb, x, outside=outside, constant=constant) )
	}
	cb = .chs.slopes.trans (nc, cx, cy, cb, constraints, transform, init.method)
	.EXTEND (f, c ("chs.integral", "kspline"),
		nc=nc, cx=cx, cy=cy, cb=cb, outside=outside, constant=constant)
}

approx.chs.derivative = function (cx, cy, cb, ...,
	constraints = chs.constraints (, ...), transform=FALSE, outside = c (NA, NA), init.method,
	apply.constraints.to=0, nth=1, trim=TRUE)
{	nc = .test.cps (cx, cy, cb)
	f = function (x)
	{	. = .THAT ()
		with (.,
			chs.eval (cx, cy, cb, x, outside=outside) )
	}
	cdefault = chs.constraints ()
	cs = if (apply.constraints.to == 0) constraints else cdefault
	cb = .chs.slopes.trans (nc, cx, cy, cb, cs, transform, init.method)
	for (i in seq_len (nth) )
	{	cs = if (apply.constraints.to == i) constraints else cdefault
		cy = cb
		cb = .chs.slopes (nc, cx, cy, NULL, init.method)
		if (trim)
		{	nc = nc - 2
			I = (2):(nc + 1)
			cx = cx [I]
			cy = cy [I]
			cb = cb [I]
		}
		cb = .apply.chs.constraints (nc, cx, cy, cb, cs)
	}
	.EXTEND (f, c ("approx.chs.derivative", "chs", "kspline"),
		nc=nc, cx=cx, cy=cy, cb=cb, outside=outside)
}

.chs.slopes.trans = function (nc, cx, cy, cb, constraints, transform, init.method)
{	if (missing (cb) )
	{	if (transform)
			stop ("slopes needed if transform is true")
		cb = .chs.slopes (nc, cx, cy, constraints, init.method)
	}
	else if (transform)
		cb = .apply.chs.constraints (nc, cx, cy, cb, constraints)
	cb
}
