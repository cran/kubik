#kubik: Cubic Hermite Splines
#Copyright (C), Abby Spurdle, 2019

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

.chs.slopes = function (nc, cx, cy, correction=TRUE)
{	sec = diff (cy) / diff (cx)
	if (nc == 2)
		cb = c (sec, sec)
	else
	{	ps = .quads.params (nc, cx, cy)
		cb = numeric (nc)
		for (i in 2:(nc - 1) )
			cb [i] = ps [i, 2] + 2 * ps [i, 3] * cx [i]
		cb [1] = sec [1]
		cb [nc] = sec [nc - 1]
		if (correction)
			cb = .chs.slopes.ext (nc, sec, cb)
	}
	cb
}

.chs.slopes.ext = function (nc, sec, cb)
{	sign.sec = sign (sec)
	sign.cb = sign (cb)
	for (i in 1:(nc - 1) )
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

chs = function (cx, cy, cb, correction=TRUE, outside = c (NA, NA) )
{	nc = .test.cps (cx, cy, cb)
	f = function (x)
	{	. = THAT ()
		chs.eval (.$nc, .$cx, .$cy, .$cb, x, .$outside)
	}
	if (missing (cb) )
		cb = .chs.slopes (nc, cx, cy, correction)
	EXTEND (f, c ("chs", "k.spline"), nc, cx, cy, cb, outside)
}

chs.derivative = function (cx, cy, cb, correction=TRUE, outside = c (NA, NA) )
{	nc = .test.cps (cx, cy, cb)
	f = function (x)
	{	. = THAT ()
		chs.derivative.eval (.$nc, .$cx, .$cy, .$cb, x, .$outside)
	}
	if (missing (cb) )
		cb = .chs.slopes (nc, cx, cy, correction)
	EXTEND (f, c ("chs.derivative", "k.spline"), nc, cx, cy, cb, outside)
}

chs.integral = function (cx, cy, cb, correction=TRUE, outside = c (NA, NA), constant=0)
{	nc = .test.cps (cx, cy, cb)
	f = function (x)
	{	. = THAT ()
		chs.integral.eval (.$nc, .$cx, .$cy, .$cb, x, .$outside, .$constant)
	}
	if (missing (cb) )
		cb = .chs.slopes (nc, cx, cy, correction)
	EXTEND (f, c ("chs.integral", "k.spline"), nc, cx, cy, cb, outside, constant)
}

approx.chs.derivative = function (cx, cy, cb, correction=TRUE, outside = c (NA, NA), nth=1)
{	nc = .test.cps (cx, cy, cb)
	f = function (x)
	{	. = THAT ()
		chs.eval (.$nc, .$cx, .$cy, .$cb, x, .$outside)
	}
	if (missing (cb) )
		cb = .chs.slopes (nc, cx, cy, correction)
	for (i in 1:nth)
	{	cy = cb
		cb = .chs.slopes (nc, cx, cb, correction)
	}
	EXTEND (f, c ("approx.chs.derivative", "chs", "k.spline"), nc, cx, cy, cb, outside)
}

chs.tangents = function (cx, cy, correction=TRUE)
{	cb = chs.slopes (cx, cy, correction)
	cbind (intercept = cy - cb * cx, cb)
}

chs.slopes = function (cx, cy, correction=TRUE)
{	nc = .test.cps (cx, cy)
	.chs.slopes (nc, cx, cy, correction)
}
