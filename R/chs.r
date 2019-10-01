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
		for (i in 1:nc)
			cb [i] = ps [i, 2] + 2 * ps [i, 3] * cx [i]
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

chs.bilinearize = function (cx, cy, cb, correction=TRUE, at,
	distance = factor * mean (diff (cx) ), factor=0.125)
{	nc = .test.cps (cx, cy, cb)
	at = sort (as.integer (at) )
	nrep = length (at)
	if (missing (cb) )
		cb = .chs.slopes (nc, cx, cy, correction)
	if (nrep == 0)
		stop ("needs one or more at values")
	if (any (diff (at) < 2) )
		stop ("needs unique non-adjacent at values")
	if (at [1] < 2 || at [nrep] > nc - 1)
		stop ("at needs to be in interval [2, nc -1]")
	if (any (distance <= 0) )
		stop ("distance <= 0")
	if (length (distance) == 1)
		distance = rep (distance, nrep)
	else if (nrep != length (distance) )
		stop ("length (distance) needs to be 1 or length (at)")
	nc2 = nc + nrep
	k = at + (0:(nrep - 1) )
	cx2 = cy2 = cb2 = rep (0, nc2)
	cx2 [- c (k, k + 1)] = cx [-at]
	cy2 [- c (k, k + 1)] = cy [-at]
	cb2 [- c (k, k + 1)] = cb [-at]
	for (i in 1:nrep)
	{	I = at [i]
		K1 = k [i]
		K2 = K1 + 1

		u1 = cx [I] - distance [i]
		u2 = cx [I] + distance [i]
		cx2 [K1] = u1
		cx2 [K2] = u2
		cy2 [K1] = cy [I - 1] + (u1 - cx [I - 1]) * cb [I - 1]
		cy2 [K2] = cy [I + 1] - (cx [I + 1] - u2) * cb [I + 1]
		cb2 [K1] = cb [I - 1]
		cb2 [K2] = cb [I + 1]
	}
	if (any (diff (cx2) <= 0) )
		stop ("\ndistance too large\n(needs unique ascending cx values, after replacement)")
	list (cx=cx2, cy=cy2, cb=cb2)
}
