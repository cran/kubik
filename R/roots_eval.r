#kubik: Cubic Hermite Splines and Related Optimization Methods
#Copyright (C), Abby Spurdle, 2020

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

root.chs.eval = function (cx, cy, cb, ..., include.implied.roots=TRUE, warning=TRUE)
{	x = roots.chs.eval (cx, cy, cb,
		include.implied.roots=include.implied.roots,
		warning=warning)
	n = length (x)
	if (n == 0)
		stop ("no roots")
	else if (n == 1)
		x
	else if (n > 1)
		stop ("multiple roots")
}

argmin.chs.eval = function (cx, cy, cb, ..., include.implied.roots=TRUE, warning=TRUE)
{	x = argmins.chs.eval (cx, cy, cb,
		include.implied.roots=include.implied.roots,
		warning=warning)
	n = length (x)
	if (n == 0)
		stop ("no minima")
	else if (n == 1)
		x
	else
	{	y = chs.eval (cx, cy, cb, x)
		I = which.min (y)
		nI = sum (y [I] == y)
		if (nI > 1)
			stop ("multiple global minima")
		x [I]
	}
}

argmax.chs.eval = function (cx, cy, cb, ..., include.implied.roots=TRUE, warning=TRUE)
{	x = argmaxs.chs.eval (cx, cy, cb,
		include.implied.roots=include.implied.roots,
		warning=warning)
	n = length (x)
	if (n == 0)
		stop ("no maxima")
	else if (n == 1)
		x
	else
	{	y = chs.eval (cx, cy, cb, x)
		I = which.max (y)
		nI = sum (y [I] == y)
		if (nI > 1)
			stop ("multiple global maxima")
		x [I]
	}
}

argflex.chs.eval = function (cx, cy, cb, ..., 
	include.implied.roots=TRUE, warning=TRUE)
{	x = argflexs.chs.eval (cx, cy, cb,
		include.implied.roots=include.implied.roots,
		warning=warning)
	n = length (x)
	if (n == 0)
		stop ("no inflection points")
	else if (n == 1)
		x
	else if (n > 1)
		stop ("multiple inflection points")
}

argmins.chs.eval = function (cx, cy, cb, ..., include.implied.roots=TRUE, warning=TRUE)
{	r = .chs.roots.derivative.eval (cx, cy, cb, include.implied.roots, warning)
	r [r [,2] == 1, 1]
}

argmaxs.chs.eval = function (cx, cy, cb, ..., include.implied.roots=TRUE, warning=TRUE)
{	r = .chs.roots.derivative.eval (cx, cy, cb, include.implied.roots, warning)
	r [r [,2] == -1, 1]
}

chs.roots.derivative.eval = function (cx, cy, cb, ..., include.implied.roots=TRUE, warning=TRUE)
{	.chs.roots.derivative.eval (cx, cy, cb,
		include.implied.roots, warning)
}
