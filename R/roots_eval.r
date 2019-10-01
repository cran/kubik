#kubik: Cubic Hermite Splines
#Copyright (C), Abby Spurdle, 2019

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

chs.root.eval = function (nc, cx, cy, cb, include.implied.roots=TRUE, warning=TRUE)
{	x = chs.roots.eval (nc, cx, cy, cb, include.implied.roots, warning)
	n = length (x)
	if (n == 0)
		stop ("no roots")
	else if (n == 1)
		x
	else if (n > 1)
		stop ("multiple roots")
}

chs.argmin.eval = function (nc, cx, cy, cb, include.implied.roots=TRUE, warning=TRUE)
{	x = chs.argmins.eval (nc, cx, cy, cb, include.implied.roots, warning)
	n = length (x)
	if (n == 0)
		stop ("no minima")
	else if (n == 1)
		x
	else
	{	y = chs.eval (nc, cx, cy, cb, x)
		I = which.min (y)
		nI = sum (y [I] == y)
		if (nI > 1)
			stop ("multiple global minima")
		x [I]
	}
}

chs.argmax.eval = function (nc, cx, cy, cb, include.implied.roots=TRUE, warning=TRUE)
{	x = chs.argmaxs.eval (nc, cx, cy, cb, include.implied.roots, warning)
	n = length (x)
	if (n == 0)
		stop ("no maxima")
	else if (n == 1)
		x
	else
	{	y = chs.eval (nc, cx, cy, cb, x)
		I = which.max (y)
		nI = sum (y [I] == y)
		if (nI > 1)
			stop ("multiple global maxima")
		x [I]
	}
}

chs.argflex.eval = function (nc, cx, cy, cb,
	include.implied.roots=TRUE, warning=TRUE, all.inflection.points=FALSE)
{	x = chs.argflexs.eval (nc, cx, cy, cb, include.implied.roots, warning, all.inflection.points)
	n = length (x)
	if (n == 0)
		stop ("no inflection points")
	else if (n == 1)
		x
	else if (n > 1)
		stop ("multiple inflection points")
}

chs.argmins.eval = function (nc, cx, cy, cb, include.implied.roots=TRUE, warning=TRUE)
{	r = .chs.roots.derivative.eval (nc, cx, cy, cb, include.implied.roots, warning)
	r [r [,2] == 1, 1]
}

chs.argmaxs.eval = function (nc, cx, cy, cb, include.implied.roots=TRUE, warning=TRUE)
{	r = .chs.roots.derivative.eval (nc, cx, cy, cb, include.implied.roots, warning)
	r [r [,2] == -1, 1]
}

chs.roots.derivative.eval = function (nc, cx, cy, cb, include.implied.roots=TRUE, warning=TRUE, all.inflection.points=FALSE)
	.chs.roots.derivative.eval (nc, cx, cy, cb, include.implied.roots, warning, !all.inflection.points)
