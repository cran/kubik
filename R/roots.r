#kubik: Cubic Hermite Splines and Related Optimization Methods
#Copyright (C), Abby Spurdle, 2020

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.root.chs = function (f, g, include.implied.roots, warning, ...)
{	. = attributes (g)
	if (inherits (g, "chs") )
	{	f (.$cx, .$cy, .$cb,
			include.implied.roots=include.implied.roots,
			warning=warning,
			...)
	}
	else
		stop ("needs chs object")
}

root.chs = function (f, ..., include.implied.roots=TRUE, warning=TRUE)
	.root.chs (root.chs.eval, f, include.implied.roots, warning)

argmin.chs = function (f, ..., include.implied.roots=TRUE, warning=TRUE)
	.root.chs (argmin.chs.eval, f, include.implied.roots, warning)

argmax.chs = function (f, ..., include.implied.roots=TRUE, warning=TRUE)
	.root.chs (argmax.chs.eval, f, include.implied.roots, warning)

argflex.chs = function (f, ..., include.implied.roots=TRUE, warning=TRUE)
	.root.chs (argflex.chs.eval, f, include.implied.roots, warning)

roots.chs = function (f, ..., include.implied.roots=TRUE, warning=TRUE)
	.root.chs (roots.chs.eval, f, include.implied.roots, warning)

argmins.chs = function (f, ..., include.implied.roots=TRUE, warning=TRUE)
	.root.chs (argmins.chs.eval, f, include.implied.roots, warning)

argmaxs.chs = function (f, ..., include.implied.roots=TRUE, warning=TRUE)
	.root.chs (argmaxs.chs.eval, f, include.implied.roots, warning)

argflexs.chs = function (f, ..., include.implied.roots=TRUE, warning=TRUE)
	.root.chs (argflexs.chs.eval, f, include.implied.roots, warning)

chs.roots.derivative = function (f, ..., include.implied.roots=TRUE, warning=TRUE)
	.root.chs (chs.roots.derivative.eval, f, include.implied.roots, warning)

solve.chs = function (a, b, ..., to.list=FALSE)
{	. = attributes (a)
	n = length (b)
	if (n == 1 && ! to.list)
		roots.chs.eval (.$cx, .$cy - b, .$cb, ...)
	else
	{	x = vector ("list", n)
		for (i in seq_len (n) )
			x [[i]] = roots.chs.eval (.$cx, .$cy - b [i], .$cb, ...)
		x
	}
}

#deprecated
chs.argmins = function (...) argmins.chs (...)
chs.argmaxs = function (...) argmaxs.chs (...)
