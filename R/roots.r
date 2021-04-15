#kubik: Cubic Hermite Splines and Related Foot Finding Methods
#Copyright (C), Abby Spurdle, 2019 to 2021

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.root.chs = function (f, g, include.implied.roots, warning, ...)
{	. = attributes (g)
	if (is (g, "CHS") )
	{	f (.$cx, .$cy, .$cb,
			include.implied.roots=include.implied.roots,
			warning=warning,
			...)
	}
	else
		stop ("needs CHS object")
}

root.chs = function (sf, ..., include.implied.roots=TRUE, warning=TRUE)
	.root.chs (root.chs.eval, sf, include.implied.roots, warning)

argmin.chs = function (sf, ..., include.implied.roots=TRUE, warning=TRUE)
	.root.chs (argmin.chs.eval, sf, include.implied.roots, warning)

argmax.chs = function (sf, ..., include.implied.roots=TRUE, warning=TRUE)
	.root.chs (argmax.chs.eval, sf, include.implied.roots, warning)

argflex.chs = function (sf, ..., include.implied.roots=TRUE, warning=TRUE)
	.root.chs (argflex.chs.eval, sf, include.implied.roots, warning)

roots.chs = function (sf, ..., include.implied.roots=TRUE, warning=TRUE)
	.root.chs (roots.chs.eval, sf, include.implied.roots, warning)

argmins.chs = function (sf, ..., include.implied.roots=TRUE, warning=TRUE)
	.root.chs (argmins.chs.eval, sf, include.implied.roots, warning)

argmaxs.chs = function (sf, ..., include.implied.roots=TRUE, warning=TRUE)
	.root.chs (argmaxs.chs.eval, sf, include.implied.roots, warning)

argflexs.chs = function (sf, ..., include.implied.roots=TRUE, warning=TRUE)
	.root.chs (argflexs.chs.eval, sf, include.implied.roots, warning)

chs.roots.derivative = function (sf, ..., include.implied.roots=TRUE, warning=TRUE)
	.root.chs (chs.roots.derivative.eval, sf, include.implied.roots, warning)

solve.CHS = function (a, b=0, ..., to.list=FALSE)
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
