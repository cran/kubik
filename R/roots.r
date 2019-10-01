#kubik: Cubic Hermite Splines
#Copyright (C), Abby Spurdle, 2019

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

.chs.root = function (f, g, include.relevant.sections, warning, ...)
{	. = attributes (g)
	if (inherits (g, "chs") )
		f (.$nc, .$cx, .$cy, .$cb, include.relevant.sections, warning, ...)
	else
		stop ("needs chs object")
}

chs.root = function (f, include.implied.roots=TRUE, warning=TRUE)
	.chs.root (chs.root.eval, f, include.implied.roots, warning)

chs.argmin = function (f, include.implied.roots=TRUE, warning=TRUE)
	.chs.root (chs.argmin.eval, f, include.implied.roots, warning)

chs.argmax = function (f, include.implied.roots=TRUE, warning=TRUE)
	.chs.root (chs.argmax.eval, f, include.implied.roots, warning)

chs.argflex = function (f, include.implied.roots=TRUE, warning=TRUE, all.inflection.points=FALSE)
	.chs.root (chs.argflex.eval, f, include.implied.roots, warning, all.inflection.points)

chs.roots = function (f, include.implied.roots=TRUE, warning=TRUE)
	.chs.root (chs.roots.eval, f, include.implied.roots, warning)

chs.argmins = function (f, include.implied.roots=TRUE, warning=TRUE)
	.chs.root (chs.argmins.eval, f, include.implied.roots, warning)

chs.argmaxs = function (f, include.implied.roots=TRUE, warning=TRUE)
	.chs.root (chs.argmaxs.eval, f, include.implied.roots, warning)

chs.argflexs = function (f, include.implied.roots=TRUE, warning=TRUE, all.inflection.points=FALSE)
	.chs.root (chs.argflexs.eval, f, include.implied.roots, warning, all.inflection.points)

chs.roots.derivative = function (f, include.implied.roots=TRUE, warning=TRUE, all.inflection.points=FALSE)
	.chs.root (chs.roots.derivative.eval, f, include.implied.roots, warning, all.inflection.points)

solve.chs = function (a, b, to.list=FALSE, ...)
{	. = attributes (a)
	n = length (b)
	if (n == 1 && ! to.list)
		chs.roots.eval (.$nc, .$cx, .$cy - b, .$cb, ...)
	else
	{	x = vector ("list", n)
		for (i in seq_len (n) )
			x [[i]] = chs.roots.eval (.$nc, .$cx, .$cy - b [i], .$cb, ...)
		x
	}
}
