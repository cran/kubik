#kubik: Cubic Hermite Splines and Related Foot Finding Methods
#Copyright (C), Abby Spurdle, 2019 to 2021

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

roots.chs.eval = function (cx, cy, cb, ..., include.implied.roots=TRUE, warning=TRUE)
{	nc = length (cx)

	ns = nc - 1
	Rc = numeric (ns)
	xs = vector ("list", ns)
	v = .interval.roots (cx [1], cx [2], cy [1], cy [2], cb [1], cb [2], TRUE, (ns == 1) )
	Rc [1] = v [[1]]
	xs [[1]] = v [[2]]
	if (ns > 1)
	{	for (i in 2:ns)
		{	include.lwr = (Rc [i - 1] >= 0)
			v = .interval.roots (cx [i], cx [i + 1], cy [i], cy [i + 1], cb [i], cb [i + 1], include.lwr, (i == ns) )
			Rc [i] = v [[1]]
			xs [[i]] = v [[2]]
		}
	}
	levs = (Rc == -1)
	if (any (levs) )
	{	if (include.implied.roots)
		{	m = .midpoints (nc, cx, levs)
			x = sort (c (m, unlist (xs) ) )
		}
		else
			x = unlist (xs)
		if (warning)
			warning ("roots coincide with level sections")
		x
	}
	else
		unlist (xs)
}

.chs.roots.derivative.eval = function (cx, cy, cb,
	include.implied.roots=TRUE, warning=TRUE, remove=FALSE)
{	nc = length (cx)

	I = rep (FALSE, nc)
	S = rep (0, nc)
	ns = nc - 1
	Rc = numeric (ns)
	xs = ss = vector ("list", ns)
	for (i in 1:ns)
	{	v = .interval.roots.derivative (cx [i], cx [i + 1], cy [i], cy [i + 1], cb [i], cb [i + 1])
		Rc [i] = v [[1]]
		xs [[i]] = v [[2]]
		ss [[i]] = v [[3]]
	}
	if (ns > 1)
	{	for (i in 2:(nc - 1) )
		{	if (cb [i] == 0 && Rc [i - 1] != -1 && Rc [i] != -1)
			{	I [i] = TRUE
				S [i] = .implied.sign (cx, cy, cb, i - 1, i)
			}
		}
	}
	if (remove)
	{	quads = .is.quadratic.suggestive (nc, cx, cy, cb)
		for (i in which (quads) )
		{	is.opt = (ss [[i]] != 0)
			xs [[i]] = xs [[i]][is.opt]
			ss [[i]] = ss [[i]][is.opt]
		}
	}
	levs = (Rc == -1)
	if (ns > 2 && any (levs) )
	{	if (include.implied.roots)
		{	v = .midpoints.1 (nc, cx, cy, cb, levs)
			r = .sort.2 (c (cx [I], v [[1]], unlist (xs) ), c (S [I], v [[2]], unlist (ss) ) )
		}
		else
			r = .sort.2 (c (cx [I], unlist (xs) ), c (S [I], unlist (ss) ) )
		if (warning)
			warning ("spline contains level sections")
		r
	}
	else
		.sort.2 (c (cx [I], unlist (xs) ), c (S [I], unlist (ss) ) )
}

argflexs.chs.eval = function (cx, cy, cb, ...,
	include.implied.roots=TRUE, warning=TRUE)
{	nc = length (cx)

	I = rep (FALSE, nc)
	ns = nc - 1
	Rc = numeric (ns)
	xs = vector ("list", ns)
	for (i in 1:ns)
	{	v = .interval.argflexs (cx [i], cx [i + 1], cy [i], cy [i + 1], cb [i], cb [i + 1])
		Rc [i] = v [[1]]
		xs [[i]] = v [[2]]
	}
	if (ns > 1)
	{	for (i in 2:(nc - 1) )
		{	if (Rc [i - 1] != -1 && Rc [i] != -1)
				I [i] = (.implied.sign (cx, cy, cb, i - 1, i) == 0)
		}
	}
	xs [.is.quadratic.suggestive (nc, cx, cy, cb)] = NULL
	levs = (Rc == -1)
	if (ns > 2 && any (levs) )
	{	if (include.implied.roots)
		{	m = .midpoints.2 (nc, cx, cy, cb, levs)
			x = sort (c (cx [I], m, unlist (xs) ) )
		}
		else
			x = sort (c (cx [I], unlist (xs) ) )
		if (warning)
			warning ("spline contains linear sections")
		x
	}
	else
		sort (c (cx [I], unlist (xs) ) )
}

.midpoints = function (nc, cx, levs)
{	ls = .levset (nc - 1, levs)
	nls = nrow (ls)
	mid.xs = numeric (nls)
	for (i in 1:nls)
		mid.xs [i] = (cx [ls [i, 1] ] + cx [ls [i, 2] + 1]) / 2
	mid.xs
}

.midpoints.1 = function (nc, cx, cy, cb, levs)
{	ls = .levset (nc - 1, levs, TRUE)
	nls = nrow (ls)
	mid.xs = mid.ss = numeric (nls)
	for (i in seq_len (nls) )
	{	mid.xs [i] = (cx [ls [i, 1] ] + cx [ls [i, 2] + 1]) / 2
		mid.ss [i] = .implied.sign (cx, cy, cb, ls [i, 1] - 1, ls [i, 2] + 1)
	}
	list (mid.xs, mid.ss)
}

.midpoints.2 = function (nc, cx, cy, cb, levs)
{	ls = .levset (nc - 1, levs, TRUE)
	nls = nrow (ls)
	mid.xs = vector ("list", nls)
	for (i in seq_len (nls) )
	{	mid.ss = .implied.sign (cx, cy, cb, ls [i, 1] - 1, ls [i, 2] + 1)
		if (mid.ss == 0)
			mid.xs [[i]] = (cx [ls [i, 1] ] + cx [ls [i, 2] + 1]) / 2
	}
	unlist (mid.xs)
}

.levset = function (ns, levs, remove.outermost.sections=FALSE)
{	k1 = diff (c (0, levs) )
	k2 = diff (c (levs, 0) )
	ls = cbind (which (k1 == 1), which (k2 == -1) )
	if (remove.outermost.sections)
	{	nls = nrow (ls)
		I = rep (TRUE, nls)
		if (ls [1, 1] == 1)
			I [1] = FALSE
		if (ls [nls, 2] == ns)
			I [nls] = FALSE
		ls = ls [I,,drop=FALSE]
	}
	ls
}

.is.quadratic.suggestive = function (nc, cx, cy, cb)
{	sec = diff (cy) / diff (cx)
	lwr.rel = sign (cb [-nc] - sec)
	upr.rel = sign (cb [-1] - sec)
	(lwr.rel == 0 | upr.rel == 0 | lwr.rel != upr.rel)
}

.implied.sign = function (cx, cy, cb, a, b)
{	dxa = cx [a + 1] - cx [a]
	dxb = cx [b + 1] - cx [b]

	p2a = .params.2nd (cy [a], cy [a + 1], dxa * cb [a], dxa * cb [a + 1])
	p2b = .params.2nd (cy [b], cy [b + 1], dxb * cb [b], dxb * cb [b + 1])
	p3a = .params.3rd (cy [a], cy [a + 1], dxa * cb [a], dxa * cb [a + 1])
	p3b = .params.3rd (cy [b], cy [b + 1], dxb * cb [b], dxb * cb [b + 1])

	lwr.sign = sign (p2a [1] + p2a [2])
	upr.sign = sign (p2b [1])
	if (lwr.sign == 0)
		lwr.sign = (-1) * sign (p3a)
	if (upr.sign == 0)
		upr.sign = sign (p3b)

	if (lwr.sign == upr.sign)
		lwr.sign
	else
		0
}

.sort.2 = function (x, s)
{	I = order (x)
	cbind (x [I], s [I])
}
