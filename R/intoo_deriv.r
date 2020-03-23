#kubik: Cubic Hermite Splines and Related Optimization Methods
#Copyright (C), Abby Spurdle, 2020

#This program is distributed without any warranty.

#This program is free software.
#You can modify it and/or redistribute it, under the terms of:
#The GNU General Public License, version 2, or (at your option) any later version.

#You should have received a copy of this license, with R.
#Also, this license should be available at:
#https://cran.r-project.org/web/licenses/GPL-2

#adapted from intoo::EXTEND
.EXTEND = function (object, class, ...)
{	class = c (class, class (object) )
	structure (object, class=class, ...)
}

#adapted from intoo::THAT
.THAT = function ()
{	this = sys.function (-1)
	attributes (this)
}
