#############################################################################
#
#  This file is a part of the R package "frbs".
#
#  Author: Lala Septem Riza
#  Co-author: Christoph Bergmeir
#  Supervisors: Francisco Herrera Triguero and Jose Manuel Benitez
#  Copyright (c) DiCITS Lab, Sci2s group, DECSAI, University of Granada.
#
#  This package is free software: you can redistribute it and/or modify it under
#  the terms of the GNU General Public License as published by the Free Software
#  Foundation, either version 2 of the License, or (at your option) any later version.
#
#  This package is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
#  A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
#############################################################################
#' This function is to transform from real-valued data into normalized data. 
#'
#' @title The data normalization
#' @param dt.ori a matrix(n x m) of the original data.
#' @param range.data a matrix(2 x n) containing the range of the data, where n is the number of variables, and
#' first and second rows are the minimum and maximum value, respectively. 
#' @param min.scale the minimum value within normalization.
#' @param max.scale the maximum value within normalization.
#' @seealso \code{\link{denorm.data}}
# @return the normalized data
# @export
norm.data <- function(dt.ori, range.data, min.scale = -1, max.scale = 1){
	row.data <- nrow(dt.ori)
	col.data <- ncol(dt.ori)
	data.norm <- matrix(nrow = row.data, ncol=col.data)
	
	# normalize all data on each column 
	for (j in 1:col.data){
		min.data <- range.data[1, j]
		max.data <- range.data[2, j]

		for (i in 1:row.data){
			data.norm[i, j] <- min.scale + (dt.ori[i, j] - min.data) * (max.scale - min.scale) / (max.data - min.data)
		}
	}
	
return(data.norm)
}