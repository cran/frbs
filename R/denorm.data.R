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
#' This function is to transform from normalized data into real-valued data. 
#'
#' @title The data de-normalization
#' @param dt.norm a matrix(n x m) of the normalized data.
#' @param range.data a matrix(2 x n) containing the range of the data, where n is the number of variables, and
#' first and second rows are the minimum and maximum value, respectively. 
#' @param min.scale the minimum value within normalization.
#' @param max.scale the maximum value within normalization.
#' @seealso \code{\link{norm.data}}
# @return the real-valued data
# @export
denorm.data <- function(dt.norm, range.data, min.scale = -1, max.scale = 1){
	row.data <- nrow(dt.norm)
	col.data <- ncol(dt.norm)
	data.denorm <- matrix(nrow = row.data, ncol=col.data)
	
	# denormalize all data on each column 
	for (j in 1:col.data){
		min.data <- range.data[1, j]
		max.data <- range.data[2, j]

		for (i in 1:row.data){
			data.denorm[i, j] <- min.data + ((dt.norm[i, j] - min.scale)*(max.data - min.data))/ (max.scale - min.scale)
		}
	}
	
return(data.denorm)
}