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
#' The purpose of this function is to generate data, which contains 
#' two input variables and one output variable, automatically for all values on a plane. 
#'
#' @title A data generator
#' @param range.input the range of the input variables, as a matrix(2 x n). 
#' @param num.grid a number representing the size of the grid on the plane.
#' @return the data
#' @examples 
#' range.input <- matrix(c(0, 100, 0, 100), nrow=2)
#' num.grid <- 10
#' data.test <- data.gen3d(range.input, num.grid)
#' @export
data.gen3d <- function(range.input, num.grid = 10){

num.multiple <- num.grid
length.data <- num.multiple ^ 2
delta.x <- matrix(nrow = ncol(range.input), ncol = 1)
data <- matrix(nrow = length.data, ncol = ncol(range.input))

for (i in 1 : ncol(range.input)){
	delta.x[i] <- (max(range.input[, i]) - min(range.input[, i])) / num.multiple
}

for (i in 1 : ncol(range.input)){
	data[1, i] <- min(range.input[, i])
}

for (i in 2 : length.data){	
	for (j in 1 : ncol(range.input)){
		if (i %% num.multiple == 1){
			data[i, 1] <- data[i - 1, 1] + delta.x[j, ]
			data[i, 2] <- data[1, 2]
		}
		else {
		data[i, 1] <- data[i - 1, 1]
		data[i, 2] <- data[i - 1, 2]  + delta.x[j, ] 		
		}
	}	
	
}


return(data)
}