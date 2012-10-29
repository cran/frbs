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
#' This function is the internal function of the MSGFS method to compute the predicted values.  
#'
#' @title MSGFS: The prediction phase
#' @param object the \code{\link{frbs-object}}.
#' @param newdata a matrix(m x n) of data for the prediction process, where m is the number of instances and 
#' n is the number of input variables.
#' @return A matrix of predicted values.
# @export
MSGFS.test <- function(object, newdata){
	mod <- object
	data.test <- newdata
	rule.gen <- mod$rule	
	res <- calc.pred.val(data.test, rule.gen)
	
	return(res)
}

# This function is the internal function of the MSGFS method to compute the predicted values.  
#
# @title MSGFS: The testing phase
# @param data.test a matrix(m x n) of data for the testing process, where m is the number of instances and 
# n is the number of input variables.
# @param rule.gen fuzzy IF-THEN rules
# @return A matrix of predicted values.
# @export
calc.pred.val <- function(data.test, rule.gen){
	output.val.pred <- matrix(nrow = nrow(data.test), ncol = 1)
	for (i in 1 : nrow(data.test)){
		deg.R <- matrix()
		deg.R <- ch.cover(matrix(data.test[i, ], nrow = 1), rule.gen)

		##defuzzification
		##get mean value of output variable on rule.gen
		output.val.mean <- matrix(rule.gen[, (ncol(rule.gen) - 1)], nrow = nrow(rule.gen))
		
		## calculate average values
		if (sum(deg.R) != 0){
			output.val.pred[i, 1] <- (t(deg.R) %*% output.val.mean) / sum(deg.R)		
		}
		else {
			output.val.pred[i, 1] <- 0.5
		}
	}
	
return(output.val.pred)
}