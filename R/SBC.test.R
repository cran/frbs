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
#' This function is the internal function of the SBC method to compute the predicted values.
#' 
#' @title SBC prediction phase
#' @param object the \code{\link{frbs-object}}. 
#' @param newdata a matrix(m x n) of data for the prediction process, where m is the number of instances and 
#' n is the number of input variables.
#' @seealso \code{\link{SBC}}
# @return A matrix of predicted value
# @export
SBC.test <- function(object, newdata){
	mod <- object
	data.test <- newdata
	cls <- mod$cls
	range.data.ori <- mod$range.data.ori
	r.a <- mod$r.a
	num.outvar <- 1
		
	cls <- norm.data(cls, range.data.ori, min.scale = -1, max.scale = 1)
	
	alpha <- 4 / (r.a)^2
	nrow.cls <- nrow(cls)
	z <- matrix(nrow = nrow(data.test))
	num.inpvar <- (ncol(cls) - num.outvar)
	nrow.dt.test <- nrow(data.test)
	miu <- matrix(nrow = nrow.cls)
	inpvar.ctr <- cls[, 1:num.inpvar]
	cum.up <- 0
	cum.down <- 0
	
	data.test <- norm.data(data.test, range.data.ori, min.scale = -1, max.scale = 1)
	
	z.ctr <- as.matrix(cls[, ncol(cls)])
		
	for (h in 1:nrow.dt.test){
		for (i in 1:nrow.cls){
			euc <- dist(rbind(data.test[h,], inpvar.ctr[i,]))	
			miu[i] <- exp(- alpha *  (euc ^ 2))
		}
		
		miu.t <- t(miu)	
		upper <- miu.t %*% z.ctr
		lower <- sum(miu)
		z[h] <- upper / lower
	}
	
	o.col <- ((ncol(range.data.ori) - num.outvar) + 1)
	res <- denorm.data(z, as.matrix(range.data.ori[, o.col]), min.scale = -1, max.scale = 1)
	
	return(res)
}