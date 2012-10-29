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
#' 
#' This function is an internal function for the prediction phase using the DENFIS method. The user should use this function not directly, but with calling \code{\link{predict}}.  
#'
#' @title DENFIS prediction function
#'
#' @param object the frbs model. See \code{\link{frbs-object}}. 
#' @param newdata a matrix(m x n) of data for the prediction process, where m is the number of instances and n is the number of input variables.
#' @seealso \code{\link{DENFIS}}
#' @return a matrix of predicted values
# @export
DENFIS.eng <- function(object, newdata){
mod <- object
data.test <- newdata

model <- mod

min.scale <- 0
max.scale <- 1

## get parameters
cluster.c <- model$cls
func.tsk <- model$func.tsk
range.data.ori <- model$range.data.ori
Dthr <- model$Dthr
d <- model$d

## normalize cluster center and data test
cluster.c <- norm.data(cluster.c, range.data.ori, min.scale, max.scale)
data.test <- norm.data(data.test, range.data.ori[, 1 : (ncol(range.data.ori) - 1)], min.scale, max.scale)

## inisialitation
num.cls <- nrow(cluster.c)
num.dt <- nrow(data.test) 
num.inpvar <- ncol(data.test)
temp <- matrix(nrow = num.cls, ncol=num.inpvar)
miu.rule <- matrix(nrow = num.dt, ncol = num.cls)
def <- matrix(nrow=num.dt, ncol = 1)

## calculate degree of MF
for (i in 1 : num.dt){
	for (j in 1 : num.cls){
		for (k in 1 : num.inpvar){
			b <- cluster.c[j, k]
			a <- b - d * Dthr
			cc <- b + d * Dthr
			x <- data.test[i, k]
			left <- (x - a)/(b - a)
			right <- (cc - x)/(cc - b)
			
			temp[j, k] <- max(min(left, right), 0)
		}
		
		miu.rule[i, j] <- prod(temp[j, ])	
	}	
}


### Calculate defuzzification
	for (j in 1 : num.dt) {
		data.k <- (data.test[j, ])
		data.m <- as.matrix(data.k)
		func.tsk.var <- as.matrix(func.tsk[, 1: (ncol(func.tsk) - 1)])
		func.tsk.cont <- t(t(func.tsk[, ncol(func.tsk)]))
		ff <- func.tsk.var %*% data.m + func.tsk.cont
		
		miu.rule.t <- as.matrix(miu.rule[j, ])
		
		cum <- sum(ff * miu.rule.t)
		
		div <- sum(miu.rule.t)
		
		if (div == 0)
			def[j] <- 0
		else {
			def[j] <- cum / div
		}
	}

range.output <- matrix(range.data.ori[, ncol(range.data.ori)], nrow = 2)
res <- denorm.data(def, range.output, min.scale, max.scale)

return (res)
}