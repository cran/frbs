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
#' The role of this function is to update parameters within the HGD method. 
#' This function is called by the main function 
#' of the HGD method, see \code{\link{HGD}}.
#'
#' @title HGD updating function
#'
#' @param data.train a matrix(m x n) of data for the training process, where m is the number of instances and 
#' n is the number of variables; the last column is the output variable.
#' @param miu.rule a matrix with the degrees of rules which is the result of the \code{\link{inference}}.
#' @param func.tsk a matrix of parameters of the function on the consequent part using the Takagi Sugeno Kang model. See \code{\link{rulebase}}.
#' @param varinp.mf a matrix of parameters of membership functions of the input variables.
#' @param step.size a real number between 0 and 1 representing the step size of the gradient descent. 
#' @param def a matrix which is obtained by the \code{\link{defuzzifier}}.
#'  
#' @seealso \code{\link{HGD}}
# @return a list of the new parameter
# @export
HGD.update <- function(data.train, miu.rule, func.tsk, varinp.mf, step.size = 0.01, def){

data <- data.train
alpha <- step.size

#TSK type

data.k <-(data[1, 1 : (ncol(data.train) - 1)])
data.m <- t(t(data.k))
	
ff <- func.tsk 
miu.rule.t <- as.matrix(miu.rule)

gal <- (def - data.train[1, ncol(data.train)])
	
func.tsk.upd <- func.tsk
for (mm in 1 : nrow(func.tsk.upd)){
	sum.miu <- sum(miu.rule.t, na.rm = TRUE)
	if (sum.miu != 0)
		func.tsk.upd[mm, 1] <- func.tsk.upd[mm, 1] - alpha * gal * (miu.rule[mm])
}


func.tsk.new <- func.tsk.upd

param.new <- list(func.tsk.new)
######################################
#=====end of code gradient descent 
#########################################

return(param.new)
}