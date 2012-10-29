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
#' This is the main function to obtain a final result as predicted values for all methods in this package. 
#' In order to get predicted values, this function is run using an \code{\link{frbs-object}}, which is typically generated using \code{\link{frbs.learn}}.
#' 
#' @title The frbs prediction stage
#'
#' @param object an \code{\link{frbs-object}}.
#' @param newdata a matrix(m x n) of data for the prediction process, where m is the number of instances and 
#' n is the number of input variables.
#' @param ... the other parameters (not used)
#' @seealso \code{\link{frbs.learn}} and \code{\link{frbs.gen}} for learning and model generation, and the internal main functions of each method for the theory:  \code{\link{WM}}, \code{\link{GFS}}, \code{\link{SBC}},
#' \code{\link{HyFIS}}, \code{\link{ANFIS}}, \code{\link{DM}}, \code{\link{DENFIS}}, \code{\link{HGD}}, \code{\link{frcs}}, and \code{\link{MSGFS}}.
#' @return The predicted values. 
#' @aliases predict
#' @examples
#' ##################################
#' ## I. Regression Problem
#' ###################################
#' data.train <- matrix(c(5.2, -8.1, 4.8, 8.8, -16.1, 4.1, 10.6, -7.8, 5.5, 10.4, -29.0, 
#'                        5.0, 1.8, -19.2, 3.4, 12.7, -18.9, 3.4, 15.6, -10.6, 4.9, 1.9, 
#'                        -25.0, 3.7, 2.2, -3.1, 3.9, 4.8, -7.8, 4.5, 7.9, -13.9, 4.8, 
#'                        5.2, -4.5, 4.9, 0.9, -11.6, 3.0, 11.8, -2.1, 4.6, 7.9, -2.0, 
#'                        4.8, 11.5, -9.0, 5.5, 10.6, -11.2, 4.5, 11.1, -6.1, 4.7, 12.8, 
#'                        -1.0, 6.6, 11.3, -3.6, 5.1, 1.0, -8.2, 3.9, 14.5, -0.5, 5.7, 
#'                        11.9, -2.0, 5.1, 8.1, -1.6, 5.2, 15.5, -0.7, 4.9, 12.4, -0.8, 
#'                        5.2, 11.1, -16.8, 5.1, 5.1, -5.1, 4.6, 4.8, -9.5, 3.9, 13.2, 
#'                        -0.7, 6.0, 9.9, -3.3, 4.9, 12.5, -13.6, 4.1, 8.9, -10.0, 
#'                        4.9, 10.8, -13.5, 5.1), ncol = 3, byrow = TRUE)
#' 
#' data.fit <- matrix(c(10.5, -0.9, 5.2, 5.8, -2.8, 5.6, 8.5, -0.2, 5.3, 13.8, -11.9,
#'                      3.7, 9.8, -1.2, 4.8, 11.0, -14.3, 4.4, 4.2, -17.0, 5.1, 6.9, 
#'                      -3.3, 5.1, 13.2, -1.9, 4.6), ncol = 3, byrow = TRUE)
#'
#' newdata <- matrix(c(10.5, -0.9, 5.8, -2.8, 8.5, -0.2, 13.8, -11.9, 9.8, -1.2, 11.0,
#'                       -14.3, 4.2, -17.0, 6.9, -3.3, 13.2, -1.9), ncol = 2, byrow = TRUE)
#'
#' range.data<-matrix(c(0.9, 15.6, -29, -0.2, 3, 6.6), ncol=3, byrow = FALSE)
#' #############################################################
#' ## I.1 The example: Implementation of Wang & Mendel
#' #############################################################
## generate model especially rule database by training process
#' method.type <- "WM"
#' 
#' control.WM <- list(num.labels = 5, type.mf = 3, type.defuz = 1, 
#'                     type.tnorm = 1, type.snorm = 1) 
#' 
#' object <- frbs.learn(data.train, range.data, method.type, control.WM)
#'
#' ## the prediction process
#' ## The following code is be used for all methods
#' res <- predict(object, newdata) 
#' 
#' @export  
#' @method predict frbs
#' @S3method predict frbs
predict.frbs <- function(object, newdata, ...) {
mod <- object
data.test <- newdata

if(!inherits(mod, "frbs")) stop("not a legitimate frbs model")
  
##############################
## Split data of frbs.learn
#############################
m.type <- mod$method.type

## check for WM approach
if (m.type == "WM") {
	
	res.comp <- frbs.eng(mod, data.test)
	res <- res.comp$predicted.val
}

## check for SBC approach
else if(m.type == "SBC"){
	
	res <- SBC.test(mod, data.test)
}

## check for GFS approach
else if(m.type == "GFS"){
	range.data.ori <- mod$range.data.ori
	range.input.ori <- range.data.ori[, 1:(ncol(range.data.ori) - 1)]
	range.output.ori <- matrix(range.data.ori[, ncol(range.data.ori)])

	data.tst.norm <- norm.data(data.test, range.input.ori, min.scale = 0, max.scale = 1)
	res.comp <- frbs.eng(mod, data.tst.norm)
	
	res.denorm <- denorm.data(res.comp$predicted.val, range.output.ori, min.scale = 0, max.scale = 1)
	res <- res.denorm
}

## check for HyFIS approach
else if(m.type == "HYFIS"){
	range.data.ori <- mod$range.data.ori
	range.input.ori <- range.data.ori[, 1:(ncol(range.data.ori) - 1)]
	range.output.ori <- matrix(range.data.ori[, ncol(range.data.ori)])

	data.tst.norm <- norm.data(data.test, range.input.ori, min.scale = 0, max.scale = 1)
	res.comp <- frbs.eng(mod, data.tst.norm)
	res.denorm <- denorm.data(res.comp$predicted.val, range.output.ori, min.scale = 0, max.scale = 1)
	res <- res.denorm
}

## check for ANFIS approach
else if(m.type == "ANFIS"){
	range.data.ori <- mod$range.data.ori
	range.input.ori <- range.data.ori[, 1:(ncol(range.data.ori) - 1)]
	range.output.ori <- matrix(range.data.ori[, ncol(range.data.ori)])

	data.tst.norm <- norm.data(data.test, range.input.ori, min.scale = 0, max.scale = 1)
	res.comp <- frbs.eng(mod, data.tst.norm)
	
	res.denorm <- denorm.data(res.comp$predicted.val, range.output.ori, min.scale = 0, max.scale = 1)
	res <- res.denorm
}

## check for frcs approach
else if(m.type == "FRCS"){
	range.data.ori <- mod$range.data.ori
	
	data.tst.norm <- norm.data(data.test, range.data.ori, min.scale = 0, max.scale = 1)
	res <- frcs.eng(mod, data.tst.norm)
}

## check for DENFIS approach
else if(m.type == "DENFIS"){
	res <- DENFIS.eng(mod, data.test)
}

## check for DM approach
else if(m.type == "DM"){
	range.data.ori <- mod$range.data.ori
	
	range.input.ori <- matrix(range.data.ori[, 1:(ncol(range.data.ori) - 1)], nrow = 2)
	range.output.ori <- matrix(range.data.ori[, ncol(range.data.ori)], nrow = 2)

	data.tst.norm <- norm.data(data.test, range.input.ori, min.scale = 0, max.scale = 1)
	
	res.comp <- frbs.eng(mod, data.tst.norm)
	
	res.denorm <- denorm.data(res.comp$predicted.val, range.output.ori, min.scale = 0, max.scale = 1)
	res <- res.denorm
}

## check for HGD approach
else if(m.type == "HGD"){
	
	range.data.ori <- mod$range.data.ori
	range.input.ori <- matrix(range.data.ori[, 1:(ncol(range.data.ori) - 1)], nrow = 2)
	range.output.ori <- matrix(range.data.ori[, ncol(range.data.ori)], nrow = 2)

	data.tst.norm <- norm.data(data.test, range.input.ori, min.scale = 0, max.scale = 1)
	
	res.comp <- frbs.eng(mod, data.tst.norm)
	
	res.denorm <- denorm.data(res.comp$predicted.val, range.output.ori, min.scale = 0, max.scale = 1)
	res <- res.denorm
}

## check for MSGFS approach
else if(m.type == "MSGFS"){
	
	range.data.ori <- mod$range.data.ori
	range.input.ori <- matrix(range.data.ori[, 1:(ncol(range.data.ori) - 1)], nrow = 2)
	range.output.ori <- matrix(range.data.ori[, ncol(range.data.ori)], nrow = 2)

	data.tst.norm <- norm.data(data.test, range.input.ori, min.scale = 0, max.scale = 1)
	
	res.comp <- MSGFS.test(mod, data.tst.norm)
	
	res.denorm <- denorm.data(res.comp, range.output.ori, min.scale = 0, max.scale = 1)
	res <- res.denorm
}

return(res)
}




