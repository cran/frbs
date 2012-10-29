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
#' This is one of the central functions of the package. This function is used to 
#' generate/learn the model from numerical data using fuzzy rule-based systems.
#' 
#' This function makes accessible all ten learning methods that are implemented 
#' in this package. All of the methods use this function as interface for the learning 
#' stage, so users do not need to call other functions in the learning phase. 
#' In order to obtain good results, users need to adjust some parameters such as the 
#' number of labels, the type of the shape of the membership function, the maximal number of iterations, 
#' the step size of the gradient descent, or other method-dependent parameters which are collected in the \code{control}
#' parameter. After creating the model using this function, it can be used to predict new data with \code{\link{predict}}.
#'
#' @title The frbs model building function
#'
#' @param data.train a matrix(m x n) of data for the training process, where m is the number of instances and 
#' n is the number of variables; the last column is the output variable.
#' @param range.data a matrix(2 x n) containing the range of the data, where n is the number of variables, and
#' first and second rows are the minimum and maximum values, respectively. 
#' @param method.type this parameter determines the learning algorithm to use. The following methods are implemented: 
#' "WM": Wang and Mendel, "SBC": subtractive clustering, "GFS": genetic fuzzy system, "HYFIS": hybrid neural fuzzy inference systems, 
#' "ANFIS": adaptive neuro-fuzzy inference systems, "FRCS": fuzzy rule-based classification systems, "DENFIS": 
#' dynamic evolving neuro-fuzzy inference systems, "HGD": heuristic and gradient descent method, "DM": fuzzy inference rules by descent method, and 
#' "MSGFS": The multi-stage genetic fuzzy systems based on iterative rule learning approach. 
#' @param control a list containing all arguments, depending on the learning algorithm to use. 
#' 
#' \bold{WM method}
#' \itemize{
#' \item num.labels: a positive integer to determine the number of labels (fuzzy terms). The default value is 7.
#' \item type.mf: the type of the membership function. 
#' 1 is for triangular, while 2, 3, 4, and 5 are for trapezoid, Gaussian, sigmoid, and generalized bell. The default value is 3.
#' \item type.defuz: the type of the defuzzification method. Here, 1 means we use the weighted average method, and 
#' 2, 3, 4, and 5 mean we use first, last, mean maxima, and modified COG, respectively. The default value is 1.
#' \item type.tnorm: the type of t-norm. 1 means standard type(min), and 2, 3, 4, and 5 mean 
#' Hamacher product, Yager class (with tao = 1), product, and bounded product, respectively. For more detail, please have a look at \code{\link{inference}}. The default value is 1.
#' \item type.snorm: the type of s-norm. 1 means standard type(max), and 2, 3, 4, and 5 mean  
#' Hamacher sum, Yager class (with tao = 1), sum, and bounded sum, respectively. For more detail, please have a look at \code{\link{inference}}. The default value is 1.
#' \item name: a name for the model. The default value is "sim-0".
#' }
#' \bold{HyFIS, ANFIS, and DM methods}
#' \itemize{
#' \item num.labels: a positive integer to determine the number of labels (fuzzy terms). The default value is 7.
#' \item max.iter: a positive integer to determine the maximal number of iterations. The default value is 100.
#' \item step.size: the step size of the gradient descent, a real number between 0 and 1. The default value is 0.01.
#' \item name: a name for the model. The default value is "sim-0".
#' }
#' \bold{SBC method}
#' \itemize{
#' \item r.a:  a positive constant which is effectively the radius defining a neighborhood. The default value is 0.5.
#' \item eps.high: an upper threshold value. The default value is 0.5.
#' \item eps.low: a lower threshold value. The default value is 0.15.
#' \item name: a name for the model. The default value is "sim-0".
#' }
#' \bold{HGD method}
#' \itemize{
#' \item num.labels: a positive integer to determine the number of labels (fuzzy terms). The default value is 7.
#' \item max.iter: a positive integer to determine the maximal number of iterations. The default value is 100.
#' \item step.size: a real number between 0 and 1. The default value is 0.01.
#' \item alpha.heuristic: a positive real number representing a heuristic value. The default value is 1.
#' \item name: a name for the model. The default value is "sim-0".
#' }
#' \bold{GFS method}
#' \itemize{
#' \item num.labels: a positive integer to determine the number of labels (fuzzy terms). The default value is 5.
#' \item popu.size: an integer number of the population size. The default value is 10.
#' \item persen_cross: a percentage of crossover. The default value is 0.9.
#' \item max.iter: a positive integer to determine the maximal number of iterations. The default value is 10.
#' \item persen_mutant: a percentage of mutation. The default value is 0.01.
#' \item classification: a boolean whether it is a classification problem or not. The default value is FALSE.
#' \item name: a name for the model. The default value is "sim-0". 
#' }
#' \bold{FRCS method}
#' \itemize{
#' \item num.labels: a positive integer to determine the number of labels (fuzzy terms). The default value is 7.
#' \item type.mf: the type of the membership function. 
#' 1 is for triangular, while 2, 3, 4, and 5 are for trapezoid, Gaussian, sigmoid, and generalized bell. The default value is 1.
#' \item name: a name for the model. The default value is "sim-0".
#' }
#' \bold{DENFIS method}
#' \itemize{
#' \item Dthr: the threshold value for the envolving clustering method (ECM), between 0 and 1. The default value is 0.1.
#' \item max.iter: a positive integer to determine the maximal number of iterations. The default value is 100.
#' \item step.size: the step size of the least squares method, between 0 and 1. The default value is 0.01.
#' \item d: a parameter for the width of the triangular membership function. The default value is 2.
#' \item name: a name for the model. The default value is "sim-0".
#' }
#' \bold{MSGFS method}
#' \itemize{
#' \item popu.size: an integer number of the population size. The default value is 50.
#' \item persen_cross: a percentage of crossover. The default value is 0.6.
#' \item max.iter: a positive integer to determine the maximal number of iterations. The default value is 100.
#' \item persen_mutant: a percentage of mutation. The default value is 0.3.
#' \item epsilon: a real number between 0 and 1 representing the boundary of covering factor. The default value is 0.05.
#' \item name: a name for the model. The default value is "sim-0". 
#' }
#'
#' @seealso \code{\link{predict}} for the prediction phase, and the following main functions of each of the methods for theoretical background and references:  \code{\link{WM}}, \code{\link{GFS}}, \code{\link{SBC}},
#' \code{\link{HyFIS}}, \code{\link{ANFIS}}, \code{\link{DM}}, \code{\link{DENFIS}}, \code{\link{HGD}}, \code{\link{frcs}}, and \code{\link{MSGFS}}.
#' @return The \code{\link{frbs-object}}. 
#' @examples
#' ##################################
#' ## I. Regression Problem
#' ## The data has two input variables and one output variable.  
#' ## We separate them into training, fitting, and testing data.
#' ## data.train, data.fit, data.test, and range.data are inputs 
#' ## for all methods except the frcs method.
#' ###################################
#' ## The simulation might take a long time dependent on the hardware we use. 
#' ## We might get better results if we take other parameters. 
#'
#' data.train <- matrix(c(5.2, -8.1, 4.8, 8.8, -16.1, 4.1, 10.6, -7.8, 5.5, 10.4, -29.0, 
#'                       5.0, 1.8, -19.2, 3.4, 12.7, -18.9, 3.4, 15.6, -10.6, 4.9, 1.9, 
#'                       -25.0, 3.7, 2.2, -3.1, 3.9, 4.8, -7.8, 4.5, 7.9, -13.9, 4.8, 
#'                       5.2, -4.5, 4.9, 0.9, -11.6, 3.0, 11.8, -2.1, 4.6, 7.9, -2.0, 
#'                       4.8, 11.5, -9.0, 5.5, 10.6, -11.2, 4.5, 11.1, -6.1, 4.7, 12.8, 
#'                       -1.0, 6.6, 11.3, -3.6, 5.1, 1.0, -8.2, 3.9, 14.5, -0.5, 5.7, 
#'                       11.9, -2.0, 5.1, 8.1, -1.6, 5.2, 15.5, -0.7, 4.9, 12.4, -0.8, 
#'                       5.2, 11.1, -16.8, 5.1, 5.1, -5.1, 4.6, 4.8, -9.5, 3.9, 13.2, 
#'                       -0.7, 6.0, 9.9, -3.3, 4.9, 12.5, -13.6, 4.1, 8.9, -10.0, 
#'                       4.9, 10.8, -13.5, 5.1), ncol = 3, byrow = TRUE)
#' 
#' data.fit <- matrix(c(10.5, -0.9, 5.2, 5.8, -2.8, 5.6, 8.5, -0.2, 5.3, 13.8, -11.9,
#'                      3.7, 9.8, -1.2, 4.8, 11.0, -14.3, 4.4, 4.2, -17.0, 5.1, 6.9, 
#'                      -3.3, 5.1, 13.2, -1.9, 4.6), ncol = 3, byrow = TRUE)
#'
#' data.test <- matrix(c(10.5, -0.9, 5.8, -2.8, 8.5, -0.2, 13.8, -11.9, 9.8, -1.2, 11.0,
#'                      -14.3, 4.2, -17.0, 6.9, -3.3, 13.2, -1.9), ncol = 2, byrow = TRUE)
#'
#' range.data<-matrix(c(0.9, 15.6, -29, -0.2, 3, 6.6), ncol=3, byrow = FALSE)
#'
#' #############################################################
#' ## I.1 The example: Implementation of Wang & Mendel
#' #############################################################
#' method.type <- "WM" 
#' 
#' ## collect control parameters into a list
#' control.WM <- list(num.labels = 3, type.mf = 3, name = "Sim-0") 
#' 
#' ## generate the model and save it as object.WM
#' object.WM <- frbs.learn(data.train, range.data, method.type, control.WM)
#'
#' #############################################################
#' ## I.2 The example: Implementation of SBC
#' #############################################################
#' method.type <- "SBC" 
#' control.SBC <- list(r.a = 0.5, eps.high = 0.5, eps.low = 0.15, name = "Sim-0")
#'
#' \dontrun{object.SBC <- frbs.learn(data.train, range.data, method.type, control.SBC)}
#'
#' #############################################################
#' ## I.3 The example: Implementation of HyFIS
#' #############################################################
#' method.type <- "HyFIS"
#' 
#' control.HyFIS <- list(num.labels = 5, max.iter = 50, step.size = 0.01, 
#'                  name = "Sim-0")
#' 
#' \dontrun{object.HyFIS <- frbs.learn(data.train, range.data, method.type, control.HyFIS)}
#'
#' #############################################################
#' ## I.4 The example: Implementation of ANFIS
#' #############################################################
#' method.type <- "ANFIS" 
#'
#' control.ANFIS <- list(num.labels = 5, max.iter = 100, step.size = 0.01, 
#'                       name = "Sim-0") 
#'
#' \dontrun{object.ANFIS <- frbs.learn(data.train, range.data, method.type, control.ANFIS)}
#'
#' #############################################################
#' ## I.5 The example: Implementation of DENFIS
#' #############################################################
#' 
#' control.DENFIS <- list(Dthr = 0.1, max.iter = 100, step.size = 0.001, d = 2, 
#'                        name = "Sim-0")
#' method.type <- "DENFIS"
#' 
#' \dontrun{object.DENFIS <- frbs.learn(data.train, range.data, method.type, control.DENFIS)}
#'
#' #############################################################
#' ## I.6 The example: Implementation of DM
#' #############################################################
#' method.type <- "DM"
#'  
#' control.DM <- list(num.labels = 5, max.iter = 100, step.size = 0.01, name = "Sim-0") 
#' \dontrun{object.DM <- frbs.learn(data.train, range.data, method.type, control.DM)}
#'
#' #############################################################
#' ## I.7 The example: Implementation of HGD
#' #############################################################
#' method.type <- "HGD" 
#'  
#' control.HGD <- list(num.labels = 5, max.iter = 100, step.size = 0.01, 
#'                alpha.heuristic = 1, name = "Sim-0") 
#' \dontrun{object.HGD <- frbs.learn(data.train, range.data, method.type, control.HGD)}
#'
#' #############################################################
#' ## I.8 The example: Implementation of GFS
#' #############################################################
#' method.type <- "GFS" 
#'  
#' control.GFS <- list(num.labels = 5, popu.size = 5, persen_cross = 0.9, 
#'                     max.iter = 5, persen_mutant = 0.1, 
#'                     classification = FALSE, name="sim-0") 
#' \dontrun{object.GFS <- frbs.learn(data.train, range.data, method.type, control.GFS)}
#'
#' #############################################################
#' ## I.9 The example: Implementation of MSGFS
#' #############################################################
#' method.type <- "MSGFS" 
#'  
#' control.MSGFS <- list(popu.size = 15, persen_cross = 0.6, 
#'                     max.iter = 20, persen_mutant = 0.3, 
#'                     epsilon = 0.05, name="sim-0") 
#' \dontrun{object.MSGFS <- frbs.learn(data.train, range.data, method.type, control.MSGFS)}
#'
#' #############################################################
#' ## II. Classification Problem 
#' #############################################################
#' ## II.1 The example: Implementation of fuzzy rule-based classification systems (frcs)
#' ## The iris dataset is shuffled and divided into training and 
#' ## testing data. Bad results in the predicted values may
#' ## result from casual imbalanced classes in the training data. 
#' 
#' data(iris)
#' irisShuffled <- iris[sample(nrow(iris)),]
#' irisShuffled[,5] <- unclass(irisShuffled[,5])
#' tra.iris <- irisShuffled[1:105,]
#' tst.iris <- irisShuffled[106:nrow(irisShuffled),1:4]
#' real.iris <- matrix(irisShuffled[106:nrow(irisShuffled),5], ncol = 1)
#'
#' range.data.input <- matrix(c(4.3, 7.9, 2.0, 4.4, 1.0, 6.9, 0.1, 2.5), nrow=2)
#' 
#' ## generate the model 
#' method.type <- "FRCS"
#' control <- list(num.labels = 7, type.mf = 1) 
#' 
#' object <- frbs.learn(tra.iris, range.data.input, method.type, control)
#' 
#' ## conduct the prediction process
#' res.test <- predict(object, tst.iris)
#' @export
frbs.learn <- function(data.train, range.data, method.type = c("WM"), control=list()){
 
## get type of method 
method.type <- toupper(method.type)

## initialize mod
mod <- NULL
    
## if we are using wang & mendel methods
if(method.type == "WM"){

	## getting all of parameters
	control <- setDefaultParametersIfMissing(control, list(num.labels = 7, type.mf = 3, type.defuz = 1, type.tnorm = 1, type.snorm = 1, name="sim-0"))
	
	## get parameters
	num.labels <- control$num.labels
	type.mf <- control$type.mf
	name <- control$name
	num.labels <- matrix(rep(num.labels, ncol(range.data)), nrow=1)
	
	## generate model
	modelSpecific <- WM(range.data, data.train, num.labels, type.mf)
	
	## collect results as model
	mod <- modelSpecific
	mod$type.model <- 1
	mod$func.tsk <- NULL
	mod$type.defuz <- control$type.defuz
	mod$type.tnorm <- control$type.tnorm
	mod$type.snorm <- control$type.snorm
	
}

## Subtractive Clustering
else if(method.type == "SBC"){
	
	## get all of parameters
	control <- setDefaultParametersIfMissing(control, list(r.a = 0.5, eps.high = 0.5, eps.low = 0.15, name ="sim-0"))
	r.a <- control$r.a
	eps.high <- control$eps.high
	eps.low <- control$eps.low
	name <- control$name
	range.data.ori <- range.data
	
	## generate model 
	modelSpecific <- SBC(data.train, range.data.ori, r.a, eps.high, eps.low)
	mod <- modelSpecific
}

## Genetic Fuzzy System
else if (method.type == "GFS"){
	
	## get all of parameters
	control <- setDefaultParametersIfMissing(control, list(num.labels = 5, popu.size = 10, persen_cross = 0.9, max.iter = 10, persen_mutant = 0.1, classification = FALSE, name="sim-0"))
	
	## getting all of parameters
	range.data.ori <- range.data
	data.train.ori <- data.train
	n.labels <- control$num.labels
	popu.size <- control$popu.size
	persen_cross <- control$persen_cross
	persen_mutant <- control$persen_mutant
	max.iter <- control$max.iter
	name <- control$name
	classification <- control$classification
	
	## normalize range of data and data training
	range.data.norm <- range.data.ori
	range.data.norm[1, ] <- 0
	range.data.norm[2, ] <- 1	
	data.train.norm <- norm.data(data.train.ori, range.data.ori, min.scale = 0, max.scale = 1)
	
	## generate labels of each variables
	num.labels <- matrix(rep(n.labels, ncol(range.data)), nrow=1)
	
	## generate model
	modelSpecific <- GFS(data.train.norm, range.data.norm, num.labels, popu.size, persen_cross, persen_mutant, max.iter, classification, range.data.ori)
	
	## collect model
	mod <- modelSpecific
	mod$type.model <- 1
	mod$func.tsk <- NULL
	mod$type.defuz <- 1
	mod$type.tnorm <- 1 
	mod$type.snorm <- 1 
}

## HyFIS
else if (method.type == "HYFIS"){

	## get all of parameters
	control <- setDefaultParametersIfMissing(control, list(num.labels = 7, max.iter = 100, step.size = 0.01, name = "sim-0"))
	n.labels <- control$num.labels
	max.iter <- control$max.iter
	step.size <- control$step.size
	name <- control$name
	range.data.ori <- range.data
	data.train.ori <- data.train

	## make labels of each variables
	num.labels <- matrix(rep(n.labels, ncol(range.data)), nrow=1)	
	
	## normalize range of data and data training
	range.data.norm <- range.data.ori
	range.data.norm[1, ] <- 0
	range.data.norm[2, ] <- 1	
	data.tra.norm <- norm.data(data.train.ori, range.data.ori, min.scale = 0, max.scale = 1)
	
	## generate model
	modelSpecific <- HyFIS(range.data.norm, data.tra.norm, num.labels, max.iter, range.data.ori, step.size)
	mod <- modelSpecific

}

#ANFIS
else if (method.type == "ANFIS"){

	## get all of parameters
	control <- setDefaultParametersIfMissing(control, list(num.labels = 7, max.iter = 100, step.size = 0.01, name="sim-0"))
	range.data.ori <- range.data
	data.train.ori <- data.train
	n.labels <- control$num.labels
	max.iter <- control$max.iter
	step.size <- control$step.size
	name <- control$name
	
	## normalize range of data and data training
	range.data.norm <- range.data.ori
	range.data.norm[1, ] <- 0
	range.data.norm[2, ] <- 1	
	data.tra.norm <- norm.data(data.train.ori, range.data.ori, min.scale = 0, max.scale = 1)
	
	## generate labels of each variables
	num.labels <- matrix(rep(n.labels, ncol(range.data.norm)), nrow=1)

	## generate model
	modelSpecific <- ANFIS(range.data.norm, data.tra.norm, num.labels, max.iter, range.data.ori, step.size)	
	mod <- modelSpecific
	
	## change rule format into TSK model 
	rule <- mod$rule	
	length.rule <- length(rule)
	temp <- matrix(rule[[1]], nrow = 1)
	for (i in 2 : length.rule){
		temp.1 <- matrix(rule[[i]], nrow = 1)
		temp <- rbind(temp, temp.1)
	}	
	m.rule <- temp
	m.rule <- m.rule[, 1 : (ncol(m.rule) - 1)]
	mod$rule <- m.rule
	
}

#Fuzzy Rule-based Classification System
else if (method.type == "FRCS"){

	## get all of parameters 
	control <- setDefaultParametersIfMissing(control, list(num.labels = 7, type.mf = 1, name="sim-0"))
	range.data.input <- range.data
	n.labels <- control$num.labels
	type.mf <- control$type.mf
	name <- control$name

	## make range of data according to class on data training
	range.data.out <- matrix(c(min(data.train[, ncol(data.train)], na.rm = TRUE) - 0.4999, max(data.train[, ncol(data.train)], na.rm = TRUE) + 0.4999), nrow = 2)
	num.class <- floor(max(range.data.out))
	
	## normalize range of data and data training
	range.data.norm <- range.data.input
	range.data.norm[1, ] <- 0
	range.data.norm[2, ] <- 1	
	range.data.ori <- range.data.input
	range.data.inout <- cbind(range.data.norm, range.data.out)	
	data.tra.norm <- norm.data(data.train[, 1 : ncol(data.train) - 1], range.data.ori, min.scale = 0, max.scale = 1)
	data.train <- cbind(data.tra.norm, matrix(data.train[, ncol(data.train)], ncol = 1))
	
	## generate model
	modelSpecific <- frcs(range.data.inout, data.train, n.labels, num.class, type.mf)
	mod <- modelSpecific
	mod$range.data.ori <- range.data.ori

}

## DENFIS 
else if (method.type == "DENFIS"){	
  
	## get all of parameters
    control <- setDefaultParametersIfMissing(control, list(Dthr = 0.1, max.iter = 100, step.size = 0.01, d = 2, name="sim-0"))
	Dthr <- control$Dthr
	max.iter <- control$max.iter
	step.size <- control$step.size
	d <- control$d
	name <- control$name
	range.data.ori <- range.data
	
	## generate model
	modelSpecific <- DENFIS(data.train, range.data.ori, Dthr, max.iter, step.size, d)
	mod <- modelSpecific

}

## DM
else if (method.type == "DM"){

	## get all of parameters
	control <- setDefaultParametersIfMissing(control, list(num.labels = 7, max.iter = 100, step.size = 0.01, name="sim-0"))	
	range.data.ori <- range.data	
	data.train.ori <- data.train
	n.labels <- control$num.labels
	max.iter <- control$max.iter
	step.size <- control$step.size
	name <- control$name

	## normalize range of data and data training
	range.data.norm <- range.data.ori
	range.data.norm[1, ] <- 0
	range.data.norm[2, ] <- 1	
	data.tra.norm <- norm.data(data.train.ori, range.data.ori, min.scale = 0, max.scale = 1)
	num.labels <- matrix(rep(n.labels, ncol(range.data.norm)), nrow=1)

	## generate model
	modelSpecific <- DM(range.data.norm, data.tra.norm, num.labels, max.iter, step.size)
	mod <- modelSpecific
	mod$range.data.ori <- range.data.ori
}

## HGD
else if (method.type == "HGD"){

	## get all of parameters
	control <- setDefaultParametersIfMissing(control, list(num.labels = 7, max.iter = 100, step.size = 0.01, alpha.heuristic = 1, name="sim-0"))
	range.data.ori <- range.data
	data.train.ori <- data.train
	n.labels <- control$num.labels
	max.iter <- control$max.iter
	step.size <- control$step.size
	name <- control$name
	alpha.heuristic <- control$alpha.heuristic
		
	## normalize range of data and data training
	range.data.norm <- range.data.ori
	range.data.norm[1, ] <- 0
	range.data.norm[2, ] <- 1	
	data.tra.norm <- norm.data(data.train.ori, range.data.ori, min.scale = 0, max.scale = 1)
	num.labels <- matrix(rep(n.labels, ncol(range.data.norm)), nrow=1)

	## generate model
	modelSpecific <- HGD(range.data.norm, data.tra.norm, num.labels, max.iter, step.size, alpha.heuristic)
	mod <- modelSpecific
	mod$range.data.ori <- range.data.ori

}

## Multi Stage Genetic Fuzzy System
else if (method.type == "MSGFS"){
	## get all of parameters
	control <- setDefaultParametersIfMissing(control, list(popu.size = 50, persen_cross = 0.6, max.iter = 100, persen_mutant = 0.3, epsilon = 0.05, name="sim-0"))
	
	## getting all of parameters
	range.data.ori <- range.data
	data.train.ori <- data.train
	popu.size <- control$popu.size
	persen_cross <- control$persen_cross
	persen_mutant <- control$persen_mutant
	max.iter <- control$max.iter
	epsilon <- control$epsilon
	name <- control$name
	
	## normalize range of data and data training
	range.data.norm <- range.data.ori
	range.data.norm[1, ] <- 0
	range.data.norm[2, ] <- 1	
	data.tra.norm <- norm.data(data.train.ori, range.data.ori, min.scale = 0, max.scale = 1)
		
	modelSpecific <- MSGFS(data.tra.norm, popu.size, persen_cross, persen_mutant, max.iter, range.data.ori, epsilon)
	mod <- modelSpecific
}

mod$method.type <- method.type

mod$name <- name

mod <- frbsObjectFactory(mod)

return(mod)
}

## checking missing parameters
setDefaultParametersIfMissing <- function(control, defaults) {
  for(i in names(defaults)) {
    if(is.null(control[[i]])) control[[i]] <- defaults[[i]]
  }
  control
}

#' This function creates objects of type \code{frbs}. Currently, its 
#' implementation is very basic and does no argument checking, as 
#' it is only used internally.
#' 
#' The members of the \code{frbs} object depend on the used learning method. The following list describes all of the members that can be present. 
#' \describe{
#' \item{range.input}{the range of the input data. Whether it is normalized to lie between 0 and 1 or not depends on the selected method.}
#' \item{range.output}{the range of the output data. Whether it is normalized to lie between 0 and 1 or not  depends on the selected method.}
#' \item{num.varinput}{the number of input variables.}
#' \item{num.fvalinput}{the number of fuzzy terms of the input variables.}
#' \item{names.varinput}{the generated names of fuzzy terms of the input variables.}
#' \item{num.fvaloutput}{the number of fuzzy terms of the output variable.}
#' \item{varout.mf}{a matrix to generate the shapes of the membership functions for the output variable. 
#' The first row represents the shape of the membership functions, the other rows contain the parameters that have been generated. 
#' Whether the values of parameters within the matrix are normalized to lie between 0 and 1 or not depends on the selected method.}
#' \item{names.varoutput}{generated names of the output variable.}
#' \item{rule}{the fuzzy IF-THEN rules; In the MSGFS case, a rule refers to the parameter values of the membership function which represents the rule.}
#' \item{varinp.mf}{a matrix to generate the shapes of the membership functions for the input variables. 
#' The first row represents the shape of the membership functions, 
#' the other rows contain the parameters that have been generated.
#' Whether the values of parameters within the matrix are normalized to lie between 0 and 1 or not depends on the selected method.}
#' \item{type.model}{the model type. Here, 1 refers to Mamdani model, and 2 refers to Takagi Sugeno Kang model on the consequence part.}
#' \item{func.tsk}{a matrix of the Takagi Sugeno Kang model consequent part of the fuzzy IF-THEN rules.}
#' \item{type.defuz}{the type of the defuzzification method.}
#' \item{type.tnorm}{the type of the t-norm method.}
#' \item{type.snorm}{the type of the s-norm method.}
#' \item{method.type}{the type of the selected method.}
#' \item{name}{the name given to the model.}
#' \item{range.data.ori}{range of the original data (before normalization).}
#' \item{cls}{cluster centers.}
#' \item{Dthr}{the boundary parameter of the DENFIS method.}
#' \item{d}{the multiplier parameters of the DENFIS method.}
#' \item{r.a}{the neighborhood factor of SBC.}
#' \item{degree.rule}{certainty degree of rules.}
#' \item{rule.data.num}{representation of the rules in matrix form.}
#' \item{grade.cert}{grade of certainty for classification problems.}
#' \item{alpha.heuristic}{a parameter for the heuristic of the HGD method.}
#' }
#' 
#' @title The object factory for frbs objects
#' @param mod a list containing all the attributes for the object
#' @return an object of type \code{frbs}
#' @aliases frbs-object
frbsObjectFactory <- function(mod){
	class(mod) <- "frbs"
	return(mod)
}
