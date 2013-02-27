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
#' This function makes accessible all fourteen learning methods that are implemented 
#' in this package. All of the methods use this function as interface for the learning 
#' stage, so users do not need to call other functions in the learning phase. 
#' In order to obtain good results, users need to adjust some parameters such as the 
#' number of labels, the type of the shape of the membership function, the maximal number of iterations, 
#' the step size of the gradient descent, or other method-dependent parameters which are collected in the \code{control}
#' parameter. After creating the model using this function, it can be used to predict new data with \code{\link{predict}}.
#'
#' @title The frbs model building function
#'
#' @param data.train a data frame or matrix (m x n) of data for the training process, 
#'        where m is the number of instances and 
#'        n is the number of variables; the last column is the output variable. It should be noted that
#'        the training data must be expressed in numbers (numerical data). 
#' @param range.data a matrix(2 x n) containing the range of the data, where n is the number of variables, and
#'        first and second rows are the minimum and maximum values, respectively. It should be noted that
#'        for "FRBCS.W", "FRBCS.CHI", "GFS.GCCL", "FH.GBML", and "SLAVE", n represents the number of input variables only
#'        (without the output variable).
#' @param method.type this parameter determines the learning algorithm to be used. 
#'        The following methods are implemented: 
#' \itemize{
#' \item "WM": Wang and Mendel's technique to handle regression tasks;
#' \item "SBC": subtractive clustering method to handle regression tasks;
#' \item "HYFIS": hybrid neural fuzzy inference systems to handle regression tasks;
#' \item "ANFIS": adaptive neuro-fuzzy inference systems to handle regression tasks;
#' \item "FRBCS.W": fuzzy rule-based classification systems with weight factor based on Ishibuchi's method 
#'                  to handle classification tasks; 
#' \item "FRBCS.CHI": fuzzy rule-based classification systems based on Chi's method to handle
#'                  classification tasks; 
#' \item "DENFIS": dynamic evolving neuro-fuzzy inference systems to handle regression tasks;
#' \item "FS.HGD": fuzzy system using heuristic and gradient descent method to handle regression tasks; 
#' \item "FIR.DM": fuzzy inference rules by descent method to handle regression tasks; 
#' \item "GFS.FR.MOGUL": genetic fuzzy systems for fuzzy rule learning based on the MOGUL methodology 
#'                    to handle regression tasks;
#' \item "GFS.THRIFT": Thrift's technique based on genetic algorithms to handle regression tasks;
#' \item "GFS.GCCL": Ishibuchi's method based on genetic cooperative-competitive learning
#'                   to handle classification tasks;
#' \item "FH.GBML": Ishibuchi's method based on hybridization of genetic cooperative-competitive learning and Pittsburgh to handle
#'                   classification tasks;
#' \item "SLAVE": structural learning algorithm on vague environment to handle classification tasks.
#' }
#' @param control a list containing all arguments, depending on the learning algorithm to use. 
#' 
#' \bold{WM method}
#' \itemize{
#' \item num.labels: a positive integer to determine the number of labels (fuzzy terms). 
#'       The default value is 7.
#' \item type.mf: the type of the membership function. 
#'       1 is for triangular, while 2, 3, 4, and 5 are for trapezoid, Gaussian, sigmoid, 
#'       and generalized bell, respectively. The default value is 3.
#' \item type.defuz: the type of the defuzzification method. Here, 1 means we use 
#'       the weighted average method, and 2, 3, 4, and 5 mean we use first, last, mean maxima, 
#'       and modified COG, respectively. The default value is 1.
#' \item type.tnorm: the type of t-norm. 1 means standard type (min), and 2, 3, 4, and 5 mean 
#'       Hamacher product, Yager class (with tao = 1), product, and bounded product, respectively. For more detail, please have a look at \code{\link{inference}}. The default value is 1.
#' \item type.snorm: the type of s-norm. 1 means standard type (max), and 2, 3, 4, and 5 mean  
#'       Hamacher sum, Yager class (with tao = 1), sum, and bounded sum, respectively. 
#'       For more detail, please have a look at \code{\link{inference}}. The default value is 1.
#' \item name: a name for the model. The default value is "sim-0".
#' }
#' \bold{HyFIS, ANFIS, and FIR.DM methods}
#' \itemize{
#' \item num.labels: a positive integer to determine the number of labels (fuzzy terms). 
#'       The default value is 7.
#' \item max.iter: a positive integer to determine the maximal number of iterations. 
#'       The default value is 100.
#' \item step.size: the step size of the gradient descent, a real number between 0 and 1. 
#'       The default value is 0.01.
#' \item name: a name for the model. The default value is "sim-0".
#' }
#' \bold{SBC method}
#' \itemize{
#' \item r.a:  a positive constant which is effectively the radius defining a neighborhood. 
#'       The default value is 0.5.
#' \item eps.high: an upper threshold value. The default value is 0.5.
#' \item eps.low: a lower threshold value. The default value is 0.15.
#' \item name: a name for the model. The default value is "sim-0".
#' }
#' \bold{FS.HGD method}
#' \itemize{
#' \item num.labels: a positive integer to determine the number of labels (fuzzy terms). 
#'       The default value is 7.
#' \item max.iter: a positive integer to determine the maximal number of iterations. 
#'       The default value is 100.
#' \item step.size: a real number between 0 and 1. The default value is 0.01.
#' \item alpha.heuristic: a positive real number representing a heuristic value. 
#'       The default value is 1.
#' \item name: a name for the model. The default value is "sim-0".
#' }
#' \bold{FRBCS.W and FRBCS.CHI method}
#' \itemize{
#' \item num.labels: a positive integer to determine the number of labels (fuzzy terms). 
#'       The default value is 7.
#' \item type.mf: the type of the membership function. 
#'       1 is for triangular, while 2, 3, 4, and 5 are for trapezoid, Gaussian, sigmoid, 
#'       and generalized bell. The default value is 1.
#' \item name: a name for the model. The default value is "sim-0".
#' }
#' \bold{DENFIS method}
#' \itemize{
#' \item Dthr: the threshold value for the envolving clustering method (ECM), between 0 and 1. 
#'       The default value is 0.1.
#' \item max.iter: a positive integer to determine the maximal number of iterations. 
#'       The default value is 100.
#' \item step.size: the step size of the least squares method, between 0 and 1. 
#'       The default value is 0.01.
#' \item d: a parameter for the width of the triangular membership function. 
#'       The default value is 2.
#' \item name: a name for the model. The default value is "sim-0".
#' }
#' \bold{GFS.FR.MOGUL method}
#' \itemize{
#' \item persen_cross: a probability of crossover. The default value is 0.6.
#' \item max.iter: a positive integer to determine the maximal number of iterations. 
#'       The default value is 10.
#' \item max.gen: a positive integer to determine the maximal number of generations of the genetic algorithm. 
#'       The default value is 10.
#' \item max.tune: a positive integer to determine the maximal number of tuning iterations.
#'       The default value is 10.
#' \item persen_mutant: a probability of mutation. The default value is 0.3.
#' \item epsilon: a real number between 0 and 1 representing the level of generalization.
#'       A high epsilon can lead to overfitting. The default value is 0.9. 
#' \item name: a name for the model. The default value is "sim-0". 
#' }
#' \bold{GFS.THRIFT method}
#' \itemize{
#' \item popu.size: the size of the population which is generated in each generation. 
#' The default value is 30.
#' \item num.labels: a matrix describing the number of fuzzy terms. The default value is 3.
#' \item persen_cross: a probability of crossover. The default value is 0.6.
#' \item persen_mutant: a probability of mutation. The default value is 0.3.
#' \item max.gen: a positive integer to determine the maximal number of generations for the genetic algorithm. 
#'       The default value is 10.
#' }
#' \bold{GFS.GCCL method}
#' \itemize{
#' \item popu.size: the size of the population which is generated in each generation. 
#' The default value is 30.
#' \item num.labels: a matrix describing the number of fuzzy terms. The default value is 3.
#' \item persen_cross: a probability of crossover. The default value is 0.6.
#' \item persen_mutant: a probability of mutation. The default value is 0.3.
#' \item max.gen: a positive integer to determine the maximal number of generations for the genetic algorithm. 
#'       The default value is 10.
#' }
#' \bold{FH.GBML method}
#' \itemize{
#' \item popu.size: the size of the population which is generated in each generation. 
#'                  The default value is 10.
#' \item max.num.rule: the maximum size of the rules.
#' \item persen_cross: a probability of crossover. The default value is 0.6.
#' \item persen_mutant: a probability of mutation. The default value is 0.3.
#' \item max.gen: a positive integer to determine the maximal number of generations for the genetic algorithm. 
#'       The default value is 10.
#' \item num.class: the number of classes.
#' \item p.dcare: a probability of "don't care" attributes. The default value is 0.5.
#' \item p.gccl: a probability of the GCCL process. The default value is 0.5.
#' }
#' \bold{SLAVE method}
#' \itemize{
#' \item persen_cross: a probability of crossover. The default value is 0.6.
#' \item persen_mutant: a probability of mutation. The default value is 0.3.
#' \item max.iter: a positive integer to determine the maximal number of iterations. 
#'       The default value is 30.
#' \item max.gen: a positive integer to determine the maximal generations of the genetic algorithm. 
#'       The default value is 30.
#' \item num.labels: the number of fuzzy terms. The default value is 3.
#' \item k.lower: a lower bound of the noise threshold with interval between 0 and 1. The default value is 0.
#' \item k.upper: an upper bound of the noise threshold with interval between 0 and 1. The default value is 1.
#' \item epsilon: a value between 0 and 1 representing the level of generalization. A high epsilon can lead to overfitting. 
#'                The default value is 0.5.
#' }
#'
#' @seealso \code{\link{predict}} for the prediction phase, and 
#' the following main functions of each of the methods for theoretical background and references:  \code{\link{WM}}, \code{\link{SBC}},
#' \code{\link{HyFIS}}, \code{\link{ANFIS}}, \code{\link{FIR.DM}}, \code{\link{DENFIS}}, 
#' \code{\link{FS.HGD}}, \code{\link{FRBCS.W}}, \code{\link{FRBCS.CHI}}, \code{\link{GFS.FR.MOGUL}},
#' \code{\link{GFS.Thrift}}, \code{\link{GFS.GCCL}}, \code{\link{FH.GBML}}, and \code{\link{SLAVE}}.
#' @return The \code{\link{frbs-object}}. 
#' @examples
#' ##################################
#' ## I. Regression Problem
#' ## Suppose data have two input variables and one output variable.  
#' ## We separate them into training, fitting, and testing data.
#' ## data.train, data.fit, data.test, and range.data are inputs 
#' ## for all regression methods.
#' ###################################
#' ## Take into account that the simulation might take a long time 
#' ## depending on the hardware you are using. The chosen parameters 
#' ## may not be optimal.
#' ## Data must be in data.frame or matrix form and the last column 
#' ## is the output variable/attribute.
#' ## The training data must be expressed in numbers (numerical data).
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
#' colnames(data.train) <- c("inp.1", "inp.2", "out.1") 
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
#' ## I.1 Example: Implementation of Wang & Mendel
#' #############################################################
#' method.type <- "WM" 
#' 
#' ## collect control parameters into a list
#' ## num.labels = 3 means we define 3 as the number of fuzzy terms
#' ## type.mf = 3 means we use Gaussian as membership function
#' control.WM <- list(num.labels = 3, type.mf = 3, name = "Sim-0") 
#' 
#' ## generate the model and save it as object.WM
#' object.WM <- frbs.learn(data.train, range.data, method.type, control.WM)
#'
#' #############################################################
#' ## I.2 Example: Implementation of SBC
#' #############################################################
#' method.type <- "SBC" 
#' control.SBC <- list(r.a = 0.5, eps.high = 0.5, eps.low = 0.15, name = "Sim-0")
#'
#' \dontrun{object.SBC <- frbs.learn(data.train, range.data, method.type, control.SBC)}
#'
#' #############################################################
#' ## I.3 Example: Implementation of HYFIS
#' #############################################################
#' method.type <- "HYFIS"
#' 
#' control.HYFIS <- list(num.labels = 5, max.iter = 50, step.size = 0.01, 
#'                  name = "Sim-0")
#' 
#' \dontrun{object.HYFIS <- frbs.learn(data.train, range.data, method.type, control.HYFIS)}
#'
#' #############################################################
#' ## I.4 Example: Implementation of ANFIS
#' #############################################################
#' method.type <- "ANFIS" 
#'
#' control.ANFIS <- list(num.labels = 5, max.iter = 100, step.size = 0.01, 
#'                       name = "Sim-0") 
#'
#' \dontrun{object.ANFIS <- frbs.learn(data.train, range.data, method.type, control.ANFIS)}
#'
#' #############################################################
#' ## I.5 Example: Implementation of DENFIS
#' #############################################################
#' 
#' control.DENFIS <- list(Dthr = 0.1, max.iter = 100, step.size = 0.001, d = 2, 
#'                        name = "Sim-0")
#' method.type <- "DENFIS"
#' 
#' \dontrun{object.DENFIS <- frbs.learn(data.train, range.data, method.type, control.DENFIS)}
#'
#' #############################################################
#' ## I.6 Example: Implementation of FIR.DM
#' #############################################################
#' method.type <- "FIR.DM"
#'  
#' control.DM <- list(num.labels = 5, max.iter = 100, step.size = 0.01, name = "Sim-0") 
#' \dontrun{object.DM <- frbs.learn(data.train, range.data, method.type, control.DM)}
#'
#' #############################################################
#' ## I.7 Example: Implementation of FS.HGD
#' #############################################################
#' method.type <- "FS.HGD" 
#'  
#' control.HGD <- list(num.labels = 5, max.iter = 100, step.size = 0.01, 
#'                alpha.heuristic = 1, name = "Sim-0") 
#' \dontrun{object.HGD <- frbs.learn(data.train, range.data, method.type, control.HGD)}
#'
#' #############################################################
#' ## I.8 Example: Implementation of GFS.FR.MOGUL
#' #############################################################
#' method.type <- "GFS.FR.MOGUL" 
#'  
#' control.GFS.FR.MOGUL <- list(persen_cross = 0.6, 
#'                     max.iter = 20, max.gen = 10, max.tune = 10, persen_mutant = 0.3, 
#'                     epsilon = 0.8, name="sim-0") 
#' \dontrun{object.GFS.FR.MOGUL <- frbs.learn(data.train, range.data, 
#'                        method.type, control.GFS.FR.MOGUL)}
#'
#' #############################################################
#' ## I.9 Example: Implementation of Thrift's method (GFS.THRIFT)
#' #############################################################
#' method.type <- "GFS.THRIFT" 
#'  
#' control.Thrift <- list(popu.size = 15, num.labels = 3, persen_cross = 1, 
#'                     max.gen = 5, persen_mutant = 1,
#'                     name="sim-0") 
#' \dontrun{object.Thrift <- frbs.learn(data.fit, range.data, method.type, control.Thrift)}
#'
#' #############################################################
#' ## II. Classification Problem 
#' #############################################################
#' ## The iris dataset is shuffled and divided into training and 
#' ## testing data. Bad results in the predicted values may result
#' ## from casual imbalanced classes in the training data.
#' ## Take into account that the simulation may take a long time 
#' ## depending on the hardware you use. 
#' ## One may get better results with other parameters. 
#' ## Data are in data.frame or matrix form and the last column is 
#' ## the output variable/attribute
#' ## The data must be expressed in numbers (numerical data).
#' 
#' data(iris)
#' irisShuffled <- iris[sample(nrow(iris)),]
#' irisShuffled[,5] <- unclass(irisShuffled[,5])
#' tra.iris <- irisShuffled[1:105,]
#' tst.iris <- irisShuffled[106:nrow(irisShuffled),1:4]
#' real.iris <- matrix(irisShuffled[106:nrow(irisShuffled),5], ncol = 1)
#'
#' ## Please take into account that the interval needed is the range of input data only.
#' range.data.input <- matrix(c(4.3, 7.9, 2.0, 4.4, 1.0, 6.9, 0.1, 2.5), nrow=2)
#' 
#' ######################################################### 
#' ## II.1 Example: Implementation of FRBCS with weighted factor based on Ishibuchi's method
#' ###############################################################
#' ## generate the model
#' method.type <- "FRBCS.W"
#' control <- list(num.labels = 7, type.mf = 1) 
#' 
#' \dontrun{object <- frbs.learn(tra.iris, range.data.input, method.type, control)}
#' 
#' ## conduct the prediction process
#' \dontrun{res.test <- predict(object, tst.iris)}
#'
#' ######################################################### 
#' ## II.2 Example: Implementation of FRBCS based on Chi's method
#' ###############################################################
#' ## generate the model
#' method.type <- "FRBCS.CHI"
#' control <- list(num.labels = 7, type.mf = 1) 
#' 
#' \dontrun{object <- frbs.learn(tra.iris, range.data.input, method.type, control)}
#' 
#' ## conduct the prediction process
#' \dontrun{res.test <- predict(object, tst.iris)}
#'
#' ######################################################### 
#' ## II.3 The example: Implementation of GFS.GCCL
#' ###############################################################
#' method.type <- "GFS.GCCL" 
#' 
#' control <- list(popu.size = 30, num.class = 3, num.labels = 5, persen_cross = 0.9, 
#'                     max.gen = 200, persen_mutant = 0.3,
#'                     name="sim-0") 
#' ## Training process
#' ## The main result of the training is a rule database which is used later for prediction.
#' \dontrun{object <- frbs.learn(tra.iris, range.data.input, method.type, control)}
#' ## Prediction process
#' \dontrun{res.test <- predict(object, tst.iris)}
#'
#' ######################################################### 
#' ## II.4 Example: Implementation of FH.GBML
#' ###############################################################
#' method.type <- "FH.GBML" 
#'	 
#'	control <- list(popu.size = 10, max.num.rule = 50, num.class = 3, 
#'				persen_cross = 0.9, max.gen = 200, persen_mutant = 0.3, p.dcare = 0.5, 
#'              p.gccl = 1, name="sim-0") 
#'	 
#'	## Training process
#'	## The main result of the training is a rule database which is used later for prediction.
#'	\dontrun{object <- frbs.learn(tra.iris, range.data.input, method.type, control)}
#'
#'	## Prediction process
#'	\dontrun{res.test <- predict(object, tst.iris)}
#'
#' ######################################################### 
#' ## II.5 The example: Implementation of SLAVE
#' ###############################################################
#' method.type <- "SLAVE" 
#'	 
#'	control <- list(num.class = 3, num.labels = 5,
#'				persen_cross = 0.9, max.iter = 50, max.gen = 30, persen_mutant = 0.3, 
#'              k.lower = 0.25, k.upper = 0.75, epsilon = 0.1, name="sim-0") 
#'	 
#'	## Training process
#'	## The main result of the training is a rule database which is used later for prediction.
#'	\dontrun{object <- frbs.learn(tra.iris, range.data.input, method.type, control)}
#'
#'	## Prediction process
#'	\dontrun{res.test <- predict(object, tst.iris)}
#' @export
frbs.learn <- function(data.train, range.data = NULL, method.type = c("WM"), control=list()){
 
## get type of method 
method.type <- toupper(method.type)

## get names of variables
colnames.var <- colnames(data.train)

## initialize mod
mod <- NULL

## condition if data.train is in data frame type
if (class(data.train) != "matrix"){
	data.train <- as.matrix(data.train)
}

## if user do not give range of data, calculate from data
if (is.null(range.data)){
	if (any(method.type == c("FRBCS.W", "FRBCS.CHI", "GFS.GCCL", "FH.GBML", "SLAVE"))){
		dt.min <- matrix(do.call(pmin, lapply(1:nrow(data.train[, -ncol(data.train), drop = FALSE]), function(i)data.train[i, -ncol(data.train), drop = FALSE])), nrow = 1)
		dt.max <- matrix(do.call(pmax, lapply(1:nrow(data.train[, -ncol(data.train), drop = FALSE]), function(i)data.train[i, -ncol(data.train), drop = FALSE])), nrow = 1)
	}
	else {
		dt.min <- matrix(do.call(pmin, lapply(1:nrow(data.train), function(i)data.train[i,])), nrow = 1)
		dt.max <- matrix(do.call(pmax, lapply(1:nrow(data.train), function(i)data.train[i,])), nrow = 1)		
	}
	range.data <- rbind(dt.min, dt.max)
}    

## 1. wang & mendel methods
if(method.type == "WM"){

	## getting all of parameters
	control <- setDefaultParametersIfMissing(control, list(num.labels = 7, type.mf = 3, type.defuz = 1, type.tnorm = 1, type.snorm = 1, name="sim-0"))
	
	## get parameters
	num.labels <- control$num.labels
	type.mf <- control$type.mf
	name <- control$name
	num.labels <- matrix(rep(num.labels, ncol(range.data)), nrow=1)

	## generate FRBS model
	modelSpecific <- WM(range.data, data.train, num.labels, type.mf)
	
	## collect results as model
	mod <- modelSpecific
	mod$type.model <- 1
	mod$func.tsk <- NULL
	mod$type.defuz <- control$type.defuz
	mod$type.tnorm <- control$type.tnorm
	mod$type.snorm <- control$type.snorm
	
}

## 2. Subtractive Clustering
else if(method.type == "SBC"){
	
	## get all of parameters
	control <- setDefaultParametersIfMissing(control, list(r.a = 0.5, eps.high = 0.5, eps.low = 0.15, name ="sim-0"))
	r.a <- control$r.a
	eps.high <- control$eps.high
	eps.low <- control$eps.low
	name <- control$name
	range.data.ori <- range.data
	
	## generate FRBS model 
	modelSpecific <- SBC(data.train, range.data.ori, r.a, eps.high, eps.low)
	mod <- modelSpecific
}

## 3. HyFIS
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
	
	## generate FRBS model
	modelSpecific <- HyFIS(range.data.norm, data.tra.norm, num.labels, max.iter, range.data.ori, step.size)
	mod <- modelSpecific

}

## 4. ANFIS
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

	## generate FRBS model
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

## 5. Fuzzy Rule-based Classification System
else if (method.type == "FRBCS.W"){

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
	data.tra.norm <- norm.data(data.train[, 1 : (ncol(data.train) - 1)], range.data.ori, min.scale = 0, max.scale = 1)
	data.train <- cbind(data.tra.norm, matrix(data.train[, ncol(data.train)], ncol = 1))
	
	## generate FRBS model
	modelSpecific <- FRBCS.W(range.data.inout, data.train, n.labels, num.class, type.mf)
	mod <- modelSpecific
	mod$range.data.ori <- range.data.ori
}

## 6. DENFIS 
else if (method.type == "DENFIS"){	
	## get all of parameters
    control <- setDefaultParametersIfMissing(control, list(Dthr = 0.1, max.iter = 100, step.size = 0.01, d = 2, name="sim-0"))
	Dthr <- control$Dthr
	max.iter <- control$max.iter
	step.size <- control$step.size
	d <- control$d
	name <- control$name
	range.data.ori <- range.data
	
	## generate FRBS model
	modelSpecific <- DENFIS(data.train, range.data.ori, Dthr, max.iter, step.size, d)
	mod <- modelSpecific
}

## 7. FIR.DM
else if (method.type == "FIR.DM"){
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

	## generate FRBS model
	modelSpecific <- FIR.DM(range.data.norm, data.tra.norm, num.labels, max.iter, step.size)
	mod <- modelSpecific
	mod$range.data.ori <- range.data.ori
}

## 8. FS.HGD
else if (method.type == "FS.HGD"){

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

	## generate FRBS model
	modelSpecific <- FS.HGD(range.data.norm, data.tra.norm, num.labels, max.iter, step.size, alpha.heuristic)
	mod <- modelSpecific
	mod$range.data.ori <- range.data.ori
}

## 9. GFS.FR.MOGUL
else if (method.type == "GFS.FR.MOGUL"){
	## get all of parameters
	control <- setDefaultParametersIfMissing(control, list(persen_cross = 0.6, max.iter = 10, 
	                 max.gen = 10, max.tune = 10, persen_mutant = 0.3, epsilon = 0.8, name="sim-0"))
	
	## getting all of parameters
	range.data.ori <- range.data
	data.train.ori <- data.train
	persen_cross <- control$persen_cross
	persen_mutant <- control$persen_mutant
	max.iter <- control$max.iter
	max.gen <- control$max.gen
	max.tune <- control$max.tune
	epsilon <- control$epsilon
	name <- control$name
	
	## normalize range of data and data training
	range.data.norm <- range.data.ori
	range.data.norm[1, ] <- 0
	range.data.norm[2, ] <- 1	
	data.tra.norm <- norm.data(data.train.ori, range.data.ori, min.scale = 0, max.scale = 1)
		
	modelSpecific <- GFS.FR.MOGUL(data.tra.norm, persen_cross, persen_mutant, 
	                             max.iter, max.gen, max.tune, range.data.ori, epsilon)
	mod <- modelSpecific
}

## 10. GFS.Thrift
else if (method.type == "GFS.THRIFT"){
	## get all of parameters
	control <- setDefaultParametersIfMissing(control, list(popu.size = 30, num.labels = 3, persen_cross = 0.6, max.gen = 10, 
	                 persen_mutant = 0.3, name="sim-0"))
	
	## getting all of parameters
	range.data.ori <- range.data
	data.train.ori <- data.train
	popu.size <- control$popu.size
	persen_cross <- control$persen_cross
	persen_mutant <- control$persen_mutant
	max.gen <- control$max.gen
	name <- control$name
	n.labels <- control$num.labels
	
	## normalize range of data and data training
	range.data.norm <- range.data.ori
	range.data.norm[1, ] <- 0
	range.data.norm[2, ] <- 1	
	data.tra.norm <- norm.data(data.train.ori, range.data.ori, min.scale = 0, max.scale = 1)
	num.labels <- matrix(rep(n.labels, ncol(range.data.norm)), nrow=1)
	
	## generate FRBS model
	modelSpecific <- GFS.Thrift(data.tra.norm, popu.size, range.data.norm, num.labels, persen_cross, persen_mutant, max.gen, range.data.ori)
	mod <- modelSpecific
}

## 11. GFS.GCCL
else if (method.type == "GFS.GCCL"){
	## get all of parameters
	control <- setDefaultParametersIfMissing(control, list(popu.size = 30, num.class = 2, num.labels = 3, persen_cross = 0.6, 
	                   max.gen = 10, persen_mutant = 0.3, name="sim-0"))
	
	## getting all of parameters
	range.data.input <- range.data
	data.train.ori <- data.train
	popu.size <- control$popu.size
	persen_cross <- control$persen_cross
	persen_mutant <- control$persen_mutant
	max.gen <- control$max.gen
	name <- control$name
	n.labels <- control$num.labels
	n.class <- control$num.class
		
	num.labels <- matrix(rep(n.labels, ncol(range.data)), nrow = 1)
	num.labels <- cbind(num.labels, n.class)
	
	## normalize range of data and data training
	range.data.norm <- range.data.input
	range.data.norm[1, ] <- 0
	range.data.norm[2, ] <- 1	
	range.data.input.ori <- range.data.input
	data.tra.norm <- norm.data(data.train[, 1 : ncol(data.train) - 1], range.data.input, min.scale = 0, max.scale = 1)
	data.train <- cbind(data.tra.norm, matrix(data.train[, ncol(data.train)], ncol = 1))
	
	## generate FRBS model
	modelSpecific <- GFS.GCCL(data.train, popu.size, range.data.norm, num.labels, persen_cross, persen_mutant, max.gen, range.data.input.ori)
	
	mod <- modelSpecific
}

## 12. Hibridization of fuzzy GBML (FH.GBML)
else if(method.type == "FH.GBML"){
	## get all of parameters
	control <- setDefaultParametersIfMissing(control, list(popu.size = 10, max.num.rule = 5, num.class = 3, persen_cross = 0.6, 
	                               max.gen = 10, persen_mutant = 0.3, p.dcare = 0.5, p.gccl = 0.5, name="sim-0"))
	
	## getting all of parameters
	range.data.input <- range.data
	data.train.ori <- data.train
	popu.size <- control$popu.size
	persen_cross <- control$persen_cross
	persen_mutant <- control$persen_mutant
	max.gen <- control$max.gen
	name <- control$name
	num.class <- control$num.class
	max.num.rule <- control$max.num.rule
	p.dcare <- control$p.dcare
	p.gccl <- control$p.gccl
		
	## normalize data.train excluded output attribute
	data.tra.norm <- norm.data(data.train.ori[, -ncol(data.train.ori), drop = FALSE], range.data.input, min.scale = 0, max.scale = 1)
	data.tra.norm <- cbind(data.tra.norm, data.train.ori[, ncol(data.train.ori), drop = FALSE])
	
	## generate FRBS model
	modelSpecific <- FH.GBML(data.tra.norm, popu.size, max.num.rule, persen_cross, persen_mutant, max.gen, 
	                         num.class, range.data.input, p.dcare, p.gccl)	
	mod <- modelSpecific
}

## 13. SLAVE
else if(method.type == "SLAVE"){
	## get all of parameters
	control <- setDefaultParametersIfMissing(control, list(num.class = 3, num.labels = 3, persen_cross = 0.6, 
				max.iter = 30, max.gen = 30, persen_mutant = 0.3, k.lower = 0, k.upper = 1, epsilon = 0.5, name="sim-0"))
	
	## getting all of parameters
	range.data.input <- range.data
	data.train.ori <- data.train
	persen_cross <- control$persen_cross
	persen_mutant <- control$persen_mutant
	max.iter <- control$max.iter
	max.gen <- control$max.gen
	name <- control$name
	num.class <- control$num.class
	num.labels <- control$num.labels
	k.lower <- control$k.lower
	k.upper <- control$k.upper
	epsilon <- control$epsilon
	
	num.labels <- matrix(rep(num.labels, ncol(range.data)), nrow = 1)
	num.labels <- cbind(num.labels, num.class)
		
	## normalize data.train excluded output attribute
	data.tra.norm <- norm.data(data.train.ori[, -ncol(data.train.ori), drop = FALSE], 
	                           range.data.input, min.scale = 0, max.scale = 1)
	data.tra.norm <- cbind(data.tra.norm, data.train.ori[, ncol(data.train.ori), drop = FALSE])
	
	## generate FRBS model
	modelSpecific <- SLAVE(data.tra.norm, persen_cross, persen_mutant, max.iter, max.gen, num.labels, range.data.input, k.lower, k.upper, epsilon)	
	mod <- modelSpecific
}

else if (method.type == "FRBCS.CHI"){
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
	modelSpecific <- FRBCS.CHI(range.data.inout, data.train, n.labels, num.class, type.mf)
	mod <- modelSpecific
	mod$range.data.ori <- range.data.ori

}
mod$method.type <- method.type
mod$name <- name

## keep colnames of training data into mod
if (!is.null(colnames.var)) {
	mod$colnames.var <- colnames.var 
}
else {
	mod$colnames.var <- paste("var", seq(1, ncol(data.train)), sep = ".")
}

mod <- frbsObjectFactory(mod)

return(mod)
}

## checking missing parameters
# @param control parameter values of each method
# @param defaults default parameter values of each method
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
#' \item{num.labels}{the number of fuzzy terms for the variables}
#' \item{varout.mf}{a matrix to generate the shapes of the membership functions for the output variable. 
#'          The first row represents the shape of the membership functions, the other rows contain the parameters that have been generated. 
#'          Whether the values of parameters within the matrix are normalized to lie between 0 and 1 or not depends on the selected method.}
#' \item{names.varoutput}{generated names of the output variable.}
#' \item{rule}{the fuzzy IF-THEN rules; In the GFS.FR.MOGUL case, a rule refers to the parameter values of the membership function 
#'          which represents the rule.}
#' \item{rule.data.num}{the fuzzy IF-THEN rules in integer format.}
#' \item{varinp.mf}{a matrix to generate the shapes of the membership functions for the input variables. 
#'           The first row represents the shape of the membership functions, 
#'           the other rows contain the non NA values representing the parameters related with their type of membership function. 
#'           For example, trapezoid, triangular, and Gaussian have four, three, and two values as their parameters, respectively. 
#'           Whether the values of parameters within the matrix are normalized to lie between 0 and 1 or not depends on the selected method.}
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
#' \item{rule.data.num}{a matrix representing the rules in integer form.}
#' \item{grade.cert}{grade of certainty for classification problems.}
#' \item{alpha.heuristic}{a parameter for the heuristic of the HGD method.}
#' \item{colnames.var}{the names of variables.}
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

#' This is the main function to obtain a final result as predicted values for all methods in this package. 
#' In order to get predicted values, this function is run using an \code{\link{frbs-object}}, which is typically generated using \code{\link{frbs.learn}}.
#' 
#' @title The frbs prediction stage
#'
#' @param object an \code{\link{frbs-object}}.
#' @param newdata a data frame or matrix (m x n) of data for the prediction process, where m is the number of instances and 
#' n is the number of input variables. It should be noted that the testing data must be expressed in numbers (numerical data).
#' @param ... the other parameters (not used)
#' @seealso \code{\link{frbs.learn}} and \code{\link{frbs.gen}} for learning and model generation, 
#' and the internal main functions of each method for the theory:  
#' \code{\link{WM}}, \code{\link{SBC}}, \code{\link{HyFIS}}, \code{\link{ANFIS}}, 
#' \code{\link{FIR.DM}}, \code{\link{DENFIS}}, \code{\link{FS.HGD}}, \code{\link{FRBCS.W}}, 
#' \code{\link{GFS.FR.MOGUL}}, \code{\link{GFS.Thrift}}, \code{\link{GFS.GCCL}}, \code{\link{FRBCS.CHI}}, 
#' \code{\link{FH.GBML}}, and \code{\link{SLAVE}}.
#' @return The predicted values. 
#' @aliases predict
#' @examples
#' ##################################
#' ## I. Regression Problem
#' ###################################
#' ## In this example, we just show how to predict using Wang and Mendel's technique but
#' ## users can do it in the same way for other methods.
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
#' ## I.1 Example: Implementation of Wang & Mendel
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
#' ## The following code can be used for all methods
#' res <- predict(object, newdata) 
#' 
#' @export  
#' @method predict frbs
#' @S3method predict frbs
predict.frbs <- function(object, newdata, ...) {
mod <- object

if(!inherits(mod, "frbs")) stop("not a legitimate frbs model")
  
if (class(newdata) != "matrix"){
	newdata <- as.matrix(newdata)
}  
##############################
## Split data of frbs.learn
#############################
m.type <- mod$method.type

## 1. WM approach
if (m.type == "WM") {
	
	res.comp <- frbs.eng(mod, newdata)
	res <- res.comp$predicted.val
}

## 2.SBC approach
else if(m.type == "SBC"){
	
	res <- SBC.test(mod, newdata)
}

## 3. HyFIS approach
else if(m.type == "HYFIS"){
	range.data.ori <- mod$range.data.ori
	range.input.ori <- range.data.ori[, 1:(ncol(range.data.ori) - 1)]
	range.output.ori <- range.data.ori[, ncol(range.data.ori), drop = FALSE]

	data.tst.norm <- norm.data(newdata, range.input.ori, min.scale = 0, max.scale = 1)
	res.comp <- frbs.eng(mod, data.tst.norm)
	res.denorm <- denorm.data(res.comp$predicted.val, range.output.ori, min.scale = 0, max.scale = 1)
	res <- res.denorm
}

## 4. ANFIS approach
else if(m.type == "ANFIS"){
	range.data.ori <- mod$range.data.ori
	range.input.ori <- range.data.ori[, 1:(ncol(range.data.ori) - 1)]
	range.output.ori <- range.data.ori[, ncol(range.data.ori), drop = FALSE]

	data.tst.norm <- norm.data(newdata, range.input.ori, min.scale = 0, max.scale = 1)
	res.comp <- frbs.eng(mod, data.tst.norm)
	
	res.denorm <- denorm.data(res.comp$predicted.val, range.output.ori, min.scale = 0, max.scale = 1)
	res <- res.denorm
}

## 5. FRBCS.W approach
else if(m.type == "FRBCS.W"){
	data.tst.norm <- norm.data(newdata, mod$range.data.ori, min.scale = 0, max.scale = 1)
	res <- FRBCS.eng(mod, data.tst.norm)
}

## 6. DENFIS approach
else if(m.type == "DENFIS"){
	res <- DENFIS.eng(mod, newdata)
}

## 7. FIR.DM approach
else if(m.type == "FIR.DM"){
	range.data.ori <- mod$range.data.ori
	
	range.input.ori <- range.data.ori[, 1:(ncol(range.data.ori) - 1), drop = FALSE]
	range.output.ori <- range.data.ori[, ncol(range.data.ori), drop = FALSE]

	data.tst.norm <- norm.data(newdata, range.input.ori, min.scale = 0, max.scale = 1)
	
	res.comp <- frbs.eng(mod, data.tst.norm)
	
	res.denorm <- denorm.data(res.comp$predicted.val, range.output.ori, min.scale = 0, max.scale = 1)
	res <- res.denorm
}

## 8. FS.HGD approach
else if(m.type == "FS.HGD"){
	
	range.data.ori <- mod$range.data.ori
	range.input.ori <- range.data.ori[, 1:(ncol(range.data.ori) - 1), drop = FALSE]
	range.output.ori <- range.data.ori[, ncol(range.data.ori), drop = FALSE]

	data.tst.norm <- norm.data(newdata, range.input.ori, min.scale = 0, max.scale = 1)
	
	res.comp <- frbs.eng(mod, data.tst.norm)
	
	res.denorm <- denorm.data(res.comp$predicted.val, range.output.ori, min.scale = 0, max.scale = 1)
	res <- res.denorm
}

## 9. GFS.FR.MOGUL approach
else if(m.type == "GFS.FR.MOGUL"){
	
	range.data.ori <- mod$range.data.ori
	range.input.ori <- range.data.ori[, 1:(ncol(range.data.ori) - 1), drop = FALSE]
	range.output.ori <- range.data.ori[, ncol(range.data.ori), drop = FALSE]

	data.tst.norm <- norm.data(newdata, range.input.ori, min.scale = 0, max.scale = 1)
	
	res.comp <- GFS.FR.MOGUL.test(mod, data.tst.norm)
	
	res.denorm <- denorm.data(res.comp, range.output.ori, min.scale = 0, max.scale = 1)
	res <- res.denorm
}

## 10. GFS.THRIFT
else if(m.type == "GFS.THRIFT"){
	
	range.data.ori <- mod$range.data.ori
	range.input.ori <- range.data.ori[, 1:(ncol(range.data.ori) - 1), drop = FALSE]
	range.output.ori <- range.data.ori[, ncol(range.data.ori), drop = FALSE]

	data.tst.norm <- norm.data(newdata, range.input.ori, min.scale = 0, max.scale = 1)

	res.comp <-	GFS.Thrift.test(mod, data.tst.norm)
	
	res.denorm <- denorm.data(res.comp, range.output.ori, min.scale = 0, max.scale = 1)
	res <- res.denorm
}

## 11. GFS.GCCL
else if(m.type == "GFS.GCCL"){
	newdata <- norm.data(newdata, mod$range.data.ori, min.scale = 0, max.scale = 1)
	res <- GFS.GCCL.eng(mod, newdata)
}

## 12. FH.GBML
else if(m.type == "FH.GBML"){
	data.tst.norm <- norm.data(newdata, mod$range.data.ori, min.scale = 0, max.scale = 1)
	res <- GFS.GCCL.eng(mod, data.tst.norm)
}

## 13. SLAVE
else if (m.type == "SLAVE"){
	data.tst.norm <- norm.data(newdata, mod$range.data.ori, min.scale = 0, max.scale = 1)
	res <- SLAVE.test(mod, data.tst.norm)
}

## 14. Chi's technique
else if(m.type == "FRBCS.CHI"){
	range.data.ori <- mod$range.data.ori
	
	data.tst.norm <- norm.data(newdata, range.data.ori, min.scale = 0, max.scale = 1)
	res <- FRBCS.eng(mod, data.tst.norm)
}


return(res)
}

#' This function enables the output of a summary of the \code{\link{frbs-object}}. 
#'
#' This function displays several components of the object. The components of 
#' one particular method can be different from components of other methods.
#' The following is a description of all components which might be printed.
#' \itemize{
#' \item The name of the model: A name given by the user representing the name of the simulation
#'                or data or model.
#' \item Model was trained using: It shows which method we have been used.
#' \item The interval of training data: It is a matrix representing the original 
#'            interval of data where the first and second rows are minimum and maximum of data, 
#'            respectively. The number of columns represents the number of variables.
#' \item The number of fuzzy terms of the input variables: Given as
#'            elements of a matrix.
#' \item The names of fuzzy terms of the input variables: These names are generated
#'           automatically by frbs expressing all fuzzy terms considered. 
#'           These names are built by two parts which are the name of variables expressed 
#'           by "v" and the name of fuzzy labels of each variables represented by "a". 
#'           For example, "v.1-a.1" means the fuzzy label "a.1" of the first variable (v.1).        
#' \item The names of fuzzy terms of the output variable: For the Mamdani model, since the frbs package only considers
#'           single output, the names of the fuzzy terms for the output variable 
#'           are simple and clear and start with "c". However, for Takagi Sugeno Kang model and
#'           fuzzy rule-based classification systems, this component is always NULL.
#' \item The parameter values of membership functions of the input variables (normalized):
#'          It is represented by a matrix (5 x n) where n depends on the number of 
#'          fuzzy terms on the input variables and the first row of the matrix describes 
#'          a type of membership function, and the rest of rows are 
#'          their parameter values. For example, label "v.1-a.2" has value 
#'          {4.0, 0.23, 0.43, 0.53, 0.73} on its column. It means that the label a.2 of variable v.1 
#'           has a parameter as follows. 
#'          4.0 on the first row shows trapezoid shape in the middle position, 
#'          while 0.23, 0.43, 0.53, and 0.73 are corner points of a trapezoid. 
#'          Furthermore, the following is the complete list of shapes of membership functions:
#'          \itemize{
#'          \item Triangular: 1 on the first row and rows 2, 3, and 4 represent corner points. 
#'          \item Trapezoid: 2, 3, or 4 on the first row means they are trapezoid in left, right and middle side, respectively,
#'                     and rows 2, 3, 4, and 5 represent corner points. But for trapezoid at left or right side the fifth row is NA. 
#'          \item Gaussian: 5 on the first row means it uses Gaussian and second and third row represent mean and variance.
#'          \item Sigmoid: 6 on the first row and two parameters (gamma and c) on second and third rows.
#'          \item Generalized bell: 7 on the first row and three parameters (a, b, c) on second, third, and fourth rows.
#'          }
#' \item The fuzzy IF-THEN rules: In this package, there are several models for representing
#'          fuzzy IF-THEN rules based on the method used. 
#'          \itemize{
#'          \item Mamdani model: they are represented as a knowledge base containing two parts: 
#'          antecedent and consequent parts which are separated by a sign "->", as for example in the
#'          following rule:
#' 
#'          \code{var.1 is v.1-a.1 and var.2 is v.2-a.2 -> var.3 is c.2}
#'          
#'          \item Takagi Sugeno Kang model: In this model, this component only represents the antecedent
#'          of rules while the consequent part will be represented by linear equations. 
#'          \item fuzzy rule-based classification systems: This model is quite similar to the Mamdani model,
#'          but the consequent part expresses pre-defined classes instead of linguistic values.
#'          \item approximate approach: Especially for GFS.FR.MOGUL, a matrix of parameters
#'          of membership functions is used to represent the fuzzy IF-THEN rules as well.  
#'          The representation of rules and membership functions is a matrix (n x (p x m)) where
#'          n is the number of rules and m is the number of variables while p is the number of corner points 
#'          of the membership function, if we are using triangular or trapezoid then p = 3 or 4, respectively. 
#'          For example, let us consider the triangular membership function and a number of variables of 3. 
#'          The representation of rules and membership functions is as follows:
#'          <<a11 a12 a13>> <<b11 b12 b13>> <<c11 c12 c13>>. 
#'          
#'          }
#' \item The linear equations on consequent parts of fuzzy IF-THEN rules: It is used in
#'         the Takagi Sugeno Kang model.
#' \item The weight of the rules or the certainty factor: For the FRBCS.W method, this shows the weight related to the rules 
#'         representing the ratio of dominance among the rules.
#' \item The cluster centers: This component is used in clustering methods representing cluster centers.
#' }
#' @title The summary function for frbs objects
#' 
#' @param object the \code{\link{frbs-object}}
#' @param ... the other parameters (not used)
#' @export  
#' @method summary frbs
#' @S3method summary frbs
summary.frbs <- function(object, ...){
	
  if(!inherits(object, "frbs")) stop("not a legitimate frbs model")
	cat("The name of model: ", object$name, "\n")
	cat("Model was trained using: ", object$method.type, "\n") 
	
	if (any(object$method.type == c("WM", "GFS.THRIFT"))){
		range.data.ori <- object$range.data.ori
		colnames(range.data.ori) <- object$colnames.var
		rownames(range.data.ori) <- c("min", "max")
		cat("The interval of training data: ", "\n")
		print(range.data.ori)
		if (object$method.type == "WM"){
			cat("Type of defuzzification technique:", "\n")
			if (object$type.defuz == 1) print("Weighted average method")
			else if (object$type.defuz == 2) print("first of maxima")
			else if (object$type.defuz == 3) print("last of maxima")
			else if (object$type.defuz == 4) print("mean of maxima")
			else print("modified COG")
			cat("Type of t-norm method:", "\n")
			if (object$type.tnorm == 1) print("Standard t-norm")
			else if (object$type.tnorm == 2) print("Hamacher product")
			else if (object$type.tnorm == 3) print("Yager class (with tao = 1)")
			else if (object$type.tnorm == 4) print("Product")
			else print("Bounded product")
			cat("Type of s-norm method:", "\n")
			if (object$type.snorm == 1) print("Standard s-norm")
			else if (object$type.snorm == 2) print("Hamacher sum")
			else if (object$type.snorm == 3) print("Yager class (with tao = 1)")
			else if (object$type.snorm == 4) print("Sum")
			else print("Bounded sum")
		}
	}
	else if (any(object$method.type == c("GFS.GCCL", "FH.GBML", "SLAVE", "FRBCS.W", "FRBCS.CHI"))){
		colnames(object$range.data.ori) <- object$colnames.var[-length(object$colnames.var)]
		rownames(object$range.data.ori) <- c("min", "max")
		cat("The interval of input data: ", "\n")
		print(object$range.data.ori)
	}	
	else {
		colnames(object$range.data.ori) <- object$colnames.var
		rownames(object$range.data.ori) <- c("min", "max")
		cat("The interval of training data: ", "\n")
		print(object$range.data.ori)
	}
		
	if (any(object$method.type == c("DENFIS", "SBC"))){
		colnames(object$cls) <- object$colnames.var
		cat("The cluster centers: ", "\n")
		print(object$cls)
	}
	else if (any(object$method.type == c("WM", "HYFIS", "GFS.THRIFT"))){
		num.labels <- cbind(object$num.fvalinput, length(object$names.varoutput))
		rule <- rep.rule(object)
		colnames(num.labels) <- object$colnames.var
		cat("The number of fuzzy terms on each variables", "\n")
		print(num.labels)
		cat("The names of fuzzy terms on the input variables: ", "\n")
		print(object$names.varinput)	
		cat("The names of fuzzy terms on the output variable: ", "\n")
		print(object$names.varoutput)
		cat("The parameter values of membership function on input variables (it might be normalized): ", "\n")
		print(object$varinp.mf)
		cat("The parameter values of membership function on the output variable (it might be normalized): ", "\n")
		print(object$varout.mf)		
		cat("The fuzzy IF-THEN rules: ", "\n")
		print(rule)
	}		
	else if (any(object$method.type == c("ANFIS", "FIR.DM", "FS.HGD", "FRBCS.W", "FRBCS.CHI"))){
		num.labels <- cbind(object$num.fvalinput, object$num.fvalinput[1,1])
		rule <- rep.rule(object)
		cat("The number of fuzzy terms on the input variables", "\n")
		num.labels.input <- num.labels[1, -ncol(num.labels), drop = FALSE]
		colnames(num.labels.input) <- object$colnames.var[-length(object$colnames.var)]
		print(num.labels.input)
		cat("The names of fuzzy terms on the input variables: ", "\n")
		print(object$names.varinput)
		cat("The names of fuzzy terms on the output variable: ", "\n")
		print(NULL)		
		cat("The parameter values of membership function on input variables (normalized): ", "\n")
		print(object$varinp.mf)
		if (any(object$method.type == c("FRBCS.W", "FRBCS.CHI"))){
			cat("The fuzzy IF-THEN rules: ", "\n")
			print(rule)
			if (any(object$method.type == c("FRBCS.W"))){
				cat("The weight of the rules", "\n")
				print(object$grade.cert[, 2, drop = FALSE])
			}
		}
		else {
			cat("The fuzzy IF-THEN rules: ", "\n")
			print(rule)
			if (ncol(object$func.tsk) > 1){	
				seq.deg <- seq(from = 1, to = (ncol(object$func.tsk) - 1), by = 1)
				coef.var <- paste("var", seq.deg, sep = ".")
				names.func <- c(coef.var, "const")	
			}
			else {
				names.func <- c("const")	
			}
			colnames(object$func.tsk) <- names.func
			cat("The linear equations on consequent parts of fuzzy IF-THEN rules: ", "\n")
			print(object$func.tsk)
		}
	}
	else if (object$method.type == "GFS.FR.MOGUL"){
		cat("The parameter values of membership functions representing the fuzzy IF-THEN rules: ", "\n")
		print(object$rule)
	}
	else if (any(object$method.type == c("GFS.GCCL", "FH.GBML"))){
		rule.data.num.inp <- object$rule[, -ncol(object$rule), drop = FALSE]
		rule.data.str <- cbind(rule.data.num.inp, object$rule[, ncol(object$rule), drop = FALSE])
		res <- generate.rule(rule.data.str, object$num.labels)
		rule <- rep.rule(object)
		names.inp.var <- c("don't_care", res$names.varinput)
		colnames(object$num.labels) <- object$colnames.var
		cat("The number of fuzzy terms on each variables", "\n")
		print(object$num.labels)
		cat("The names of fuzzy terms on the input variables: ", "\n")
		print(names.inp.var)		
		cat("The parameter values of membership function on input variables: ", "\n")
		print(object$varinp.mf)
		cat("The fuzzy IF-THEN rules: ", "\n")
		print(rule)
		cat("The certainty factor:", "\n")
		print(object$grade.cert)

	}
	else if (object$method.type == "SLAVE"){			
		res <- generate.rule(object$rule, object$num.labels)
		rule <- rep.rule(object)
		colnames(object$num.labels) <- object$colnames.var
		cat("The number of fuzzy terms on each variables", "\n")
		print(object$num.labels)
		cat("The names of fuzzy terms on the input variables: ", "\n")
		print(res$names.varinput)	
		cat("The names of fuzzy terms on the output variable: ", "\n")
		print(res$names.varoutput)	
		cat("The parameter values of membership function on input variables: ", "\n")
		print(object$varinp.mf)
		cat("The fuzzy IF-THEN rules: ", "\n")
		print(rule)		
	}
	
  invisible(object)	
}

#' This function can be used to plot the shapes of the membership functions.
#'
#' @title The plotting function
#' 
#' @param object an \code{\link{frbs-object}} or a list of parameters to plot membership functions when we build the frbs model without learning. 
#'        There are several parameters that must be inserted in params as follows.
#'        \itemize{
#'        \item var.mf: a matrix of membership function of input and output variables. Please see \code{\link{fuzzifier}}.
#'        \item range.data: a matrix(2 x n) containing the range of the data, where n is the number of variables, and
#'        first and second rows are the minimum and maximum values, respectively. 
#'        \item num.labels: the number of fuzzy terms of the input and output variables. 
#' 
#'              For example: \code{num.labels <- matrix(c(3, 3, 3), nrow = 1)}
#'
#'              It means we have 3 fuzzy values/labels for two input variables and one output variable.
#'        \item names.variables: a list of names of variables. 
#'
#'              For example: names.variables <- c("input1", "input2", "output1")
#'        }
#' @examples
#' ## The following examples contain two different cases which are
#' ## using an frbs-object and the manual way. 
#' ## 
#' ## 1. Plotting using frbs object.
#' data(iris)
#' irisShuffled <- iris[sample(nrow(iris)),]
#' irisShuffled[,5] <- unclass(irisShuffled[,5])
#' tra.iris <- irisShuffled[1:105,]
#' tst.iris <- irisShuffled[106:nrow(irisShuffled),1:4]
#' real.iris <- matrix(irisShuffled[106:nrow(irisShuffled),5], ncol = 1)
#'
#' ## Please take into account that the interval needed is the range of input data only.
#' range.data.input <- matrix(c(4.3, 7.9, 2.0, 4.4, 1.0, 6.9, 0.1, 2.5), nrow=2)
#' 
#' ## generate the model
#' method.type <- "FRBCS.W"
#' control <- list(num.labels = 7, type.mf = 1) 
#' \dontrun{object <- frbs.learn(tra.iris, range.data.input, method.type, control)} 
#' 
#' ## plot the frbs object
#' \dontrun{plotMF(object)}
#'
#' ## 2. Plotting using params.
#' ## Define shape and parameters of membership functions of input variables.
#' ## Please see the fuzzifier function of how to contruct the matrix.
#' varinp.mf <- matrix(c(2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA,
#'                       2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA,
#'                       2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA,
#'                       2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA),
#'                       nrow = 5, byrow = FALSE)
#' ## Define the shapes and parameters of the membership functions of the output variables.
#' varout.mf <- matrix(c(2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA),
#'                       nrow = 5, byrow = FALSE)
#' var.mf <- cbind(varinp.mf, varout.mf)
#' range.data <- matrix(c(0,100, 0, 100, 0, 100, 0, 100, 0, 100), nrow=2)
#' num.labels <- matrix(c(3,3,3,3,3), nrow = 1)
#' names.variables <- c("input1", "input2", "input3", "input4", "output1")
#' ## plot the membership function
#' \dontrun{plotMF(object = list(var.mf = var.mf, range.data = range.data, 
#'           num.labels = num.labels, names.variables = names.variables))}
#' @export
plotMF <- function(object) {
  
  if (inherits(object, "frbs")) { 
	method.type <- object$method.type
  }
  else if (class(object) == "list") {
	method.type <- c("MANUAL")
  }
  else {
	stop("please input a frbs object or a list of parameters")
  }
  
  if (any(method.type == c("WM", "HYFIS", "ANFIS", "FS.HGD", "GFS.THRIFT", "FIR.DM", "FRBCS.W", "FRBCS.CHI", "GFS.GCCL", "FH.GBML", "SLAVE", "MANUAL"))){
	  if (any(method.type == c("WM", "HYFIS", "GFS.THRIFT"))){
		  range.input <- object$range.input
		  range.output <- object$range.output
		  num.varinput <- object$num.varinput
		  num.fvalinput <- object$num.fvalinput
		  varinp.mf <- object$varinp.mf
		  varout.mf <- object$varout.mf		   
		  range.data <- cbind(range.input, range.output) 
		  num.varinput <- num.varinput + 1
		  num.fvalinput <- cbind(num.fvalinput, num.fvalinput[1]) 
		  var.mf <- cbind(varinp.mf, varout.mf) 
	  }
	 
	 else if (any(method.type == c("ANFIS", "FS.HGD", "FIR.DM", "FRBCS.W", "FRBCS.CHI"))){
		  range.input <- object$range.input
		  range.output <- object$range.output
		  num.fvalinput <- object$num.fvalinput	  
		  range.data <- object$range.input
		  num.varinput <- object$num.varinput
		  var.mf <- object$varinp.mf
	  }
	  else if (any(method.type == c("GFS.GCCL", "FH.GBML", "SLAVE"))){
		  range.data <- object$range.data.ori
		  range.data[1, ] <- 0
		  range.data[2, ] <- 1
		  
		  var.mf <- object$varinp.mf
		  num.varinput <- ncol(range.data)
		  num.fvalinput <- object$num.labels[1, -ncol(object$num.labels), drop = FALSE]
	  }
	  
	  else if (method.type == c("MANUAL")){
		  if (is.null(object$num.labels) || is.null(object$range.data)) {
			  stop("please input the matrix of num.labels and range.data")
		  }
		  else {
			  num.varinput <- ncol(object$num.labels)
			  var.mf <- object$var.mf
			  range.data <- object$range.data
			  num.fvalinput <- object$num.labels
			  if (!is.null(object$names.variables)){ 
				names.variables <- object$names.variables
			  } 
			  else {
				names.variables <- paste("var", i, sep = ".")
			  }
		  }
	  }
	  
	  ##get number of column of var.mf
	  col.var.mf <- ncol(var.mf)
	  
	  ## counter is used to make continue index j
	  counter <- 1
	  
	  ## k is used to index number fuzzy value on each variable on varinput.
	  k <- 1
	  
	  ## make row plot 
	  op <- par(mfrow = c(ceiling(num.varinput/2), 2))
	  
	  ## set the names of var.mf
	  #colnames(var.mf) <- (names.varinput)
	  
	  ## loop as many as number of input variable
	  for (i in 1 : num.varinput){
		j <- counter
		
		## Initialize plot
		MF.degree <- function(x){
		  y <- x - x
		  return (y)
		}
		if (!is.null(object$colnames.var)) { 
			names <- object$colnames.var[i]
		}
		else {
			names <- names.variables[i]
		}
		
		
		curve(MF.degree, range.data[1, i], range.data[2, i], ylim = c(0, 1), col = "white")
		title(main = names)

		## loop as many as number of column of var.mf
		for (j in counter : col.var.mf){
		  
		  ## make boundary point
		  oo <- range.data[1, i]
		  aa <- var.mf[2, j]
		  bb <- var.mf[3, j]
		  cc <- var.mf[4, j]
		  dd <- var.mf[5, j]
		  mm <- range.data[2, i]
		  
		  ##condition for triangular type
		  if (var.mf[1, j] == 1){				
			
			## make a function for plotting, args (x) is sequence data
			f.y1 <- function(x){
			  
			  ## range of input data
			  p.0 <- x[x >= oo & x <= aa]
			  p.1 <- x[x > aa & x <= bb]
			  p.2 <- x[x > bb & x <= cc]
			  p.3 <- x[x > cc & x <= mm]
			  
			  ## build functions 
			  y0 <- (p.0 - p.0)
			  y1 <- (p.1 - aa) / (bb - aa)			
			  y2 <- (p.2 - cc) / (bb - cc)			
			  y4 <- (p.3 - p.3)
			  y <- c(y0, y1, y2, y4)
			  
			  return (y)
			}
			
			curve(f.y1, oo, mm, add = TRUE, col = "violet")
		  }
		  ## condition for trapezoid in left side
		  else if (var.mf[1, j] == 2){
			
			f.y2 <- function(x){
			  
			  p.1 <- x[x <= bb]
			  p.2 <- x[x > bb & x <= cc]
			  p.3 <- x[x > cc & x <= mm]
			  
			  y1 <- (p.1 - p.1 + 1)			
			  y2 <- (p.2 - cc) / (bb - cc)			
			  y3 <- (p.3 - p.3)
			  
			  y <- c(y1, y2, y3)
			  
			  return (y)
			}
			
			curve(f.y2, oo, mm, add=TRUE, col = "blue")
			
		  }
		  ## condition for trapezoid in right side
		  else if (var.mf[1, j] == 3){
			f.y3 <- function(x){
			  
			  p.0 <- x[x >= oo & x <= aa]
			  p.1 <- x[x > aa & x <= bb]
			  p.2 <- x[x > bb & x <= mm]
			  
			  y0 <- (p.0 - p.0)
			  y1 <- (p.1 - aa) / (bb - aa)
			  y2 <- (p.2 - p.2 + 1)			
			  y <- c(y0, y1, y2)
			  
			  return (y)
			}
			
			curve(f.y3, oo, mm, add = TRUE, col = "green")
			
		  }
		  ## condition for trapezoid in the middle
		  else if (var.mf[1, j] == 4){
			
			f.y4 <- function(x){
			  
			  p.0 <- x[x >= oo & x <= aa]
			  p.1 <- x[x > aa & x <= bb]
			  p.2 <- x[x > bb & x <= cc]
			  p.3 <- x[x > cc & x <= dd]
			  p.4 <- x[x > dd & x <= mm]
			  
			  y0 <- (p.0 - p.0)
			  y1 <- (p.1 - aa) / (bb - aa)
			  y2 <- (p.2 - p.2 + 1)			
			  y3 <- (p.3 - dd) / (cc - dd)
			  y4 <- (p.4 - p.4)
			  
			  y <- c(y0, y1, y2, y3, y4)
			  
			  return (y)
			}
			
			curve(f.y4, oo, mm, add = TRUE, col = "red")
			
		  }
		  ## condition for gaussian shape
		  else if (var.mf[1, j] == 5){
			
			f.y5 <- function(x){
			  y <- exp(- 0.5 * (x - aa) ^ 2 / bb ^ 2)
			  return (y)
			}
			## plot the functions
			curve(f.y5, oo, mm, add = TRUE, col = "gray")
		  }
		  ## condition for sigmoid/logistic
		  else if (var.mf[1, j] == 6){
			
			f.y6 <- function(x){
			  y <- 1 / (1 + exp(- aa * (x - bb)))
			  return (y)
			}
			#plot the functions
			curve(f.y6, oo, mm, add = TRUE, col = "black")
		  }
		  ## condition for generalized bell
		  else if (var.mf[1, j] == 7){
			
			f.y7 <- function(x){
			  y <- 1 / (1 + abs((x - cc)/aa) ^ (2 * bb))
			  return (y)
			}
			
			## plot the functions
			curve(f.y7, oo, mm, add = TRUE, col = "black")
		  }
		  
		  counter <- j + 1 
		  k <- k + 1
		  if (k > num.fvalinput[1, i]){
			k <- 1
			break
		  }
		}
	  }

	  par(op)
  }
 else {
	print("The plot is not supported by the used method, please have a look the documentation")
 }
}

# This function can be used to make representations of rules.
#
# @title The rule representing function
# 
# @param object an \code{\link{frbs-object}}.
# @export
rep.rule <- function(object){
	if(!inherits(object, "frbs")) stop("not a legitimate frbs model")
	colnames.var <- object$colnames.var
	
	## make description on rule
	if (any(object$method.type == c("WM", "HYFIS", "GFS.THRIFT", "ANFIS", "FIR.DM", "FS.HGD",
					"FRBCS.W", "FRBCS.CHI", "GFS.GCCL", "FH.GBML", "SLAVE"))){
		if (!is.null(object$num.varinput)){
			num.varinput <- object$num.varinput
		}
		else {
			num.varinput <- ncol(object$rule.data.num) - 1
		}
		
		if (!is.null(object$num.labels)){
			num.labels <- object$num.labels
		}
		else {
			num.labels <- cbind(object$num.fvalinput, object$num.fvalinput[1,1])
		}
		
		if (!is.null(object$rule.data.num)){
			res <- generate.rule(object$rule.data.num, num.labels)
			rule <- res$rule
		}
		else {
			rule <- object$rule
		}
		
		new.rule <- matrix(nrow = nrow(rule), ncol = (ncol(rule) + 2 * (num.varinput + 1)))
		k <- 1
		for (j in 1 : num.varinput){				
			new.rule[, k] <- colnames.var[j] 
			new.rule[, k + 1] <- "is"
			new.rule[, k + 2] <- rule[, 2 * j - 1]
		
			if (j < num.varinput){
				new.rule[, k + 3] <- "and"
			}
			else {
				new.rule[, k + 3] <- "->"
			}
			k <- k + 4
		}
		new.rule[, (ncol(new.rule) - 2)] <- colnames.var[num.varinput + 1] 
		new.rule[, (ncol(new.rule) - 1)] <- "is"
		new.rule[, ncol(new.rule)] <- rule[, ncol(rule)]
		
		## TSK model
		if (any(object$method.type == c("ANFIS", "FIR.DM", "FS.HGD"))){
			rule <- new.rule[, 1 : (ncol(new.rule) - 3), drop = FALSE]
		}
		
		## FRBCS model
		else if (any(object$method.type == c("FRBCS.W", "FRBCS.CHI"))){
			rule <- new.rule[, 1 : (ncol(new.rule) - 3), drop = FALSE]
			rule <- cbind(rule, colnames.var[num.varinput + 1], "is", object$class)
		}
		
		## GFS.GCCL and FH.GBML methods
		else if (any(object$method.type == c("GFS.GCCL", "FH.GBML", "SLAVE"))){
			rule <- new.rule[, 1 : (ncol(new.rule) - 3), drop = FALSE]
			rule <- cbind(rule, colnames.var[num.varinput + 1], "is", object$rule.data.num[, ncol(object$rule.data.num)])
		}
		
		## Mamdani model
		else {
			rule <- new.rule
		}
	}
	
	return (rule)
} 

#' The purpose of this function is to generate a FRBS model from user-given 
#' input without a learning process.
#'
#' It can be used if rules have already been obtained manually, without employing the 
#' learning process. 
#' In the examples shown, we generate a fuzzy model using \code{frbs.gen} and generate the
#' fuzzy rule-based systems step by step manually.   
#'
#' @title The frbs model generator
#' @param range.input a matrix(2 x n) containing the range of the input data. 
#' @param range.output a matrix(2 x n) containing the range of the output data. 
#' @param num.fvalinput a matrix with the number of fuzzy terms of each input variable.
#' 
#' For example: \code{num.fvalinput <- matrix(c(3,2), nrow = 1)}
#' 
#' means that there are two variables where the first variable has three fuzzy terms and the second one has two fuzzy terms.
#' @param varinp.mf a matrix for constructing the shapes of the membership functions. See \code{\link{fuzzifier}}.
#' @param names.varinput a list for giving names to the fuzzy terms. See \code{\link{rulebase}}.
#' @param num.fvaloutput the number of fuzzy terms of the output variable. 
#'
#' For example: \code{num.fvaloutput <- matrix(3, nrow = 1)}
#' 
#' means there are 3 fuzzy terms for the first variable (in this case, there is only one variable).
#' @param varout.mf a matrix for constructing the membership functions of the output variable. 
#' The form is the same as for the \code{varinp.mf} parameter. Please see \code{\link{fuzzifier}}.
#' @param names.varoutput a list for giving names of the fuzzy terms. The form is the same as 
#' for the \code{names.varinput} parameter. Please see \code{\link{rulebase}}.
#' @param rule a list of fuzzy IF-THEN rules. Please see \code{\link{rulebase}}.
#' @param type.model the type of the model. Please see \code{\link{defuzzifier}}.
#' @param type.defuz the type of the defuzzification method. Please see \code{\link{defuzzifier}}.
#' @param type.tnorm the type of the t-norm method. Please see \code{\link{inference}}.
#' @param type.snorm the type of the s-norm method. Please see \code{\link{inference}}.
#' @param func.tsk a matrix of parameters of the function on the consequent part using the Takagi Sugeno Kang model. Please see \code{\link{rulebase}}.
#' @param colnames.var a list of names of input and output variables.
#' @param method.type the type of the selected method. Please see \code{\link{frbs.learn}}.
#' @param name a name of the simulation.
#' @return The \code{\link{frbs-object}}. 
#' @examples 
#' 
#' ## This example shows how to use frbs without 
#' ## learning process.
#'
#' ## Define shape and parameters of membership functions of input variables.
#' ## Please see fuzzifier function to contruct the matrix.
#' varinp.mf <- matrix(c(2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA,
#'                       2, 0, 35, 75, NA, 3, 35, 75, 100, NA,
#'                       2, 0, 20, 40, NA, 1, 20, 50, 80, NA, 3, 60, 80, 100, NA,
#'                       2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA),
#'                       nrow = 5, byrow = FALSE)
#'
#' ## Define number of fuzzy terms of input variables.
#' ## Suppose, we have 3, 2, 3, and 3 numbers of fuzzy terms 
#' ## for first, second, third and fourth variables, respectively.
#' num.fvalinput <- matrix(c(3, 2, 3, 3), nrow=1)
#' 
#' ## Give the names of the fuzzy terms of each input variable.
#' ## It should be noted that the names of the fuzzy terms must be unique,
#' ## so we put a number for making it unique.
#' varinput.1 <- c("a1", "a2", "a3")
#' varinput.2 <- c("b1", "b2")
#' varinput.3 <- c("c1", "c2", "c3")
#' varinput.4 <- c("d1", "d2", "d3")
#' names.varinput <- c(varinput.1, varinput.2, varinput.3, varinput.4)
#'
#' ## Set interval of data.
#' range.input <- matrix(c(0,100, 0, 100, 0, 100, 0, 100), nrow=2)
#' range.output <- matrix(c(0,100), nrow=2)
#'
#' ## Define number of fuzzy terms of output variable.
#' ## In this case, we set the number of fuzzy terms to 3.
#' num.fvaloutput <- matrix(c(3), nrow=1)
#'
#' ## Give the names of the fuzzy terms of the output variable.
#' ## Note: the names of the fuzzy terms must be unique.
#' varoutput.1 <- c("e1", "e2", "e3")
#' names.varoutput <- c(varoutput.1)
#'
#' ## Define the shapes and parameters of the membership functions of the output variables.
#' varout.mf <- matrix(c(2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA),
#'                       nrow = 5, byrow = FALSE)
#'
#' ## Set type of model which is 1 or 2 for Mamdani or Takagi Sugeno Kang model, respectively.
#' ## In this case, we choose Mamdani model.
#' type.model <- 1
#' ## Set weighted average method to be used as defuzzification method.
#' type.defuz <- 1
#' ## We are using standard t-norm and s-norm.
#' type.tnorm <- 1
#' type.snorm <- 1
#'
#' ## Since we don't generate the fuzzy model by learning from data, 
#' ## we have to set Wang and Mendel's technique as type of method.
#' method.type <- "WM"
#' ## Give the name of simulation.
#' name <- "Sim-0"
#' ## Define the fuzzy IF-THEN rules; 
#' ## there are two kinds of model: Mamdani and Takagi Sugeno Kang model
#' ## if we use the Mamdani model then the consequent part is a linguistic term,
#' ## but if we use Takagi Sugeno Kang then we build a matrix representing 
#' ## linear equations in the consequent part.
#' ## In this example we are using the Mamdani model 
#' ## (see the type.model parameter). 
#' ## We make sure that each rule has a "->" sign. 
#' rule <- matrix(c("a1","and","b1","and","c1","and","d1","->","e1",
#'                  "a2","and","b2","and","c2","and","d2", "->", "e2", 
#'                  "a3","and","b2","and","c2","and","d1", "->", "e3"), 
#'                  nrow=3, byrow=TRUE) 
#'
#' ## Define function of TSK if we use it or 
#' ## set NULL if we use the Mamdani model.
#' func.tsk<-matrix(c(1, 1, 5, 2, 1, 3, 1, 0.5, 0.1, 2, 1, 3, 2, 2, 2), nrow=3, byrow=TRUE)
#' ## Provide new data for testing. 
#' newdata<- matrix(c(25, 40, 35, 15, 45, 75, 78, 70), nrow= 2, byrow = TRUE)
#' ## the names of variables
#' colnames.var <- c("input1", "input2", "input3", "input4", "output1")
#' ######################
#' ## 1. The following codes show how to generate a fuzzy model using the frbs.gen function
#' ######################
#' ## Generate a fuzzy model with frbs.gen.
#' object <- frbs.gen(range.input, range.output, num.fvalinput, names.varinput, 
#'                 num.fvaloutput, varout.mf, names.varoutput, rule, varinp.mf,
#'                 type.model, type.defuz, type.tnorm, type.snorm, func.tsk, 
#'                 colnames.var, method.type, name)
#' 
#' ## We can plot the membership function
#' plotMF(object)
#'
#' ## Predicting using new data.
#' res <- predict(object, newdata)
#'
#' ######################
#' ## 2. Using the same data as in the previous example, this example performs 
#' ## step by step of the generation of a fuzzy rule-based system
#' ######################
#' ## Check input data given by user.
#' rule <- rulebase(type.model, rule, func.tsk)
#' 
#' ## Fuzzification Module:
#' ## In this function, we convert crisp values into fuzzy values 
#' ## based on the data and the parameters of the membership function.
#' ## The output: a matrix representing the degree of the membership of the data
#' num.varinput <- ncol(num.fvalinput)
#' MF <- fuzzifier(newdata, num.varinput, num.fvalinput, varinp.mf)
#' 
#' ## Inference Module:
#' ## In this function, we will calculate the confidence factor on the antecedent for each rule
#' ## considering t-norm and s-norm.
#' miu.rule <- inference(MF, rule, names.varinput, type.tnorm, type.snorm)
#'
#' ## Defuzzification Module
#' ## In this function, we calculate and convert the fuzzy values back into crisp values. 
#' result <- defuzzifier(newdata, rule, range.output, names.varoutput,
#'                   varout.mf, miu.rule, type.defuz, type.model, func.tsk)
#' 
#' @export
frbs.gen <- function (range.input, range.output, num.fvalinput, names.varinput, num.fvaloutput, varout.mf, names.varoutput, rule, varinp.mf,
                type.model = 1, type.defuz = 1, type.tnorm = 1, type.snorm = 1, func.tsk = NULL, colnames.var = NULL, method.type = "WM", name = "Sim-0"){

	if (any(method.type == c("WM", "HYFIS", "ANFIS", "FIR.DM", "FS.HGD"))) {
		
		if (any(method.type == c("ANFIS", "FIR.DM", "FS.HGD")) && is.null(func.tsk)) {
			stop("Generating using this method, the consequent part should be given by linear equations as Takagi Sugeno Model") 			
		}
		else {
			#get number of input variable
			num.varinput <- ncol(num.fvalinput)
		
			#get number of output variable
			num.varoutput <- ncol(num.fvaloutput)

			## keep colnames of training data into mod
			if (is.null(colnames.var)) {
				colnames.var <- paste("var", seq(1, (ncol(range.input) + 1)), sep = ".")
			}

			mod <- list(range.input = range.input, range.output = range.output, num.varinput = num.varinput, num.fvalinput = num.fvalinput,
         			names.varinput = names.varinput, num.varoutput = num.varoutput, num.fvaloutput = num.fvaloutput, varout.mf = varout.mf,
 					names.varoutput = names.varoutput, rule = rule, varinp.mf = varinp.mf, type.model = type.model, type.defuz = type.defuz,
					type.tnorm = type.tnorm, type.snorm = type.snorm, func.tsk = func.tsk, method.type = method.type, name = name, 
					colnames.var = colnames.var)
			## build into frbs class
			mod <- frbsObjectFactory(mod)
		}
	}
	else {
		stop ("the built model is not supported by this package, please read the manual")
	}	
return(mod)

}


#' This function is to transform from normalized data into real-valued data. 
#'
#' @title The data de-normalization
#' @param dt.norm a matrix(n x m) of the normalized data.
#' @param range.data a matrix(2 x n) containing the range of the data, where n is the number of variables, and
#' first and second rows are the minimum and maximum value, respectively. 
#' @param min.scale the minimum value within normalization.
#' @param max.scale the maximum value within normalization.
#' @seealso \code{\link{norm.data}}
#' @return the real-valued data
#' @export
denorm.data <- function(dt.norm, range.data, min.scale = 0, max.scale = 1){
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

#' This function is to transform from real-valued data into normalized data. 
#'
#' @title The data normalization
#' @param dt.ori a matrix(n x m) of the original data.
#' @param range.data a matrix(2 x n) containing the range of the data, where n is the number of variables, and
#' first and second rows are the minimum and maximum value, respectively. 
#' @param min.scale the minimum value within normalization.
#' @param max.scale the maximum value within normalization.
#' @seealso \code{\link{denorm.data}}
#' @return the normalized data
#' @export
norm.data <- function(dt.ori, range.data, min.scale = 0, max.scale = 1){
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

