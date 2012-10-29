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
#' The purpose of this function is to generate the frbs model from user-given input. 
#' It can be used if rules have already been obtained manually, without employing the 
#' learning process, e.g. from expert knowledge.  
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
#' @param method.type the type of the selected method. Please see \code{\link{frbs.learn}}.
#' @param name a name of the simulation.
#' @return The \code{\link{frbs-object}}. 
#' @examples 
#' 
#' ##########
#' ## This example shows how to use frbs without 
#' ## learning process if we have already rules.
#' #########
#' ## define shape of membership functions of input variables
#' varinp.mf <- matrix(c(2,1,3,2,3,2,3,2,3,0,30,60,0,40,0,40,0,40,20,50,80,
#'                       30,80,30,80,30, 80,40,70, 100,60,100,60,100,60,100,
#'                       0,0,100,0,100,0,100,0,100), nrow=5, byrow=TRUE)
#'
#' ## define number of fuzzy terms of input variables
#' num.fvalinput <- matrix(c(3, 2, 2, 2), nrow=1)
#'
#' ## give the names of the fuzzy terms of each input variable
#' ## Note: the names of the fuzzy terms must be unique.
#' varinput.1 <- c("a1", "a2", "a3")
#' varinput.2 <- c("b1", "b2")
#' varinput.3 <- c("c1", "c2")
#' varinput.4 <- c("d1", "d2")
#' names.varinput <- c(varinput.1, varinput.2, varinput.3, varinput.4)
#' range.input <- matrix(c(0,100, 0, 100, 0, 100, 0, 100), nrow=2)
#' range.output <- matrix(c(0,100), nrow=2)
#'
#' ## define number of fuzzy terms of output variable
#' num.fvaloutput <- matrix(c(3), nrow=1)
#'
#' ## give the names of the fuzzy terms of the output variable
#' ## Note: the names of the fuzzy terms must be unique.
#' varoutput.1 <- c("e1", "e2", "e3")
#' names.varoutput <- c(varoutput.1)
#'
#' ## define the shapes of the membership functions of the input variables
#' varout.mf <- matrix(c(2,1,3,0,30,60,20,50,80,40,70,100,0,0,100),
#'                       nrow=5, byrow=TRUE)
#'
#' type.model <- 1
#' type.defuz <- 1
#' type.tnorm <- 1
#' type.snorm <- 1
#' method.type <- "WM"
#' name <- "Sim-0"
#' ## define fuzzy IF-THEN rules, 
#' ## Please make sure that each rule has "->" sign. 
#' ## If we use the TSK model, we need to define the 
#' ## fuzzy terms of the consequent part.
#' rule <- matrix(c("a1","and","b1","and","c1","and","d1","->","e1",
#'                  "a2","or","b2","and","c2","and","d2", "->", "e2", 
#'                  "a3","or","b2","and","c2","and","d1", "->", "e3"), 
#'                  nrow=3, byrow=TRUE) 
#'
#' ## define function of TSK if we use it or 
#' ## set NULL if we use the Mamdani model
#' func.tsk<-matrix(c(1, 1, 5, 2, 1, 3, 1, 2, 2), nrow=3, ncol=3, byrow=TRUE)
#'
#' ## generate model by frbs.gen that will be used as input for predict
#' object <- frbs.gen(range.input, range.output, num.fvalinput, names.varinput, 
#'                 num.fvaloutput, varout.mf, names.varoutput, rule, varinp.mf,
#'                 type.model, type.defuz, type.tnorm, type.snorm, func.tsk, 
#'                 method.type, name)
#' 
#' 
#' ## testing data 
#' newdata <- matrix(c(10, 10, 10, 10, 10, 20, 12, 23, 30, 30, 23, 23), nrow =3)
#'
#' ## make the predictions 
#' res <- predict(object, newdata)
#' 
#' @export
frbs.gen <- function (range.input, range.output, num.fvalinput, names.varinput, num.fvaloutput, varout.mf, names.varoutput, rule, varinp.mf, type.model = 1, type.defuz = 1, type.tnorm = 1, type.snorm = 1, func.tsk = NULL, method.type = "WM", name = "Sim-0"){

	if (any(method.type == c("WM", "HYFIS", "FRCS", "ANFIS", "DM", "HGD", "GFS"))) {
		
		if (any(method.type == c("ANFIS", "FRCS", "DM", "HGD")) && is.null(func.tsk)) {
			stop("Generating using this method, the consequent part should be given by linear equations as Takagi Sugeno Model") 			
		}
		else {
			#get number of input variable
			num.varinput <- ncol(num.fvalinput)
		
			#get number of output variable
			num.varoutput <- ncol(num.fvaloutput)
			mod <- list(range.input = range.input, range.output = range.output, num.varinput = num.varinput, num.fvalinput = num.fvalinput, names.varinput = names.varinput, num.varoutput = num.varoutput, num.fvaloutput = num.fvaloutput, varout.mf = varout.mf, names.varoutput = names.varoutput, rule = rule, varinp.mf = varinp.mf, type.model = type.model, type.defuz = type.defuz, type.tnorm = type.tnorm, type.snorm = type.snorm, method.type = method.type, name = name)
			## build into frbs class
			mod <- frbsObjectFactory(mod)
		}
	}
	else {
		stop ("the built model is not supported by this package, please read the manual")
	}	
return(mod)

}