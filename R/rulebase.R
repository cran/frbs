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
#' This function checks the consistency of a rule definition (given by the user).
#' The rulebase consists of several fuzzy IF-THEN rules. 
#'  
#' For rules of the Mamdani model, there are 2 parts in each rule, the antecedent and the consequent part, which are separated by "->". 
#'
#' For example:  \code{r1 <- c("a1","and","b1","->", "c1")}
#'
#' Here, ("a1", "and", "b1") is the antecedent, with "a1" and "b1" being fuzzy terms, and ("c1") is the consequent part. 
#' 
#' A fuzzy IF-THEN rule base with several rules is defined in the following way:
#' 
#' \code{r1 <- c("not a1","and","b1", "->", "c1")}
#'
#' \code{r2 <- c("a2","or","b2", "->", "c2")}
#'
#' \code{r3 <- c("a3","or","b2", "->", "c3")}
#'
#' \code{rule <- list(r1,r2,r3)}
#'
#' For rules of the Takagi Sugeno Kang model, the rules are at first defined without the consequent part, e.g.: 
#'
#' \code{r1 <- c("a1",1,"b1","->")}
#'
#' \code{r2 <- c("a2",2,"b2", "->")}
#'
#' \code{r3 <- c("a3","2","b2", "->")}
#'
#' \code{rule <- list(r1,r2,r3)}
#'
#' The consequences are defined then as a matrix \code{fun_tsk}, which contains the linear equations of the consequences of the rules.
#' The dimension of this matrix is [<number_of_rules>, <number_of_variables> + 1]. The matrix has one extra column for the constants. 
#' If there is no constant, a zero is put.
#'
#' So, for example, if we have 3 rules and 2 fuzzy variables (A, B), the matrix \code{fun_tsk} has dim(3,3), as in: 
#' 
#' \code{func.tsk <- matrix(c(1, 1, 5, 2, 1, 3, 1, 2, 2), nrow=3, ncol=3, byrow = TRUE)}
#'
#' Furthermore, we can represent linguistic hedge within the rules. The kinds of hedges used are
#'
#' \itemize{
#' \item "extremely" reduces the truth value. For example, "extremely a1" means membership function a1 = miu(a1)^3 
#'
#' \item "very" reduces the truth value. For example, "very a1" means membership function a1 = miu(a1)^2
#'
#' \item "somewhat" increases the truth value. For example, "somewhat a1" means membership function 
#' a1 = miu(a1)^0.5
#'
#' \item "slightly" increases the truth value. For example, "slightly a1" means membership function 
#' a1 = miu(a1)^0.33
#' }
#' 
#' An example of fuzzy IF-THEN rules using linguistic hedge is: 
#' 
#' \code{r1 <- c("very a1","and","b1","->","c1")}
#'
#' \code{r2 <- c("a2",2,"b2", "->", "c2")}
#'
#' \code{r3 <- c("a3","2","slightly b2", "->", "c3")}
#'
#' \code{rule <- list(r1,r2,r3)} 
#' 
#' Furthermore, the following is an example in order to give names to the fuzzy terms in the input and output variables.
#'
#' \code{varinput.1 <- c("a1", "a2", "a3")}
#'
#' \code{varinput.2 <- c("b1", "b2")}
#'
#' \code{names.varinput <- c(varinput.1, varinput.2)}
#' 
#' \code{names.varoutput <- c("c1", "c2", "c3")}
#' 
#' Note that the names of the fuzzy terms must be unique and
#' if we are using the learning methods, the fuzzy IF-THEN rules will be generated automatically 
#' as the outputs of \code{\link{frbs.learn}}.
#'
#' @title The rule checking function
#' @param type.model a value determining the type of model to use. Here, 1 and 2 mean Mamdani and Takagi Sugeno Kang model, respectively. 
#' @param rule a matrix or list of rules.
#' @param func.tsk a matrix representing the consequent parts of rules in Takagi Sugeno Kang formulation.
#' @seealso \code{\link{defuzzifier}}, \code{\link{inference}}, and \code{\link{fuzzifier}}
#' @return fuzzy IF-THEN rule base
#' @export 
rulebase <- function(type.model, rule, func.tsk = NULL){

if (class(rule) == "matrix"){
	rule.list <- list() 
	for (i in 1 : nrow(rule)){
		rule.list[i] <- list(rule[i, ])
	}
	rule <- rule.list	
}


##condition for Mamdani model
if (type.model == 1) {
	##Checking the sign of separating between antecedent and consequence 
	for (i in 1 : length(rule)){
		equal <- unlist(rule[i])
		equal <- as.matrix(equal)
		check <- any(equal == "->")
		if (check == FALSE )
			stop("rule must contain \"->\" as separating between antecendent and consequence") 			
	}
}	

#condition for Takagi Sugeno Kang model
if (type.model == 2) {
	
	if (is.null(func.tsk))
		stop("You are using Takagi Sugeno Kang model, so the consequent part must be linear equations, please insert their values")
	##Checking the sign of separating between antecedent and consequence 
	for (i in 1 : length(rule)){
		equal <- unlist(rule[i])
		equal <- as.matrix(equal)
		check <- any(equal == "->")
		if (check == FALSE )
			stop("rule must contain \"->\" as separating between antecendent and consequence") 				
	}
}

return(rule)
}

