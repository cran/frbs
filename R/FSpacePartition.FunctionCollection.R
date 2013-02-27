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
#' Fuzzification refers to the process of transforming a crisp set into fuzzy terms. 
#'
#' In this function, there are five shapes of membership functions implemented, 
#' namely triangular, trapezoid, Gaussian, sigmoid, and generalized bell.
#' 
#' @title Transform from crisp set into fuzzy terms
#'
#' @param data a matrix of data containing numerical elements.
#' @param num.varinput number of input variables.
#' @param num.fvalinput the number of labels of the input variables.
#' @param varinp.mf a matrix containing the parameters to form the membership functions. 
#' The dimension of the matrix is (5, n) where n is 
#' the number of fuzzy terms/labels and the number of variables.
#' The rows of the matrix represent:
#' The first row is the type of membership function, where 1 means triangular, 
#' 2 means trapezoid 1a (left side), 
#' 3 means trapezoid 1b (right side), 4 means trapezoid 2 (in the middle), 5 means Gaussian, 
#' 6 means sigmoid, 7 means generalized bell. The second until fifth row indicate 
#' the critical points to construct the functions. 
#' \itemize{
#' \item triangular has three parameters (a, b, c), where b is the center point of the triangular,
#' and a and c are the left and right points, respectively.
#' \item trapezoid has four parameters (a, b, c, d).
#' \item Gaussian has two parameters (mean and variance).
#' \item sigmoid has two parameters (gamma and c).
#' \item generalized bell has three parameters (a, b, c).
#' }
#' 
#' For example:
#' 
#' \code{varinp.mf <- matrix(c(2,1,3,2,3,0,30,60,0,40,20,50,80,}
#'
#' \code{30,80,40,70,100,60,100,0,0,100,0,100), nrow=5, byrow=TRUE)}
#'
#' @seealso \code{\link{defuzzifier}}, \code{\link{rulebase}}, and \code{\link{inference}}
#' @return A matrix of the degree of each fuzzy term based on the shape of 
#' the membership functions 
#' @export
fuzzifier <- function(data, num.varinput, num.fvalinput, varinp.mf){

##count number of column of data
ncol.data <- ncol(data)

##count number of column of varinp.mf (matrix used to build membership function) 
ncol.var <- ncol(varinp.mf)

##Inisialitation matrix of Membership function
MF <- matrix(nrow = nrow(data), ncol = ncol.var)

##check 
if (ncol.data != num.varinput)
	stop("data is not the same as the number of variable")
if (ncol.var != sum(num.fvalinput))
	stop("the parameter of membership function is not the same with variable")

##h is index equal to number of data
##i is index for numbering variable
##j is index equal to number of varinp.mf column
##ii is used for counting how many iteration has been done in the following loop
##jj is used for keeping the iteration continueing to next index in varinp.mf


##iterate as along number of data
for (h in 1 : nrow(data)){
	jj <- 1
	##iterate for each crisp value on each data 
	for (i in 1: ncol(data)){

		##counter for break
		ii <- 1
		
		##loop all column on varinp.mf
		for (j in jj : ncol(varinp.mf)){		
			
			##
			##checking for type 1: Triangular, if varinp.mf[1,] == 1 
			##parameter=(a,b,c), where a < b < c
			##a=varinp.mf[2,]
			##b=varinp.mf[3,]
			##c=varinp.mf[4,]
			if (varinp.mf[1, j] == 1){
				if (data[h, i] <= varinp.mf[2, j])
				temp <- 0
				else if (data[h, i] <= varinp.mf[3, j])
				temp <- (data[h, i] - varinp.mf[2, j]) / (varinp.mf[3, j] - varinp.mf[2, j])
				else if (data[h, i] < varinp.mf[4, j])
				temp <- (data[h, i] - varinp.mf[4, j]) / (varinp.mf[3, j] - varinp.mf[4, j])
				else
				temp <- 0
			}

			##checking for type 2: Trapezoid_1a, if varinp.mf[1,] ==2
			##Trapezoid_1a is the edge on the left: vertical
			##parameter=(a,b,c)
			##a=varinp.mf[2,]
			##b=varinp.mf[3,]
			##c=varinp.mf[4,]
			else if (varinp.mf[1, j] == 2){
				if (data[h, i] <= varinp.mf[3, j])
				temp <- 1
				else if (data[h, i] <= varinp.mf[4, j])
				temp <- (data[h, i] - varinp.mf[4, j]) / (varinp.mf[3, j] - varinp.mf[4, j])
				else
				temp <- 0				
			}

			##checking for type 3: Trapezoid_1b, if varinp.mf[1,] == 3
			##Trapezoid_1b is the edge on the right: vertical
			##parameter=(a,b,c)
			##a=varinp.mf[2,]
			##b=varinp.mf[3,]
			##c=varinp.mf[4,]
			else if (varinp.mf[1, j] == 3){
				if (data[h, i] <= varinp.mf[2, j])
				temp <- 0
				else if (data[h, i] < varinp.mf[3, j])
				temp <- (data[h, i] - varinp.mf[2, j]) / (varinp.mf[3, j] - varinp.mf[2, j])
				else 
				temp <- 1
			}

			##checking for type 4: Trapezoid_2
			##parameter=(a,b,c,d)
			##a=varinp.mf[2,]
			##b=varinp.mf[3,]
			##c=varinp.mf[4,]
			##d=varinp.mf[5,]
			else if (varinp.mf[1, j] == 4){
				if (data[h,i] <= varinp.mf[2, j] || data[h,i] > varinp.mf[5, j])
				temp <- 0
				else if (data[h, i] > varinp.mf[3, j] && data[h, i] <= varinp.mf[4, j])
				temp <- 1
				else if (data[h, i] > varinp.mf[2, j] && data[h, i] <= varinp.mf[3, j])
				temp <- (data[h, i] - varinp.mf[2, j]) / (varinp.mf[3, j] - varinp.mf[2, j])
				else if (data[h,i] > varinp.mf[4, j] && data[h,i] <= varinp.mf[5,j])
				temp <- (data[h, i] - varinp.mf[5, j]) / (varinp.mf[4, j] - varinp.mf[5, j])
			}

			##checking for type 5: Gaussian
			##parameter=(mean a, standard deviation b)
			##a=varinp.mf[2,]
			##b=varinp.mf[3,]
			else if (varinp.mf[1, j] == 5){
				temp <- exp(- 0.5 * (data[h, i] - varinp.mf[2, j])^2 / varinp.mf[3, j]^2)
			}

			##checking for type 6: Sigmoid/logistic
			##parameter=(gamma,c)
			##gamma=varinp.mf[2,]
			##c=varinp.mf[3,]
			else if (varinp.mf[1, j] == 6) {
				temp <- 1/(1 + exp(- varinp.mf[2, j] * (data[h, i] - varinp.mf[3, j])))
			}

			##checking for type 7: Generalized Bell
			##parameter=(a,b,c)
			##a=varinp.mf[2,]
			##b=varinp.mf[3,]
			##c=varinp.mf[4,]
			else if (varinp.mf[1, j] == 7) {
				temp <- 1/(1 + abs((data[h, i] - varinp.mf[4, j])/varinp.mf[2, ]) ^ (2 * varinp.mf[3, ]))   
			}

		##save membership function on MF for each data		
		MF[h, j] <- temp
		
		ii <- ii + 1
		jj <- jj + 1
		##this checking is used for control the number of fuzzy value for each variable
			if (ii > num.fvalinput[1, i])
			break
		
		}
	}

}

return(MF)
}

#' This function checks the consistency of a rule definition (given by the user).
#' The rulebase consists of several fuzzy IF-THEN rules. 
#'  
#' For rules of the Mamdani model, there are 2 parts in each rule, the antecedent and 
#' the consequent part, which are separated by "->". 
#'
#' For example:  \code{r1 <- c("a1","and","b1","->", "c1")}
#'
#' Here, ("a1", "and", "b1") is the antecedent, with "a1" and "b1" being fuzzy terms, 
#' and ("c1") is the consequent part. 
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

#' Inference refers to the process of fuzzy reasoning. 
#' 
#' There are two methods of inference for fuzzy systems based on linguistic rules: 
#' The Mamdani and Takagi Sugeno Kang model. 
#' 
#' \bold{The Mamdani model:}
#' A fuzzy system with, e.g., two inputs x1 and x2 (antecedents) and a single output y (consequent)
#' is described by the following fuzzy IF-THEN rule:
#'
#' \code{IF x1 is A1 and x2 is A2 THEN y is B}
#'
#' where A1 and A2 are the fuzzy sets representing the antecent pairs and
#' B is the fuzzy set representing the consequent.
#'
#' \bold{The Takagi Sugeno Kang model:}
#' Suppose we have two inputs x1 and x2 and output y, then the fuzzy IF-THEN rule is as follows:
#'
#' \code{IF x1 is A1 and x2 is A2 THEN y is y = f(x1, x2)}
#'
#' where y = f(x1, x2) is a crisp function in the consequent part which is usually a polynomial function,
#' and A1 and A2 are the fuzzy sets representing the antecent pairs.
#' 
#' Futhermore, this function has the following capabilities: 
#' \itemize{ 
#' \item It supports unary operators (not) and binary operators (AND and OR).
#' \item It provides linguistic hedge (extremely, very, somewhat, slightly).
#' \item there are several methods for the t-norm and s-norm.
#' }
#' @title The process of fuzzy reasoning
#'
#' @param MF a matrix of the degrees of membership functions which is a result of the \code{\link{fuzzifier}}.
#' @param rule a matrix or list of fuzzy IF-THEN rules. See \code{\link{rulebase}}.
#' @param names.varinput a list of names of the input variables. 
#' @param type.tnorm a value between 1 and 5 which represents the type of t-norm to be used: 
#' \itemize{
#' \item \code{1} means standard t-norm: min(x1, x2).
#' \item \code{2} means Hamacher product: (x1 * x2)/(x1 + x2 - x1 * x2).
#' \item \code{3} means Yager class (with tao = 1): 1- min(1, ((1 - x1) + (1 - x2))).
#' \item \code{4} means product: (x1 * x2).
#' \item \code{5} means bounded product: max(0, x1 + x2 - 1).
#' }
#' @param type.snorm a value between 1 and 5 which represents the type of s-norm to be used: 
#' \itemize{
#' \item \code{1} means standard s-norm: max(x1, x2).
#' \item \code{2} means Hamacher sum: (x1 + x2 - 2x1 * x2) / 1 - x1 * x2.
#' \item \code{3} means Yager class (with tao = 1): min(1, (x1 + x2)).
#' \item \code{4} means the sum: (x1 + x2 - x1 * x2).
#' \item \code{5} means the bounded sum: min(1, x1 + x2).
#' }
#' @seealso \code{\link{defuzzifier}}, \code{\link{rulebase}}, and \code{\link{fuzzifier}}.
#' @return a matrix of the degrees of the rules. 
#' @export
inference<-function(MF, rule, names.varinput, type.tnorm, type.snorm){

##calculate number of data
nMF <- nrow(MF)

##calculate number of rule
nrule <- length(rule)
	
##allocate memory for membership function on antecedent
miu.rule <- matrix(nrow = nMF, ncol = nrule)

##give names on each column for detecting which rule will be fired
colnames(MF) <- c(names.varinput)

##Iteration for n data
for(k in 1 : nMF){
	i <- 1
	##Iteration for each rule
	for(i in 1 : nrule){
		##change list of rule into matrix
		temp <- unlist(rule[i])
		
		##detect location of "->" sign as separation between antecedent and consequence
		loc <- which(temp == "->")
				
			##Inisialization for calculating MF of antecendet
			temp.fvalue <- strsplit(temp[1], split=" ")
						
			Mtemp.fvalue <- unlist(temp.fvalue)
			
			length.MtempValue <- length(Mtemp.fvalue)
			
			##check that the fuzzy value has unary operator ("not" and linguistic hedge)
			if (length.MtempValue > 1){
				if (Mtemp.fvalue[1] == "not"){
					if (Mtemp.fvalue[2] == "extremely") {
						val.antecedent.temp <- (MF[k, Mtemp.fvalue[3]])^3
						val.antecedent <- 1 - val.antecedent.temp
					}
					else if (Mtemp.fvalue[2] == "very"){
						val.antecedent.temp <- (MF[k, Mtemp.fvalue[3]])^2
						val.antecedent <- 1 - val.antecedent.temp
					}
					else if (Mtemp.fvalue[2] == "somewhat"){
						val.antecedent.temp <- (MF[k, Mtemp.fvalue[3]])^0.5
						val.antecedent <- 1 - val.antecedent.temp
					}
					else if (Mtemp.fvalue[2] == "slightly"){
						val.antecedent.temp <- (MF[k, Mtemp.fvalue[3]])^0.333
						val.antecedent <- 1 - val.antecedent.temp
					}
					else 
						val.antecedent <- 1 - MF[k, Mtemp.fvalue[2]]
				}
				else {
					if (Mtemp.fvalue[1] == "extremely")
					val.antecedent <- (MF[k, Mtemp.fvalue[2]])^3
					else if (Mtemp.fvalue[1] == "very")
					val.antecedent <- (MF[k, Mtemp.fvalue[2]])^2
					else if (Mtemp.fvalue[1] == "somewhat")
					val.antecedent <- (MF[k, Mtemp.fvalue[2]])^0.5
					else if (Mtemp.fvalue[1] == "slightly")
					val.antecedent <- (MF[k, Mtemp.fvalue[2]])^0.333
					else if (Mtemp.fvalue[1] == "not")
					val.antecedent <- 1 - (MF[k, Mtemp.fvalue[2]])
				}			
			} 
			
			else{
				if (Mtemp.fvalue[1] == "don't_care"){
					val.antecedent <- 1
				}
				else{
					val.antecedent <- MF[k, temp[1]]
				}
			}	
			## iterate to calculate degree of MF until find "->" sign (equals to loc)
			seqq <- seq(from = 1, to = loc, by = 2)
			for (j in seqq) {
			
			temp.split <- strsplit(temp[j + 2], split = " ")
			Mtemp.fvalue <- unlist(temp.split)
			length.MtempValue <- length(Mtemp.fvalue)
			
			if (length.MtempValue > 1) {
				if (Mtemp.fvalue[1] == "not"){
					if (Mtemp.fvalue[2] == "extremely") {
						val.antecedent.temp <- (MF[k, Mtemp.fvalue[3]])^3
						val.antecedent.b <- 1 - val.antecedent.temp
					}
					else if (Mtemp.fvalue[2] == "very"){
						val.antecedent.temp <- (MF[k, Mtemp.fvalue[3]])^2
						val.antecedent.b <- 1 - val.antecedent.temp
					}
					else if (Mtemp.fvalue[2] == "somewhat"){
						val.antecedent.temp <- (MF[k, Mtemp.fvalue[3]])^0.5
						val.antecedent.b <- 1 - val.antecedent.temp
					}
					else if (Mtemp.fvalue[2] == "slightly"){
						val.antecedent.temp <- (MF[k, Mtemp.fvalue[3]])^0.333
						val.antecedent.b <- 1 - val.antecedent.temp
					}
					else 
						val.antecedent.b <- 1 - MF[k, Mtemp.fvalue[2]]
				}
				else {
					if (Mtemp.fvalue[1] == "extremely")
					val.antecedent.b <- (MF[k, Mtemp.fvalue[2]])^3
					else if (Mtemp.fvalue[1] == "very")
					val.antecedent.b <- (MF[k, Mtemp.fvalue[2]])^2
					else if (Mtemp.fvalue[1] == "somewhat")
					val.antecedent.b <- (MF[k, Mtemp.fvalue[2]])^0.5
					else if (Mtemp.fvalue[1] == "slightly")
					val.antecedent.b <- (MF[k, Mtemp.fvalue[2]])^0.333
					else if (Mtemp.fvalue[1] == "not")
					val.antecedent.b <- 1 - (MF[k, Mtemp.fvalue[2]])
				}		
			} 
			else{
				if (Mtemp.fvalue[1] == "don't_care"){
					val.antecedent.b <- 1
				}
				else{
					val.antecedent.b <- MF[k, temp[j + 2]]
				}
				
			}
			##condition for conjunction operator (AND)
			if ((temp[j + 1] == "1") || (temp[j + 1] == "and")){
				##condition for type.tnorm used is standard type(min)
				if (type.tnorm == 1){								
					if (val.antecedent.b < val.antecedent)
					val.antecedent <- val.antecedent.b
				}
				##condition for type.tnorm used is Hamacher product
				else if (type.tnorm == 2) {									
					val.antecedent <- (val.antecedent * val.antecedent.b) / (val.antecedent + val.antecedent.b - val.antecedent * val.antecedent.b)
				}
				
				##condition for type.tnorm used is yager class (with tao = 1)
				else if (type.tnorm == 3) {										
					temp.val.ante <- (1 - val.antecedent) + (1 - val.antecedent.b)
					if (temp.val.ante <= 1){
						val.antecedent <- temp.val.ante
					}
					else
					val.antecedent <- 1
				}
				
				##condition for type.tnorm used is product
				else if (type.tnorm == 4) {				
					val.antecedent <- val.antecedent * val.antecedent.b
				}
				
				##condition for type.tnorm used is bounded product
				else if (type.tnorm == 5){
					temp.val.ante <- (val.antecedent * val.antecedent.b - 1)
					if (temp.val.ante > 0){
						val.antecedent <- temp.val.ante
					}
					else
					val.antecedent <- 0
				}
			}

			##condition for disjunction operator (OR)
			else if ((temp[j + 1] == "2") || (temp[j + 1] == "or")){
				
					##condition for type.snorm used is standard type (max)
					if (type.snorm == 1){						
						if (val.antecedent.b > val.antecedent)
						val.antecedent <- val.antecedent.b
					}
					
					##condition for type.snorm used is Hamacher sum
					else if (type.snorm == 2) {							
						val.antecedent <- (val.antecedent + val.antecedent.b - 2 * val.antecedent * val.antecedent.b) / (1 - val.antecedent * val.antecedent.b)
					}
					
					##condition for type.snorm used is yager class (with tao = 1)
					else if (type.snorm == 3){							
						temp.val.ante <- (val.antecedent + val.antecedent.b)
						if (temp.val.ante <= 1){
							val.antecedent <- temp.val.ante
						}
						else
							val.antecedent <- 1						
					}
					
					##condition for type.snorm used is sum
					else if (type.snorm == 4){
						val.antecedent <- (val.antecedent + val.antecedent.b - val.antecedent * val.antecedent.b)
					}
					
					##condition for type.snorm used is bounded sum
					else if (type.snorm == 5){
						temp.val.ante <- (val.antecedent + val.antecedent.b)
						if (temp.val.ante <= 1){
							val.antecedent <- temp.val.ante
						}
						else
							val.antecedent <- 1					
					}				
			}
			
			j <- j + 2			
			##Checking 
			if (j == loc - 1)
				break	
			}
		
		##save value of MF on each rule			
		miu.rule[k, i] <- c(val.antecedent)		
	}	
}
## result 
## number of row is based on number of data.
## number of column is based on number of rule.
return(miu.rule)
}

#' Defuzzification is a transformation that extracts the crisp values from the fuzzy terms. 
#' 
#' In this function, there exist two kinds of models which are based on the Mamdani model and 
#' the Takagi Sugeno Kang model. 
#' For the Mamdani model there are five methods for defuzzifying a fuzzy term A of a universe 
#' of discourse Z. 
#' They are as follows:
#' \enumerate{
#' \item weighted average method
#' \item first of maxima
#' \item last of maxima
#' \item mean of maxima
#' \item modified COG
#' }
#' 
#' @title Defuzzifier to transform from fuzzy terms to crisp values
#'
#' @param data a matrix(m x n) of data, where m is the number of instances and 
#' n is the number of variables.
#' @param rule a list or matrix of fuzzy IF-THEN rules, as discussed in \code{\link{rulebase}}.
#' @param range.output a matrix(2 x n) containing the range of the output data. 
#' @param names.varoutput a list for giving names to the fuzzy terms. See \code{\link{rulebase}}.
#' @param varout.mf a matrix constructing the membership function of the output variable. 
#'        See \code{\link{fuzzifier}}.
#' @param miu.rule the results of the inference module. See \code{\link{inference}}.
#' @param type.defuz the type of defuzzification to be used, where 1 means weighted average method, 
#'        and 2, 3, 4 and 5 mean first, last, mean maxima and modified COG, respectively.
#' @param type.model the type of the model that will be used in the simulation. 
#'         Here, 1 or 2 means we use Mamdani or Takagi Sugeno Kang model, respectively.
#' @param func.tsk a matrix used to build the linear equation for the consequent part 
#'         if we are using Takagi Sugeno Kang. See also \code{\link{rulebase}}.
#' @seealso \code{\link{fuzzifier}}, \code{\link{rulebase}}, and \code{\link{inference}}
#' @return A matrix of crisp values
#' @export
defuzzifier <- function(data, rule, range.output, names.varoutput = NULL, varout.mf = NULL, miu.rule, type.defuz = 1, type.model = 1, func.tsk = NULL){

## copy rule
rule.temp <- matrix(unlist(rule), nrow = length(rule), byrow= TRUE)

## Inisialitation
def <- matrix(0, nrow=nrow(data), ncol = 1)
def.temp <- matrix(0, nrow=nrow(data), ncol = 1)

## Mamdani type
if (type.model == 1){
	
	## check names.varoutput
	if (is.null(names.varoutput)){
		stop("please define the names of output variable")
	}

	## check parameters of membership functions on output variables
	if (is.null(varout.mf)){
		stop("please define the parameters of membership functions on the output variable")
	}

	## Inisialitation
	cum <- matrix(0, nrow=nrow(data), ncol = ncol(varout.mf)) 
	div <- matrix(0, nrow=nrow(data), ncol = ncol(varout.mf))

	for (k in 1 : nrow(data)){
		##Check not zero on miu.rule
		chck <- which(miu.rule[k, ] != 0)
		l.chck <- length(chck)
		cum <- matrix(0, nrow=nrow(data), ncol = l.chck) 
		div <- matrix(0, nrow=nrow(data), ncol = l.chck)
		
		## initialize
		temp <- 0
		temp.d <- 0
		temp1 <- 0
		temp2 <- 0
		
		#if there are some no zero element on miu rule
		if (l.chck != 0) {
			indx <- matrix(nrow = l.chck, ncol = 2)
			temp.indx <- matrix(nrow = 1, ncol = 2)
			## along number of not zero on miu.rule, check and save string related on names.varoutput and its value
			for (ii in 1 : l.chck){
				#get the names of output value on each rule
				
				strg <- c(rule.temp[chck[ii], ncol(rule.temp)])
				aaa <- which(strg == names.varoutput)
				if (length(aaa) != 0) {
				indx[ii, 1] <- aaa
				indx[ii, 2] <- miu.rule[k, chck[ii]]
				}
			}
			
			## check duplication on indx, choose with has a max of degree
			if (nrow(indx) >= 2){
				for (jj in 1 : (nrow(indx) - 1)){
					temp.indx[1, ] <- indx[jj, ]
					for (kk in (jj + 1) : nrow(indx)){
						ck <- length(which(temp.indx[1, 1] == indx[kk, 1]))
						if (ck != 0){
							if (temp.indx[1, 2] >= indx[kk, 2])
								indx[kk, ] <- NA
							else 
								indx[jj, ] <- NA
						}
					}
				}				
			}
			
			# erase NA value on indx
			indx <- na.omit(indx)
			
			#defuzzification procedure for  Centroid
			if (type.defuz == 1) {								
				for (i in 1 : nrow(indx)){										
					#calculate modified centroid
					# update for gaussian
					
					if (varout.mf[1, indx[i, 1]] == 5 || varout.mf[1, indx[i, 1]] == 6 || varout.mf[1, indx[i, 1]] == 7) {
						## center point
						av.point <- varout.mf[2, indx[i, 1]]
						
						## indx is fired miu.rule/rule
						temp1 <- indx[i, 2] * av.point 
						temp2 <- indx[i, 2] 
												
						temp <- temp + temp1
						temp.d <- temp.d + temp2
					}
					
					else {						
						av.point <- varout.mf[3, indx[i, 1]] 	
						
						# check first which one greater between indx[i, 2] with the value of MF (for centroid)
						temp1 <- indx[i, 2] * av.point 
						
						temp2 <- indx[i, 2]
						temp <- temp + temp1
						temp.d <- temp.d + temp2
					}
					cum[k, i] <- temp
					div[k, i] <- temp.d   					
				}
							
				if (sum(div[k, ]) == 0){
					def[k, 1] <- (min(range.output, na.rm=TRUE) + max(range.output, na.rm = TRUE))/2
				}
				else{
					def.temp[k, 1] <- sum(cum[k, ]) / sum(div[k, ])
					if (def.temp[k, 1] <= min(range.output, na.rm=TRUE)){
						def[k, 1] <- (min(range.output, na.rm=TRUE) + max(range.output, na.rm = TRUE))/2
					}
					else if (def.temp[k, 1] >= max(range.output, na.rm=TRUE)){
						def[k, 1] <- (min(range.output, na.rm=TRUE) + max(range.output, na.rm = TRUE))/2
					}
					else {
						def[k, 1] <- def.temp[k, 1]
					}
				}
			}
			
			## procedure for type.defuz == 2 (fisrt of Maxima)
			else if (type.defuz == 2){
				max.temp <- max(indx[, 2], na.rm = TRUE)
				max.indx <- which(indx[, 2] == max.temp)			
				
				aa <- varout.mf[2, indx[max.indx[1], 1]]
				bb <- varout.mf[3, indx[max.indx[1], 1]]
				cc <- varout.mf[4, indx[max.indx[1], 1]]
				dd <- varout.mf[5, indx[max.indx[1], 1]]
									
				# check shape of membership function
				if (varout.mf[1, indx[max.indx[1], 1]] == 1){
					seqq <- seq(from = varout.mf[2, indx[max.indx[1], 1]], to = varout.mf[4, indx[max.indx[1], 1]], by = (varout.mf[4, indx[max.indx[1], 1]] - varout.mf[2, indx[max.indx[1], 1]]) / 10)
											
					for (j in 1:length(seqq)){
							if (seqq[j] < aa){
								temp.miu <- 1
							}
							else if (seqq[j] >= aa & seqq[j] < bb){
								temp.miu <- (seqq[j] - aa) / (bb - aa)
							}
							else if (seqq[j] >= bb & seqq[j] < cc) {
								temp.miu <- (seqq[j] - cc) / (bb - cc)
							}
							else 
								temp.miu <- 1
							
							if (temp.miu >= indx[max.indx[1], 2]) {
								def[k, 1] <- seqq[j]
								break
							}
					}					
			
				}
				
				else if (varout.mf[1, indx[max.indx[1], 1]] == 2){
					seqq <- seq(from = varout.mf[2, indx[max.indx[1], 1]], to = varout.mf[4, indx[max.indx[1], 1]], by = (varout.mf[4, indx[max.indx[1], 1]] - varout.mf[2, indx[max.indx[1], 1]]) / 10)
					
					for (j in 1:length(seqq)){
						if (seqq[j] < bb){
							temp.miu <- 1
						}
						else if (seqq[j] >= bb & seqq[j] < cc){
							temp.miu <- (seqq[j] - cc) / (bb - cc)
						}
						else if (seqq[j] > cc) {
							temp.miu <- 1
						}
														
						if (temp.miu >= indx[max.indx[1], 2]) {
							def[k, 1] <- seqq[j]
							break
						}
					}										
				}
				else if (varout.mf[1, indx[max.indx[1], 1]] == 3){
					seqq <- seq(from = varout.mf[2, indx[max.indx[1], 1]], to = varout.mf[4, indx[max.indx[1], 1]], by = (varout.mf[4, indx[max.indx[1], 1]] - varout.mf[2, indx[max.indx[1], 1]]) / 10)
					for (j in 1:length(seqq)){
						if (seqq[j] < aa){
							temp.miu <- 0
						}
						else if (seqq[j] >= aa & seqq[j] < bb){
							temp.miu <- (seqq[j] - aa) / (bb - aa)
						}
						else if (seqq[j] > cc) {
							temp.miu <- 1
						}
														
						if (temp.miu >= indx[max.indx[1], 2]) {
							def[k, 1] <- seqq[j]
							break
						}
					}
				}
				else if (varout.mf[1, indx[max.indx[1], 1]] == 4) {
					seqq <- seq(from = varout.mf[2, indx[max.indx[1], 1]], to = varout.mf[4, indx[max.indx[1], 1]], by = (varout.mf[4, indx[max.indx[1], 1]] - varout.mf[2, indx[max.indx[1], 1]]) / 10)
					
					for (j in 1:length(seqq)){
						if (seqq[j] < aa){
							temp.miu <- 0
						}
						else if (seqq[j] >= aa & seqq[j] < bb){
							temp.miu <- (seqq[j] - aa) / (bb - aa)
						}
						else if (seqq[j] >= bb & seqq[j] < cc) {
							temp.miu <- 1
						}
						else if (seqq[j] >= cc & seqq[j] < dd) {
							temp.miu <- (seqq[j] - dd) / (cc - dd)
						}
						else {
							temp.miu <- 0
						}
						
						if (temp.miu >= indx[max.indx[1], 2]) {
							def[k, 1] <- seqq[j]
							break
						}
					}
				}
				
				else if (varout.mf[1, indx[max.indx[1], 1]] == 5) {
					seqq <- seq(from = min(range.output) , to = max(range.output), by = (max(range.output) - min(range.output)) / 100)
					
					for (j in 1:length(seqq)){
						
						temp.miu <- exp(- 0.5 * (seqq[j] - aa) ^ 2 / bb ^ 2)
													
						if (temp.miu >= indx[max.indx[1], 2]) {
							def[k, 1] <- seqq[j]
							break
						}
						else {
							def[k, 1] <- (min(range.output, na.rm=TRUE) + max(range.output, na.rm = TRUE))/2
						}
					}
				}
				
				else if (varout.mf[1, indx[max.indx[1], 1]] == 6) {
					seqq <- seq(from = min(range.output) , to = max(range.output), by = (max(range.output) - min(range.output)) / 10)
					for (j in 1:length(seqq)){
						
						temp.miu <- 1 / (1 + exp(- aa * (seqq[j] - bb)))
													
						if (temp.miu >= indx[max.indx[1], 2]) {
							def[k, 1] <- seqq[j]
							break
						}
						else {
							def[k, 1] <- (min(range.output, na.rm=TRUE) + max(range.output, na.rm = TRUE))/2
						}
					}
				}
				else if (varout.mf[1, indx[max.indx[1], 1]] == 7) {
					seqq <- seq(from = min(range.output) , to = max(range.output), by = (max(range.output) - min(range.output)) / 10)
					for (j in 1:length(seqq)){
						
						temp.miu <- 1 / (1 + abs((seqq[j] - cc)/aa) ^ (2 * bb))
													
						if (temp.miu >= indx[max.indx[1], 2]) {
							def[k, 1] <- seqq[j]
							break
						}
						else {
							def[k, 1] <- (min(range.output, na.rm=TRUE) + max(range.output, na.rm = TRUE))/2
						}
					}
				}
			}	
			
			## procedure for type.defuz == 3 (last of maxima)
			else if (type.defuz == 3) {
				max.temp <- max(indx[, 2], na.rm = TRUE)
				max.indx <- which(indx[, 2] == max.temp)			
				
				aa <- varout.mf[2, indx[max.indx[1], 1]]
				bb <- varout.mf[3, indx[max.indx[1], 1]]
				cc <- varout.mf[4, indx[max.indx[1], 1]]
				dd <- varout.mf[5, indx[max.indx[1], 1]]
									
				# check shape of membership function
				if (varout.mf[1, indx[max.indx[1], 1]] == 1){
					seqq <- seq(from = varout.mf[2, indx[max.indx[1], 1]], to = varout.mf[4, indx[max.indx[1], 1]], by = (varout.mf[4, indx[max.indx[1], 1]] - varout.mf[2, indx[max.indx[1], 1]]) / 10)
											
					for (j in 1:length(seqq)){
							if (seqq[j] < aa){
								temp.miu <- 1
							}
							else if (seqq[j] >= aa & seqq[j] < bb){
								temp.miu <- (seqq[j] - aa) / (bb - aa)
							}
							else if (seqq[j] >= bb & seqq[j] < cc) {
								temp.miu <- (seqq[j] - cc) / (bb - cc)
							}
							else 
								temp.miu <- 1
							
							if (temp.miu >= indx[max.indx[1], 2]) {
								def[k, 1] <- seqq[j]
							}
					}					
			
				}
				
				else if (varout.mf[1, indx[max.indx[1], 1]] == 2){
					seqq <- seq(from = varout.mf[2, indx[max.indx[1], 1]], to = varout.mf[4, indx[max.indx[1], 1]], by = (varout.mf[4, indx[max.indx[1], 1]] - varout.mf[2, indx[max.indx[1], 1]]) / 10)
					
					for (j in 1:length(seqq)){
						if (seqq[j] < bb){
							temp.miu <- 1
						}
						else if (seqq[j] >= bb & seqq[j] < cc){
							temp.miu <- (seqq[j] - cc) / (bb - cc)
						}
						else if (seqq[j] > cc) {
							temp.miu <- 1
						}
														
						if (temp.miu >= indx[max.indx[1], 2]) {
							def[k, 1] <- seqq[j]
						}
					}										
				}
				else if (varout.mf[1, indx[max.indx[1], 1]] == 3){
					seqq <- seq(from = varout.mf[2, indx[max.indx[1], 1]], to = varout.mf[4, indx[max.indx[1], 1]], by = (varout.mf[4, indx[max.indx[1], 1]] - varout.mf[2, indx[max.indx[1], 1]]) / 10)
					for (j in 1:length(seqq)){
						if (seqq[j] < aa){
							temp.miu <- 0
						}
						else if (seqq[j] >= aa & seqq[j] < bb){
							temp.miu <- (seqq[j] - aa) / (bb - aa)
						}
						else if (seqq[j] > cc) {
							temp.miu <- 1
						}
														
						if (temp.miu >= indx[max.indx[1], 2]) {
							def[k, 1] <- seqq[j]
						}
					}
				}
				else if (varout.mf[1, indx[max.indx[1], 1]] == 4) {
					seqq <- seq(from = varout.mf[2, indx[max.indx[1], 1]], to = varout.mf[4, indx[max.indx[1], 1]], by = (varout.mf[4, indx[max.indx[1], 1]] - varout.mf[2, indx[max.indx[1], 1]]) / 10)
					
					for (j in 1:length(seqq)){
						if (seqq[j] < aa){
							temp.miu <- 0
						}
						else if (seqq[j] >= aa & seqq[j] < bb){
							temp.miu <- (seqq[j] - aa) / (bb - aa)
						}
						else if (seqq[j] >= bb & seqq[j] < cc) {
							temp.miu <- 1
						}
						else if (seqq[j] >= cc & seqq[j] < dd) {
							temp.miu <- (seqq[j] - dd) / (cc - dd)
						}
						else {
							temp.miu <- 0
						}
						
						if (temp.miu >= indx[max.indx[1], 2]) {
							def[k, 1] <- seqq[j]
						}
					}
				}
				
				else if (varout.mf[1, indx[max.indx[1], 1]] == 5) {
					seqq <- seq(from = min(range.output) , to = max(range.output), by = (max(range.output) - min(range.output)) / 10)
					
					for (j in 1:length(seqq)){
						
						temp.miu <- exp(- 0.5 * (seqq[j] - aa) ^ 2 / bb ^ 2)
													
						if (temp.miu >= indx[max.indx[1], 2]) {
							def[k, 1] <- seqq[j]
						}
						else {
							def[k, 1] <- (min(range.output, na.rm=TRUE) + max(range.output, na.rm = TRUE))/2
						}
					}
				}
				
				else if (varout.mf[1, indx[max.indx[1], 1]] == 6) {
					seqq <- seq(from = min(range.output) , to = max(range.output), by = (max(range.output) - min(range.output)) / 10)
					for (j in 1:length(seqq)){
						
						temp.miu <- 1 / (1 + exp(- aa * (seqq[j] - bb)))
													
						if (temp.miu >= indx[max.indx[1], 2]) {
							def[k, 1] <- seqq[j]
						}
						else {
							def[k, 1] <- (min(range.output, na.rm=TRUE) + max(range.output, na.rm = TRUE))/2
						}
					}
				}
				else if (varout.mf[1, indx[max.indx[1], 1]] == 7) {
					seqq <- seq(from = min(range.output) , to = max(range.output), by = (max(range.output) - min(range.output)) / 10)
					for (j in 1:length(seqq)){
						
						temp.miu <- 1 / (1 + abs((seqq[j] - cc)/aa) ^ (2 * bb))
													
						if (temp.miu >= indx[max.indx[1], 2]) {
							def[k, 1] <- seqq[j]
						}
						else {
							def[k, 1] <- (min(range.output, na.rm=TRUE) + max(range.output, na.rm = TRUE))/2
						}
					}
				}	
			}
			
			## procedure for type.defuz == 4 (mean of maxima)
			else if (type.defuz == 4) {
				max.temp <- max(indx[, 2], na.rm = TRUE)
				max.indx <- which(indx[, 2] == max.temp)			
							
				def[k, 1] <- 0.5 * (max(varout.mf[2:5, indx[max.indx[1], 1]], na.rm = TRUE) + (varout.mf[2, indx[max.indx[1], 1]]))
			}
			
			else if (type.defuz == 5) {
				for (i in 1 : nrow(indx)){										
					#calculate modified centroid
					# update for gaussian
					if (varout.mf[1, indx[i, 1]] == 5 || varout.mf[1, indx[i, 1]] == 6 || varout.mf[1, indx[i, 1]] == 7) {
						## center point
						
						av.point <- varout.mf[2, indx[i, 1]]
						
						## indx is fired miu.rule/rule
						temp1 <- indx[i, 2] * av.point * varout.mf[3, indx[i, 1]]
						temp2 <- indx[i, 2] * varout.mf[3, indx[i, 1]]

						temp <- temp + temp1
						temp.d <- temp.d + temp2
					}
					
					else {				
						average <- (varout.mf[2, indx[i, 1]] + max(varout.mf[2 : 5, indx[i, 1]], na.rm = TRUE)) / 2
						av.point <- varout.mf[3, indx[i, 1]]
						# check first which one greater between indx[i, 2] with the value of MF (for centroid)
						temp1 <- indx[i, 2] * av.point #* abs((average - max(varout.mf[2 : 5, indx[i, 1]], na.rm = TRUE))) / 2
						temp2 <- indx[i, 2] # * (abs((average - max(varout.mf[2 : 5, indx[i, 1]], na.rm = TRUE))) / 2)
						temp <- temp + temp1
						temp.d <- temp.d + temp2
					}
					cum[k, i] <- temp
					div[k, i] <- temp.d   					
				}
				
				if (sum(div[k, ]) == 0){
					def[k, 1] <- (min(range.output, na.rm=TRUE) + max(range.output, na.rm = TRUE))/2
				}
				else{
					def.temp[k, 1] <- sum(cum[k, ]) / sum(div[k, ])
					if (def.temp[k, 1] <= min(range.output, na.rm=TRUE)){
						def[k, 1] <- (min(range.output, na.rm=TRUE) + max(range.output, na.rm = TRUE))/2
					}
					else if (def.temp[k, 1] >= max(range.output, na.rm=TRUE)){
						def[k, 1] <- (min(range.output, na.rm=TRUE) + max(range.output, na.rm = TRUE))/2
					}
					else {
						def[k, 1] <- def.temp[k, 1]
					}
				}
			}
			
		}
		else {
			def[k, 1] <- (min(range.output, na.rm=TRUE) + max(range.output, na.rm = TRUE))/2
		}
	}	

}


#TSK type
else if (type.model == 2){
	
	## check the linear equation on consequent part
	if (is.null(func.tsk)){
		stop("please define the linear equation as the consequent part on fuzzy IF-THEN rules")
	}
		
	for (k in 1 : nrow(data)){
		data.k <-(data[k, ])
		data.m <- as.matrix(data.k)
		
		if (ncol(func.tsk) > 1){
			func.tsk.var <- func.tsk[, 1: ncol(func.tsk) - 1]
			func.tsk.cont <- t(t(func.tsk[, ncol(func.tsk)]))
		
			ff <- func.tsk.var %*% data.m + func.tsk.cont
		}
		else if (ncol(func.tsk) == 1){
			ff <- func.tsk 
		}
		
		miu.rule.t <- as.matrix(miu.rule[k, ])
		cum <- sum(miu.rule.t * ff)	
		div <- sum(miu.rule.t)
		
		def[k, 1] <- cum / div
		if (div == 0)
			def[k, 1] <- 0
		else
		{
			def[k, 1] <- cum / div
			if (def[k, 1] > max(range.output, na.rm=TRUE))
				def[k, 1] <- max(range.output, na.rm=TRUE)
			else if (def[k, 1] < min(range.output, na.rm=TRUE))
				def[k, 1] <- min(range.output, na.rm=TRUE)
		}
	}
}

return(def)
}

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
