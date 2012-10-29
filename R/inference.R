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
# @export
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
				
			##Inisialitation for calculating MF of antecendet
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
			
			else
			val.antecedent <- MF[k, temp[1]]

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
			else
				val.antecedent.b <- MF[k, temp[j + 2]]
			
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

