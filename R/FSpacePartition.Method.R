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
#' This is the internal function that implements the model proposed by L.X. Wang and J.M. 
#' Mendel. It is used to solve regression task. Users do not need to call it directly,
#' but just use \code{\link{frbs.learn}} and \code{\link{predict}}
#' 
#' The fuzzy rule-based system for learning from L.X. Wang and J.M. 
#' Mendel's paper is implemented in this function. For the learning process, there are three stages: 
#' Step 1 divides the input and output spaces of the given numerical data into fuzzy regions. 
#' Step 2 generates fuzzy IF-THEN rules from the training data. 
#' Step 3 determines a degree of each rule. 
#' In the prediction phase, there are four steps: fuzzification, checking the rules, inference, and defuzzification.
#' 
#' @title WM model building 
#'
#' @param data.train a matrix(m x n) of normalized data for the training process, where m is the number of instances and 
#'        n is the number of variables; the last column is the output variable. 
#'        Note the data must be normalized between 0 and 1. 
#' @param num.labels a matrix(1 x n), whose elements represent the number of labels (fuzzy terms); 
#' n is the number of variables.
#' @param type.mf the type of the membership function. See \code{\link{frbs.learn}}.
#' @param type.tnorm a value which represents the type of t-norm. See \code{\link{inference}}.
#' @param type.implication.func a value representing type of implication function. Let us consider a rule, a -> b,  
#' \itemize{
#' \item \code{DIENES_RESHER} means (b > 1 - a? b : 1 - a).
#' \item \code{LUKASIEWICZ} means (b < a ? 1 - a + b : 1).
#' \item \code{ZADEH} means (a < 0.5 || 1 - a > b ? 1 - a : (a < b ? a : b)).
#' \item \code{GOGUEN} means (a < b ? 1 : b / a).
#' \item \code{GODEL} means (a <= b ? 1 : b).
#' \item \code{SHARP} means (a <= b ? 1 : 0).
#' \item \code{MIZUMOTO} means (1 - a + a * b).
#' \item \code{DUBOIS_PRADE} means (b == 0 ? 1 - a : (a == 1 ? b : 1)).
#' \item \code{MIN} means (a < b ? a : b)
#' }
#' @param classification a boolean representing whether it is a classification problem or not.
#' @param range.data a matrix representing interval of data.
#' @seealso \code{\link{frbs.learn}}, \code{\link{predict}} and \code{\link{frbs.eng}}.
# @return an object of type \code{frbs}. Please have a look at An \code{\link{frbs-object}} for looking its complete components.
#' @references 
#' L. X. Wang and J.M. Mendel, "Generating fuzzy rule by learning from examples", IEEE Trans. Syst., Man, and Cybern.,
#' vol. 22, no. 6, pp. 1414 - 1427 (1992).
# @export
WM <- function(data.train, num.labels, type.mf = "GAUSSIAN", type.tnorm = "PRODUCT", type.implication.func = "ZADEH", classification = FALSE, range.data = NULL) {
	
	if (!is.null(range.data)) {
		 range.data = range.data
	}
	else if (classification == FALSE && is.null(range.data)) {
		range.data <- matrix(nrow = 2, ncol = ncol(data.train))
		range.data[1, ] <- 0
		range.data[2, ] <- 1
	}
	else {
		## make range of data according to class on data training
	    range.data.out <- matrix(c(0.5001, num.labels[1, ncol(num.labels)] + 0.4999), nrow = 2)
		range.data.inp <- matrix(nrow = 2, ncol = (ncol(data.train) - 1))
		range.data.inp[1, ] <- 0
		range.data.inp[2, ] <- 1
		range.data <- cbind(range.data.inp, range.data.out)
	}
	
	## build the names of fuzzy values
	fuzzyTerm <- create.fuzzyTerm(num.labels)
	names.fvalinput <- fuzzyTerm$names.fvalinput
	names.fvaloutput <- fuzzyTerm$names.fvaloutput
	names.fvalues <- c(names.fvalinput, names.fvaloutput)
	
	#get number of row and column
	nrow.data.train <- nrow(data.train)
	ncol.data.train <- ncol(data.train)
	
	##initialize
	seg <- matrix(ncol = 1)
	center.value <- matrix(ncol = 1)
	var.mf <- matrix(0, nrow = 5, ncol = sum(num.labels))
	
	## use data.train as data
	data <- data.train
	
	## get number of input variables as number of column of data
	num.varinput <- ncol.data.train


	####### Wang and Mendel's Steps begin ####################
	###Step 1: Divide the Input and Output Spaces Into Fuzzy Regions
	
	## initialize
	jj <- 0
	
	## loop for all variable
	for (i in 1 : ncol.data.train){
		## initialize
		seg <- num.labels[1, i]
		kk <- 1
			
		## Make the depth on each fuzzy value, assumed it's similar on all region. 
		delta.point <- (range.data[2, i] - range.data[1, i]) / (seg + 1)
		
			##contruct matrix of parameter of membership function (var.mf) on each variable (max(var) is ncol.data.train
			for (j in 1 : num.labels[1, i]){
			
				##counter to continue index of column var.mf					
				jj <- jj + 1
				
				## type.mf <- 1 means using trapezoid and triangular
				if (type.mf == 1 || type.mf == "TRIANGLE") {
				
					delta.tri <- (range.data[2, i] - range.data[1, i]) / (seg - 1)
					## on the left side
					if (kk %% num.labels[1, i] == 1){   
						var.mf[1, jj] <- 1	
						var.mf[2, jj] <- range.data[1, i]
						var.mf[3, jj] <- var.mf[2, jj]
						var.mf[4, jj] <- var.mf[3, jj] + delta.tri
						var.mf[5, jj] <- NA						
					}
					
					## on the right side
					else if (kk %% num.labels[1, i] == 0){
						var.mf[1, jj] <- 1	
						var.mf[4, jj] <- range.data[2, i]
						var.mf[3, jj] <- var.mf[4, jj] 
						var.mf[2, jj] <- var.mf[3, jj] - delta.tri
						var.mf[5, jj] <- NA						
					}
					
					## on the middle
					else{
						var.mf[1, jj] <- 1	
						var.mf[3, jj] <- range.data[1, i] + (j - 1) * delta.tri
						var.mf[2, jj] <-  var.mf[3, jj] - delta.tri
						var.mf[4, jj] <- var.mf[3, jj] + delta.tri
						var.mf[5, jj] <- NA						
					}	
				}
				## type.mf == 2 means we use trapezoid
				else if (type.mf == 2 || type.mf == "TRAPEZOID") {
					delta.tra <- (range.data[2, i] - range.data[1, i]) / (seg + 2)
					
					## on the left side
					if (kk %% num.labels[1, i] == 1){   
						var.mf[1, jj] <- 2	
						var.mf[2, jj] <- range.data[1, i]
						var.mf[3, jj] <- var.mf[2, jj] + delta.tra
						var.mf[4, jj] <- var.mf[3, jj] + delta.tra
						var.mf[5, jj] <- NA												
					}
					
					## on the right side
					else if (kk %% num.labels[1, i] == 0){
						var.mf[1, jj] <- 3	
						var.mf[4, jj] <- range.data[2, i]
						var.mf[3, jj] <- var.mf[4, jj] - delta.tra
						var.mf[2, jj] <- var.mf[3, jj] - delta.tra
						var.mf[5, jj] <- NA												
					}
					
					## on the middle
					else{
						var.mf[1, jj] <- 4	
						var.mf[2, jj] <- range.data[1, i] + (j - 1) * 1.15 * delta.tra
						var.mf[3, jj] <- var.mf[2, jj] + delta.tra
						var.mf[4, jj] <- var.mf[3, jj] + 0.5 * delta.tra					
						var.mf[5, jj] <- var.mf[4, jj] + delta.tra											
					}
				}
				##Type 5: Gaussian
				##parameter=(mean a, standard deviation b)
				##a=var.mf[2,]
				##b=var.mf[3,]
				else if (type.mf == 3 || type.mf == "GAUSSIAN") {
					delta.gau <- (range.data[2, i] - range.data[1, i]) / (seg - 1)
					##On the left side
					if (kk %% num.labels[1, i] == 1){   
						var.mf[1, jj] <- 5	
						var.mf[2, jj] <- range.data[1, i]
						var.mf[3, jj] <- 0.35 * delta.gau
						var.mf[4, jj] <- NA
						var.mf[5, jj] <- NA
					}
					##On the right side
					else if (kk %% num.labels[1, i] == 0){
						var.mf[1, jj] <- 5	
						var.mf[2, jj] <- range.data[2, i]
						var.mf[3, jj] <- 0.35 * delta.gau
						var.mf[4, jj] <- NA
						var.mf[5, jj] <- NA
					}
					## On the middle
					else {
					var.mf[1, jj] <- 5	
					var.mf[2, jj] <- range.data[1, i] + (j - 1) * delta.gau
					var.mf[3, jj] <- 0.35 * delta.gau
					var.mf[4, jj] <- NA
					var.mf[5, jj] <- NA
					}
				}
				##Type 6: Sigmoid/logistic
				##parameter=(gamma,c)
				##gamma=var.mf[2,]
				##c=var.mf[3,]
				else if (type.mf == 4 || type.mf == "SIGMOID") {
					delta.sig <- (range.data[2, i] - range.data[1, i]) / (seg + 1)
					##On the left side
					if (kk %% num.labels[1, i] == 1){   
						var.mf[1, jj] <- 6	
						var.mf[2, jj] <- range.data[1, i]
						var.mf[3, jj] <- delta.sig
						var.mf[4, jj] <- NA
						var.mf[5, jj] <- NA
					}
					##On the right side
					else if (kk %% num.labels[1, i] == 0){
						var.mf[1, jj] <- 6	
						var.mf[2, jj] <- range.data[2, i]
						var.mf[3, jj] <- delta.sig
						var.mf[4, jj] <- NA
						var.mf[5, jj] <- NA
					}
					## On the middle
					else {
					var.mf[1, jj] <- 6	
					var.mf[2, jj] <- range.data[1, i] + (j - 1) * delta.sig
					var.mf[3, jj] <- delta.sig
					var.mf[4, jj] <- NA
					var.mf[5, jj] <- NA
					}
				}
				##checking for type 7: Generalized Bell
				##parameter=(a,b,c)
				##a=var.mf[2,]
				##b=var.mf[3,]
				##c=var.mf[4,]
				else if (type.mf == 5 || type.mf == "BELL") {
					##On the left side
					if (kk %% num.labels[1, i] == 1){   
						var.mf[1, jj] <- 7	
						var.mf[2, jj] <- 0.6 * delta.point
						var.mf[3, jj] <- 0.6 * delta.point
						var.mf[4, jj] <- range.data[1, i]
						var.mf[5, jj] <- NA
					}
					##On the right side
					else if (kk %% num.labels[1, i] == 0){
						var.mf[1, jj] <- 7	
						var.mf[2, jj] <- 0.6 * delta.point
						var.mf[3, jj] <- 0.6 * delta.point
						var.mf[4, jj] <- range.data[2, i]
						var.mf[5, jj] <- NA
					}
					## On the middle
					else {
					
					var.mf[1, jj] <- 7	
					var.mf[2, jj] <- 0.6 * delta.point
					var.mf[3, jj] <- 0.6 * delta.point
					var.mf[4, jj] <- range.data[1, i] + delta.point * j
					var.mf[5, jj] <- NA
					}
				}						
				kk <- kk + 1
			}
	}

	if (classification == TRUE){
		range.data.out <- range.data[, ncol(range.data), drop = FALSE]
		num.fvaloutput <- num.labels[1, ncol(num.labels)]
		seg <- num.labels[1, ncol(num.labels)]
		delta.tra <- (range.data.out[2, 1] - range.data.out[1, 1]) / seg
		ncol.var.mf <- ncol(var.mf)
	
		k <- 1
		for (j in (ncol.var.mf - num.fvaloutput + 1) : ncol.var.mf){
			##On the left side	
			var.mf[1, j] <- 4	
			var.mf[2, j] <- range.data.out[1, 1] + (k - 1) * delta.tra 
			var.mf[3, j] <- var.mf[2, j] 
			var.mf[4, j] <- var.mf[3, j] + delta.tra
			var.mf[5, j] <- var.mf[4, j] 					
			k <- k + 1
		}
	}
	
	## Step 2: Generate Fuzzy Rules from Given Data Pairs.
	## Step 2a: Determine the degree of data pairs.
	## MF is matrix membership function. The dimension of MF is n x m, where n is number of data and m is num.labels * input variable (==ncol(var.mf))
		
	## get degree of membership by fuzzification
	MF <- fuzzifier(data, num.varinput, num.labels, var.mf)	
	
	colnames(MF) <- c(names.fvalues)
		
	####get max value of degree on each variable to get one rule.
	
	MF.max <- matrix(0, nrow = nrow(MF), ncol = ncol(MF))
	k <- 1
	for (i in 1 : length(num.labels)){
		start <- k
		end <- start + num.labels[1, i] - 1
		MF.temp <- MF[, start : end]
		
		for (m in 1 : nrow(MF)){
			max.MFTemp <- max(MF.temp[m, ], na.rm = TRUE)
			max.loc <- which.max(MF.temp[m, ])
			MF.max[m, k + max.loc - 1] <- max.MFTemp		
		}
		k <- k + num.labels[1, i]	
	}
	
	colnames(MF.max) <- c(names.fvalues)
	rule.matrix <- MF.max
		
	## Step III
	## determine the degree of the rule
	degree.rule <- matrix(nrow = nrow(rule.matrix), ncol =1)
	degree.ante <- matrix(nrow = nrow(rule.matrix), ncol =1)

	## calculate t-norm in antecedent part and implication function with consequent part
	for (i in 1:nrow(rule.matrix)){
		temp.ant.rule.degree <- c(1)
		max.ant.indx <- c(ncol(rule.matrix) - num.labels[1, ncol(num.labels)])
		for (j in 1 : max.ant.indx){
			if (rule.matrix[i, j] != 0){
				if (type.tnorm == "PRODUCT"){
					temp.ant.rule.degree <- temp.ant.rule.degree * rule.matrix[i, j]					
				}			
				else if (type.tnorm == "MIN"){
					temp.ant.rule.degree <- min(temp.ant.rule.degree, rule.matrix[i, j])
				}
				else if (type.tnorm == "HAMACHER"){
					temp.ant.rule.degree <- (temp.ant.rule.degree * rule.matrix[i, j]) / (temp.ant.rule.degree + rule.matrix[i, j] - temp.ant.rule.degree * rule.matrix[i, j])
				}
				else if (type.tnorm == "YAGER"){
					temp.ant.rule.degree <- 1 - min(1, ((1 - temp.ant.rule.degree) + (1 - rule.matrix[i, j])))
				}
				else if (type.tnorm == "BOUNDED"){
					temp.ant.rule.degree <- max(0, (temp.ant.rule.degree + rule.matrix[i, j] - 1))
				}
			}
		}
		degree.ante[i] <- temp.ant.rule.degree
		
		for (k in (max.ant.indx + 1) : ncol(rule.matrix)){
			if (rule.matrix[i, k] != 0){
				temp.rule.degree <- calc.implFunc(degree.ante[i], rule.matrix[i, k], type.implication.func)
			}
		}
		degree.rule[i] <- temp.rule.degree
	}

	rule.matrix[rule.matrix > 0] <- 1
	rule.matrix.bool <- rule.matrix

	## find the same rows on matrix rule
	temp <- matrix(nrow = 1, ncol = ncol(rule.matrix.bool))
	kkk <- 1
	rule.red.mat <- rule.matrix
	degree.rule.temp <- degree.rule
	
	for (i in 1 : (nrow(rule.matrix.bool) - 1)){
		temp.rule <- matrix(rule.matrix.bool[i, ], nrow = 1)
		ii <- i + 1
		for (j in ii : nrow(rule.matrix.bool)){
			sim.check = dist(rbind(temp.rule, matrix(rule.matrix.bool[j, ], nrow = 1)))
			if (sim.check == 0){
				if (degree.rule[i] >= degree.rule[j]){
					rule.red.mat[j, ] <- NA
					degree.ante[j] <- NA
					degree.rule.temp[j] <- NA
				}
				else {
					rule.red.mat[i, ] <- NA
					degree.ante[i] <- NA
					degree.rule.temp[i] <- NA
				}
			}
		}
	}	
	
	degree.ante <- na.omit(degree.ante)
	rule.red.mat <- na.omit(rule.red.mat)
	degree.rule <- na.omit(degree.rule.temp)
	
	## Build Rule database (in list of string) from matrix rule.matrix
	## rule.data.num is the numbers representing the sequence of string of variable names 	
	rule.data.num <- matrix(nrow = nrow(rule.red.mat), ncol = ncol(num.labels))	
	for (i in 1 : nrow(rule.red.mat)){
		k <- 0
		for (j in 1 : ncol(rule.red.mat)){
			if (rule.red.mat[i, j] != 0){
				k <- k + 1
				rule.data.num[i, k] <- j 
				if (k >= ncol(num.labels))
				break
			}	
		}
	}

	## check is.na on rule.data.num and degree.ante
	for (i in 1 : nrow(rule.data.num)){
		if(any(is.na(rule.data.num[i, ])== TRUE)){
			degree.ante[i, ] <- NA 
			degree.rule[i, ] <- NA
		}	
	}	

	rule.data.num <- na.omit(rule.data.num)
	degree.ante <- na.omit(degree.ante)
	degree.rule <- na.omit(degree.rule)

	
	## build the rule into list of string
	res <- generate.rule(rule.data.num, num.labels)
	rule <- res$rule
	
	#############################################################
	### Collect the Input data for "frbs for testing"
	#############################################################
		
	## indx is index of output variable on num.labels
	indx <- length(num.labels)
	length.namesvar <- length(names.fvalues)
		
	## cut membership function of output variable on var.mf
	ll <- (ncol(var.mf) - num.labels[1, indx])
	varinp.mf <- var.mf[, 1 : ll, drop = FALSE]
	colnames(varinp.mf) <- names.fvalinput
	## get varout.mf
	varout.mf <- var.mf[, (ll + 1) : ncol(var.mf), drop = FALSE]
	colnames(varout.mf) <- names.fvaloutput
	
	## clean rule
	rule <- na.omit(rule)
	
	## range of original data
	range.data.ori <- range.data
	
	mod <- list(num.labels = num.labels, varout.mf = varout.mf, rule = rule, varinp.mf = varinp.mf, degree.ante = degree.ante, rule.data.num = rule.data.num, 
			   degree.rule = degree.rule, range.data.ori = range.data.ori, type.mf = type.mf, type.tnorm = type.tnorm, type.implication.func = type.implication.func)

	return (mod)
}  

#' This is the internal function that implements the fuzzy rule-based classification 
#' system using Chi's technique (FRBCS.CHI). It is used to solve classification tasks. 
#' Users do not need to call it directly,
#' but just use \code{\link{frbs.learn}} and \code{\link{predict}}. This method is
#' suitable only for classification problems.
#' 
#' This method was proposed by Zheru Chi, Hong Yan, and Tuan Pham that extends
#' Wang and Mendel's method. The method consists of the following five steps:
#' \itemize{
#' \item Step 1: Fuzzify the input space.
#' \item Step 2: Generate fuzzy rules from given training data pairs.
#' \item Step 3: Assign a degree to each rule.
#' \item Step 4: Create a combined rule bank.
#' \item Step 5: Determine the mapping by using a defuzzification method.
#' }
#'
#' @title FRBCS.CHI model building 
#'
#' @param range.data a matrix (2 x n) containing the range of the normalized data, where n is the number of variables, and
#' first and second rows are the minimum and maximum values, respectively. 
#' @param data.train a matrix (m x n) of normalized data for the training process, where m is the number of instances and 
#' n is the number of variables; the last column is the output variable. Note the data must be normalized between 0 and 1. 
#' @param num.labels a matrix (1 x n), whose elements represent the number of labels (fuzzy terms); 
#' n is the number of variables.
#' @param num.class an integer number representing the number of labels (fuzzy terms).
#' @param type.mf the type of the shape of the membership functions. See \code{\link{fuzzifier}}.
#' @param type.tnorm the type of t-norm. See \code{\link{inference}}.
#' @param type.snorm the type of s-norm. See \code{\link{inference}}.
#' @param type.implication.func the type of implication function. See \code{\link{WM}}.
#' @seealso \code{\link{FRBCS.eng}}, \code{\link{frbs.learn}}, and \code{\link{predict}}
# @return a list of the model data. Please have a look at \code{\link{frbs.learn}} for looking its complete components. 
#' @references
#' Z. Chi, H. Yan, T. Pham, "Fuzzy algorithms with applications to image processing 
#' and pattern recognition", World Scientific, Singapore (1996).
# @export
FRBCS.CHI <- function(range.data, data.train, num.labels, num.class, type.mf = "TRIANGLE", type.tnorm = "MIN", type.snorm = "MAX", type.implication.func = "ZADEH") {

	num.labels[1, ncol(num.labels)] <- num.class
	
	## generate initial model using WM
	mod <- WM(data.train, num.labels, type.mf, type.tnorm = type.tnorm, type.implication.func = type.implication.func, classification = TRUE)
	
	names.fvaloutput <- as.matrix(colnames(mod$varout.mf))
	rule <- mod$rule
	rule.constant <- mod$rule
	degree.ante <- mod$degree.ante
	 for (i in 1 : nrow(names.fvaloutput)){
		 for (j in 1 : nrow(rule.constant)){
			 if (rule[j, ncol(rule)] == names.fvaloutput[i])
				 rule.constant[j, ncol(rule.constant)] <- i
		 }
	 }
	
	classes <- matrix(strtoi(rule.constant[, ncol(rule.constant)])) 
	
	grade.cert <- cbind(classes, degree.ante)
	
	mod$class <- classes
	mod$rule[, ncol(mod$rule)] <- classes
	mod$grade.cert <- grade.cert
	mod$varout.mf <- NULL
	mod$type.tnorm <- type.tnorm
	mod$type.snorm <- type.snorm
	mod$type.model <- "FRBCS"
	mod$type.mf <- type.mf
	mod$type.implication.func
	return (mod)	
}

#' This is the internal function that implements the fuzzy rule-based classification 
#' system with weight factor (FRBCS.W). It is used to solve classification tasks. 
#' Users do not need to call it directly,
#' but just use \code{\link{frbs.learn}} and \code{\link{predict}}. This method is
#' suitable only for classification problems.
#' 
#' This method is adopted from Hisao Ishibuchi and Tomoharu Nakashima's paper. 
#' Each fuzzy IF-THEN rule consists of antecedent linguistic values and a single consequent class with certainty grades 
#' (weights). The antecedent part is determined by a grid-type fuzzy partition from 
#' the training data. The consequent class is defined as the dominant class in 
#' the fuzzy subspace corresponding to the antecedent part of each fuzzy IF-THEN rule and 
#' the certainty grade is calculated from the ratio among the consequent class. 
#' A class of the new instance is determined by the consequent class of the rule with 
#' the maximal product of the compatibility grade and the certainty grade. 
#'
#' @title FRBCS.W model building 
#'
#' @param range.data a matrix (2 x n) containing the range of the normalized data, where n is the number of variables, and
#' first and second rows are the minimum and maximum values, respectively. 
#' @param data.train a matrix (m x n) of normalized data for the training process, where m is the number of instances and 
#' n is the number of variables; the last column is the output variable. Note the data must be normalized between 0 and 1. 
#' @param num.labels a matrix (1 x n), whose elements represent the number of labels (fuzzy terms); 
#' n is the number of variables.
#' @param num.class an integer number representing the number of labels (fuzzy terms).
#' @param type.mf the type of the shape of the membership functions.
#' @param type.tnorm the type of t-norm. See \code{\link{inference}}.
#' @param type.snorm the type of s-norm. See \code{\link{inference}}.
#' @param type.implication.func the type of implication function. See \code{\link{WM}}.
#' @seealso \code{\link{FRBCS.eng}}, \code{\link{frbs.learn}}, and \code{\link{predict}}
# @return a list of the model data. Please have a look at \code{\link{frbs.learn}} for looking its complete components. 
#' @references
#' H. Ishibuchi and T. Nakashima, "Effect of rule weights in fuzzy rule-based classification systems", 
#' IEEE Transactions on Fuzzy Systems, vol. 1, pp. 59 - 64 (2001).
# @export
FRBCS.W <- function(range.data, data.train, num.labels, num.class, type.mf, type.tnorm = "MIN", type.snorm = "MAX", type.implication.func = "ZADEH") {

	num.labels[1, ncol(num.labels)] <- num.class
	
	## generate initial model using WM
	mod <- WM(data.train, num.labels, type.mf, type.tnorm = type.tnorm, type.implication.func = type.implication.func, classification = TRUE)
	
	names.fvaloutput <- as.matrix(colnames(mod$varout.mf))
	rule <- mod$rule
	rule.constant <- mod$rule
	degree.ante <- mod$degree.ante
	 for (i in 1 : nrow(names.fvaloutput)){
		 for (j in 1 : nrow(rule.constant)){
			 if (rule[j, ncol(rule)] == names.fvaloutput[i])
				 rule.constant[j, ncol(rule.constant)] <- i
		 }
	 }
	
	classes <- matrix(strtoi(rule.constant[, ncol(rule.constant)])) 
	grade.cert.temp <- cbind(classes, degree.ante)
	class.fact <- factor(grade.cert.temp[, 1], exclude = "")
	beta.class <- as.double(as.matrix(aggregate(grade.cert.temp[, 2], by = list(class.fact), sum)))
	beta.class <- matrix(beta.class, nrow = num.class)	
	
	CF <- matrix()
	 
	 for(i in 1 : nrow(beta.class)){
		CF[i] <- 1 + (beta.class[i, 2] - (sum(beta.class[, 2]) - (beta.class[i, 2] / (nrow(beta.class) - 1))))/sum(beta.class[, 2])
	 }
	grade.cert <- grade.cert.temp 
	for(i in 1 : nrow(grade.cert.temp)){
		grade.cert[i, 2] <- CF[grade.cert[i,1]] 
	}
	mod$class <- classes
	mod$rule[, ncol(mod$rule)] <- classes
	mod$grade.cert <- grade.cert
	mod$varout.mf <- NULL
	mod$type.tnorm <- type.tnorm
	mod$type.snorm <- type.snorm
	mod$type.model <- "FRBCS"
	mod$type.mf <- type.mf
	mod$type.implication.func <- type.implication.func
	
	return (mod)
}  
