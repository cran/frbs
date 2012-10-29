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
#' Mendel. Users do not need to call it directly,
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
#' @param range.data a matrix(2 x n) containing the range of the data, where n is the number of variables, and
#' first and second rows are the minimum and maximum values, respectively. 
#' @param data.train a matrix(m x n) of data for the training process, where m is the number of instances and 
#' n is the number of variables; the last column is the output variable.
#' @param num.labels a matrix(1 x n), whose elements represent the number of labels (fuzzy terms); 
#' n is the number of variables.
#' @param type.mf the type of the membership function. See \code{\link{frbs.learn}}.
#' @param classification a boolean representing whether it is a classification problem or not.
#' @seealso \code{\link{frbs.learn}}, \code{\link{predict}} and \code{\link{frbs.eng}}.
# @return an object of type \code{frbs}. Please have a look at An \code{\link{frbs-object}} for looking its complete components.
#' @references 
#' L. X. Wang, "Fuzzy systems are universal approximators," 
#' in Proc. IEEE Int. Conf. Fuzzy Systems, San Diego, CA (1992). 
#'
#' L. X. Wang and J.M. Mendel, "Generating fuzzy rule by learning from examples", in Proc. 6th Int. Symp.
#' Intelligent Control (Washington, DC), 1991, pp. 263-268; also IEEE Trans. Syst., Man, Cybern.,
#' vol. 22, No. 6, (1992).
#'
#' L. X. Wang and J. M. Mendel, "Fuzzy basis function, universal approximation, and orthogonal least squares learning," 
#' IEEE Int. Conf. Neural Network, vol. 3 no. 5, pp. 807 - 814 (1992).
# @export
WM <- function(range.data, data.train, num.labels = 5, type.mf = 3, classification = FALSE) {
	
	#get names of variable
	num.inp.var <- sum(num.labels[1, 1 : (ncol(num.labels) - 1)])
	seq.inp.num <- seq(from = 1, to = num.inp.var, by = 1)
	names.inp.var <- paste("a", seq.inp.num, sep = ".")

	num.out.var <- num.labels[1, ncol(num.labels)]
	seq.out.num <- seq(from = 1, to = num.out.var, by = 1)
	names.out.var <- paste("c", seq.out.num, sep = ".")

	names.variable <- c(names.inp.var, names.out.var)
	
	#get number of row and column
	nrow.data.train <- nrow(data.train)
	ncol.data.train <- ncol(data.train)
	
	##initialize
	min.data.train <- matrix(ncol = 1)
	max.data.train <- matrix(ncol = 1)
	seg <- matrix(ncol = 1)
	center.value <- matrix(ncol = 1)
	var.mf <- matrix(0, nrow = 5, ncol = sum(num.labels))
	num.fvalinput <- matrix(nrow = 1, ncol = ncol.data.train)
	
	## use data.train as data
	data <- data.train
	
	## get number of input variables as number of column of data
	num.varinput <- ncol.data.train

	## get matrix of number of fuzzy value on each variabel
	## for example (3,2) means there 2 variable where the first variable has 3 fuzzy value and the scond one has 2 fuzzy value
	num.fvalinput <- num.labels
	
	####### Wang and Mendel's Steps begin ####################
	###Step 1: Divide the Input and Output Spaces Into Fuzzy Regions
	
	## initialize
	jj <- 0
	
	## loop for all variable
	for (i in 1 : ncol.data.train){
		## initialize
		min.data.train <- min(data.train[, i], na.rm=TRUE)
		max.data.train <- max(data.train[, i], na.rm=TRUE)
		seg <- num.labels[1, i]
		kk <- 1
			
		## Make the depth on each fuzzy value, assumed it's similar on all region. 
		delta.point <- (range.data[2, i] - range.data[1, i]) / (seg + 1)
		
			##contruct matrix of parameter of membership function (var.mf) on each variable (max(var) is ncol.data.train
			for (j in 1 : num.labels[1, i]){
			
				##counter to continue index of column var.mf					
				jj <- jj + 1
				
				## type.mf <- 1 means using trapezoid and triangular
				if (type.mf == 1) {
				
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
				else if (type.mf == 2) {
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
				else if (type.mf == 3) {
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
				else if (type.mf == 4) {
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
				else if (type.mf == 5) {
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
		range.data.out <- matrix(range.data[, ncol(range.data)], nrow = 2)
		seg <- num.labels[1, ncol(num.labels)]
		delta.tra <- (range.data.out[2, 1] - range.data.out[1, 1]) / seg
		ncol.var.mf <- ncol(var.mf)
	
		k <- 1
		for (j in (ncol.var.mf - num.out.var + 1) : ncol.var.mf){
			##On the left side	
			var.mf[1, j] <- 4	
			var.mf[2, j] <- range.data.out[1, 1] + (k - 1) * delta.tra 
			var.mf[3, j] <- var.mf[2, j] #+ delta.tra
			var.mf[4, j] <- var.mf[3, j] + delta.tra
			var.mf[5, j] <- var.mf[4, j] 					
			k <- k + 1
		}
	}
	

	## Step 2: Generate Fuzzy Rules from Given Data Pairs.
	## Step 2a: Determine the degree of data pairs.
	## MF is matrix membership function. The dimension of MF is n x m, where n is number of data and m is num.labels * input variable (==ncol(var.mf))
		
	## get degree of membership by fuzzification
	MF <- fuzzifier(data, num.varinput, num.fvalinput, var.mf)	
	
	colnames(MF) <- c(names.variable)
		
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
	
	colnames(MF.max) <- c(names.variable)
	rule.matrix <- MF.max
		
	## Step III
	## determine the degree of the rule
	degree.rule <- matrix(nrow = nrow(rule.matrix), ncol =1)
	degree.ante <- matrix(nrow = nrow(rule.matrix), ncol =1)
	
	for (i in 1:nrow(rule.matrix)){
		prod.d <- 1
		for (j in 1:ncol(rule.matrix)){
			if (rule.matrix[i, j] != 0){
			prod.d <- prod.d * rule.matrix[i, j]
			}
		}
		degree.rule[i] <- prod.d
	}

	# calculate degree of antecedent on each rule
	for (i in 1:nrow(rule.matrix)){
		prod.a <- 1
		for (j in 1 : (ncol(rule.matrix) - num.labels[1, ncol(num.labels)])){
			if (rule.matrix[i, j] != 0){
			prod.a <- prod.a * rule.matrix[i, j]
			}
		}
		degree.ante[i] <- prod.a
	}
		
	## convert element on rule.matrix into binary for detecting the same rows with each others.
	rule.matrix.bool <- rule.matrix
	check.cum <- matrix(nrow=nrow(rule.matrix))
	for (i in 1 : nrow(rule.matrix)){
		for (j in 1 : ncol(rule.matrix)){
			if (rule.matrix[i, j] != 0)
				rule.matrix.bool[i, j] = 1
			else
				rule.matrix.bool[i, j] = 0
		}		
		check.cum[i] <- sum(rule.matrix.bool[i, ])		
	}
		

	## find the same rows on matrix rule
	temp <- matrix(nrow = 1, ncol = ncol(rule.matrix.bool))
	kkk <- 1
	rule.red.mat <- rule.matrix
	degree.rule.temp <- degree.rule
	
	for (i in 1 : (nrow(rule.matrix.bool) - 1)){
		temp[1, ] <- rule.matrix.bool[i, ]		
		for (j in (i + 1) : nrow(rule.matrix.bool)) {
			chk <- which(temp[1, ] == rule.matrix.bool[j, ])
			if (length(chk) == length(temp)){
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
				kkk <- kkk + 1
			}
	
			if (j >= nrow(rule.matrix.bool)){
			break
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
	

	remove.row <- matrix()
	
	## check is.na on rule.data.num and degree.ante
	for (i in 1 : nrow(rule.data.num)){
		if(any(is.na(rule.data.num[i, ])== TRUE)){
			remove.row[i] <- i
		}	
	}
	
	for(i in 1:length(remove.row)){
		degree.ante[remove.row[i], ] <- NA 
		degree.rule[remove.row[i], ] <- NA
	}

	rule.data.num <- na.omit(rule.data.num)
	degree.ante <- na.omit(degree.ante)
	degree.rule <- na.omit(degree.rule)

	
	## build the rule into list of string
	rule <- matrix(nrow = nrow(rule.data.num), ncol = 2 * ncol(rule.data.num) - 1)
	
	for (i in 1 : nrow(rule.data.num)){
		k <- 0
		for (j in 1 : ncol(rule.data.num)){
			k <- k + 1
			if (j == ncol(rule.data.num) - 1){
				rule[i, k] <- c(names.variable[rule.data.num[i, j]])
				rule[i, k + 1] <- c("->")
				k <- k + 1
			}
			else if (j == ncol(rule.data.num)){
				rule[i, k] <- c(names.variable[rule.data.num[i, j]])
			}
			else{
				rule[i, k] <- c(names.variable[rule.data.num[i, j]])
				rule[i, k + 1] <- c("and")
				k <- k + 1
			}			
		}
	
	}
	
	## check a duplication on rule
	temp.rule <- matrix(nrow = 1, ncol = ncol(rule))
	for (i in 1 : (nrow(rule) - 1)){
		temp.rule[1, ] <- rule[i, ]
		
		for (j in (i + 1) : nrow(rule)) {
			chk <- which(temp.rule[1, ] == rule[j, ])
			if (length(chk) == length(temp.rule)){
				rule[i, ] <- NA
				degree.ante[i] <- NA
				rule.data.num[i, ] <- NA
				degree.rule[i] <- NA
			}
		}
	}
	
	rule <- na.omit(rule)
	degree.ante <- na.omit(degree.ante)
	rule.data.num <- na.omit(rule.data.num)
	degree.rule <- na.omit(degree.rule)

	
	#############################################################
	### Collect the Input data for "frbs for testing"
	#############################################################
	#get range.data
	range.input <- range.data[, 1 : ncol(range.data) - 1]
	
	##get range.output
	range.output.temp <- t(range.data[, ncol(range.data)])
	range.output <- t(range.output.temp)
	
	##  Change number of variable as input variable
	num.varinput <- num.varinput - 1
	
	## get num.fvalinput
	jj <- ncol(num.labels) - 1
	num.fvalinput <- matrix(num.labels[1, 1 : jj], nrow = 1)
	
	## indx is index of output variable on num.labels
	indx <- length(num.labels)
	length.namesvar <- length(names.variable)
		
	## convert names.variable into names.varinput
	names.varinput <- names.variable[1 : (length.namesvar - num.labels[1, indx])]

	## get number of output variable
	num.varoutput <- 1
	
	## cut membership function of output variable on var.mf
	ll <- (ncol(var.mf) - num.labels[1, indx])
	varinp.mf <- matrix(var.mf[, 1 : ll], nrow = 5)
	
	## get varout.mf
	varout.mf <- matrix(var.mf[, (ll + 1) : ncol(var.mf)], nrow = 5)
	
	## get names.varoutput
	mm <- (length.namesvar - num.labels[1, indx]) + 1
	names.varoutput <- names.variable[mm : length.namesvar]

	## clean rule
	rule <- na.omit(rule)
	
	mod <- list(range.input = range.input, range.output = range.output, num.varinput = num.varinput, num.fvalinput = num.fvalinput, names.varinput = names.varinput, 
	           varout.mf = varout.mf, names.varoutput = names.varoutput, rule = rule, varinp.mf = varinp.mf, degree.ante = degree.ante, rule.data.num = rule.data.num, 
			   degree.rule = degree.rule)
			   
	return (mod)
}  




