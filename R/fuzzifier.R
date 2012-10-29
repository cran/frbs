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
#' The dimension of the matrix is (5, n) where n is the number of fuzzy terms/labels and the number of variables.
#' The rows of the matrix represent:
#' The first row is the type of membership function, where 1 means triangular, 2 means trapezoid 1a (left side), 
#' 3 means trapezoid 1b (right side), 4 means trapezoid 2 (in the middle), 5 means Gaussian, 6 means sigmoid, 7 means 
#' generalized bell. The second until fifth row indicate the critical points to construct the functions. 
#' \itemize{
#' \item triangular has three parameters (a, b, c), where b is the center point of the triangular, and a and c are the left and right points, respectively.
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
#' @return A matrix of the degree of each fuzzy term based on the shape of the membership functions 
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

