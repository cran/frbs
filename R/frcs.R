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
#' This is the internal function that implements the fuzzy rule-based classification 
#' system (frcs) model. Users do not need to call it directly,
#' but just use \code{\link{frbs.learn}} and \code{\link{predict}}. This method is
#' suitable only for classification problems.
#' 
#' This method is adopted from Hisao Ishibuchi and Tomoharu Nakashima's paper and 
#' Oscar Cordon, Maria Jose del Jesus and Francisco Herrera's paper. 
#' The difference of a usual FRBS and this method is on the consequent part. The 
#' consequent part of frcs is a class with a certainty grade (i.e. rule weight). 
#' This method uses the techniques of Wang and Mendel for determining the antecedent 
#' part. Whereas on the consequent part, the value is associated to a class of the training data. 
#' After getting the fuzzy IF-THEN rules, we calculate a certainty grade (CF) for each rule. 
#' The prediction phase of this model is conducted by the internal function \code{\link{frcs.eng}}.
#'
#' @title frcs model building 
#'
#' @param range.data a matrix(2 x n) containing the range of the normalized data, where n is the number of variables, and
#' first and second rows are the minimum and maximum values, respectively. 
#' @param data.train a matrix(m x n) of data for the training process, where m is the number of instances and 
#' n is the number of variables; the last column is the output variable.
#' @param label.inp a matrix(1 x n) whose elements represent the number of labels (fuzzy terms), where n is the number of variables.
#' @param num.class an integer number representing the number of labels (fuzzy terms).
#' @param type.mf the type of the shape of the membership functions.
#' @seealso \code{\link{frcs.eng}}, \code{\link{frbs.learn}}, and \code{\link{predict}}
# @return a list of the model data. Please have a look at \code{\link{frbs.learn}} for looking its complete components. 
#' @references
#' Hisao Ishibuchi and Takashi Yamamoto, "Rule weight specification in fuzzy rule-based classification systems," 
#' IEEE Transactions on Fuzzy Systems, Vol. 13, No. 4 (2005).
#' 
#' Hisao Ishibuchi and Tomoharu Nakashima, "Effect of rule weights in fuzzy rule-based classification systems", 
#' IEEE Transactions on Fuzzy Systems, Vol. 9, No. 4 (2001).
#'
#' Oscar Cordon, Maria Jose del Jesus and Francisco Herrera, "A proposal on reasoning methods in fuzzy rule-based classification systems",
#' International Journal of Approximate Reasoning 20, pp. 21 - 45 (1999).
#'
#' Tomoharu N., Yasuyuki Y., Hisao Isibuchi, "Learning fuzzy IF-THEN rules for pattern classification with weighted training patterns," 
#' Proceedings of Joint 4th Conference of the European Society for Fuzzy Logic and Technology and the 11th Rencontres Francophones 
#' sur la Logique Floue et ses Applications, (2005).
# @export
frcs <- function(range.data, data.train, label.inp, num.class, type.mf) {

	## create labels on each variables
	num.labels.inp <- matrix(rep(label.inp, (ncol(range.data) - 1)), nrow=1)
	num.labels.out <- matrix(rep(num.class, 1), nrow=1)	
	num.labels <- cbind(num.labels.inp, num.labels.out)
	
	## generate initial model using WM
	mod <- WM(range.data, data.train, num.labels, type.mf, classification = TRUE)
	
	##getting names.varoutput and rule
	names.varoutput <- as.matrix(mod$names.varoutput)
	rule <- mod$rule
	rule.constant <- mod$rule
	degree.ante <- mod$degree.ante
	 for (i in 1 : nrow(names.varoutput)){
		 for (j in 1 : nrow(rule.constant)){
			 if (rule[j, ncol(rule)] == names.varoutput[i])
				 rule.constant[j, ncol(rule.constant)] <- i
		 }
	 }
	
	func.tsk <- matrix(strtoi(rule.constant[, ncol(rule.constant)])) 
	grade.cert.temp <- cbind(func.tsk, degree.ante)
	
	class.fact <- factor(grade.cert.temp[, 1], exclude = "")
	beta.class <- as.real(as.matrix(aggregate(grade.cert.temp[, 2], by = list(class.fact), sum)))
	beta.class <- matrix(beta.class, nrow = num.class)	
	
	CF <- matrix()
	 
	 for(i in 1 : nrow(beta.class)){
		CF[i] <- 1 + (beta.class[i, 2] - (sum(beta.class[, 2]) - (beta.class[i, 2] / (nrow(beta.class) - 1))))/sum(beta.class[, 2])
	 }
	 
	grade.cert <- grade.cert.temp 
	for(i in 1 : nrow(grade.cert.temp)){
		grade.cert[i, 2] <- CF[grade.cert[i,1]] 
	}


	mod$func.tsk <- func.tsk
	mod$grade.cert <- grade.cert

	return (mod)
}  




