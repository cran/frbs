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
#' This function enables the output of a summary of the \code{\link{frbs-object}}. 
#'
#' @title The summary function for frbs objects
#' 
#' @param object the \code{\link{frbs-object}}
#' @param ... the other parameters (not used)
#' @export  
#' @method summary frbs
#' @S3method summary frbs
summary.frbs <- function(object, ...){
	
  if(!inherits(object, "frbs")) stop("not a legitimate frbs model")
  
	cat("Model name: ", object$name, "\n")
	cat("Model was trained with: ", object$method.type, "\n") 
	
	if (object$method.type == "WM"){
		range.data.ori <- cbind(object$range.input, object$range.output)
		print("Range of data (original): ")
		print(range.data.ori)
	}
	else {
		print("Range of data (original): ")
		print(object$range.data.ori)
	}

	if (any(object$method.type == c("GFS", "HYFIS", "ANFIS", "DM", "FRCS", "HGD"))){
		print("Range of the input variables (normalized): ")
		print(object$range.input)
		print("Range of the output variable (normalized): ")
		print(object$range.output)
	}
		
	if (any(object$method.type == c("DENFIS", "SBC"))){
		print("Cluster center: ")
		print(object$cls)
	}

	else if (any(object$method.type == c("WM", "GFS", "HYFIS", "ANFIS", "DM", "FRCS", "HGD"))){
		print("Fuzzy terms on input variables: ")
		print(object$names.varinput)
		print("The fuzzy IF-THEN rules: ")
		print(object$rule)
		print("Parameters of membership function on input variables: ")
		print(object$varinp.mf)
	}
		
	if (any(object$method.type == c("WM", "HYFIS", "GFS"))) {
		print("Fuzzy terms on output variables: ")
		print(object$names.varoutput)
		print("Parameters of membership function on the output variable: ")
		print(object$varout.mf)
	}

	if (any(object$method.type == c("ANFIS", "DM", "HGD", "FRCS"))){
		print("The linear equations on consequent parts of fuzzy IF-THEN rules: ")
		print(object$func.tsk)
	}
  
	if (object$method.type == "MSGFS"){
		print("The parameter values of membership functions representing the fuzzy IF-THEN rules: ")
		print(object$rule)
	}
	
  invisible(object)	
}