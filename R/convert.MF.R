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
# This function is to convert population which its elements are integer of membership function parameters
# into the matrix of membership function parameters.
#
# @title the converting function
# @param MF.param.popu the population on GA representing parameters of the membership function.
# @param varinp.mf a matrix of parameters of membership function on the input variables.
# @param varout.mf a matrix of parameters of membership function on the output variable.
# @param num.label number of label used as the fuzzy terms.
# @param num.varinput a number of input variables.
# @return a matrix representing a population of parameters of membership function.
# @export
convert.MF <- function(MF.param.popu, varinp.mf, varout.mf, num.label, num.varinput) {
	k <- 1
	j <- 1
	for (i in seq(1, num.varinput * num.label * 2, by = 2)){
		varinp.mf[2, k] <- MF.param.popu[1, i + 1]
		varinp.mf[3, k] <- MF.param.popu[1, i]
		k <- k + 1
	}
	
	for (i in seq(num.varinput * num.label * 2 + 1, ncol(MF.param.popu), by = 2)){
		varout.mf[2, j] <- MF.param.popu[1, i + 1]
		varout.mf[3, j] <- MF.param.popu[1, i]
		j <- j + 1
	}
	
	var.mf <- list(varinp.mf = varinp.mf, varout.mf = varout.mf)

	return(var.mf)
}