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
# This function is to convert population which its elements are integer into rule database.  
#
# @title the rule converting function
# @param popu the population on GA representing the fuzzy IF-THEN rules in integer data type.
# @param n.inputvar the number of input variables.
# @param num.label the number of label used as the fuzzy terms.
# @param names.variable the names of the fuzzy terms of variabels.
# @return a rule database
# @export
convert.rule <- function(popu, n.inputvar, num.label, names.variable) {

n.outputvar <- 1 
num.var <- n.inputvar + n.outputvar
num.rule <- ncol(popu) / num.var

rule.temp <- popu[1, 1 : num.var]
for (i in 2 : num.rule){
	start.rule <- (i - 1) * num.var + 1
	end.rule <- start.rule + num.var - 1
	temp <- popu[1, start.rule : end.rule]
	
	rule.temp <- rbind(rule.temp, temp)	
}

popu <- rule.temp

rule <- matrix(nrow = nrow(popu), ncol = 2 * ncol(popu) - 1)

for (i in 1 : nrow(popu)){
		k <- 0
		for (j in 1 : ncol(popu)){
			k <- k + 1
			if (j == ncol(popu) - 1){
				rule[i, k] <- c(names.variable[popu[i, j]])
				rule[i, k + 1] <- c("->")
				k <- k + 1
			}
			else if (j == ncol(popu)){
				rule[i, k] <- c(names.variable[popu[i, j]])
			}
			else{
				rule[i, k] <- c(names.variable[popu[i, j]])
				rule[i, k + 1] <- c("and")
				k <- k + 1
			}			
		}
	
	}

return(rule)
}
