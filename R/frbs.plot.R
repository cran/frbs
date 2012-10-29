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
#' This function can be used to plot the shapes of the membership functions.
#'
#' @title The plotting function
#' 
#' @param object an \code{\link{frbs-object}}.
# @S3method plot frbs
#' @export
plotMF <- function(object){
  
  res <- object
  
  if(!inherits(object, "frbs")) stop("not a legitimate frbs model")

  if (any(object$method.type == c("WM", "HYFIS", "GFS", "ANFIS", "HGD", "DM", "FRCS"))){
	  ### argument: range.data, num.varinput, num.fvalinput, var.mf, names.varinput
	  
	  if (any(object$method.type == c("WM", "HYFIS", "GFS", "FRCS"))){
		  range.input <- object$range.input
		  range.output <- object$range.output
		  num.varinput <- object$num.varinput
		  num.fvalinput <- object$num.fvalinput
		  varinp.mf <- object$varinp.mf
		  varout.mf <- object$varout.mf
		  names.varinput <- object$names.varinput
		  names.varoutput <- object$names.varoutput
		   
		  range.data <- cbind(range.input, range.output) 
		  num.varinput <- num.varinput + 1
		  num.fvalinput <- cbind(num.fvalinput, num.fvalinput[1]) 
		  var.mf <- cbind(varinp.mf, varout.mf) 
		  names.varinput <- c(names.varinput, names.varoutput)
	  }
	 
	 else if (any(object$method.type == c("ANFIS", "HGD", "DM"))){
		### argument: range.data, num.varinput, num.fvalinput, var.mf, names.varinput
		  range.input <- object$range.input
		  range.output <- object$range.output
		  num.fvalinput <- object$num.fvalinput
		  names.varinput <- object$names.varinput
		  
		  range.data <- object$range.input
		  num.varinput <- object$num.varinput
		  var.mf <- object$varinp.mf
	  }
	  ##get number of column of var.mf
	  col.var.mf <- ncol(var.mf)
	  
	  ## counter is used to make continue index j
	  counter <- 1
	  
	  ## k is used to index number fuzzy value on each variable on varinput.
	  k <- 1
	  
	  ## make row plot 
	  op <- par(mfrow = c(num.varinput, 1), mar = c(3, 2, 2, 1) + 0.1)
	  
	  ## set the names of var.mf
	  colnames(var.mf) <- (names.varinput)
	  
	  ## loop as many as number of input variable
	  for (i in 1 : num.varinput){
		j <- counter
		
		## Initialize plot
		MF <- function(x){
		  y <- x - x
		  return (y)
		}
		names <- paste("The plot of Membership Function of Variable", i, sep= " ")
		curve(MF, range.data[1, i], range.data[2, i], ylim = c(0, 1), col = "white", ylab = "degree", main = names)
		
		## loop as many as number of column of var.mf
		for (j in counter : col.var.mf){
		  
		  ## make boundary point
		  oo <- range.data[1, i]
		  aa <- var.mf[2, j]
		  bb <- var.mf[3, j]
		  cc <- var.mf[4, j]
		  dd <- var.mf[5, j]
		  mm <- range.data[2, i]
		  
		  ##condition for triangular type
		  if (var.mf[1, j] == 1){				
			
			## make a function for plotting, args (x) is sequence data
			f.y1 <- function(x){
			  
			  ## range of input data
			  p.0 <- x[x >= oo & x <= aa]
			  p.1 <- x[x > aa & x <= bb]
			  p.2 <- x[x > bb & x <= cc]
			  p.3 <- x[x > cc & x <= mm]
			  
			  ## build functions 
			  y0 <- (p.0 - p.0)
			  y1 <- (p.1 - aa) / (bb - aa)			
			  y2 <- (p.2 - cc) / (bb - cc)			
			  y4 <- (p.3 - p.3)
			  y <- c(y0, y1, y2, y4)
			  
			  return (y)
			}
			
			curve(f.y1, oo, mm, add = TRUE, col = "violet")
		  }
		  ## condition for trapezoid in left side
		  else if (var.mf[1, j] == 2){
			
			f.y2 <- function(x){
			  
			  p.1 <- x[x <= bb]
			  p.2 <- x[x > bb & x <= cc]
			  p.3 <- x[x > cc & x <= mm]
			  
			  y1 <- (p.1 - p.1 + 1)			
			  y2 <- (p.2 - cc) / (bb - cc)			
			  y3 <- (p.3 - p.3)
			  
			  y <- c(y1, y2, y3)
			  
			  return (y)
			}
			
			curve(f.y2, oo, mm, add=TRUE, col = "blue")
			
		  }
		  ## condition for trapezoid in right side
		  else if (var.mf[1, j] == 3){
			f.y3 <- function(x){
			  
			  p.0 <- x[x >= oo & x <= aa]
			  p.1 <- x[x > aa & x <= bb]
			  p.2 <- x[x > bb & x <= mm]
			  
			  y0 <- (p.0 - p.0)
			  y1 <- (p.1 - aa) / (bb - aa)
			  y2 <- (p.2 - p.2 + 1)			
			  y <- c(y0, y1, y2)
			  
			  return (y)
			}
			
			curve(f.y3, oo, mm, add = TRUE, col = "green")
			
		  }
		  ## condition for trapezoid in the middle
		  else if (var.mf[1, j] == 4){
			
			f.y4 <- function(x){
			  
			  p.0 <- x[x >= oo & x <= aa]
			  p.1 <- x[x > aa & x <= bb]
			  p.2 <- x[x > bb & x <= cc]
			  p.3 <- x[x > cc & x <= dd]
			  p.4 <- x[x > dd & x <= mm]
			  
			  y0 <- (p.0 - p.0)
			  y1 <- (p.1 - aa) / (bb - aa)
			  y2 <- (p.2 - p.2 + 1)			
			  y3 <- (p.3 - dd) / (cc - dd)
			  y4 <- (p.4 - p.4)
			  
			  y <- c(y0, y1, y2, y3, y4)
			  
			  return (y)
			}
			
			curve(f.y4, oo, mm, add = TRUE, col = "red")
			
		  }
		  ## condition for gaussian shape
		  else if (var.mf[1, j] == 5){
			
			f.y5 <- function(x){
			  y <- exp(- 0.5 * (x - aa) ^ 2 / bb ^ 2)
			  return (y)
			}
			## plot the functions
			curve(f.y5, oo, mm, add = TRUE, col = "gray")
		  }
		  ## condition for sigmoid/logistic
		  else if (var.mf[1, j] == 6){
			
			f.y6 <- function(x){
			  y <- 1 / (1 + exp(- aa * (x - bb)))
			  return (y)
			}
			#plot the functions
			curve(f.y6, oo, mm, add = TRUE, col = "black")
		  }
		  ## condition for generalized bell
		  else if (var.mf[1, j] == 7){
			
			f.y7 <- function(x){
			  y <- 1 / (1 + abs((x - cc)/aa) ^ (2 * bb))
			  return (y)
			}
			
			## plot the functions
			curve(f.y7, oo, mm, add = TRUE, col = "black")
		  }
		  
		  counter <- j + 1 
		  k <- k + 1
		  if (k > num.fvalinput[1, i]){
			k <- 1
			break
		  }
		}
	  }

	  par(op)
  }
 else {
	print("The plot is not supported by the used method, please have a look the documentation")
 }
}

