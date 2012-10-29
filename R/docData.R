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
#' The package includes embedded versions of the Mackey-Glass chaotic time series and the Gas Furnance dataset.
#'
#' \bold{Mackey-Glass chaotic time series}
#' 
#' The Mackey-Glass chaotic time series is defined by the following delayed differential equation:
#' 
#' dx(t) / dt = (a * x(t - tau) / (1 + x(t - tau) ^ 10)) - b * x(t)
#' 
#' For this dataset, we generated 1000 samples, with input parameters as follows:
#' \itemize{
#' \item a = 0.2
#' \item b = 0.1
#' \item tau = 17
#' \item x0 = 1.2
#' \item dt = 1
#' }
#' 
#' The dataset is embedded in the following way: 
#'
#' input variables: x(t - 18), x(t - 12), x(t - 6), x(t)
#'
#' output variable: x(t + 6)
#'
#' \bold{Gas Furnance dataset}
#' 
#' The Gas Furnance dataset is taken from Box and Jenkins. It consists of 292 consecutive 
#' values of methane at time (t - 4), and the CO2 produced in a furnance at time (t - 1) as input 
#' variables, with the produced CO2 at time (t) as an output variable. So, each training data 
#' point consists of [u(t - 4), y(t - 1), y(t)], where u is methane and y is CO2.
#'
#' @title Data set of the package
#' @name frbsData
#' @docType data
#' @references 
#' A. Lapedes and R. Farber, "Nonlinear signal processing using neural networks: Prediction and system modelling, " LA-UR-87-2662 (1987).
#'
#' Breiman, L., "Bagging Predictors," Machine Learning, 24(3), 123 - 140, Kluwer Academic Publishres (1996).
#'
#' Box, G. E. P., & Jenkins, G. M. "Time Series Analysis, forecasting and control, San Fransisco, CA: Holden Day (1970).
#' 
#' Mackey, M., & Glass, L., "Oscillation and chaos in physiological control systems, " Science, vol. 197, pp. 287 - 289 (1977).
#' 
#' @keywords data
NULL