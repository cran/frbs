\name{plotMF}
\alias{plotMF}
\title{The plotting function}
\usage{
  plotMF(object)
}
\arguments{
  \item{object}{an \code{\link{frbs-object}} or a list of
  parameters to plot membership functions when we build the
  frbs model without learning.  For plotting using the
  list, there are several parameters that must be inserted
  in params as follows.  \itemize{ \item \code{var.mf}: a
  matrix of membership function of input and output
  variables. Please see \code{\link{fuzzifier}}.  \item
  \code{range.data.ori}: a matrix (\eqn{2 \times n})
  containing the range of the data, where \eqn{n} is the
  number of variables, and first and second rows are the
  minimum and maximum values, respectively.  \item
  \code{num.labels}: the number of linguistic terms of the
  input and output variables.

  For example: \code{num.labels <- matrix(c(3, 3, 3), nrow
  = 1)}

  It means we have 3 linguistic values/labels for two input
  variables and one output variable.  \item
  \code{names.variables}: a list of names of variables.

  For example: \code{names.variables <- c("input1",
  "input2", "output1")} }}
}
\description{
  This function can be used to plot the shapes of the
  membership functions.
}
\examples{
## The following examples contain two different cases which are
## using an frbs-object and the manual way.
##
## 1. Plotting using frbs object.
data(iris)
irisShuffled <- iris[sample(nrow(iris)),]
irisShuffled[,5] <- unclass(irisShuffled[,5])
tra.iris <- irisShuffled[1:105,]
tst.iris <- irisShuffled[106:nrow(irisShuffled),1:4]
real.iris <- matrix(irisShuffled[106:nrow(irisShuffled),5], ncol = 1)

## Please take into account that the interval needed is the range of input data only.
range.data.input <- matrix(c(4.3, 7.9, 2.0, 4.4, 1.0, 6.9, 0.1, 2.5), nrow=2)

## generate the model
method.type <- "FRBCS.W"
control <- list(num.labels = 7, type.mf = 1)
\dontrun{object <- frbs.learn(tra.iris, range.data.input, method.type, control)}

## plot the frbs object
\dontrun{plotMF(object)}

## 2. Plotting using params.
## Define shape and parameters of membership functions of input variables.
## Please see the fuzzifier function of how to contruct the matrix.
varinp.mf <- matrix(c(2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA,
                      2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA,
                      2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA,
                      2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA),
                      nrow = 5, byrow = FALSE)

## Define the shapes and parameters of the membership functions of the output variables.
varout.mf <- matrix(c(2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA),
                      nrow = 5, byrow = FALSE)
var.mf <- cbind(varinp.mf, varout.mf)
range.data <- matrix(c(0,100, 0, 100, 0, 100, 0, 100, 0, 100), nrow=2)
num.labels <- matrix(c(3,3,3,3,3), nrow = 1)
names.variables <- c("input1", "input2", "input3", "input4", "output1")

## plot the membership function.
\dontrun{plotMF(object = list(var.mf = var.mf, range.data.ori = range.data,
          num.labels = num.labels, names.variables = names.variables))}
}

