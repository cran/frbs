\name{frbs.gen}
\alias{frbs.gen}
\title{The frbs model generator}
\usage{
  frbs.gen(range.input, range.output, num.fvalinput,
    names.varinput, num.fvaloutput, varout.mf,
    names.varoutput, rule, varinp.mf, type.model = 1,
    type.defuz = 1, type.tnorm = 1, type.snorm = 1,
    func.tsk = NULL, colnames.var = NULL,
    method.type = "WM", name = "Sim-0")
}
\arguments{
  \item{range.input}{a matrix(2 x n) containing the range
  of the input data.}

  \item{range.output}{a matrix(2 x n) containing the range
  of the output data.}

  \item{num.fvalinput}{a matrix with the number of fuzzy
  terms of each input variable.

  For example: \code{num.fvalinput <- matrix(c(3,2), nrow =
  1)}

  means that there are two variables where the first
  variable has three fuzzy terms and the second one has two
  fuzzy terms.}

  \item{varinp.mf}{a matrix for constructing the shapes of
  the membership functions. See \code{\link{fuzzifier}}.}

  \item{names.varinput}{a list for giving names to the
  fuzzy terms. See \code{\link{rulebase}}.}

  \item{num.fvaloutput}{the number of fuzzy terms of the
  output variable.

  For example: \code{num.fvaloutput <- matrix(3, nrow = 1)}

  means there are 3 fuzzy terms for the first variable (in
  this case, there is only one variable).}

  \item{varout.mf}{a matrix for constructing the membership
  functions of the output variable. The form is the same as
  for the \code{varinp.mf} parameter. Please see
  \code{\link{fuzzifier}}.}

  \item{names.varoutput}{a list for giving names of the
  fuzzy terms. The form is the same as for the
  \code{names.varinput} parameter. Please see
  \code{\link{rulebase}}.}

  \item{rule}{a list of fuzzy IF-THEN rules. Please see
  \code{\link{rulebase}}.}

  \item{type.model}{the type of the model. Please see
  \code{\link{defuzzifier}}.}

  \item{type.defuz}{the type of the defuzzification method.
  Please see \code{\link{defuzzifier}}.}

  \item{type.tnorm}{the type of the t-norm method. Please
  see \code{\link{inference}}.}

  \item{type.snorm}{the type of the s-norm method. Please
  see \code{\link{inference}}.}

  \item{func.tsk}{a matrix of parameters of the function on
  the consequent part using the Takagi Sugeno Kang model.
  Please see \code{\link{rulebase}}.}

  \item{colnames.var}{a list of names of input and output
  variables.}

  \item{method.type}{the type of the selected method.
  Please see \code{\link{frbs.learn}}.}

  \item{name}{a name of the simulation.}
}
\value{
  The \code{\link{frbs-object}}.
}
\description{
  The purpose of this function is to generate a FRBS model
  from user-given input without a learning process.
}
\details{
  It can be used if rules have already been obtained
  manually, without employing the learning process. In the
  examples shown, we generate a fuzzy model using
  \code{frbs.gen} and generate the fuzzy rule-based systems
  step by step manually.
}
\examples{
## This example shows how to use frbs without
## learning process.

## Define shape and parameters of membership functions of input variables.
## Please see fuzzifier function to contruct the matrix.
varinp.mf <- matrix(c(2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA,
                      2, 0, 35, 75, NA, 3, 35, 75, 100, NA,
                      2, 0, 20, 40, NA, 1, 20, 50, 80, NA, 3, 60, 80, 100, NA,
                      2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA),
                      nrow = 5, byrow = FALSE)

## Define number of fuzzy terms of input variables.
## Suppose, we have 3, 2, 3, and 3 numbers of fuzzy terms
## for first, second, third and fourth variables, respectively.
num.fvalinput <- matrix(c(3, 2, 3, 3), nrow=1)

## Give the names of the fuzzy terms of each input variable.
## It should be noted that the names of the fuzzy terms must be unique,
## so we put a number for making it unique.
varinput.1 <- c("a1", "a2", "a3")
varinput.2 <- c("b1", "b2")
varinput.3 <- c("c1", "c2", "c3")
varinput.4 <- c("d1", "d2", "d3")
names.varinput <- c(varinput.1, varinput.2, varinput.3, varinput.4)

## Set interval of data.
range.input <- matrix(c(0,100, 0, 100, 0, 100, 0, 100), nrow=2)
range.output <- matrix(c(0,100), nrow=2)

## Define number of fuzzy terms of output variable.
## In this case, we set the number of fuzzy terms to 3.
num.fvaloutput <- matrix(c(3), nrow=1)

## Give the names of the fuzzy terms of the output variable.
## Note: the names of the fuzzy terms must be unique.
varoutput.1 <- c("e1", "e2", "e3")
names.varoutput <- c(varoutput.1)

## Define the shapes and parameters of the membership functions of the output variables.
varout.mf <- matrix(c(2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA),
                      nrow = 5, byrow = FALSE)

## Set type of model which is 1 or 2 for Mamdani or Takagi Sugeno Kang model, respectively.
## In this case, we choose Mamdani model.
type.model <- 1
## Set weighted average method to be used as defuzzification method.
type.defuz <- 1
## We are using standard t-norm and s-norm.
type.tnorm <- 1
type.snorm <- 1

## Since we don't generate the fuzzy model by learning from data,
## we have to set Wang and Mendel's technique as type of method.
method.type <- "WM"
## Give the name of simulation.
name <- "Sim-0"
## Define the fuzzy IF-THEN rules;
## there are two kinds of model: Mamdani and Takagi Sugeno Kang model
## if we use the Mamdani model then the consequent part is a linguistic term,
## but if we use Takagi Sugeno Kang then we build a matrix representing
## linear equations in the consequent part.
## In this example we are using the Mamdani model
## (see the type.model parameter).
## We make sure that each rule has a "->" sign.
rule <- matrix(c("a1","and","b1","and","c1","and","d1","->","e1",
                 "a2","and","b2","and","c2","and","d2", "->", "e2",
                 "a3","and","b2","and","c2","and","d1", "->", "e3"),
                 nrow=3, byrow=TRUE)

## Define function of TSK if we use it or
## set NULL if we use the Mamdani model.
func.tsk<-matrix(c(1, 1, 5, 2, 1, 3, 1, 0.5, 0.1, 2, 1, 3, 2, 2, 2), nrow=3, byrow=TRUE)
## Provide new data for testing.
newdata<- matrix(c(25, 40, 35, 15, 45, 75, 78, 70), nrow= 2, byrow = TRUE)
## the names of variables
colnames.var <- c("input1", "input2", "input3", "input4", "output1")
######################
## 1. The following codes show how to generate a fuzzy model using the frbs.gen function
######################
## Generate a fuzzy model with frbs.gen.
object <- frbs.gen(range.input, range.output, num.fvalinput, names.varinput,
                num.fvaloutput, varout.mf, names.varoutput, rule, varinp.mf,
                type.model, type.defuz, type.tnorm, type.snorm, func.tsk,
                colnames.var, method.type, name)

## We can plot the membership function
plotMF(object)

## Predicting using new data.
res <- predict(object, newdata)

######################
## 2. Using the same data as in the previous example, this example performs
## step by step of the generation of a fuzzy rule-based system
######################
## Check input data given by user.
rule <- rulebase(type.model, rule, func.tsk)

## Fuzzification Module:
## In this function, we convert crisp values into fuzzy values
## based on the data and the parameters of the membership function.
## The output: a matrix representing the degree of the membership of the data
num.varinput <- ncol(num.fvalinput)
MF <- fuzzifier(newdata, num.varinput, num.fvalinput, varinp.mf)

## Inference Module:
## In this function, we will calculate the confidence factor on the antecedent for each rule
## considering t-norm and s-norm.
miu.rule <- inference(MF, rule, names.varinput, type.tnorm, type.snorm)

## Defuzzification Module
## In this function, we calculate and convert the fuzzy values back into crisp values.
result <- defuzzifier(newdata, rule, range.output, names.varoutput,
                  varout.mf, miu.rule, type.defuz, type.model, func.tsk)
}

