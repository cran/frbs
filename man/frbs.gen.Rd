\name{frbs.gen}
\alias{frbs.gen}
\title{The frbs model generator}
\usage{
  frbs.gen(range.data, num.fvalinput, names.varinput,
    num.fvaloutput = NULL, varout.mf = NULL,
    names.varoutput = NULL, rule, varinp.mf,
    type.model = "MAMDANI", type.defuz = "WAM",
    type.tnorm = "MIN", type.snorm = "MAX",
    func.tsk = NULL, colnames.var = NULL,
    type.implication.func = "ZADEH", name = "Sim-0")
}
\arguments{
  \item{range.data}{a matrix (\eqn{2 \times n}) containing
  the range of the data, where \eqn{n} is the number of
  variables, and first and second rows are the minimum and
  maximum values, respectively.}

  \item{num.fvalinput}{a matrix representing the number of
  linguistic terms of each input variables.

  For example: \code{num.fvalinput <- matrix(c(3,2), nrow =
  1)}

  means that there are two variables where the first
  variable has three linguistic terms and the second one
  has two linguistic terms.}

  \item{varinp.mf}{a matrix for constructing the shapes of
  the membership functions. See how to construct it in
  \code{\link{fuzzifier}}.}

  \item{names.varinput}{a list giving names to the
  linguistic terms for input variables. See
  \code{\link{rulebase}}.}

  \item{num.fvaloutput}{the number of linguistic terms of
  the output variable. This parameter is required for
  Mamdani model only.

  For example: \code{num.fvaloutput <- matrix(3, nrow = 1)}

  means there are 3 linguistic terms for the output
  variable.}

  \item{varout.mf}{a matrix for constructing the membership
  functions of the output variable. The form is the same as
  for the \code{varinp.mf} parameter. This parameter is
  required for Mamdani model only. See
  \code{\link{fuzzifier}}.}

  \item{names.varoutput}{a list giving names of the
  linguistic terms for the output variable. The form is the
  same as for the \code{names.varinput} parameter. This
  parameter is required for Mamdani model only. See
  \code{\link{rulebase}}.}

  \item{rule}{a list of fuzzy IF-THEN rules. There are some
  types of rule structures, for example: Mamdani, Takagi
  Sugeno Kang, and fuzzy rule-based classification systems
  (FRBCS). It allows to involve linguistic hedge, negation
  operator, etc. For more detail, see
  \code{\link{rulebase}}.}

  \item{type.model}{the type of the model. There are three
  types as follows. \itemize{ \item \code{MAMDANI} means we
  are using Mamdani model. \item \code{TSK} means we are
  using Takagi Sugeno Kang model. \item \code{FRBCS} means
  we are using fuzzy rule-based classification systems
  (FRBCS). }}

  \item{type.defuz}{the type of the defuzzification method.
  It is used in Mamdani model only.  See
  \code{\link{defuzzifier}}.}

  \item{type.tnorm}{the type of the t-norm method. See
  \code{\link{inference}}.}

  \item{type.snorm}{the type of the s-norm method. See
  \code{\link{inference}}.}

  \item{func.tsk}{a matrix of parameters of the function on
  the consequent part using the Takagi Sugeno Kang model.
  This parameter must be defined when we are using Takagi
  Sugeno Kang. See \code{\link{rulebase}}.}

  \item{colnames.var}{a list of names of input and output
  variables.}

  \item{type.implication.func}{a type of implication
  function. See \code{\link{WM}}.}

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
  step by step manually. Additionally, the examples show
  several scenarios as follows. \itemize{ \item Using
  \code{frbs.gen} for constructing Mamdani model on a
  regression task. \item Using \code{frbs.gen} for
  constructing Takagi Sugeno Kang model on a regression
  task. \item Constructing Mamdani model by executing
  internal functions such as \code{rulebase},
  \code{fuzzifier}, \code{inference}, and
  \code{defuzzifier} for Mamdani model. \item Using
  \code{frbs.gen} for constructing fuzzy rule-based
  classification systems (FRBCS) model. }
}
\examples{
#################################################
## 1. The following codes show how to generate a fuzzy model
## using the frbs.gen function for regression tasks.
## The following are three scenarios:
## 1a. Using Mamdani model
## 1b. Using Takagi Sugeno Kang model
## 1c. Using Mamdani model and internal functions: fuzzifier, etc.
## Note:
## In the examples, let us consider four input variabels and one output variable.
## Some variables could be shared together for other examples.
#################################################

## Define shape and parameters of membership functions of input variables.
## Please see fuzzifier function to construct the matrix.
varinp.mf <- matrix(c(2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA,
                      2, 0, 35, 75, NA, 3, 35, 75, 100, NA,
                      2, 0, 20, 40, NA, 1, 20, 50, 80, NA, 3, 60, 80, 100, NA,
                      2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA),
                      nrow = 5, byrow = FALSE)

## Define number of linguistic terms of input variables.
## Suppose, we have 3, 2, 3, and 3 numbers of linguistic terms
## for first, second, third and fourth variables, respectively.
num.fvalinput <- matrix(c(3, 2, 3, 3), nrow=1)

## Give the names of the linguistic terms of each input variables.
## It should be noted that the names of the linguistic terms must be unique,
## so we put a number for making it unique.
varinput.1 <- c("a1", "a2", "a3")
varinput.2 <- c("b1", "b2")
varinput.3 <- c("c1", "c2", "c3")
varinput.4 <- c("d1", "d2", "d3")
names.varinput <- c(varinput.1, varinput.2, varinput.3, varinput.4)

## Set interval of data.
range.data <- matrix(c(0,100, 0, 100, 0, 100, 0, 100, 0, 100), nrow=2)

## Define inference parameters.
type.defuz <- "WAM"
type.tnorm <- "MIN"
type.snorm <- "MAX"
type.implication.func <- "ZADEH"

## Give the name of simulation.
name <- "Sim-0"

## Provide new data for testing.
newdata<- matrix(c(25, 40, 35, 15, 45, 75, 78, 70), nrow= 2, byrow = TRUE)
## the names of variables
colnames.var <- c("input1", "input2", "input3", "input4", "output1")

###################################################################
## 1a. Using Mamdani Model
####################################################################
## Define number of linguistic terms of output variable.
## In this case, we set the number of linguistic terms to 3.
num.fvaloutput <- matrix(c(3), nrow=1)

## Give the names of the linguistic terms of the output variable.
## Note: the names of the linguistic terms must be unique.
varoutput.1 <- c("e1", "e2", "e3")
names.varoutput <- c(varoutput.1)

## Define the shapes and parameters of the membership functions of the output variables.
varout.mf <- matrix(c(2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA),
                      nrow = 5, byrow = FALSE)

## Set type of model which is "MAMDANI" or "TSK" for Mamdani or
## Takagi Sugeno Kang model, respectively.
## In this case, we choose Mamdani model.
type.model <- "MAMDANI"

## Define the fuzzy IF-THEN rules;
## there are two kinds of model: Mamdani and Takagi Sugeno Kang model
## if we use the Mamdani model then the consequent part is a linguistic term,
## but if we use Takagi Sugeno Kang then we build a matrix representing
## linear equations in the consequent part.
## In this example we are using the Mamdani model
## (see the type.model parameter).
## Note:
## "a1", "and", "b1, "->", "e1" means that
## "IF inputvar.1 is a1 and inputvar.2 is b1 THEN outputvar.1 is e1"
## Make sure that each rule has a "->" sign.
rule <- matrix(c("a1","and","b1","and","c1","and","d1","->","e1",
                 "a2","and","b2","and","c2","and","d2", "->", "e2",
                 "a3","and","b2","and","c2","and","d1", "->", "e3"),
                 nrow=3, byrow=TRUE)

## Generate a fuzzy model with frbs.gen.
object <- frbs.gen(range.data, num.fvalinput, names.varinput,
                 num.fvaloutput, varout.mf, names.varoutput, rule,
                 varinp.mf, type.model, type.defuz, type.tnorm,
                 type.snorm, func.tsk = NULL, colnames.var, type.implication.func, name)

## Plot the membership function.
plotMF(object)

## Predicting using new data.
res <- predict(object, newdata)$predicted.val

#####################################################################
## 1b. Using Takagi Sugeno Kang (TSK) Model
#####################################################################
## Define "TSK" for Takagi Sugeno Kang model
type.model <- "TSK"

## Define linear equation for consequent parts.
func.tsk <- matrix(c(1, 1, 5, 2, 1, 3, 1, 0.5, 0.1, 2, 1, 3, 2, 2, 2), nrow=3, byrow=TRUE)

## Define the fuzzy IF-THEN rules;
## For TSK model, it isn't necessary to put linguistic term in consequent parts.
## Make sure that each rule has a "->" sign.
rule <- matrix(c("a1","and","b1","and","c1","and","d1","->",
                 "a2","and","b2","and","c2","and","d2", "->",
                 "a3","and","b2","and","c2","and","d1", "->"),
                 nrow=3, byrow=TRUE)

## Generate a fuzzy model with frbs.gen.
## It should be noted that for TSK model, we do not need to input:
## num.fvaloutput, varout.mf, names.varoutput, type.defuz.
object <- frbs.gen(range.data, num.fvalinput, names.varinput,
             num.fvaloutput = NULL, varout.mf = NULL, names.varoutput = NULL, rule,
				varinp.mf, type.model, type.defuz = NULL, type.tnorm, type.snorm,
             func.tsk, colnames.var, type.implication.func, name)

## Plot the membership function.
plotMF(object)

## Predicting using new data.
res <- predict(object, newdata)$predicted.val

######################
## 1c. Using the same data as in the previous example, this example performs
## step by step of the generation of a fuzzy rule-based system
######################
## Using Mamdani model.
type.model <- "MAMDANI"

## Construct rules.
rule <- matrix(c("a1","and","b1","and","c1","and","d1","->","e1",
                 "a2","and","b2","and","c2","and","d2", "->", "e2",
                 "a3","and","b2","and","c2","and","d1", "->", "e3"),
                 nrow=3, byrow=TRUE)

## Check input data given by user.
rule <- rulebase(type.model, rule, func.tsk = NULL)

## Fuzzification Module:
## In this function, we convert crisp into linguistic values/terms
## based on the data and the parameters of the membership function.
## The output: a matrix representing the degree of the membership of the data
num.varinput <- ncol(num.fvalinput)
MF <- fuzzifier(newdata, num.varinput, num.fvalinput, varinp.mf)

## Inference Module:
## In this function, we will calculate the confidence factor on the antecedent for each rule
## considering t-norm and s-norm.
miu.rule <- inference(MF, rule, names.varinput, type.tnorm, type.snorm)

## Defuzzification Module.
## In this function, we calculate and convert the linguistic values back into crisp values.
range.output <- range.data[, ncol(range.data), drop = FALSE]
result <- defuzzifier(newdata, rule, range.output, names.varoutput,
                  varout.mf, miu.rule, type.defuz, type.model, func.tsk = NULL)


#################################################
## 2. The following codes show how to generate a fuzzy model
## using the frbs.gen function for classification tasks using Mamdani model.
#################################################
## define range of data.
## Note. we only define range of input data.
range.data.input <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1), nrow=2)

## Define shape and parameters of membership functions of input variables.
## Please see fuzzifier function to construct the matrix.
## In this case, we are using TRIANGLE for membership functions.
varinp.mf <- matrix(c(1, 0, 0, 0.5, NA, 1, 0, 0.5, 1, NA, 1, 0.5, 1, 1, NA,
                      1, 0, 0, 0.5, NA, 1, 0, 0.5, 1, NA, 1, 0.5, 1, 1, NA,
                      1, 0, 0, 0.5, NA, 1, 0, 0.5, 1, NA, 1, 0.5, 1, 1, NA,
                      1, 0, 0, 0.5, NA, 1, 0, 0.5, 1, NA, 1, 0.5, 1, 1, NA),
                      nrow = 5, byrow = FALSE)

## Define number of linguistic terms of input variables.
## Suppose, we have 3, 3, 3, and 3 numbers of linguistic terms
## for first up to fourth variables, respectively.
num.fvalinput <- matrix(c(3, 3, 3, 3), nrow=1)

## Give the names of the linguistic terms of each input variable.
## It should be noted that the names of the linguistic terms must be unique,
## so we put a number for making it unique.
varinput.1 <- c("v.1_a.1", "v.1_a.2", "v.1_a.3")
varinput.2 <- c("v.2_a.1", "v.2_a.2", "v.2_a.3")
varinput.3 <- c("v.3_a.1", "v.3_a.2", "v.3_a.3")
varinput.4 <- c("v.4_a.1", "v.4_a.2", "v.4_a.3")
names.varinput <- c(varinput.1, varinput.2, varinput.3, varinput.4)

## Provide inference parameters.
type.tnorm <- "MIN"
type.snorm <- "MAX"
type.implication.func <- "ZADEH"
type.model <- "FRBCS"

## Give the name of simulation.
name <- "Sim-0"

## Provide new data for testing.
newdata<- matrix(c(0.45, 0.5, 0.89, 0.44, 0.51, 0.99, 0.1, 0.98, 0.51,
                 0.56, 0.55, 0.5), nrow= 3, byrow = TRUE)

## the names of variables
colnames.var <- c("input1", "input2", "input3", "input4", "output1")

## Construct rules.
## Take into account in consequent part, which expresses classes.
rule <- matrix(c("v.1_a.2","and","v.2_a.2","and","v.3_a.3","and","v.4_a.2","->","3",
         "v.1_a.2","and","v.2_a.3","and","v.3_a.1","and","v.4_a.3","->","1",
         "v.1_a.2","and","v.2_a.2","and","v.3_a.2","and","v.4_a.2","->","2"),
         nrow=3, byrow=TRUE)

## Generate frbs object.
object <- frbs.gen(range.data = range.data.input, num.fvalinput,
             names.varinput, num.fvaloutput = NULL, varout.mf = NULL,
             names.varoutput = NULL, rule, varinp.mf, type.model,
             type.defuz = NULL, type.tnorm, type.snorm, func.tsk = NULL,
             colnames.var, type.implication.func, name)

## Plot the shape of membership functions.
plotMF(object)

## Predicting using new data.
res <- predict(object, newdata)
}

