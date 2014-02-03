\name{frbs.learn}
\alias{frbs.learn}
\title{The frbs model building function}
\usage{
  frbs.learn(data.train, range.data = NULL,
    method.type = c("WM"), control = list())
}
\arguments{
  \item{data.train}{a data frame or matrix (\eqn{m \times
  n}) of data for the training process, where \eqn{m} is
  the number of instances and \eqn{n} is the number of
  variables; the last column is the output variable. It
  should be noted that the training data must be expressed
  in numbers (numerical data). And, especially for
  classification tasks, the last column representing class
  names/symbols isn't allowed to have values 0 (zero). In
  the other words, the categorical values 0 should be
  replaced with other values.}

  \item{range.data}{a matrix (\eqn{2 \times n}) containing
  the range of the data, where \eqn{n} is the number of
  variables, and first and second rows are the minimum and
  maximum values, respectively. It should be noted that for
  \code{"FRBCS.W"}, \code{"FRBCS.CHI"}, \code{"GFS.GCCL"},
  \code{"FH.GBML"}, and \code{"SLAVE"}, \eqn{n} represents
  the number of input variables only (without the output
  variable). It will be assigned as min/max of training
  data if it is omitted.}

  \item{method.type}{this parameter determines the learning
  algorithm to be used.  The following methods are
  implemented: \itemize{ \item \code{"WM"}: Wang and
  Mendel's technique to handle regression tasks. See
  \code{\link{WM}}; \item \code{"SBC"}: subtractive
  clustering method to handle regression tasks. See
  \code{\link{SBC}}; \item \code{"HYFIS"}: hybrid neural
  fuzzy inference systems to handle regression tasks. See
  \code{\link{HyFIS}}; \item \code{"ANFIS"}: adaptive
  neuro-fuzzy inference systems to handle regression tasks.
  See \code{\link{ANFIS}}; \item \code{"FRBCS.W"}: fuzzy
  rule-based classification systems with weight factor
  based on Ishibuchi's method to handle classification
  tasks. See \code{\link{FRBCS.W}}; \item
  \code{"FRBCS.CHI"}: fuzzy rule-based classification
  systems based on Chi's method to handle classification
  tasks. See \code{\link{FRBCS.CHI}}; \item
  \code{"DENFIS"}: dynamic evolving neuro-fuzzy inference
  systems to handle regression tasks. See
  \code{\link{DENFIS}}; \item \code{"FS.HGD"}: fuzzy system
  using heuristic and gradient descent method to handle
  regression tasks. See \code{\link{FS.HGD}}; \item
  \code{"FIR.DM"}: fuzzy inference rules by descent method
  to handle regression tasks. See \code{\link{FIR.DM}};
  \item \code{"GFS.FR.MOGUL"}: genetic fuzzy systems for
  fuzzy rule learning based on the MOGUL methodology to
  handle regression tasks. See \code{\link{GFS.FR.MOGUL}};
  \item \code{"GFS.THRIFT"}: Thrift's technique based on
  genetic algorithms to handle regression tasks. See
  \code{\link{GFS.Thrift}}; \item \code{"GFS.GCCL"}:
  Ishibuchi's method based on genetic
  cooperative-competitive learning to handle classification
  tasks. See \code{\link{GFS.GCCL}}; \item
  \code{"FH.GBML"}: Ishibuchi's method based on
  hybridization of genetic cooperative-competitive learning
  and Pittsburgh to handle classification tasks. See
  \code{\link{FH.GBML}}; \item \code{"SLAVE"}: structural
  learning algorithm on vague environment to handle
  classification tasks. See \code{\link{SLAVE}}; \item
  \code{"GFS.LT.RS"}: genetic algorithm for lateral tuning
  and rule selection. See \code{\link{GFS.LT.RS}} }}

  \item{control}{a list containing all arguments, depending
  on the learning algorithm to use. The following list are
  parameters required for each methods, whereas their
  descriptions will be explained later on. \itemize{ \item
  \code{WM}:

  \code{list(num.labels, type.mf, type.tnorm, type.snorm,
  type.defuz,}

  \code{type.implication.func, name)}

  \item \code{HYFIS}:

  \code{list(num.labels, max.iter, step.size, type.tnorm,}

  \code{type.snorm, type.defuz, type.implication.func,
  name)}

  \item \code{ANFIS} and \code{FIR.DM}:

  \code{list(num.labels, max.iter, step.size,}

  \code{type.tnorm, type.snorm, type.implication.func ,
  name)}

  \item \code{SBC}:

  \code{list(r.a, eps.high, eps.low, name)}

  \item \code{FS.HGD}:

  \code{list(num.labels, max.iter, step.size,
  alpha.heuristic,}

  \code{type.tnorm, type.snorm, type.implication.func,
  name)}

  \item \code{FRBCS.W} and \code{FRBCS.CHI}:

  \code{list(num.labels, type.mf, type.tnorm,}

  \code{type.snorm, type.implication.func, name)}

  \item \code{DENFIS} method:

  \code{list(Dthr, max.iter, step.size, d, name)}

  \item \code{GFS.FR.MOGUL}:

  \code{list(persen_cross, max.iter, max.gen, max.tune,}

  \code{persen_mutant, epsilon, name)}

  \item \code{GFS.THRIFT} method:

  \code{list(popu.size, num.labels, persen_cross,}

  \code{max.gen, persen_mutant, type.tnorm, type.snorm,
  type.defuz,}

  \code{type.implication.func, name)}

  \item \code{GFS.GCCL}:

  \code{list(popu.size, num.class, num.labels,
  persen_cross,}

  \code{max.gen, persen_mutant, name)}

  \item \code{FH.GBML}:

  \code{list(popu.size, max.num.rule, num.class,
  persen_cross,}

  \code{max.gen, persen_mutant, p.dcare, p.gccl, name)}

  \item \code{SLAVE}:

  \code{list(num.class, num.labels, persen_cross,
  max.iter,}

  \code{max.gen, persen_mutant, k.lower, k.upper, epsilon,
  name)}

  \item \code{GFS.LT.RS}:

  \code{list(popu.size, num.labels, persen_mutant,
  max.gen,}

  \code{mode.tuning, type.tnorm, type.snorm,
  type.implication.func,}

  \code{type.defuz, rule.selection, name)}

  }

  \bold{Description of the \code{control} Parameters}
  \itemize{ \item \code{num.labels}: a positive integer to
  determine the number of labels (linguistic terms).  The
  default value is 7. \item \code{type.mf}: the following
  type of the membership function. The default value is
  \code{GAUSSIAN}. For more detail, see \code{fuzzifier}.
  \itemize{ \item \code{TRIANGLE}: it refers triangular
  shape.  \item \code{TRAPEZOID}: it refers trapezoid
  shape.  \item \code{GAUSSIAN}: it refers gaussian shape.
  \item \code{SIGMOID}: it refers sigmoid.  \item
  \code{BELL}: it refers generalized bell.  } \item
  \code{type.defuz}: the type of the defuzzification method
  as follows. The default value is \code{WAM}. For more
  detail, see \code{defuzzifier}.  \itemize{ \item
  \code{WAM}: the weighted average method.  \item
  \code{FIRST.MAX}: the first maxima.  \item
  \code{LAST.MAX}: the last maxima.  \item \code{MEAN.MAX}:
  the mean maxima.  \item \code{COG}: the modified center
  of gravity (COG).  } \item \code{type.tnorm}: the type of
  conjunction operator (t-norm). The following are options
  of t-norm available. For more detail, please have a look
  at \code{\link{inference}}.  The default value is
  \code{MIN}.  \itemize{ \item \code{MIN} means standard
  type (minimum).  \item \code{HAMACHER} means Hamacher
  product.  \item \code{YAGER} means Yager class (with tao
  = 1).  \item \code{PRODUCT} means product.  \item
  \code{BOUNDED} mean bounded product.  } \item
  \code{type.snorm}: the type of disjunction operator
  (s-norm). The following are options of s-norm available.
  For more detail, please have a look at
  \code{\link{inference}}.  The default value is
  \code{MAX}.  \itemize{ \item \code{MAX} means standard
  type (maximum).  \item \code{HAMACHER} means Hamacher
  sum.  \item \code{YAGER} means Yager class (with tao =
  1).  \item \code{SUM} means sum.  \item \code{BOUNDED}
  mean bounded sum.  } \item \code{type.implication.func}:
  the type of implication function. The following are
  options of implication function available:
  \code{DIENES_RESHER}, \code{LUKASIEWICZ}, \code{ZADEH},
  \code{GOGUEN}, \code{GODEL}, \code{SHARP},
  \code{MIZUMOTO}, \code{DUBOIS_PRADE}, and \code{MIN}.
  For more detail, please have a look at \code{\link{WM}}.
  The default value is \code{ZADEH}. \item \code{name}: a
  name for the model. The default value is \code{"sim-0"}.
  \item \code{max.iter}: a positive integer to determine
  the maximal number of iterations.  The default value is
  10. \item \code{step.size}: the step size of the gradient
  descent, a real number between 0 and 1.  The default
  value is 0.01. \item \code{r.a}: a positive constant
  which is effectively the radius defining a neighborhood.
  The default value is 0.5. \item \code{eps.high}: an upper
  threshold value. The default value is 0.5. \item
  \code{eps.low}: a lower threshold value. The default
  value is 0.15. \item \code{alpha.heuristic}: a positive
  real number representing a heuristic value.  The default
  value is 1. \item \code{Dthr}: the threshold value for
  the envolving clustering method (ECM), between 0 and 1.
  The default value is 0.1. \item \code{d}: a parameter for
  the width of the triangular membership function.  The
  default value is 2. \item \code{persen_cross}: a
  probability of crossover. The default value is 0.6. \item
  \code{max.gen}: a positive integer to determine the
  maximal number of generations of the genetic algorithm.
  The default value is 10. \item \code{max.tune}: a
  positive integer to determine the maximal number of
  tuning iterations.  The default value is 10. \item
  \code{persen_mutant}: a probability of mutation. The
  default value is 0.3. \item \code{epsilon}: a real number
  between 0 and 1 representing the level of generalization.
  A high epsilon can lead to overfitting. The default value
  is 0.9. \item \code{popu.size}: the size of the
  population which is generated in each generation. The
  default value is 10. \item \code{max.num.rule}: the
  maximum size of the rules. The default value is 5. \item
  \code{num.class}: the number of classes. \item
  \code{p.dcare}: a probability of "don't care" attributes.
  The default value is 0.5. \item \code{p.gccl}: a
  probability of the GCCL process. The default value is
  0.5. \item \code{k.lower}: a lower bound of the noise
  threshold with interval between 0 and 1. The default
  value is 0. \item \code{k.upper}: an upper bound of the
  noise threshold with interval between 0 and 1. The
  default value is 1. \item \code{mode.tuning}: a type of
  lateral tuning which are \code{"LOCAL"} or
  \code{"GLOBAL"}. The default value is \code{"GLOBAL"}.
  \item \code{rule.selection}:a boolean value representing
  whether performs rule selection or not.  The default
  value is \code{"TRUE"}. }}
}
\value{
  The \code{\link{frbs-object}}.
}
\description{
  This is one of the central functions of the package. This
  function is used to generate/learn the model from
  numerical data using fuzzy rule-based systems.
}
\details{
  This function makes accessible all learning methods that
  are implemented in this package. All of the methods use
  this function as interface for the learning stage, so
  users do not need to call other functions in the learning
  phase. In order to obtain good results, users need to
  adjust some parameters such as the number of labels, the
  type of the shape of the membership function, the maximal
  number of iterations, the step size of the gradient
  descent, or other method-dependent parameters which are
  collected in the \code{control} parameter. After creating
  the model using this function, it can be used to predict
  new data with \code{\link{predict}}.
}
\examples{
##################################
## I. Regression Problem
## Suppose data have two input variables and one output variable.
## We separate them into training, fitting, and testing data.
## data.train, data.fit, data.test, and range.data are inputs
## for all regression methods.
###################################
## Take into account that the simulation might take a long time
## depending on the hardware you are using. The chosen parameters
## may not be optimal.
## Data must be in data.frame or matrix form and the last column
## is the output variable/attribute.
## The training data must be expressed in numbers (numerical data).
data.train <- matrix(c(5.2, -8.1, 4.8, 8.8, -16.1, 4.1, 10.6, -7.8, 5.5, 10.4, -29.0,
                      5.0, 1.8, -19.2, 3.4, 12.7, -18.9, 3.4, 15.6, -10.6, 4.9, 1.9,
                      -25.0, 3.7, 2.2, -3.1, 3.9, 4.8, -7.8, 4.5, 7.9, -13.9, 4.8,
                      5.2, -4.5, 4.9, 0.9, -11.6, 3.0, 11.8, -2.1, 4.6, 7.9, -2.0,
                      4.8, 11.5, -9.0, 5.5, 10.6, -11.2, 4.5, 11.1, -6.1, 4.7, 12.8,
                      -1.0, 6.6, 11.3, -3.6, 5.1, 1.0, -8.2, 3.9, 14.5, -0.5, 5.7,
                      11.9, -2.0, 5.1, 8.1, -1.6, 5.2, 15.5, -0.7, 4.9, 12.4, -0.8,
                      5.2, 11.1, -16.8, 5.1, 5.1, -5.1, 4.6, 4.8, -9.5, 3.9, 13.2,
                      -0.7, 6.0, 9.9, -3.3, 4.9, 12.5, -13.6, 4.1, 8.9, -10.0,
                      4.9, 10.8, -13.5, 5.1), ncol = 3, byrow = TRUE)
colnames(data.train) <- c("inp.1", "inp.2", "out.1")

data.fit <- data.train[, -ncol(data.train)]

data.test <- matrix(c(10.5, -0.9, 5.8, -2.8, 8.5, -0.6, 13.8, -11.9, 9.8, -1.2, 11.0,
                     -14.3, 4.2, -17.0, 6.9, -3.3, 13.2, -1.9), ncol = 2, byrow = TRUE)

range.data <- matrix(apply(data.train, 2, range), nrow = 2)

#############################################################
## I.1 Example: Implementation of Wang & Mendel
#############################################################
method.type <- "WM"

## collect control parameters into a list
## num.labels = 3 means we define 3 as the number of linguistic terms
control.WM <- list(num.labels = 3, type.mf = "GAUSSIAN", type.tnorm = "MIN", type.snorm = "MAX",
type.defuz = "WAM", type.implication.func = "ZADEH", name = "Sim-0")

## generate the model and save it as object.WM
object.WM <- frbs.learn(data.train, range.data, method.type, control.WM)

#############################################################
## I.2 Example: Implementation of SBC
#############################################################
\dontrun{method.type <- "SBC"
control.SBC <- list(r.a = 0.5, eps.high = 0.5, eps.low = 0.15, name = "Sim-0")

object.SBC <- frbs.learn(data.train, range.data, method.type, control.SBC)

#############################################################
## I.3 Example: Implementation of HYFIS
#############################################################
method.type <- "HYFIS"

control.HYFIS <- list(num.labels = 5, max.iter = 50, step.size = 0.01, type.tnorm = "MIN",
                      type.snorm = "MAX", type.defuz = "COG",
                      type.implication.func = "ZADEH", name = "Sim-0")

object.HYFIS <- frbs.learn(data.train, range.data, method.type, control.HYFIS)

#############################################################
## I.4 Example: Implementation of ANFIS
#############################################################
method.type <- "ANFIS"

control.ANFIS <- list(num.labels = 5, max.iter = 10, step.size = 0.01, type.tnorm = "MIN",
                     type.snorm = "MAX", type.implication.func = "ZADEH", name = "Sim-0")

object.ANFIS <- frbs.learn(data.train, range.data, method.type, control.ANFIS)

#############################################################
## I.5 Example: Implementation of DENFIS
#############################################################

control.DENFIS <- list(Dthr = 0.1, max.iter = 10, step.size = 0.001, d = 2,
                       name = "Sim-0")
method.type <- "DENFIS"

object.DENFIS <- frbs.learn(data.train, range.data, method.type, control.DENFIS)

#############################################################
## I.6 Example: Implementation of FIR.DM
#############################################################
method.type <- "FIR.DM"

control.DM <- list(num.labels = 5, max.iter = 10, step.size = 0.01, type.tnorm = "MIN",
                     type.snorm = "MAX", type.implication.func = "ZADEH", name = "Sim-0")
object.DM <- frbs.learn(data.train, range.data, method.type, control.DM)

#############################################################
## I.7 Example: Implementation of FS.HGD
#############################################################
method.type <- "FS.HGD"

control.HGD <- list(num.labels = 5, max.iter = 10, step.size = 0.01,
               alpha.heuristic = 1, type.tnorm = "MIN", type.snorm = "MAX",
               type.implication.func = "ZADEH", name = "Sim-0")
object.HGD <- frbs.learn(data.train, range.data, method.type, control.HGD)

#############################################################
## I.8 Example: Implementation of GFS.FR.MOGUL
#############################################################
method.type <- "GFS.FR.MOGUL"

control.GFS.FR.MOGUL <- list(persen_cross = 0.6,
                    max.iter = 5, max.gen = 2, max.tune = 2, persen_mutant = 0.3,
                    epsilon = 0.8, name="sim-0")
object.GFS.FR.MOGUL <- frbs.learn(data.train, range.data,
                       method.type, control.GFS.FR.MOGUL)

#############################################################
## I.9 Example: Implementation of Thrift's method (GFS.THRIFT)
#############################################################
method.type <- "GFS.THRIFT"

control.Thrift <- list(popu.size = 6, num.labels = 3, persen_cross = 1,
                      max.gen = 5, persen_mutant = 1, type.tnorm = "MIN",
                      type.snorm = "MAX", type.defuz = "COG",
                      type.implication.func = "ZADEH", name="sim-0")
object.Thrift <- frbs.learn(data.train, range.data, method.type, control.Thrift)

##############################################################
## I.10 Example: Implementation of
##      genetic for lateral tuning and rule selection (GFS.LT.RS)
#############################################################
## Set the method and its parameters
method.type <- "GFS.LT.RS"

control.lt.rs <- list(popu.size = 5, num.labels = 5, persen_mutant = 0.3,
	               max.gen = 10, mode.tuning = "LOCAL", type.tnorm = "MIN", type.snorm = "MAX",
                type.implication.func = "ZADEH", type.defuz = "WAM",
                rule.selection = TRUE, name="sim-0")

## Generate fuzzy model
object.lt.rs <- frbs.learn(data.train, range.data, method.type, control.lt.rs)}

#############################################################
## II. Classification Problems
#############################################################
## The iris dataset is shuffled and divided into training and
## testing data. Bad results in the predicted values may result
## from casual imbalanced classes in the training data.
## Take into account that the simulation may take a long time
## depending on the hardware you use.
## One may get better results with other parameters.
## Data are in data.frame or matrix form and the last column is
## the output variable/attribute
## The data must be expressed in numbers (numerical data).

data(iris)
irisShuffled <- iris[sample(nrow(iris)),]
irisShuffled[,5] <- unclass(irisShuffled[,5])
tra.iris <- irisShuffled[1:105,]
tst.iris <- irisShuffled[106:nrow(irisShuffled),1:4]
real.iris <- matrix(irisShuffled[106:nrow(irisShuffled),5], ncol = 1)

## Please take into account that the interval needed is the range of input data only.
range.data.input <- matrix(apply(iris[, -ncol(iris)], 2, range), nrow = 2)

#########################################################
## II.1 Example: Implementation of FRBCS with weighted factor based on Ishibuchi's method
###############################################################
## generate the model
method.type <- "FRBCS.W"
control <- list(num.labels = 3, type.mf = "TRIANGLE", type.tnorm = "MIN",
               type.snorm = "MAX", type.implication.func = "ZADEH", name = "sim-0")

object <- frbs.learn(tra.iris, range.data.input, method.type, control)

## conduct the prediction process
res.test <- predict(object, tst.iris)

#########################################################
## II.2 Example: Implementation of FRBCS based on Chi's method
###############################################################
## generate the model
\dontrun{method.type <- "FRBCS.CHI"
control <- list(num.labels = 7, type.mf = "TRIANGLE", type.tnorm = "MIN",
               type.snorm = "MAX", type.implication.func = "ZADEH", name = "sim-0")

object <- frbs.learn(tra.iris, range.data.input, method.type, control)

## conduct the prediction process
res.test <- predict(object, tst.iris)

#########################################################
## II.3 The example: Implementation of GFS.GCCL
###############################################################
method.type <- "GFS.GCCL"

control <- list(popu.size = 5, num.class = 3, num.labels = 5, persen_cross = 0.9,
                    max.gen = 2, persen_mutant = 0.3,
                    name="sim-0")
## Training process
## The main result of the training is a rule database which is used later for prediction.
object <- frbs.learn(tra.iris, range.data.input, method.type, control)

## Prediction process
res.test <- predict(object, tst.iris)

#########################################################
## II.4 Example: Implementation of FH.GBML
###############################################################
method.type <- "FH.GBML"

	control <- list(popu.size = 5, max.num.rule = 5, num.class = 3,
				persen_cross = 0.9, max.gen = 2, persen_mutant = 0.3, p.dcare = 0.5,
             p.gccl = 1, name="sim-0")

	## Training process
	## The main result of the training is a rule database which is used later for prediction.
	object <- frbs.learn(tra.iris, range.data.input, method.type, control)

	## Prediction process
	res.test <- predict(object, tst.iris)

#########################################################
## II.5 The example: Implementation of SLAVE
###############################################################
method.type <- "SLAVE"

	control <- list(num.class = 3, num.labels = 5,
				persen_cross = 0.9, max.iter = 5, max.gen = 3, persen_mutant = 0.3,
             k.lower = 0.25, k.upper = 0.75, epsilon = 0.1, name="sim-0")

	## Training process
	## The main result of the training is a rule database which is used later for prediction.
	object <- frbs.learn(tra.iris, range.data.input, method.type, control)

	## Prediction process
	res.test <- predict(object, tst.iris)}
}
\seealso{
  \code{\link{predict}} for the prediction phase, and the
  following main functions of each of the methods for
  theoretical background and references: \code{\link{WM}},
  \code{\link{SBC}}, \code{\link{HyFIS}},
  \code{\link{ANFIS}}, \code{\link{FIR.DM}},
  \code{\link{DENFIS}}, \code{\link{FS.HGD}},
  \code{\link{FRBCS.W}}, \code{\link{FRBCS.CHI}},
  \code{\link{GFS.FR.MOGUL}}, \code{\link{GFS.Thrift}},
  \code{\link{GFS.GCCL}}, \code{\link{FH.GBML}},
  \code{\link{GFS.LT.RS}}, and \code{\link{SLAVE}}.
}

