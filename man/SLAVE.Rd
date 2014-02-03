\name{SLAVE}
\alias{SLAVE}
\title{SLAVE model building}
\usage{
  SLAVE(data.train, persen_cross = 0.6,
    persen_mutant = 0.3, max.iter = 10, max.gen = 10,
    num.labels, range.data.input, k.lower = 0.25,
    k.upper = 0.75, epsilon = 0.1)
}
\arguments{
  \item{data.train}{a matrix (\eqn{m \times n}) of
  normalized data for the training process, where \eqn{m}
  is the number of instances and \eqn{n} is the number of
  variables; The last column is the output variable.  Note
  the data must be normalized between 0 and 1.}

  \item{persen_cross}{a real number between 0 and 1
  representing the probability of crossover.}

  \item{persen_mutant}{a real number between 0 and 1
  representing the probability of mutation.}

  \item{max.iter}{the maximal number of iterations.}

  \item{max.gen}{the maximal number of generations for the
  genetic algorithm.}

  \item{num.labels}{a number of the linguistic terms.}

  \item{range.data.input}{a matrix containing the ranges of
  the normalized input data.}

  \item{k.lower}{a lower bound of the noise threshold.}

  \item{k.upper}{an upper bound of the noise threshold.}

  \item{epsilon}{a value between 0 and 1 representing the
  covering factor.}
}
\description{
  This is the internal function that implements the
  structural learning algorithm on vague environment
  (SLAVE). It is used to handle classification tasks. Users
  do not need to call it directly, but just use
  \code{\link{frbs.learn}} and \code{\link{predict}}.
}
\details{
  This method is adopted from A. Gonzalez and R. Perez's
  paper which is applied for classification problems. SLAVE
  is based on the iterative rule learning approach which
  means that we get only one fuzzy rule in each execution
  of the genetic algorithm. In order to eliminate the
  irrelevant variables in a rule, SLAVE has a structure
  composed of two parts: the first part is to represent the
  relevance of variables and the second one is to define
  values of the parameters. The following steps are
  conducted in order to obtain fuzzy rules: \itemize{ \item
  Step 1: Use the genetic algorithm process to obtain ONE
  RULE for the system. \item Step 2: Collect the rule into
  the final set of rules. \item Step 3: Check and penalize
  this rule. \item Step 4: If the stopping criteria is
  satisfied, the system returns the set of rules as
  solution. Otherwise, back to Step 1. } This method uses
  binary codes as representation of the population and
  applies the basic genetic operators, i.e., selection,
  crossover, and mutation on it. And, the best rule is
  obtained by calculating the degree of consistency and
  completeness.
}
\references{
  A. Gonzalez, R. Perez, "Selection of relevant features in
  a fuzzy genetic learning algorithm", IEEE Transactions on
  Systems, Man, and Cybernetics, Part B: Cybernetics, vol.
  31, no. 3, pp.  417 - 425 (2001).
}

