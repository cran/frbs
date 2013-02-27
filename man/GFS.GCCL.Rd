\name{GFS.GCCL}
\alias{GFS.GCCL}
\title{GFS.GCCL model building}
\usage{
  GFS.GCCL(data.train, popu.size, range.data.input,
    num.labels, persen_cross, persen_mutant, max.gen,
    range.data.ori)
}
\arguments{
  \item{data.train}{a matrix(m x n) of data for the
  training process, where m is the number of instances and
  n is the number of variables; the last column is the
  output variable.}

  \item{popu.size}{the size of the population which is
  generated in each generation.}

  \item{range.data.input}{a matrix containing the ranges of
  the normalized input data.}

  \item{num.labels}{a matrix describing the number of fuzzy
  terms.}

  \item{persen_cross}{a real number between 0 and 1
  representing the probability of crossover.}

  \item{persen_mutant}{a real number between 0 and 1
  representing the probability of mutation.}

  \item{max.gen}{the maximal number of generations for the
  genetic algorithm.}

  \item{range.data.ori}{a matrix containing the ranges of
  the input data.}
}
\description{
  This is the internal function that implements the
  Ishibuchi's method based on genetic
  cooperative-competitive learning (GFS.GCCL). It is used
  to handle classification tasks. Users do not need to call
  it directly, but just use \code{\link{frbs.learn}} and
  \code{\link{predict}}.
}
\details{
  This method is based on Ishibuchi's method. In this
  method, a chromosome describes each linguistic IF-THEN
  rule using integer as its representation of the
  antecedent part. In the consequent part of the fuzzy
  rules, the heuristic method is applied to automatically
  generate the class. The evaluation is calculated for each
  rule which means that the performance is not based on the
  entire rule set. The outline of the method is as follows.
  \itemize{ \item Step 1: Generate an initial population of
  fuzzy IF-THEN rules. \item Step 2: Evaluate each fuzzy
  IF-THEN rule in the current population. \item Step 3:
  Generate new fuzzy IF-THEN rules by genetic operators.
  \item Step 4: Replace a part of the current population
  with the newly generated rules. \item Step 5: Terminate
  the algorithm if a stopping condition is satisfied,
  otherwise return to Step 2. } Additionally, to handle
  high dimensional data, this method uses "don't care"
  attributes on the antecedent fuzzy set.
}
\references{
  H. Ishibuchi, et.al, "Performance Evaluation of fuzzy
  classifier systems for multidimensional pattern
  classification problems", IEEE trans. on Systems, Man,
  and Cybernetics - Part B: Sybernetics, vol. 29. no. 5,
  pp. 601 - 618 (1999).
}

