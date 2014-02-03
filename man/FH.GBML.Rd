\name{FH.GBML}
\alias{FH.GBML}
\title{FH.GBML model building}
\usage{
  FH.GBML(data.train, popu.size = 10, max.num.rule = 5,
    persen_cross = 0.6, persen_mutant = 0.3, max.gen = 10,
    num.class, range.data.input, p.dcare = 0.5,
    p.gccl = 0.5)
}
\arguments{
  \item{data.train}{a matrix (\eqn{m \times n}) of
  normalized data for the training process, where \eqn{m}
  is the number of instances and \eqn{n} is the number of
  variables; the last column is the output variable.  Note
  the data must be normalized between 0 and 1.}

  \item{popu.size}{the size of the population which is
  generated in each generation.}

  \item{max.num.rule}{the maximum number of rules.}

  \item{persen_cross}{a real number between 0 and 1
  determining the probability of crossover.}

  \item{persen_mutant}{a real number between 0 and 1
  determining the probability of mutation.}

  \item{max.gen}{the maximal number of generations for the
  genetic algorithms.}

  \item{num.class}{a number of the classes.}

  \item{range.data.input}{a matrix containing the ranges of
  the normalized input data.}

  \item{p.dcare}{a probability of "don't care" attributes
  occurred.}

  \item{p.gccl}{a probability of GCCL process occurred.}
}
\description{
  This is the internal function that implements the
  Ishibuchi's method based on hybridization of genetic
  cooperative-competitive learning (GCCL) and Pittsburgh
  (FH.GBML). It is used to solve classification tasks.
  Users do not need to call it directly, but just use
  \code{\link{frbs.learn}} and \code{\link{predict}}.
}
\details{
  This method is based on Ishibuchi's method using the
  hybridization of GCCL and the Pittsburgh approach for
  genetic fuzzy systems. The algorithm of this method is as
  follows. \itemize{ \item Step 1: Generate population
  where each individual in the population is a fuzzy rule
  set. \item Step 2: Calculate the fitness value of each
  rule set in the current population. \item Step 3:
  Generate new rule sets by the selection, crossover, and
  mutation in the same manner as the Pittsburgh-style
  algorithm. Then, apply iterations of the GCCL to each of
  the generated rule sets with a probability. \item Step 4:
  Add the best rule set in the current population to newly
  generated rule sets to form the next population. \item
  Step 5: Return to Step 2 if the prespecified stopping
  condition is not satisfied. }
}
\references{
  H. Ishibuchi, T. Yamamoto, and T. Nakashima,
  "Hybridization of fuzzy GBML approaches for pattern
  classification problems," IEEE Trans. on Systems, Man,
  and Cybernetics-Part B: Cybernetics, vol. 35, no. 2, pp.
  359 - 365 (2005).
}

