\name{GFS}
\alias{GFS}
\title{GFS model building}
\usage{
  GFS(data.train, range.data, num.labels, popu.size,
    persen_cross, persen_mutant, max.iter,
    classification = FALSE, range.data.ori)
}
\arguments{
  \item{data.train}{a matrix(m x n) of data for the
  training process, where m is the number of instances and
  n is the number of variables; the last column is the
  output variable.}

  \item{range.data}{a matrix(2 x n) containing the range of
  the normalized data, where n is the number of variables,
  and first and second rows are the minimum and maximum
  value, respectively.}

  \item{num.labels}{a matrix(1 x n), whose elements
  represent the number of labels (fuzzy terms); n is the
  number of variables.}

  \item{popu.size}{the size of the population which is
  generated in each generation.}

  \item{persen_cross}{a real number between 0 and 1
  representing the probability of crossover.}

  \item{persen_mutant}{a real number between 0 and 1
  representing the probability of mutation.}

  \item{max.iter}{the maximal number of iterations.}

  \item{classification}{a boolean representing whether it
  is a classification problem or not.}

  \item{range.data.ori}{a matrix(2 x n) containing the
  range of the original data, where n is the number of
  variables, and first and second rows are the minimum and
  maximum value, respectively.}
}
\description{
  This is the internal function that implements the genetic
  fuzzy systems (GFS) model. Users do not need to call it
  directly, but just use \code{\link{frbs.learn}} and
  \code{\link{predict}}
}
\details{
  This function is adopted from Jindrich Liska and Stephen
  S. Melsheimer's paper. GFS implements a genetic algorithm
  in order to determine the structure of the fuzzy IF-THEN
  rules and the membership function parameters. Genetic
  algorithms (GA) are a stochastic optimization technique
  that mimics natural selection. We conduct selection,
  uniform crossover, mutation, and the creep operator for
  the GA. One important thing when using a GA is to choose
  the representation of the solutions within the
  individuals. We split every individual to become two
  matrices which are a matrix of membership function
  parameters for all variables and a matrix of the
  structure of the fuzzy IF-THEN rules. Initially, one
  individual is generated using the techniques of Wang and
  Mendel, and the rest is generated at random. So, this
  method is always better or at least the same as the Wang
  and Mendel method. The fitness function is computed as
  the root mean square error for each individual.
}
\references{
  Davis L., "Adaptive operator probabilities in genetic
  algorithm," ICGA'89 (1989).

  Holland J. H., "Adaptation in natural and artificial
  systems," University of Michigan Press, Ann Arbor (1975).

  Jindrich Liska and Stephen S. Melsheimer, "Complete
  design of fuzzy logic systems using genetic algorithm,"
  Fuzzy Systems, IEEE World Congress on Computational
  Intelligence (1994).

  Jindrich Liska and Stephen S. Melsheimer, "Design of
  fuzzy logic systems for nonlinear model identification,"
  American Control Conference (1994).
}

