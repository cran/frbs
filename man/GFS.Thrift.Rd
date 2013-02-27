\name{GFS.Thrift}
\alias{GFS.Thrift}
\title{GFS.Thrift model building}
\usage{
  GFS.Thrift(data.train, popu.size, range.data, num.labels,
    persen_cross, persen_mutant, max.gen, range.data.ori)
}
\arguments{
  \item{data.train}{a matrix(m x n) of data for the
  training process, where m is the number of instances and
  n is the number of variables; the last column is the
  output variable.}

  \item{popu.size}{the size of the population which is
  generated in each generation.}

  \item{range.data}{a matrix containing the ranges of the
  normalized data.}

  \item{num.labels}{a matrix describing the number of fuzzy
  terms.}

  \item{persen_cross}{a real number between 0 and 1
  representing the probability of crossover.}

  \item{persen_mutant}{a real number between 0 and 1
  representing the probability of mutation.}

  \item{max.gen}{the maximal number of generations for the
  genetic algorithm.}

  \item{range.data.ori}{a matrix containing the ranges of
  the original data.}
}
\description{
  This is the internal function that implements the
  Thrift's technique based on a genetic algorithm. It is
  used to tackle regression tasks. Users do not need to
  call it directly, but just use \code{\link{frbs.learn}}
  and \code{\link{predict}}.
}
\details{
  This method was developed by Thrift using Mamdani's model
  as fuzzy IF-THEN rules. In this method, we consider a
  table as a genotype with alleles that are fuzzy set
  indicators over the output domain. The phenotype is
  produced by the behavior produced by the fuzzification,
  max-* composition, and defuzzification operations. A
  chromosome (genotype) is formed from the decision table
  by going rowwise and producing a string of numbers from
  the code set. Standard crossover and mutation operators
  can act on these string.
}
\references{
  P. Thrift, "Fuzzy logic synthesis with genetic
  algorithms", In Proceedings of the Fourth International
  Conference on Genetic Algorithms (ICGA91), San Diego
  (United States of America), pp. 509 - 513 (1991).
}
\seealso{
  \code{\link{GFS.Thrift.test}}, \code{\link{frbs.learn}},
  and \code{\link{predict}}
}

