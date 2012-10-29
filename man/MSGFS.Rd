\name{MSGFS}
\alias{MSGFS}
\title{MSGFS model building}
\usage{
  MSGFS(data.train, popu.size, persen_cross, persen_mutant,
    max.iter, range.data.ori, epsilon)
}
\arguments{
  \item{data.train}{a matrix(m x n) of data for the
  training process, where m is the number of instances and
  n is the number of variables; the last column is the
  output variable.}

  \item{popu.size}{the size of the population which is
  generated in each generation.}

  \item{persen_cross}{a real number between 0 and 1
  representing the probability of crossover.}

  \item{persen_mutant}{a real number between 0 and 1
  representing the probability of mutation.}

  \item{max.iter}{the maximal number of iterations.}

  \item{range.data.ori}{a matrix containing the ranges of
  the original data.}

  \item{epsilon}{a real number between 0 and 1 representing
  the boundary of covering factor.}
}
\description{
  This is the internal function that implements the
  multi-stage genetic fuzzy systems (MSGFS) based on
  iterative rule learning approach. Users do not need to
  call it directly, but just use \code{\link{frbs.learn}}
  and \code{\link{predict}}.
}
\details{
  This method was proposed by Herrera et al. MSGFS
  implements a genetic algorithm in order to determine the
  structure of the fuzzy IF-THEN rules and the membership
  function parameters. There are two general types of fuzzy
  IF-THEN rules, namely the descriptive and the
  approximative/free semantic approaches. A descriptive
  approach means that the linguistic labels represent a
  real-world semantic; the linguistic labels are uniformly
  defined for all rules. In contrast, in the approximative
  approach there isn't any associated linguistic label.
  This method is based on the latter one. We model a fuzzy
  IF-THEN rule on a chromosome which consists of the
  parameter values of the membership function. So, every
  rule has its own membership function values. A population
  contains many such generated chromosomes, based on the
  iterative rule learning approach (IRL). IRL means that
  the chromosomes will be generated one by one, taking into
  account the fitness value and covering factor, until
  there are sufficient chromosomes in the population. After
  having obtained the population, the genetic algorithm is
  started, using the genetic operators selection, mutation,
  and crossover.
}
\references{
  Antonio Gonzalez and Francisco Herrera, "Multi-stage
  genetic fuzzy systems based on the iterative rule
  learning approach," Mathware & Soft Computing 4, pp. 233
  - 249 (1997).

  F. Herrera, M. Lozano, J.L. Verdegay, "Tuning fuzzy logic
  controllers by genetic algorithms," Internat. J. Approx.
  Reasoning 12, pp. 299 - 315 (1995).

  F. Herrera, M. Lozano, J.L. Verdegay, "Generating rules
  from examples using genetic algorithms," In: B. Bounchon,
  R. Yager, L. Zadeh (Eds.), Fuzzy Logic and Soft
  Computing, Word Scientific, pp. 11 - 20 (1995).

  F. Herrera, M. Lozano, and J.L. Verdegay, "A learning
  process for fuzzy control rules using genetic
  algorithms", Fuzzy Sets and Systems, 100, pp. 143 - 158
  (1998).

  O. Cordon, M.J. del Jesus, F. Herrera, M. Lozano, "MOGUL:
  A methodology to obtain genetic fuzzy rule-based systems
  under the iterative rule learning approach,"
  International Journal of Intelligent Systems, vol. 14,
  pp. 1123 - 1153 (1999).
}
\seealso{
  \code{\link{MSGFS.test}}, \code{\link{frbs.learn}}, and
  \code{\link{predict}}
}

