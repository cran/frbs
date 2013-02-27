\name{GFS.FR.MOGUL}
\alias{GFS.FR.MOGUL}
\title{GFS.FR.MOGUL model building}
\usage{
  GFS.FR.MOGUL(data.train, persen_cross, persen_mutant,
    max.iter, max.gen, max.tune, range.data.ori, epsilon)
}
\arguments{
  \item{data.train}{a matrix(m x n) of data for the
  training process, where m is the number of instances and
  n is the number of variables; the last column is the
  output variable.}

  \item{persen_cross}{a real number between 0 and 1
  determining the probability of crossover.}

  \item{persen_mutant}{a real number between 0 and 1
  determining the probability of mutation.}

  \item{max.iter}{the maximal number of iterations.}

  \item{max.gen}{the maximal number of generations of the
  genetic algorithm.}

  \item{max.tune}{the maximal number of tuning iterations.}

  \item{range.data.ori}{a matrix containing the ranges of
  the original data.}

  \item{epsilon}{a real number between 0 and 1 determining
  the boundary of covering factor.}
}
\description{
  This is the internal function that implements genetic
  fuzzy systems for fuzzy rule learning based on the MOGUL
  methodology (GFS.FR.MOGUL). It is used to solve
  regression tasks. Users do not need to call it directly,
  but just use \code{\link{frbs.learn}} and
  \code{\link{predict}}.
}
\details{
  This method was proposed by Herrera et al. GFS.FR.MOGUL
  implements a genetic algorithm determining the structure
  of the fuzzy IF-THEN rules and the membership function
  parameters. There are two general types of fuzzy IF-THEN
  rules, namely the descriptive and the approximative/free
  semantic approaches.  A descriptive approach means that
  the linguistic labels represent a real-world semantic;
  the linguistic labels are uniformly defined for all
  rules. In contrast, in the approximative approach there
  isn't any associated linguistic label. This method is
  based on the latter one. We model a fuzzy IF-THEN rule on
  a chromosome which consists of the parameter values of
  the membership function. So, every rule has its own
  membership function values. A population contains many
  such generated chromosomes, based on the iterative rule
  learning approach (IRL). IRL means that the chromosomes
  will be generated one by one, taking into account the
  fitness value and covering factor, until there are
  sufficient chromosomes in the population. After having
  obtained the population, the genetic algorithm is
  started, using the genetic operators selection, mutation,
  and crossover.
}
\references{
  F. Herrera, M. Lozano, and J.L. Verdegay, "A learning
  process for fuzzy control rules using genetic
  algorithms", Fuzzy Sets and Systems, vol. 100, pp. 143 -
  158 (1998).

  O. Cordon, M.J. del Jesus, F. Herrera, M. Lozano, "MOGUL:
  A methodology to obtain genetic fuzzy rule-based systems
  under the iterative rule learning approach",
  International Journal of Intelligent Systems, vol. 14,
  pp. 1123 - 1153 (1999).
}
\seealso{
  \code{\link{GFS.FR.MOGUL.test}},
  \code{\link{frbs.learn}}, and \code{\link{predict}}
}

