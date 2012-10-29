\docType{package}
\name{frbs-package}
\alias{frbs}
\alias{frbs-package}
\title{Getting started with the frbs package}
\description{
  Fuzzy rule-based systems (FRBS) are based on the fuzzy
  concept proposed by Zadeh in 1965, which represents the
  reasoning of human experts in production rules (a set of
  IF-THEN rules) to handle real-life problems from domains
  such as control, prediction and inference, data mining,
  bioinformatics data processing, robotics, and speech
  recognition. FRBSs are also known as fuzzy inference
  systems and fuzzy models. When applied to specific tasks,
  they may also be known under specific names such as fuzzy
  associative memories or fuzzy controllers. In this
  package, we consider systems with multiple inputs and a
  single output (MISO), with real-valued data.
}
\details{
  FRBSs are a competitive alternative to other classic
  models and algorithms in order to solve classification
  and regression problems. Generally, an FRBS consists of
  four functional parts: \itemize{ \item a fuzzification
  interface which transforms the crisp inputs into degrees
  of membership function of the fuzzy term of a linguistic
  variable. See \code{\link{fuzzifier}}. \item a knowledge
  base composed of a database and a rulebase. While the
  database includes the fuzzy set definitions, the rulebase
  contains the fuzzy IF-THEN rules. We will represent the
  knowledge as a set of rules. Each one has the following
  structure.

  \code{IF premise (antecedent), THEN conclusion
  (consequent)}. See \code{\link{rulebase}}. \item an
  inference engine which performs the inference operations
  on the fuzzy IF-THEN rules. There are two kinds of
  inference for fuzzy systems based on linguistic rules:
  The Mamdani and the Takagi Sugeno Kang model. See
  \code{\link{inference}}. \item a defuzzification process
  to obtain the crisp values. There are several methods for
  defuzzification such as the weighted average, centroid,
  etc. See \code{\link{defuzzifier}}. }

  Since it may be difficult to obtain information from
  human experts in the form required, an alternative and
  effective way to acquire the knowledge is to generate the
  fuzzy IF-THEN rules automatically from the input data. In
  general, when modeling an FRBS, there are two important
  processes which should be conducted, namely structure
  identification and parameter estimation. Structure
  identification is a process to find appropriate fuzzy
  IF-THEN rules and to determine the overall number of
  rules. Parameter estimation is applied to tune the
  parameters on the consequent and/or antecedent parts of
  the fuzzy IF-THEN rules. Many approaches have been
  proposed in order to perform this modeling such as a
  table-lookup scheme, heuristic procedures, neuro-fuzzy
  techniques, clustering methods, genetic algorithms, least
  squares methods, gradient descent, etc. In this package,
  the following approaches to generate fuzzy IF-THEN rules
  have been implemented: \itemize{ \item Wang and Mendel's
  Method (WM): \code{\link{WM}} and \code{\link{frbs.eng}}.
  \item The genetic fuzzy system (GFS): \code{\link{GFS}}
  and \code{\link{frbs.eng}}. \item The hybrid neural fuzzy
  inference system (HyFIS): \code{\link{HyFIS}} and
  \code{\link{frbs.eng}}. \item The adaptive-network-based
  fuzzy inference system (ANFIS): \code{\link{ANFIS}} and
  \code{\link{frbs.eng}}. \item The fuzzy rule-based
  classification system (frcs): \code{\link{frcs}} and
  \code{\link{frcs.eng}}. \item The dynamic evolving
  neural-fuzzy inference system (DENFIS):
  \code{\link{DENFIS}} and \code{\link{DENFIS.eng}}. \item
  The subtractive clustering and fuzzy c-means (SBC):
  \code{\link{SBC}} and \code{\link{SBC.test}}. \item The
  fuzzy inference rules by descent method (DM):
  \code{\link{DM}} and \code{\link{frbs.eng}}. \item The
  FRBS using heuristics and gradient descent method (HGD):
  \code{\link{HGD}} and \code{\link{frbs.eng}}. \item The
  multi-stage genetic fuzzy systems (MSGFS) based on
  iterative rule learning approach: \code{\link{MSGFS}} and
  \code{\link{MSGFS.test}}. } The functions documented in
  the manual for the single methods are all called
  internally by \code{\link{frbs.learn}}, which is the
  central function of the package. However, in the
  documentation of each of the internal learning functions,
  we give some theoretical background and references to the
  original literature.

  \bold{Usage of the package:}

  If you have problems using the package, find a bug, or
  have suggestions, please contact the package maintainer
  by email, instead of writing to the general R lists or to
  other internet forums and mailing lists.

  The main functions of the package are the following:
  \itemize{ \item The function \code{\link{frbs.learn}}
  allows to generate the model by creating fuzzy IF-THEN
  rules or cluster centers from training data. The
  different algorithms mentioned above are all accessible
  through this function. The outcome of the function is an
  \code{\link{frbs-object}}. \item Even though the main
  purpose of this package is to generate the FRBS models
  automatically, we provide the function
  \code{\link{frbs.gen}}, which can be used to build a
  model manually without using a learning method. \item The
  purpose of the function \code{\link{predict}} is to
  obtain predicted values according to the testing data and
  the model (analogous to the \code{predict} function that
  is implemented in many other R packages). \item There
  exist functions \code{\link{summary.frbs}} and
  \code{\link{plotMF}} to show a summary about an
  \code{\link{frbs-object}}, and to plot the shapes of the
  membership functions. }

  To get started with the package, the user can have a look
  at the examples included in the documentation of the
  functions \code{\link{frbs.learn}} for generating models
  and \code{\link{predict}} for the prediction phase.

  Also, there are many demos that ship with the package. To
  get a list of them, type:

  \code{demo()}

  Then, to start a demo, type
  \code{demo(<demo_name_here>)}. All the demos are present
  as R scripts in the package sources in the "demo"
  subdirectory.

  Currently, there are the following demos available:

  Regression using the Gas Furnance dataset:

  \code{demo(WM.GasFur)}, \code{demo(SBC.GasFur)},
  \code{demo(ANFIS.GasFur)},

  \code{demo(HGD.GasFur)}, \code{demo(DENFIS.GasFur)},
  \code{demo(HyFIS.GasFur)},

  \code{demo(GFS.GasFur)}, \code{demo(DM.GasFur)},
  \code{demo(MSGFS.GasFur)}

  Regression using the Mackey-Glass dataset:

  \code{demo(WM.MG1000)}, \code{demo(SBC.MG1000)},
  \code{demo(ANFIS.MG1000)},

  \code{demo(HGD.MG1000)}, \code{demo(DENFIS.MG1000)},
  \code{demo(HyFIS.MG1000)},

  \code{demo(GFS.MG1000)}, \code{demo(DM.MG1000)},
  \code{demo(MSGFS.MG1000)}

  Classification using the Iris dataset:

  \code{demo(WM.Iris)}, \code{demo(SBC.Iris)},
  \code{demo(ANFIS.Iris)},

  \code{demo(HGD.Iris)}, \code{demo(DENFIS.Iris)},
  \code{demo(HyFIS.Iris)},

  \code{demo(GFS.Iris)}, \code{demo(DM.Iris)},
  \code{demo(MSGFS.Iris)}, \code{demo(frcs.Iris)}

  The Gas Furnance data and Mackey-Glass data are included
  in the package, please see \code{\link{frbsData}}. The
  Iris data is the standard Iris dataset that ships with R.

  Also have a look at the package webpage
  \url{http://sci2s.ugr.es/dicits/software/FRBS}, where we
  provide a more extensive introduction as well as
  additional explanations of the procedures.
}
\author{
  Lala Septem Riza \email{lala.s.riza@decsai.ugr.es},

  Christoph Bergmeir \email{c.bergmeir@decsai.ugr.es},

  Francisco Herrera Triguero \email{herrera@decsai.ugr.es},

  and Jose Manuel Benitez Sanchez
  \email{j.m.benitez@decsai.ugr.es}

  DiCITS Lab, SCI2S group, DECSAI, University of Granada.

  \url{http://dicits.ugr.es}, \url{http://sci2s.ugr.es}
}
\references{
  C.C. Lee, "Fuzzy Logic in control systems: Fuzzy Logic
  controller part I", IEEE Trans. Syst., Man, Cybern., vol.
  20, no.2, pp. 404 - 418 (1990).

  C.C. Lee, "Fuzzy Logic in control systems: Fuzzy Logic
  controller part II",IEEE Trans. Syst., Man, Cybern., vol.
  20, no.2, pp. 419 - 435 (1990).

  D. Dubois and H. Prade, "Fuzzy Sets and Systems: Theory
  and Applications," New York: Academic (1980).

  L.A. Zadeh, "Fuzzy sets", Information and Control, vol.
  8, pp. 338 - 353 (1965).

  Mamdani, E. H., & Assilian, S., "An experiment in
  linguistic synthesis with a fuzzy logic controller,"
  International Journal of Man Machine Studies, 7(1), pp. 1
  - 13 (1975).

  M. Sugeno and G. T. Kang, "Structure identification of
  fuzzy model," Fuzzy Sets Syst., vol. 28, pp. 15 - 33
  (1988).

  Takagi, T., Sugeno, M., "Fuzzy identification of systems
  and its application to modelling and control", IEEE
  Transactions on Systems, Man and Cybernetics 15(1), pp.
  116 - 132 (1985).

  W. Pedrycz, "Fuzzy Control and Fuzzy Systems," New York:
  Wiley (1989).
}
\seealso{
  \code{\link{frbs.learn}} and \code{\link{predict}} for
  the learning and prediction stage, respectively.
}

