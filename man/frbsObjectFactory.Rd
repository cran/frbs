\name{frbsObjectFactory}
\alias{frbs-object}
\alias{frbsObjectFactory}
\title{The object factory for frbs objects}
\usage{
  frbsObjectFactory(mod)
}
\arguments{
  \item{mod}{a list containing all the attributes for the
  object}
}
\value{
  an object of type \code{frbs}
}
\description{
  This function creates objects of type \code{frbs}.
  Currently, its implementation is very basic and does no
  argument checking, as it is only used internally.
}
\details{
  The members of the \code{frbs} object depend on the used
  learning method. The following list describes all of the
  members that can be present. \describe{
  \item{range.input}{the range of the input data. Whether
  it is normalized to lie between 0 and 1 or not depends on
  the selected method.} \item{range.output}{the range of
  the output data. Whether it is normalized to lie between
  0 and 1 or not depends on the selected method.}
  \item{num.varinput}{the number of input variables.}
  \item{num.fvalinput}{the number of fuzzy terms of the
  input variables.} \item{names.varinput}{the generated
  names of fuzzy terms of the input variables.}
  \item{num.fvaloutput}{the number of fuzzy terms of the
  output variable.} \item{num.labels}{the number of fuzzy
  terms for the variables} \item{varout.mf}{a matrix to
  generate the shapes of the membership functions for the
  output variable.  The first row represents the shape of
  the membership functions, the other rows contain the
  parameters that have been generated.  Whether the values
  of parameters within the matrix are normalized to lie
  between 0 and 1 or not depends on the selected method.}
  \item{names.varoutput}{generated names of the output
  variable.} \item{rule}{the fuzzy IF-THEN rules; In the
  GFS.FR.MOGUL case, a rule refers to the parameter values
  of the membership function which represents the rule.}
  \item{rule.data.num}{the fuzzy IF-THEN rules in integer
  format.} \item{varinp.mf}{a matrix to generate the shapes
  of the membership functions for the input variables.  The
  first row represents the shape of the membership
  functions, the other rows contain the non NA values
  representing the parameters related with their type of
  membership function.  For example, trapezoid, triangular,
  and Gaussian have four, three, and two values as their
  parameters, respectively.  Whether the values of
  parameters within the matrix are normalized to lie
  between 0 and 1 or not depends on the selected method.}
  \item{type.model}{the model type. Here, 1 refers to
  Mamdani model, and 2 refers to Takagi Sugeno Kang model
  on the consequence part.} \item{func.tsk}{a matrix of the
  Takagi Sugeno Kang model consequent part of the fuzzy
  IF-THEN rules.} \item{type.defuz}{the type of the
  defuzzification method.} \item{type.tnorm}{the type of
  the t-norm method.} \item{type.snorm}{the type of the
  s-norm method.} \item{method.type}{the type of the
  selected method.} \item{name}{the name given to the
  model.} \item{range.data.ori}{range of the original data
  (before normalization).} \item{cls}{cluster centers.}
  \item{Dthr}{the boundary parameter of the DENFIS method.}
  \item{d}{the multiplier parameters of the DENFIS method.}
  \item{r.a}{the neighborhood factor of SBC.}
  \item{degree.rule}{certainty degree of rules.}
  \item{rule.data.num}{a matrix representing the rules in
  integer form.} \item{grade.cert}{grade of certainty for
  classification problems.} \item{alpha.heuristic}{a
  parameter for the heuristic of the HGD method.}
  \item{colnames.var}{the names of variables.} }
}

