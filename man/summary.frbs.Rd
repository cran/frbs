\name{summary.frbs}
\alias{summary.frbs}
\title{The summary function for frbs objects}
\usage{
  \method{summary}{frbs} (object, ...)
}
\arguments{
  \item{object}{the \code{\link{frbs-object}}}

  \item{...}{the other parameters (not used)}
}
\description{
  This function enables the output of a summary of the
  \code{\link{frbs-object}}.
}
\details{
  This function displays several components of the object.
  The components of one particular method can be different
  from components of other methods. The following is a
  description of all components which might be printed.
  \itemize{ \item The name of the model: A name given by
  the user representing the name of the simulation or data
  or model. \item Model was trained using: It shows which
  method we have been used. \item The interval of training
  data: It is a matrix representing the original interval
  of data where the first and second rows are minimum and
  maximum of data, respectively. The number of columns
  represents the number of variables. \item The number of
  fuzzy terms of the input variables: Given as elements of
  a matrix. \item The names of fuzzy terms of the input
  variables: These names are generated automatically by
  frbs expressing all fuzzy terms considered.  These names
  are built by two parts which are the name of variables
  expressed by "v" and the name of fuzzy labels of each
  variables represented by "a".  For example, "v.1-a.1"
  means the fuzzy label "a.1" of the first variable (v.1).
  \item The names of fuzzy terms of the output variable:
  For the Mamdani model, since the frbs package only
  considers single output, the names of the fuzzy terms for
  the output variable are simple and clear and start with
  "c". However, for Takagi Sugeno Kang model and fuzzy
  rule-based classification systems, this component is
  always NULL. \item The parameter values of membership
  functions of the input variables (normalized): It is
  represented by a matrix (5 x n) where n depends on the
  number of fuzzy terms on the input variables and the
  first row of the matrix describes a type of membership
  function, and the rest of rows are their parameter
  values. For example, label "v.1-a.2" has value {4.0,
  0.23, 0.43, 0.53, 0.73} on its column. It means that the
  label a.2 of variable v.1 has a parameter as follows.
  4.0 on the first row shows trapezoid shape in the middle
  position, while 0.23, 0.43, 0.53, and 0.73 are corner
  points of a trapezoid.  Furthermore, the following is the
  complete list of shapes of membership functions:
  \itemize{ \item Triangular: 1 on the first row and rows
  2, 3, and 4 represent corner points.  \item Trapezoid: 2,
  3, or 4 on the first row means they are trapezoid in
  left, right and middle side, respectively, and rows 2, 3,
  4, and 5 represent corner points. But for trapezoid at
  left or right side the fifth row is NA.  \item Gaussian:
  5 on the first row means it uses Gaussian and second and
  third row represent mean and variance.  \item Sigmoid: 6
  on the first row and two parameters (gamma and c) on
  second and third rows.  \item Generalized bell: 7 on the
  first row and three parameters (a, b, c) on second,
  third, and fourth rows.  } \item The fuzzy IF-THEN rules:
  In this package, there are several models for
  representing fuzzy IF-THEN rules based on the method
  used.  \itemize{ \item Mamdani model: they are
  represented as a knowledge base containing two parts:
  antecedent and consequent parts which are separated by a
  sign "->", as for example in the following rule:

  \code{var.1 is v.1-a.1 and var.2 is v.2-a.2 -> var.3 is
  c.2}

  \item Takagi Sugeno Kang model: In this model, this
  component only represents the antecedent of rules while
  the consequent part will be represented by linear
  equations.  \item fuzzy rule-based classification
  systems: This model is quite similar to the Mamdani
  model, but the consequent part expresses pre-defined
  classes instead of linguistic values.  \item approximate
  approach: Especially for GFS.FR.MOGUL, a matrix of
  parameters of membership functions is used to represent
  the fuzzy IF-THEN rules as well.  The representation of
  rules and membership functions is a matrix (n x (p x m))
  where n is the number of rules and m is the number of
  variables while p is the number of corner points of the
  membership function, if we are using triangular or
  trapezoid then p = 3 or 4, respectively.  For example,
  let us consider the triangular membership function and a
  number of variables of 3.  The representation of rules
  and membership functions is as follows: <<a11 a12 a13>>
  <<b11 b12 b13>> <<c11 c12 c13>>.

  } \item The linear equations on consequent parts of fuzzy
  IF-THEN rules: It is used in the Takagi Sugeno Kang
  model. \item The weight of the rules or the certainty
  factor: For the FRBCS.W method, this shows the weight
  related to the rules representing the ratio of dominance
  among the rules. \item The cluster centers: This
  component is used in clustering methods representing
  cluster centers. }
}

