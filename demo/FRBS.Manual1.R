 ## This example shows how to use frbs without 
 ## learning process.

 ## Define shape and parameters of membership functions of input variables.
 ## Please see fuzzifier function to contruct the matrix.
 varinp.mf <- matrix(c(2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA,
                       2, 0, 35, 75, NA, 3, 35, 75, 100, NA,
                       2, 0, 20, 40, NA, 1, 20, 50, 80, NA, 3, 60, 80, 100, NA,
                       2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA),
                       nrow = 5, byrow = FALSE)

 ## Define number of fuzzy terms of input variables.
 ## Suppose, we have 3, 2, 3, and 3 number of fuzzy terms 
 ## for first, second, third and fourth variables, respectively.
 num.fvalinput <- matrix(c(3, 2, 3, 3), nrow=1)
 
 ## Give the names of the fuzzy terms of each input variable.
 ## It should be noted that the names of the fuzzy terms must be unique,
 ## so user might put some number for making to be unique.
 varinput.1 <- c("a1", "a2", "a3")
 varinput.2 <- c("b1", "b2")
 varinput.3 <- c("c1", "c2", "c3")
 varinput.4 <- c("d1", "d2", "d3")
 names.varinput <- c(varinput.1, varinput.2, varinput.3, varinput.4)

 ## Set interval of data.
 range.input <- matrix(c(0,100, 0, 100, 0, 100, 0, 100), nrow=2)
 range.output <- matrix(c(0,100), nrow=2)

 ## Define number of fuzzy terms of output variable.
 ## In this case, we set the number of fuzzy terms is 3.
 num.fvaloutput <- matrix(c(3), nrow=1)

 ## Give the names of the fuzzy terms of the output variable.
 ## Note: the names of the fuzzy terms must be unique.
 varoutput.1 <- c("e1", "e2", "e3")
 names.varoutput <- c(varoutput.1)

 ## Define the shapes and parameters of the membership functions of the output variables.
 varout.mf <- matrix(c(2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA),
                       nrow = 5, byrow = FALSE)

 ## Set type of model which is 1 or 2 for Mamdani or Takagi Sugeno Kang model, respectively.
 ## In this case, we choose Mamdani model.
 type.model <- 1
 ## Set weighted average method to be used as defuzzification method.
 type.defuz <- 1
 ## We are using standard t-norm and s-norm.
 type.tnorm <- 1
 type.snorm <- 1

 ## Since we don't generate fuzzy model by learning from data, 
 ## we have to set Wang and Mendel's technique as type of method.
 method.type <- "WM"
 ## Give the name of simulation.
 name <- "Sim-0"
 ## Define fuzzy IF-THEN rules, 
 ## there are two kind of model: Mamdani and Takagi Sugeno Kang model
 ## if we use Mamdani model then the consequent part is in linguistic term,
 ## but if we set Takagi Sugeno Kang then we build a matrix representing 
 ## linear equations at consequent part.
 ## In this example we are using Mamdani model 
 ## (It should refer to the previous value of type.model parameter). 
 ## Please make sure that each rule has "->" sign. 
 rule <- matrix(c("a1","and","b1","and","c1","and","d1","->","e1",
                  "a2","and","b2","and","c2","and","d2", "->", "e2", 
                  "a3","and","b2","and","c2","and","d1", "->", "e3"), 
                  nrow=3, byrow=TRUE) 

 ## Define function of TSK if we use it or 
 ## set NULL if we use the Mamdani model.
 func.tsk<-matrix(c(1, 1, 5, 2, 1, 3, 1, 0.5, 0.1, 2, 1, 3, 2, 2, 2), nrow=3, byrow=TRUE)
 ## Provide the new data as testing. 
 newdata<- matrix(c(25, 40, 35, 15, 45, 75, 78, 70), nrow= 2, byrow = TRUE)
 ## the names of variables
 colnames.var <- c("input1", "input2", "input3", "input4", "output1")
 ######################
 ## 1. The following codes show that we generate fuzzy model using frbs.gen function
 ######################
 ## Generate fuzzy model by frbs.gen.
 object <- frbs.gen(range.input, range.output, num.fvalinput, names.varinput, 
                 num.fvaloutput, varout.mf, names.varoutput, rule, varinp.mf,
                 type.model, type.defuz, type.tnorm, type.snorm, func.tsk, 
                 colnames.var, method.type, name)
 
 ## We can display the membership function in graphical mode
 plotMF(object)

 ## Predicting phase using new data.
 res <- predict(object, newdata)