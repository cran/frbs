### R code from vignette source 'lala2015pmml.Rtex'
### Encoding: UTF-8

###################################################
### code chunk number 1: lala2015pmml.Rtex:72-79
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
r = getOption("repos") # hard code the UK repo for CRAN
r["CRAN"] = "http://cran.uk.r-project.org"
options(repos = r)
rm(r)
set.seed(2)
library(frbs)


###################################################
### code chunk number 2: lala2015pmml.Rtex:860-869
###################################################
fun <- function(input.xy){
  z <- 1 / (input.xy[1]^4 + input.xy[2]^4
  -2 * input.xy[1]^2 - 2 * input.xy[2]^2 + 3)
}
input.xy <- expand.grid(seq(-2, 2, 0.14),
  seq(-2, 2, by = 0.14))
z <- apply(input.xy, 1, fun)
data <- cbind(input.xy, z)
colnames(data)<- c("X", "Y", "Z")


###################################################
### code chunk number 3: lala2015pmml.Rtex:872-879
###################################################
data <- data[sample(nrow(data)), ]
cut.indx <- round(0.8 * nrow(data))
data.tra <- data[1 : cut.indx, ]
data.tst <- data[(cut.indx + 1) : nrow(data), 
  1 : 2]
real.val <- data[(cut.indx + 1) : nrow(data), 
  3, drop = FALSE]


###################################################
### code chunk number 4: lala2015pmml.Rtex:882-883
###################################################
range.data <- apply(data, 2, range)


###################################################
### code chunk number 5: lala2015pmml.Rtex:887-890
###################################################
method.type <- "WM"
control <- list(num.labels = 5, type.mf = "GAUSSIAN", type.defuz = "WAM", 
 type.tnorm = "MIN", type.implication.func = "LUKASIEWICZ", name = "fourhill") 


###################################################
### code chunk number 6: lala2015pmml.Rtex:893-894
###################################################
mod.reg <- frbs.learn(data.tra, range.data, method.type, control)


###################################################
### code chunk number 7: lala2015pmml.Rtex:902-903
###################################################
write.frbsPMML(frbsPMML(mod.reg), "modRegress")


###################################################
### code chunk number 8: lala2015pmml.Rtex:906-907
###################################################
frbsPMML(mod.reg)


###################################################
### code chunk number 9: lala2015pmml.Rtex:1012-1014
###################################################
objReg <- read.frbsPMML("modRegress.frbsPMML")
res.test <- predict(objReg, data.tst)


###################################################
### code chunk number 10: lala2015pmml.Rtex:1017-1019
###################################################
err.MSE <- mean((real.val - res.test)^2)
print(err.MSE) 


###################################################
### code chunk number 11: lala2015pmml.Rtex:1032-1033
###################################################
data(iris)


###################################################
### code chunk number 12: lala2015pmml.Rtex:1036-1038
###################################################
set.seed(2)
irisShuffled <- iris[sample(nrow(iris)), ]


###################################################
### code chunk number 13: lala2015pmml.Rtex:1041-1048
###################################################
irisShuffled[,5] <- unclass(
   irisShuffled[, 5])
tra.iris <- irisShuffled[1 : 105, ]
tst.iris <- irisShuffled[106 : 
   nrow(irisShuffled), 1:4]
real.iris <- matrix(irisShuffled
  [106 : nrow(irisShuffled), 5], ncol = 1)


###################################################
### code chunk number 14: lala2015pmml.Rtex:1051-1053
###################################################
range.data.input <- apply(iris[,-ncol(iris)],
  2, range)


###################################################
### code chunk number 15: lala2015pmml.Rtex:1058-1061
###################################################
method.type <- "GFS.GCCL"
control <- list(popu.size = 30, num.class = 3, num.labels = 3, 
 persen_cross = 0.9, max.gen = 200, persen_mutant = 0.3, name="sim-Iris") 


###################################################
### code chunk number 16: lala2015pmml.Rtex:1064-1065
###################################################
mod.class <- frbs.learn(tra.iris, range.data.input, method.type, control)


###################################################
### code chunk number 17: lala2015pmml.Rtex:1072-1073
###################################################
write.frbsPMML(frbsPMML(mod.class), "modClass")


###################################################
### code chunk number 18: lala2015pmml.Rtex:1076-1077
###################################################
frbsPMML(mod.class)


###################################################
### code chunk number 19: lala2015pmml.Rtex:1173-1175
###################################################
objectClass <- read.frbsPMML("modClass.frbsPMML")
res.test <- predict(objectClass, tst.iris)


###################################################
### code chunk number 20: lala2015pmml.Rtex:1178-1180
###################################################
err = 100 * sum(real.iris != res.test)/nrow(real.iris)
print(err) 


