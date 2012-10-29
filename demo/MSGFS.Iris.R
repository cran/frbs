library(frbs)

## Start input data
data(iris)
set.seed(2)
irisShuffled <- iris[sample(nrow(iris)),]
irisShuffled[,5] <- unclass(irisShuffled[,5])

tra.iris <- irisShuffled[1:105,]
tst.iris <- irisShuffled[106:nrow(irisShuffled),1:4]

real.iris <- matrix(irisShuffled[106:nrow(irisShuffled),5], ncol = 1)


range.data<-matrix(c(4.3, 7.9, 2.0, 4.4, 1.0, 6.9, 0.1, 2.5, 0.501, 3.499), nrow=2)
 
## generate model especially rule database by training process
method.type <- "MSGFS"

 
control <- list(popu.size = 80, persen_cross = 0.6, max.iter = 150, persen_mutant = 0.3, epsilon = 0.05, name="sim-0")


object <- frbs.learn(tra.iris, range.data, method.type, control)

## conduct the testing process
res <- predict(object, tst.iris)

y.pred <- round(res)

y.real <- real.iris
	
bench <- cbind(y.pred, y.real)

#### Measure error for classification
counter <- 0

for (i in 1 : nrow(bench)){
	if (bench[i, 1] != bench[i, 2]){
		counter = counter + 1
	}
}

err <- counter / nrow(bench) * 100
print("The result: ")
print(y.pred)

print("MSGFS: percentage Error on Iris: ")
print(err) 
