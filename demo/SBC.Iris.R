library(frbs)

data(iris)

irisShuffled <- iris[sample(nrow(iris)),]
irisShuffled[,5] <- unclass(irisShuffled[,5])

tra.iris <- irisShuffled[1:105,]
tst.iris <- irisShuffled[106:nrow(irisShuffled),1:4]

real.iris <- matrix(irisShuffled[106:nrow(irisShuffled),5], ncol = 1)

range.data<-matrix(c(4.3, 7.9, 2.0, 4.4, 1.0, 6.9, 0.1, 2.5, 0.501, 3.499), nrow=2)

##########################
### 2. Training Process
method.type <- "SBC"

control <- list(r.a = 0.5, eps.high = 0.5, eps.low = 0.15, name = "iris")

object <- frbs.learn(tra.iris, range.data, method.type, control)


res.test <- predict(object, tst.iris)

##########################
### Error Measurement
##########################
y.pred <- round(res.test)

y.real <- real.iris

bench <- cbind(y.pred, y.real)

#### Measure error for classification
counter1 <- 0

for (i in 1 : nrow(bench)){
	if (bench[i, 1] != bench[i, 2]){
		counter1 = counter1 + 1
	}
}

err <- counter1 / nrow(bench) * 100
print("The result: ")
print(y.pred)

print("SBC: percentage Error on Iris: ")
print(err) 


