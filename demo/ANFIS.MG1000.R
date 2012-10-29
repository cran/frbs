library(frbs)

data(frbsData)
data.train <- frbsData$MackeyGlass1000.dt[1: 500, ]
data.fit <- data.train[, 1 : 4]
data.tst <- frbsData$MackeyGlass1000.dt[501 : 1000, 1 : 4]
real.val <- matrix(frbsData$MackeyGlass1000.dt[501 : 1000, 5], ncol = 1)

#####################
### 1. Input Process
range.data <- matrix(c(0.43462, 1.3105, 0.43462, 1.3105, 0.43462, 1.3105, 0.43462, 1.3105, 0.43462, 1.3105), nrow=2)

method.type <- "ANFIS"
 
control <- list(num.labels = 7, max.iter = 500, step.size = 0.01, name = "MG1000")

##########################
### 2. Training Process
### The main result on this process is rule database which is used on testing process.
object <- frbs.learn(data.train, range.data, method.type, control)

##########################
### 3a. This process is a part of fitting the model using data training. 
res.fit <- predict(object, data.fit)

##########################
### 3b. Testing process
res.test <- predict(object, data.tst)

##########################
### 4a. Error Measurement
##########################
y.pred <- res.test
y.real <- real.val

bench <- cbind(y.pred, y.real)
colnames(bench) <- c("pred. val.", "real. val.")
print("Comparison ANFIS Vs Real Value on Mackey Glass Data Set")
print(bench)


residuals <- (y.real - y.pred)
MSE <- mean(residuals^2)
RMSE <- sqrt(mean(residuals^2))
SMAPE <- mean(abs(residuals)/(abs(y.real) + abs(y.pred))/2)*100

err <- c(MSE, RMSE, SMAPE)
names(err) <- c("MSE", "RMSE", "SMAPE")

print("ANFIS: Error Measurement: ")
print(err) 


##########################
### 4b. Visualisation and compare between simulation and real data
### In this step, we use the same data on training and testing phase
op <- par(mfrow = c(2, 1))
x1 <- seq(from = 1, to = nrow(res.fit))
result.fit <- cbind(data.train[, 5], res.fit)
plot(x1, result.fit[, 1], col="red", main = "Fitting phase: the training data(red) Vs Sim. result(blue)", type = "l", ylab = "MG")
lines(x1, result.fit[, 2], col="blue")



###########################
### 4c. The same as before, this process is for visualisation and compare between real data and prediction on testing.

result.test <- cbind(real.val, res.test)
x2 <- seq(from = 1, to = nrow(result.test))
plot(x2, result.test[, 1], col="red", main = "Predicting phase: the Real Data(red) Vs Sim. result(blue)", type = "l", ylab = "MG")
lines(x2, result.test[, 2], col="blue", type = "l")

par(op)

