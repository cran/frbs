library(frbs)

data(frbsData)

data.train <- frbsData$GasFurnance.dt[1 : 204, ]
data.fit <- data.train[, 1 : 2]
data.tst <- frbsData$GasFurnance.dt[205 : 292, 1 : 2]
real.val <- matrix(frbsData$GasFurnance.dt[205 : 292, 3], ncol = 1)

#####################
### 1. Input Process
range.data<-matrix(c(-2.716, 2.834, 45.6, 60.5, 45.6, 60.5), nrow=2)

method.type <- "ANFIS"
 
control <- list(num.labels = 5, max.iter = 100, step.size = 0.01, type.mf = 3, name = "GasFur")

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

#summary(res)
##########################
### 4a. Error Measurement
##########################
y.pred <- res.test
y.real <- real.val

bench <- cbind(y.pred, y.real)
colnames(bench) <- c("pred. val.", "real. val.")
print("Comparison ANFIS Vs Real Value on Gas Furnace Data Set")
print(bench)


residuals <- (y.real - y.pred)
MSE <- mean(residuals^2)
RMSE <- sqrt(mean(residuals^2))
SMAPE <- mean(abs(residuals)/(abs(y.real) + abs(y.pred))/2)*100

err <- c(MSE, RMSE, SMAPE)
names(err) <- c("MSE", "RMSE", "SMAPE")

print("Error Measurement: ")
print(err) 


##########################
### 4b. Visualisation and compare between simulation and real data
### In this step, we use the same data on training and testing phase
op <- par(mfrow = c(2, 1))
x1 <- seq(from = 1, to = nrow(res.fit))
result.fit <- cbind(data.train[, 3], res.fit)
plot(x1, result.fit[, 1], col="red", main = "Fitting phase: the training data(red) Vs Sim. result(blue)", type = "l", ylab = "CO2")
lines(x1, result.fit[, 2], col="blue")



###########################
### 4c. The same as before, this process is for visualisation and compare between real data and prediction on testing.

result.test <- cbind(real.val, res.test)
x2 <- seq(from = 1, to = nrow(result.test))
plot(x2, result.test[, 1], col="red", main = "Predicting phase: the Real Data(red) Vs Sim. result(blue)", type = "l", ylab = "CO2")
lines(x2, result.test[, 2], col="blue", type = "l")

par(op)

