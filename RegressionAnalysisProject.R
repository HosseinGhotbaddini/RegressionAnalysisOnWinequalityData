# Regression Analysis Project — regression and classification on winequality-white dataset.

library(ggplot2) 
library(dplyr)
library(broom)
library(readxl)
library(readr)
library(varhandle)
library(ggthemes)
library(Rmisc)
library(class)
library(MLmetrics)
library(leaps)
library(glmnet)
library(ISLR)
library(gam)
library(corrgram)
library(rmdformats)
library(corrplot)
library(MASS)
library(olsrr)
library(factoextra)
library(car)
library(ggfortify)
library(plotmo)


set.seed(100)


setwd("/Users/QotboddiniH/Documents/R/")
df <- read.table("winequality-white.csv", sep = ";", header = T)
dim(df)
head(df)
summary(df)
sum(is.na(df))

par(mfcol = c(1,1))
corrplot(cor(df))
ggplot(data = df, mapping = aes(x = density, y = alcohol)) + geom_point(color = "red", alpha = 0.4)
ggplot(data = df, mapping = aes(x = density, y = residual.sugar)) + geom_point(color = "red", alpha = 0.1)
ggplot(data = df, mapping = aes(x = fixed.acidity, y = pH)) + geom_point(color = "red", alpha = 0.1)


par(mfcol = c(3,4))
for (i in 1:12) {
  plot(density(df[,i], cut = F), xlab = "value", main = colnames(df)[i])
}


pca <- prcomp(df[, -12], scale = TRUE)
fviz_eig(pca, , ncp = 11)

model <- lm(quality ~ . , data = df)
#ols_step_all_possible(model)
ols_step_backward_p(model)
ols_step_forward_p(model)
df <- df[, -c(3,5,7)]

#regression data
df.rg <- df

#classification data
df.cl <- df
df.cl$quality <- floor(df.cl$quality / 6)


model <- lm(quality ~ . , data = df.rg)
vif(model)
df.rg <- df.rg[, -5]
model <- lm(quality ~ . , data = df.rg)
vif(model)

model <- glm(quality ~ ., family = binomial, data = df.cl)
vif(model)
df.cl <- df.cl[, -5]
model <- glm(quality ~ ., family = binomial, data = df.cl)
vif(model)

# Remove Outliers
par(mfcol = c(2,4))
for (i in 1:7) {
  boxplot(df.rg[[i]])
  mtext(names(df.rg)[i], cex = 0.8, side = 1, line = 2)
}

outliers = c()
for (i in 1:7) {
  stats <- boxplot.stats(df.rg[[i]])$stats
  bottom_outlier_rows = which(df.rg[[i]] < stats[1])
  top_outlier_rows = which(df.rg[[i]] > stats[5])
  outliers <- c(outliers , top_outlier_rows[ !top_outlier_rows %in% outliers ] )
  outliers <- c(outliers , bottom_outlier_rows[ !bottom_outlier_rows %in% outliers ] )
}

#rg
model <- lm(quality ~ ., data = df.rg)
cooksd <- cooks.distance(model)
plot(cooksd, pch = "*", cex = 2, main = "Influential Obs by Cooks distance")
abline(h = 4*mean(cooksd, na.rm = T), col = "red")
coutliers <- as.numeric(rownames(df.rg[cooksd > 4 * mean(cooksd, na.rm=T), ]))
outliers <- c(outliers , coutliers[ !coutliers %in% outliers ] )

#cl
model <- glm(quality ~ ., family = binomial, data = df.cl)
cooksd <- cooks.distance(model)
plot(cooksd, pch = "*", cex = 2, main = "Influential Obs by Cooks distance")
abline(h = 4*mean(cooksd, na.rm = T), col = "red")
coutliers <- as.numeric(rownames(df.cl[cooksd > 4 * mean(cooksd, na.rm=T), ]))
outliers <- c(outliers , coutliers[ !coutliers %in% outliers ] )

df.rg <- df.rg[-outliers, ]
df.cl <- df.cl[-outliers, ]


par(mfcol = c(3,3))
for ( i in 1:8 ) {
  truehist(df.rg[[i]], xlab = names(df.rg)[i], col = 'lightgreen', main = paste("Average =", signif(mean(df.rg[[i]]),3)), nbins = 10)
}
truehist(df.cl[[8]], xlab = names(df.cl)[8], col = 'lightgreen', main = paste("Average =", signif(mean(df.cl[[8]]),3)), nbins = 10)
# now we are almost clean





# Regression
df.rg <- df.rg[,c(8,1,2,3,4,5,6,7)]

for (i in 2:7) {
  for (j in (i+1):8) {
    df.rg <- cbind(df.rg, df.rg[i]*df.rg[j])
    colnames(df.rg)[ncol(df.rg)] <- paste0(colnames(df.rg)[i], ":", colnames(df.rg)[j])
  }
}
for (i in 2:8) {
  df.rg <- cbind(df.rg, log(df.rg[i]))
  colnames(df.rg)[ncol(df.rg)] <- paste0("log(", colnames(df.rg)[i], ")")
  df.rg <- cbind(df.rg, df.rg[i]*df.rg[i])
  colnames(df.rg)[ncol(df.rg)] <- paste0("(", colnames(df.rg)[i], ")^2")
}
df.rg <- df.rg[is.finite(rowSums(df.rg)),]


size <- dim(df.rg)[1]
sample <- sample(x = (1:size), size = (size / 5), replace = FALSE)
#test <- df.rg[sample,]
#train <- df.rg[-sample,]

# Backward Selection
is.done <- F
while (!is.done) {
  lmf <- lm(quality ~ .-quality, data = df.rg[-sample,])
  pval <- data.frame(summary(lmf)$coefficients[,4])
  
  removable <- c()
  for (i in 2:nrow(pval)) {
    if (pval[i,] > 0.05)
      removable <- c(removable, i)
  }
  
  if (is.null(removable)) {
    is.done <- T
  } else {
    df.rg <- df.rg[-removable]
  }

  #print(summary(lmf))
}

#lmf <- lm(quality ~ .-quality, data = df.rg[-sample,])
#print(summary(lmf))
#print(vif(lmf))

# VIF Selection
is.done <- F
while (!is.done) {
  lmf <- lm(quality ~ .-quality, data = df.rg[-sample,])
  vif <- as.data.frame(vif(lmf))

  maxVif <- 1
  for (i in 2:nrow(vif)) {
    if (vif[i,] > vif[maxVif,])
      maxVif <- i
  }
  
  if (vif[maxVif,] < 10) {
    is.done <- T
  } else {
    df.rg <- df.rg[-(maxVif+1)]
  }
  
  #print(maxVif)
  #print(vif)
}

lmf <- lm(quality ~ .-quality, data = df.rg[-sample,])
print(summary(lmf))
print(vif(lmf))




# Lasso
x = model.matrix(quality~ ., df.rg[-sample,])[, -1]
y = df.rg[-sample, "quality"]

lassoRes = glmnet(scale(x), y, alpha = 1)
par(mfcol = c(1,1))
plot_glmnet(lassoRes)


cvLassoRes = cv.glmnet(scale(x), y, alpha = 1, lambda = 10^((-200:20)/80))
plot(cvLassoRes)
cvLassoRes$lambda.min
cvLassoRes$lambda.1se

predict(lassoRes, type = "coefficients", s = cvLassoRes$lambda.1se)


# Ridge
ridgeResScaled = glmnet(scale(x), y, alpha = 0)
plot_glmnet(ridgeResScaled)


cvRidgeResScaled = cv.glmnet(scale(x), y, alpha = 0, lambda = 10^((-50:60)/20))
plot(cvRidgeResScaled)
cvRidgeResScaled$lambda.min
cvRidgeResScaled$lambda.1se

predict(ridgeResScaled, type = "coefficients", s = cvRidgeResScaled$lambda.1se)

finalModel <- glmnet(scale(x), y, alpha = 0, lambda = cvRidgeResScaled$lambda.1se)

test = df.rg[sample,]
pred <- predict(finalModel, scale(as.matrix(test[,-1])),  type="response")
mean((pred - test$quality)*(pred - test$quality))

# Compare with Naive Model
naive = mean(df.rg[-sample,]$quality)
naivePred = seq(naive, naive, length.out = dim(test)[1])
mean((test$quality - naivePred)*(test$quality - naivePred))





# Classification
df.cl <- df.cl[,c(8,1,2,3,4,5,6,7)]

size <- dim(df.cl)[1]
sample <- sample(x = (1:size), size = (size / 5), replace = FALSE)

test <- df.cl[sample,]
train <- df.cl[-sample,]
histogram(train$quality)

# Acurracy is a good criterion for this data

# Logistic
logistic <- glm(quality ~ ., family = binomial, data = train)
logistic.func <- function(th) {
  logistic.pred <- ifelse(predict(logistic, test[,-1], type="response") > th, 1, 0)
  return (mean(logistic.pred == test$quality))
}
#table(logistic.pred, test$admit)

thvals = seq(from = 0.1, to = 0.9, by = 0.1)
acurracies.lg = sapply(X = thvals, FUN = logistic.func)
acurracies.lg

# bishtarin deqat be ezaaye:
max(acurracies.lg)
for (i in 1:9) {
  if(acurracies.lg[i] == max(acurracies.lg)) {
    print(i / 10)
    bestth <- (i / 10)
  }
}

# Acurracy
acurracies.lg[bestth * 10]





# KNN
knn.func <- function(kval) {
  knn.pred <- knn(train[,-1], test[,-1], train$quality, k = kval)
  return (mean(knn.pred == test$quality))
}
#table(knn.pred, test$admit)

kvals = seq(from = 1, to = 15, by = 1)
acurracies = sapply(X = kvals, FUN = knn.func)
acurracies

# bishtarin deqat be ezaaye:
max(acurracies)
for (i in 1:15) {
  if(acurracies[i] == max(acurracies)) {
    print(i)
    bestk <- i
  }
}

knn.pred <- knn(train[,-1], test[,-1], train$quality, k = bestk)
table(knn.pred, test$quality)

# Acurracy
acurracies[bestk]




# LDA
lda.fit = lda(quality~., data = train)
plot(lda.fit)

lda.pred = predict(lda.fit, test[,-1])

predictions = lda.pred$class
accuracy = mean(test$quality == predictions)

# Acurracy
accuracy




# QDA
qda.fit = qda(quality~.,data = train)
qda.pred = predict(qda.fit, test[,-1])

predictions = qda.pred$class
accuracy = mean(test$quality == predictions)

# Acurracy
accuracy




# Compare with Naive Model
naivePred = floor(runif(dim(test)[1], min = 0, max = 2))
# Acurracy
mean(test$quality == naivePred)


# So the best model is Logistic
pca <- prcomp(test[,-1], scale = TRUE, center = TRUE)

par(mfcol = c(1,2))
autoplot(pca, data = test, colour = "quality")

logistic <- glm(quality ~ ., family = binomial, data = train)
logistic.pred <- ifelse(predict(logistic, test[,-1], type="response") > bestth, 1, 0)
test[,1] = logistic.pred
autoplot(pca, data = test, colour = "quality")



