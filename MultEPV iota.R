# Install Rtools version 3.3 first. # In the wizard, tick the change path box.
# install_version("VGAM", "0.9-8") # New versions are not available for this version of R.
library(VGAM)
library(mlogit)
library(glmnet)
library(remotes)

sim.path <- "path/to/simulation/directory/"
source(paste(sim.path, "Coerce.R", sep = "")) # github.com/VMTdeJong/Multinomial-Predictive-Performance/tree/master/Simulation%20Study
source(paste(sim.path, "Polcal.R", sep = ""))
source(paste(sim.path, "PDI.R"   , sep = ""))

source("iota functions.R")
source("SampleIota.R")

#################################################################################################################################
################################################### The data ####################################################################
#################################################################################################################################
iota.full <- read.table("path/case/study/data", header = T)

R <- 6 # Number of predictors
J <- 3 # Number of outcome categories

# Sizes of smallest category
n3  <- 3  * (J - 1) * R # for EPVm = 3
n5  <- 5  * (J - 1) * R # for EPVm = 5
n10 <- 10 * (J - 1) * R # for EPVm = 10

# Sample 
iota.split <- SampleIota(iota.full, n10)
iota.val   <- iota.split$ex
iota.10    <- iota.split$inc
iota.5     <- SampleIota(iota.10,   n5)$inc
iota.3     <- SampleIota(iota.5,    n3)$inc

# Relative outcome frequencies are the same for all development and validation sets:
round(table(iota.full$outcome) / table(iota.val$outcome), 2)
round(table(iota.full$outcome) / table(iota.10$outcome) , 2)
round(table(iota.full$outcome) / table(iota.5$outcome)  , 2)
round(table(iota.full$outcome) / table(iota.3$outcome)  , 2)

# EPVm is fixed at 10, 5 and 3.
min(table(iota.10$outcome)) / (J - 1) / R
min(table(iota.5$outcome))  / (J - 1) / R
min(table(iota.3$outcome))  / (J - 1) / R

dev.samples <- list("epv3" = iota.3, "epv5" = iota.5, "epv10" = iota.10)
print(lapply(dev.samples, nrow))
print(nrow(iota.val))
#################################################################################################################################
################################################## Evaluation ###################################################################
#################################################################################################################################
for (s in seq_len(length(dev.samples))) {
  dev <- dev.samples[[s]]
  
  # Some necessary data preparation
  # dev.ml <- mlogit.data(dev, choice = "outcome", shape = "wide", alt.var = "out.type", chid.var = "id")  # original code, not compatible with new R version
  dev.ml <- mlogit.data(dev, choice = "outcome", shape = "wide", chid.var = "id")  # Possible fix?
  
  dev.y <- as.matrix(dev[["outcome"]], ncol = 1)
  dev.x <- as.matrix(dev[c("age", "soldmax", "papflow", "wallreg", "shadows", "ascites")])
  
  # Run the cross-validation procedure, to obtain optimal lambdas for lasso and ridge.
  cv.lasso <- cv.glmnet(dev.x, dev.y, alpha = 1, nfolds = 10, family = "multinomial") 
  cv.ridge <- cv.glmnet(dev.x, dev.y, alpha = 0, nfolds = 10, family = "multinomial") 
  
  which(cv.lasso$lambda.min == cv.lasso$lambda)/length(cv.lasso$lambda) # (One of) the last lambdas is selected
  which(cv.ridge$lambda.min == cv.ridge$lambda)/length(cv.ridge$lambda) # Thus, we need to extend the vector and run again
  
  lasso.lambdas <- GLMnetLambdas(cv.lasso$lambda, 100)
  ridge.lambdas <- GLMnetLambdas(cv.ridge$lambda, 100)
  
  cv.lasso <- cv.glmnet(dev.x, dev.y, alpha = 1, nfolds = 10, family = "multinomial", lambda = lasso.lambdas) 
  cv.ridge <- cv.glmnet(dev.x, dev.y, alpha = 0, nfolds = 10, family = "multinomial", lambda = ridge.lambdas) 
  
  which(cv.lasso$lambda.min == cv.lasso$lambda)/length(cv.lasso$lambda) # Now none of the lambdas on the edges are chosen.
  which(cv.ridge$lambda.min == cv.ridge$lambda)/length(cv.ridge$lambda)
  
  # Rewrite the penalized multinomial logit into the reference category multinomial logit
  lasso.b <- GLMnetAsRefCat(GLMnetExtractCoefs(cv.lasso))[-1, ]
  ridge.b <- GLMnetAsRefCat(GLMnetExtractCoefs(cv.ridge))[-1, ]
  
  # And run ml
  (summary(ml <- mlogit(outcome ~ 0 | age + soldmax + papflow + wallreg + shadows + ascites, data = dev.ml)))
  ml.b <- data.frame(matrix(coef(ml), nrow = 2))
  
  ### Quantify coefficients and shrinkage
  out.names <- c("Borderline", "Invasive")
  colnames(ml.b) <- colnames(lasso.b) <- colnames(ridge.b) <- c("intercept", colnames(dev.x))
  rownames(ml.b)    <- paste("ML", out.names, sep = " ")
  rownames(lasso.b) <- paste("Lasso", out.names, sep = " ")
  rownames(ridge.b) <- paste("Ridge", out.names, sep = " ")
  
  dev.coefs.temp <- cbind(data.frame(round(t(ml.b), digits = 2)),
                          CombIndStrings(round(t(lasso.b), digits = 2), data.frame(round(t((ml.b - lasso.b) / ml.b * 100), digits = 0)), 
                                         sep = " (", append.by = "%)"),
                          CombIndStrings(round(t(ridge.b), digits = 2), data.frame(round(t((ml.b - ridge.b) / ml.b * 100), digits = 0)), 
                                         sep = " (", append.by = "%)")
  )
  if (s <= 1) dev.coefs <- dev.coefs.temp else dev.coefs <- rbind(dev.coefs, dev.coefs.temp)
  print(dev.coefs)
  # write.csv(dev.coefs, file = "results/iota_coefs.csv")
  
  ### Predict new cases
  val.y  <- as.matrix(iota.val[["outcome"]], ncol = 1)
  val.x  <- as.matrix(iota.val[c("age", "soldmax", "papflow", "wallreg", "shadows", "ascites")])
  # val.ml <- mlogit.data(iota.val, choice = "outcome", shape = "wide", alt.var = "out.type", chid.var = "id") # original code, not compatible with new R version
  val.ml <- mlogit.data(iota.val, choice = "outcome", shape = "wide", chid.var = "id") # Possible fix?
  
  p.ml    <- predict(ml, newdata = val.ml, type = "probs")
  p.lasso <- predict(cv.lasso, newx = val.x, type = "response", s = cv.lasso$lambda.min)[ , , 1]
  p.ridge <- predict(cv.ridge, newx = val.x, type = "response", s = cv.ridge$lambda.min)[ , , 1]
  
  # Recode outcome into disease indicators.
  val.y.i <- matrix(0, nrow = nrow(p.ml), ncol = ncol(p.ml))
  for (i in 1:ncol(p.ml)) val.y.i[val.y == colnames(p.ml)[i], i] <- 1
  
  ### Quantify performance
  perf.ml <- Performance(val.y, val.y.i, p.ml)
  row.names(perf.ml) <- "ML"
  
  perf.lasso <- Performance(val.y, val.y.i, p.lasso)
  row.names(perf.lasso) <- "Lasso"
  
  perf.ridge <- Performance(val.y, val.y.i, p.ridge)
  row.names(perf.ridge) <- "Ridge"
  
  perf.temp <- rbind(perf.ml, perf.lasso, perf.ridge)
  if (s <= 1) perf <- perf.temp else perf <- rbind(perf, perf.temp)
  
  # write.csv(perf, file = "results/iota_perf.csv" )
  
  print(perf)
}
