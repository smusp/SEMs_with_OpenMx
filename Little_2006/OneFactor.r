
#### Methods of Scaling and Identification

## Some easier examples.
## Demonstrates three methods of scaling in a one-factor model:
## 1. Reference-Group Method - Constrain latent variable's variance and mean;
## 2. Marker-Variable Method - Constrain one loading and that indicator's intercept;
## 3. Effects-Scaling Method - Constrain sums of loadings and intercepts.

## Compare results with lavaan's results

## One Factor:
#  "Positive Affect" factor for 7th grade

# Data in Appendix A of:
# Little, T., Slegers, D., & Card, N. (2006). A non-arbitrary method of
# identifying and scaling latent variables in SEM and MACS models.
# Structural Equation Modeling, 13(1), 59-72.

## Load package
library(OpenMx)

## Get data
# Vectors of correlations (row-by-row), standard deviations, and means, and sample size.
vcor <- c(
   1.00000,
   0.75854,  1.00000,
   0.76214,  0.78705,  1.00000)

vmean <- c(3.13552, 2.99061, 3.06945)
vsd <- c(0.66770, 0.68506, 0.70672)
n <- 380

# Variable names
names <- c("pos1", "pos2", "pos3")

# Get full correlation matrix
mcor <- matrix( , 3, 3)                          # Empty matrix
mcor[upper.tri(mcor, diag = TRUE)] <- vcor       # Fill the upper triangle
mcor <- pmax(mcor, t(mcor), na.rm = TRUE)        # Fill the lower triangle

# Get co/variance matrix
mcov <- outer(vsd, vsd) * mcor

# Name the rows and columns
dimnames(mcov) <- list(names, names); mcov

names(vmean) <- names   # OpenMx requires the means be named

# Get data into OpenMx format
data <- mxData(observed = mcov, type = "cov", means = vmean, numObs = n)


### Method 1: Reference-Group Method
## Constrain latent variance to 1
## Constrain latent mean to 0
## See OneFactor1.svg in the `images` folder for the model diagram.

## Collect the bits and pieces needed by OpenMx
# Factor loadings
loadings <- mxPath(from = "POS", to = names, arrows = 1,
   free = TRUE, values = 0.5,
   labels = c("lambda1", "lambda2", "lambda3"))

# Factor variance - Constrain variance to 1
varFac <- mxPath(from = "POS", arrows = 2,
   free = FALSE, values = 1,
   labels = "phi")

# Factor mean - Constrain mean to 0
means <- mxPath(from = "one", to = "POS", arrows = 1,
   free = FALSE, values = 0,
   labels = "kappa")

# Residual variances
varRes <- mxPath(from = names, arrows = 2,
   free = TRUE, values = 1,
   labels = c("theta1", "theta2", "theta3"))

# Intercepts
intercepts <- mxPath(from = "one", to = names, arrows = 1,
   free = TRUE, values = 1,
   labels = c("tau1", "tau2", "tau3"))

## Setup the model
model1 <- mxModel("One Factor Model", type = "RAM",
   manifestVars = names, latentVars = "POS",
   data, loadings, varFac, means, varRes, intercepts)

## Run the model and get summary
fit1 <- mxRun(model1)
summary(fit1)

# These models are just-identified.
# Number of variables is 3;
# Therefore, number of pieces of information in co/variance matrix: (3 X 4) / 2 = 6
# plus 3 means = 9 pieces of information.
# Number of parameters:
#   3 loadings
#   3 redisual variances
#   3 intercepts
#   Total of 9 parameters

# Therefore, degrees of freedom is zero,
# chi square is 0,
# and other fit indices are either 1 or 0.
# There is a small discrepancy - chi sq is not quite 0.
# Note: the log likelihoods for the model and the saturated model differ.
# OpenMx needs to estimate reference (saturated and independence) models.

summary(fit1, refModels = mxRefModels(fit1, run = TRUE))
coef(fit1)

########################
## Check with lavaan
library(lavaan)

m1 <- "
  # Loadings
  POS =~ NA*lambda1*pos1 + lambda2*pos2 + lambda3*pos3

  # Latent variance - Constrained to 1
  POS ~~ 1*phi*POS

  # Latent mean - Constrained to 0
  POS ~ 0*kappa*1

  # Residual variances
  pos1 ~~ theta1*pos1
  pos2 ~~ theta2*pos2
  pos3 ~~ theta3*pos3

  # Intercepts 
  pos1 ~ tau1*1
  pos2 ~ tau2*1
  pos3 ~ tau3*1
"

m1_fit <- sem(m1, sample.cov = mcov, sample.nobs = n, sample.mean = vmean)
OpenMx <- coef(fit1)
lavaan <- coef(m1_fit)

# Get coefs in same order
lavaan <- lavaan[match(names(OpenMx), names(lavaan))]
cbind(OpenMx, lavaan)
########################


### Method 2: Marker-Variable Method
## Constrain third loading to 1
## Constrain third intercept to 0
## See OneFactor2.svg in the `images` folder for the model diagram.

## Collect the bits and pieces needed by OpenMx
# Factor loadings - Constrain 3rd loading to 1
loadings <- mxPath(from = "POS", to = names, arrows = 1,
   free = c(TRUE, TRUE, FALSE), values = c(0.5, 0.5, 1),
   labels = c("lambda1", "lambda2", "lambda3"))

# Factor variance
varFac <- mxPath(from = "POS", arrows = 2,
   free = TRUE, values = 1,
   labels = "phi")

# Factor mean
means <- mxPath(from = "one", to = "POS", arrows = 1,
   free = TRUE, values = 1,
   labels = "kappa")

# Residual variances
varRes <- mxPath(from = names, arrows = 2,
   free = TRUE, values = 1,
   labels = c("theta1", "theta2", "theta3"))

# Intercepts - Constrain 3rd intercept to 0
intercepts <- mxPath(from = "one", to = names, arrows = 1,
   free = c(TRUE, TRUE, FALSE), values = c(1, 1, 0),
   labels = c("tau1", "tau2", "tau3"))

## Setup the model
model2 <- mxModel("One Factor Model", type = "RAM",
   manifestVars = names, latentVars = "POS",
   data, loadings, varFac, means, varRes, intercepts)

## Run the model and get summary
fit2 <- mxRun(model2)
summary(fit2, refModels = mxRefModels(fit2, run = TRUE))
coef(fit2)

########################
## Check with lavaan
library(lavaan)

m2 <- "
  # Loadings - Constrain 3rd loading to 1
  POS =~ NA*lambda1*pos1 + lambda2*pos2 + 1*lambda3*pos3

  # Latent variance
  POS ~~ phi*POS

  # Latent mean
  POS ~ kappa*1

  # Residual variances
  pos1 ~~ theta1*pos1
  pos2 ~~ theta2*pos2
  pos3 ~~ theta3*pos3

  # Intercepts - Constrain 3rd intercept to 0
  pos1 ~ tau1*1
  pos2 ~ tau2*1
  pos3 ~ 0*tau3*1
"

m2_fit <- sem(m2, sample.cov = mcov, sample.nobs = n, sample.mean = vmean)
OpenMx <- coef(fit2)
lavaan <- coef(m2_fit)

# Get coefs in same order
lavaan <- lavaan[match(names(OpenMx), names(lavaan))]      
cbind(OpenMx, lavaan)
########################


### Method 3: Effects-Scaling Method
## See model diagram in OneFactor3.svg
## Constrain sum of loadings to equal number of loadings
## Constrain sum of intercepts to equal 0
## See OneFactor3.svg in the `images` folder for the model diagram.

## Collect the bits and pieces needed by OpenMx
# Factor loadings
loadings <- mxPath(from = "POS", to = names, arrows = 1,
   free = TRUE, values = 0.5,
   labels = c("lambda1", "lambda2", "lambda3"))

# Factor variance
varFac <- mxPath(from = "POS", arrows = 2,
   free = TRUE, values = 1,
   labels = "phi")

# Factor mean
means <- mxPath(from = "one", to = "POS", arrows = 1,
   free = TRUE, values = 1,
   labels = "kappa")

# Residual variances
varRes <- mxPath(from = names, arrows = 2,
   free = TRUE, values = 1,
   labels = c("theta1", "theta2", "theta3"))

# Intercepts
intercepts <- mxPath(from = "one", to = names, arrows = 1,
   free = TRUE, values = 1,
   labels = c("tau1", "tau2", "tau3"))

# Constraints
conLoad <- mxConstraint(lambda1 + lambda2 + lambda3 == 3)
conInter <- mxConstraint(tau1 + tau2 + tau3 == 0)

## Setup the model
model3 <- mxModel("One Factor Model", type = "RAM",
   manifestVars = names, latentVars = "POS",
   data, loadings, means, varFac, varRes, intercepts,
   conLoad, conInter)

## Run the model and get summary
fit3 <- mxRun(model3)
summary(fit3, refModels = mxRefModels(fit3, run = TRUE))
coef(fit3)

########################
## Check with lavaan
library(lavaan)

m3 <- "
  # Loadings
  POS =~ NA*lambda1*pos1 + lambda2*pos2 + lambda3*pos3

  # Latent variance
  POS ~~ phi*POS

  # Latent mean
  POS ~ kappa*1

  # Residual variances
  pos1 ~~ theta1*pos1
  pos2 ~~ theta2*pos2
  pos3 ~~ theta3*pos3

  # Intercepts 
  pos1 ~ tau1*1
  pos2 ~ tau2*1
  pos3 ~ tau3*1

  # Constraints
  lambda1 + lambda2 + lambda3 == 3
  tau1 + tau2 + tau3 == 0
"

m3_fit <- sem(m3, sample.cov = mcov, sample.nobs = n, sample.mean = vmean)
OpenMx <- coef(fit3)
lavaan <- coef(m3_fit)

# Get coefs in same order
lavaan <- lavaan[match(names(OpenMx), names(lavaan))]
cbind(OpenMx, lavaan)
########################
