
#### Methods of Scaling and Identification


## Demonstrates three methods of scaling in one-factor, two-group model:
## 1. Reference-Group Method - Constrain latent variable's variance and mean;
## 2. Marker-Variable Method - Constrain one loading and that indicator's intercept;
## 3. Effects-Scaling Method - Constrain sums of loadings and intercepts.


## Following Little et al's lead, assume strong metric invariance:
## corresponding loadings and intercepts constrained to equality across groups


## Compare results from OpenMx with lavaan's results


## One Factor, Two Groups
#  "Positive Affect" factor for 7th and 8th grades
#  Data in Appendix A of:
#
# Little, T., Slegers, D., & Card, N. (2006). A non-arbitrary method of
# identifying and scaling latent variables in SEM and MACS models.
# Structural Equation Modeling, 13(1), 59-72.


## Load package
library(OpenMx)


## Get data 
# Vectors of correlations (row-by-row), standard deviations, and means, and sample size.
# 7th grade
vcor7 <- c(
   1.00000,
   0.75854,  1.00000,
   0.76214,  0.78705,  1.00000)

vmean7 <- c(3.13552, 2.99061, 3.06945)
vsd7 <- c(0.66770, 0.68506, 0.70672)
n7 <- 380

# 8th grade
vcor8 <- c(
   1.00000,
   0.81366,  1.00000,
   0.84980,  0.83523,  1.00000)

vmean8 <- c(3.07338, 2.84716, 2.97882)
vsd8 <- c(0.70299, 0.71780, 0.76208)
n8 <- 379

# Variable names
names <- c("pos1", "pos2", "pos3")


# Get full correlation matrix for each Grade
mcor7 = matrix( , 3, 3)                            # Empty matrix
mcor7[upper.tri(mcor7, diag = TRUE)] <- vcor7      # Fill the upper triangle
mcor7 = pmax(mcor7, t(mcor7), na.rm = TRUE)        # Fill the lower triangle

mcor8 = matrix( , 3, 3)                            # Empty matrix
mcor8[upper.tri(mcor8, diag = TRUE)] <- vcor8      # Fill the upper triangle
mcor8 = pmax(mcor8, t(mcor8), na.rm = TRUE)        # Fill the lower triangle

# Get (co)variance matrix
mcov7 <- outer(vsd7, vsd7) * mcor7
mcov8 <- outer(vsd8, vsd8) * mcor8

# Name the rows and columns
dimnames(mcov7) <- list(names, names)
dimnames(mcov8) <- list(names, names)
mcov7; mcov8

names(vmean7) = names   # OpenMx requires the means be named
names(vmean8) = names

# Put data into lists - used in lavaan analysis
mcov = list("Grade 7" = mcov7, "Grade 8" = mcov8)
vmean = list(vmean7, vmean8)
n = list(n7, n8)

# Get data into OpenMx format
data7 <- mxData(observed = mcov7, type = "cov", means = vmean7, numObs = n7)
data8 <- mxData(observed = mcov8, type = "cov", means = vmean8, numObs = n8)



### Method 1 
## See model diagram in OneFactorTwoGroups1.svg
#  Labels are Little el al's Greek labels

## Constrain latent variance to 1
## Constrain latent mean to 0
## These constraints apply to Grade 7 only.


## Collect the bits and pieces needed by OpenMx


# Factor loadings
loadings <- mxPath(from = "POS", to = names, arrows = 1,
   free = TRUE, values = 0.5,
   labels = c("lambda1", "lambda2", "lambda3"))

# Factor variances - Constrain Grade 7 variance to 1
varFac7 <- mxPath(from = "POS", arrows = 2,
   free = FALSE, values = 1, labels = "phi7")

varFac8 <- mxPath(from = "POS", arrows = 2,
   free = TRUE, values = 1, labels = "phi8")

# Factor mean - Constrain Grade 7 mean to 0
means7 <- mxPath(from = "one", to = "POS", arrows = 1,
   free = FALSE, values = 0, labels = "kappa7")

means8 <- mxPath(from = "one", to = "POS", arrows = 1,
   free = TRUE, values = 1, labels = "kappa8") 

# Residual variances
varRes7 <- mxPath(from = names, arrows = 2,
   free = TRUE, values = 1,
   labels = c("theta71", "theta72", "theta73"))

varRes8 <- mxPath(from = names, arrows = 2,
   free = TRUE, values = 1, 
   labels = c("theta81", "theta82", "theta83"))

# Intercepts
intercepts <- mxPath(from = "one", to = names, arrows = 1,
   free = TRUE, values = 1,
   labels = c("tau1", "tau2", "tau3"))


## Setup models for each Grade
modGr7 <- mxModel("Grade7", type = "RAM",
   manifestVars = names, latentVars = "POS",
   data7, loadings, varFac7, means7, varRes7, intercepts)

modGr8 <- mxModel("Grade8", type = "RAM",
   manifestVars = names, latentVars = "POS",
   data8, loadings, varFac8, means8, varRes8, intercepts)


## Combine the two models
fun <- mxFitFunctionMultigroup(c("Grade7.fitfunction", "Grade8.fitfunction"))
model1 <- mxModel("One Factor Two Group Model", modGr7, modGr8, fun)


## Run the model and get summary
fit1 <- mxRun(model1)
summary(fit1)

# Number of variables: 3
# Number of pices of information in (co)variance matrix: (3 X 4) / 2 = 6
# plus 3 means = 9 pieces of information for each group;
# that is, 18 for the model

# Number of parameters:
#   3 loadings (constrained to equality across groups)
#   1 latent mean
#   1 latent variance
#   6 residual variances (3 per group)
#   3 intercepts (constrained to equality across groups)
#   Total of 14 parameters

# Therefore, degrees of freedom = 4,
# and chi sq and fit indices are calculated.
# But make sure OpenMx estimates the reference models
# upon which to base calculations for chi sq and fit indices.

summary(fit1, refModels = mxRefModels(fit1, run = TRUE))
coef(fit1)


########################
## Check with lavaan
library(lavaan)

m1 <- "
  # Loadings
  POS =~ c(NA,NA)*c(lambda1, lambda1)*pos1 + c(lambda2, lambda2)*pos2 + c(lambda3, lambda3)*pos3

  # Latent variance - Consteain Grade 7 variance to 1
  POS ~~ c(1, NA)*c(phi7, phi8)*POS

  # Latent means - Constrain Grade 7 mean to 0
  POS ~ c(0, NA)*c(kappa7, kappa8)*1

  # Residual variances
  pos1 ~~ c(theta71, theta81)*pos1
  pos2 ~~ c(theta72, theta82)*pos2
  pos3 ~~ c(theta73, theta83)*pos3
  
  # Intercepts 
  pos1 ~ c(tau1, tau1)*1
  pos2 ~ c(tau2, tau2)*1
  pos3 ~ c(tau3, tau3)*1
"

m1_fit = sem(m1, sample.cov = mcov, sample.nobs = n, sample.mean = vmean)
summary(m1_fit)
OpenMx = coef(fit1)
lavaan = coef(m1_fit)

# Get coefs in same order
lavaan <- lavaan[match(names(OpenMx), names(lavaan))]
t(rbind(OpenMx, lavaan))
########################



### Method 2
## See model diagram in OneFactorTwoGroups2.svg
## Constrain third loading in POS to 1
## Constrain third intercept to 0
## With strong measurement invariance, 
## these constraints apply to both 

# Factor loadings - Constrain 3rd loading to 1
loadings <- mxPath(from = "POS", to = names, arrows = 1,
   free = c(TRUE, TRUE, FALSE), values = c(0.5, 0.5, 1),
   labels = c("lambda1", "lambda2", "lambda3"))

# Factor variances
varFac7 <- mxPath(from = "POS", arrows = 2,
   free = TRUE, values = 1, labels = "phi7")

varFac8 <- mxPath(from = "POS", arrows = 2,
   free = TRUE, values = 1, labels = "phi8")

# Factor means
means7 <- mxPath(from = "one", to = "POS", arrows = 1,
   free = TRUE, values = 1, labels = "kappa7")

means8 <- mxPath(from = "one", to = "POS", arrows = 1,
   free = TRUE, values = 1, labels = "kappa8")
   
# Residual variances
varRes7 <- mxPath(from = names, arrows = 2,
   free = TRUE, values = 1, 
   labels = c("theta71", "theta72", "theta73"))

varRes8 <- mxPath(from = names, arrows = 2,
   free = TRUE, values = 1, 
   labels = c("theta81", "theta82", "theta83")) 

# Intercepts - Constrain 3rd intercept to 0
intercepts <- mxPath(from = "one", to = names, arrows = 1,
   free = c(TRUE, TRUE, FALSE), values = c(1, 1, 0),
   labels = c("tau1", "tau2", "tau3"))


## Setup models for each Grade
modGr7 <- mxModel("Grade7", type = "RAM",
   manifestVars = names, latentVars = "POS",
   data7, loadings, varFac7, means7, varRes7, intercepts)

modGr8 <- mxModel("Grade8", type = "RAM",
   manifestVars = names, latentVars = "POS",
   data8, loadings, varFac8, means8, varRes8, intercepts)


## Combine the two models using "mxFitFunctionMultigroup()"
fun <- mxFitFunctionMultigroup(c("Grade7.fitfunction", "Grade8.fitfunction"))
model2 <- mxModel("One Factor Two Group Model", modGr7, modGr8, fun)


## Run the model and get summary
fit2 <- mxRun(model2)
summary(fit2, refModels = mxRefModels(fit2, run = TRUE))
coef(fit2)


########################
## Check with lavaan
library(lavaan)

m2 <- "
  # Loadings - Constrain 3rd loading to 1 in both Grades
  POS =~ c(NA,NA)*c(lambda1, lambda1)*pos1 + c(lambda2, lambda2)*pos2 + c(1,1)*c(lambda3, lambda3)*pos3

  # Latent variances
  POS ~~ c(phi7, phi8)*POS

  # Latent means
  POS ~ c(kappa7, kappa8)*1

  # Residual variances
  pos1 ~~ c(theta71, theta81)*pos1
  pos2 ~~ c(theta72, theta82)*pos2
  pos3 ~~ c(theta73, theta83)*pos3
  
  # Intercepts - Constrain 3rd intercept to 0 in both Grades
  pos1 ~ c(tau1, tau1)*1
  pos2 ~ c(tau2, tau2)*1
  pos3 ~ c(0,0)*c(tau3, tau3)*1
"

m2_fit = sem(m2, sample.cov = mcov, sample.nobs = n, sample.mean = vmean)
summary(m2_fit)
OpenMx = coef(fit2)
lavaan = coef(m2_fit)

# Get coefs in same order
lavaan <- lavaan[match(names(OpenMx), names(lavaan))]      
t(rbind(OpenMx, lavaan))
########################



### Method 3
## See model diagram in OneFactorTwoGroups3.svg
## Constrain sum of loadings to equal number of loadings
## Constrain sum of intercepts to equal 0
## With strong measurement invariance, 
## these constraints apply to both Grades.

# Factor loadings
loadings <- mxPath(from = "POS", to = names, arrows = 1,
   free = TRUE, values = 0.5,
   labels = c("lambda1", "lambda2", "lambda3"))

# Factor variances
varFac7 <- mxPath(from = "POS", arrows = 2,
   free = TRUE, values = 1, labels = "phi7")

varFac8 <- mxPath(from = "POS", arrows = 2,
   free = TRUE, values = 1, labels = "phi8")

# Factor means
means7 <- mxPath(from = "one", to = "POS", arrows = 1,
   free = TRUE, values = 0, labels = "kappa7")

means8 <- mxPath(from = "one", to = "POS", arrows = 1,
   free = TRUE, values = 1, labels = "kappa8")

# Residual variances
varRes7 <- mxPath(from = names, arrows = 2,
   free = TRUE, values = 1, 
   labels = c("theta71", "theta72", "theta73"))

varRes8 <- mxPath(from = names, arrows = 2,
   free = TRUE, values = 1,
   labels = c("theta81", "theta82", "theta83"))

# Intercepts
intercepts <- mxPath(from = "one", to = names, arrows = 1,
   free = TRUE, values = 1,
   labels = c("tau1", "tau2", "tau3"))


## Setup models for each Grade
modGr7 <- mxModel("Grade7", type = "RAM",
   manifestVars = names, latentVars = "POS",
   data7, loadings, varFac7, means7, varRes7, intercepts)

modGr8 <- mxModel("Grade8", type = "RAM",
   manifestVars = names, latentVars = "POS",
   data8, loadings, varFac8, means8, varRes8, intercepts)


## Constraints
conLoad <- mxConstraint(lambda1 + lambda2 + lambda3 == 3)
conInter <- mxConstraint(tau1 + tau2 + tau3 == 0)


## Combine the two models
fun <- mxFitFunctionMultigroup(c("Grade7.fitfunction", "Grade8.fitfunction"))
model3 <- mxModel("One Factor Two Group Model", modGr7, modGr8,
   conLoadings, conIntercepts, fun)

# Note: Constraints are added to final model, not to each of the Grade7 and Grade 8 models;
# otherwise, OpenMx will count them as 4 constraints, accounting for 4 degrees of freedom,
# instead of 2.


## Run the model and get summary
fit3 <- mxRun(model3)
summary(fit3, refModels = mxRefModels(fit3, run = TRUE))

# Counting degrees of freedom.
# Number of variables: 3
# Number of pieces of information in the covariance matrix: 6
# plus 3 means = 9 per group;
# that is, 18 for the model

# Number of parameters:
# 3 loadings (constrained to equality across the groups)
# 3 intercepts (constrained to equality across the groups)
# 6 residual variances (3 per group)
# 2 latent means (1 per group)
# 2 latent variances (1 per group)
# Total of 16
# but take away 2 for the 2 constraints = 14,
# resulting in 18 - 14 = 4 degrees of freedom.
#
# If the constraints were added to each Grade's model,
# OpenMx would have counted them as 4 constraints (even thought they are identical),
# resulting in 6 degrees of freedom for the model.


########################
## Check with lavaan
library(lavaan)

m3 <- "
  # Loadings
  POS =~ c(NA,NA)*c(lambda1, lambda1)*pos1 + c(lambda2, lambda2)*pos2 + c(lambda3, lambda3)*pos3

  # Latent variances
  POS ~~ c(phi7, phi8)*POS

  # Latent means 
  POS ~ c(kappa7, kappa8)*1

  # Residual variances
  pos1 ~~ c(theta71, theta81)*pos1
  pos2 ~~ c(theta72, theta82)*pos2
  pos3 ~~ c(theta73, theta83)*pos3
  
  # Intercepts 
  pos1 ~ c(tau1, tau1)*1
  pos2 ~ c(tau2, tau2)*1
  pos3 ~ c(tau3, tau3)*1
  
  # Constraints
  lambda1 + lambda2 + lambda3 == 3
  tau1 + tau2 + tau3 == 0
"

m3_fit = sem(m3, sample.cov = mcov, sample.nobs = n, sample.mean = vmean)
summary(m3_fit)
OpenMx = coef(fit3)
lavaan = coef(m3_fit)

# Get coefs in same order
lavaan <- lavaan[match(names(OpenMx), names(lavaan))] 
t(rbind(OpenMx, lavaan))
########################
