
#### Methods of Scaling and Identification


## Demonstrates three methods of scaling in two-factor model:
## 1. Reference-Group Method - Constrain latent variables' variances and means;
## 2. Marker-Variable Method - Constrain one loading and that indicator's 
##    intercept in both factors;
## 3. Effects-Scaling Method - Constrain sums of loadings and intercepts
##    for both factors.


## Compare results from OpenMx with lavaan's results


## Two Factors
#  "Positive Affect" and "Negative Affect" factors for 7th grade
#  Data in Appendix A of:
#
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
   0.76214,  0.78705,  1.00000,
   0.02766,  0.00973, -0.05762,  1.00000,
  -0.06112, -0.06105, -0.14060,  0.78501,  1.00000,
  -0.02222, -0.05180, -0.10250,  0.81616,  0.81076,  1.00000)

vmean <- c(3.13552, 2.99061, 3.06945, 1.70069, 1.52705, 1.54483)
vsd <- c(0.66770, 0.68506, 0.70672, 0.71418, 0.66320, 0.65276)
n <- 380

# Variable names
names <- c("pos1", "pos2", "pos3", "neg1", "neg2", "neg3")

# Get full correlation matrix
mcor = matrix( , 6, 6)                           # Empty matrix
mcor[upper.tri(mcor, diag = TRUE)] <- vcor       # Fill the upper triangle
mcor = pmax(mcor, t(mcor), na.rm = TRUE)         # Fill the lower triangle

# Get (co)variance matrix 
mcov <- outer(vsd, vsd) * mcor

# Name the rows and columns
dimnames(mcov) <- list(names, names)
mcov

names(vmean) = names   # OpenMx requires the means be named

# Get data into OpenMx format
data <- mxData(observed = mcov, type = "cov", means = vmean, numObs = n)



### Method 1
## See model diagram in TwoFactors1.svg
# Labels are Little el al's Greek labels 

## Constrain latent variances to 1
## Constrain latent means to 0


## Collect the bits and pieces needed by OpenMx

# Factor loadings
loadings1 <- mxPath(from = "POS", to = c("pos1", "pos2", "pos3"), arrows = 1,
   free = TRUE, values = 0.5,
   labels = c("lambda1", "lambda2", "lambda3"))

loadings2 <- mxPath(from = "NEG", to = c("neg1", "neg2", "neg3"), arrows = 1,
   free = TRUE, values = 0.5
   labels = c("lambda4", "lambda5", "lambda6"))

# Factor variances and covariance - constrain variances to 1
varFac <- mxPath(from = c("POS", "NEG"), arrows = 2, connect = "unique.pairs",
   free = c(FALSE, TRUE, FALSE), values = 1, 
   labels = c("phi11", "phi12", "phi22"))

# Factor means - constrain means to 0
means <- mxPath(from = "one", to = c("POS", "NEG"), arrows = 1,
   free = FALSE, values = 0,
   labels = c("kappa1", "kappa2"))

# Residual variances
varRes <- mxPath(from = names, arrows = 2,
   free = TRUE, values = 1, 
   labels = c("theta1", "theta2", "theta3", "theta4", "theta5", "theta6"))

# Intercepts
intercepts <- mxPath(from = "one", to = names, arrows = 1,
   free = TRUE, values = 1,
   labels = c("tau1", "tau2", "tau3", "tau4", "tau5", "tau6"))

## Setup the model
model1 <- mxModel("Two Factor Model", type = "RAM",
   manifestVars = names, latentVars = c("POS", "NEG"),
   data, loadings1, loadings2, varFac, means, varRes, intercepts)


## Run the model and get summary
fit1 <- mxRun(model1)
summary(fit1)   

# Number of variables is 6;
# Number of pieces of information in covariance matrix: (6 X 7) / 2 = 21
# plus 6 means = 27 pieces of information
# Number of parameters: 
#   6 loadings (3 per factor)
#   1 covariance between factors
#   6 residual variances (3 per factor)
#   6 intercepts (3 per factor)
#   Total of 19 parameters
 
# Therefore, degrees of freedom = 8,
# and chi sq and fit indices are calculated. 
# But make sure OpenMx estimates reference models (saturated and independence)
# upon which to base calculations for chi sq and fit indices.

summary(fit1, refModels = mxRefModels(fit1, run = TRUE)) 
coef(fit1)


########################
## Check with lavaan
library(lavaan)

m1 <- "
  # Loadings
  POS =~ NA*lambda1*pos1 + lambda2*pos2 + lambda3*pos3
  NEG =~ NA*lambda4*neg1 + lambda5*neg2 + lambda6*neg3

  # Latent variances and covariance - constrain variances to 1
  POS ~~ 1*phi11*POS
  NEG ~~ 1*phi22*NEG
  POS ~~ phi12*NEG

  # Latent means - constrain means to 0
  POS ~ 0*kappa1*1
  NEG ~ 0*kappa2*1

  # Residual variances
  pos1 ~~ theta1*pos1
  pos2 ~~ theta2*pos2
  pos3 ~~ theta3*pos3
  neg1 ~~ theta4*neg1
  neg2 ~~ theta5*neg2
  neg3 ~~ theta6*neg3

  # Intercepts 
  pos1 ~ tau1*1
  pos2 ~ tau2*1
  pos3 ~ tau3*1
  neg1 ~ tau4*1
  neg2 ~ tau5*1
  neg3 ~ tau6*1
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
## See model diagram in TwoFactors2.svg
## Constrain 3rd loading for POS and 1st loading for NEG to 1
## Constrain 3rd intercept for POS and 1st intercept for NEG to 0

# Factor loadings - Constrain 3rd loading for POS & 1st loading for NEG to 1
loadings1 <- mxPath(from = "POS", to = c("pos1", "pos2", "pos3"), arrows = c(0.5, 0.5, 1),
   free = c(TRUE, TRUE, FALSE), values = 1,
   labels = c("lambda1", "lambda2", "lambda3"))

loadings2 <- mxPath(from = "NEG", to = c("neg1", "neg2", "neg3"), arrows = c(1, 0.5, 0.5),
   free = c(FALSE, TRUE, TRUE), values = 1,
   labels = c("lambda4", "lambda5", "lambda6"))

# Factor variances and covariance
varFac <- mxPath(from = c("POS", "NEG"), arrows = 2, connect = "unique.pairs",
   free = TRUE, values = c(1, 0.5, 1), 
   labels = c("phi11", "phi12", "phi22"))

# Factor means
means <- mxPath(from = "one", to = c("POS", "NEG"), arrows = 1,
   free = TRUE, values = 1,
   labels = c("kappa1", "kappa2"))

# Residual variances
varRes <- mxPath(from = names, arrows = 2,
   free = TRUE, values = 1,
   labels = c("theta1", "theta2", "theta3", "theta4", "theta5", "theta6"))

# Intercepts - constrain 3rd intercept for POS & 1st intercept for NEG to 0
intercepts <- mxPath(from = "one", to = names, arrows = 1,
   free = c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE), 
   values = c(1, 1, 0, 0, 1, 1),  
   labels = c("tau1", "tau2", "tau3", "tau4", "tau5", "tau6"))

## Setup the model
model2 <- mxModel("Two Factor Model", type = "RAM",
   manifestVars = names, latentVars = c("POS", "NEG"),
   data, loadings1, loadings2, varFac, means, varRes, intercepts)


## Run the model and get summary
fit2 <- mxRun(model2)
summary(fit2, refModels = mxRefModels(fit2, run = TRUE))
coef(fit2)


########################
## Check with lavaan
library(lavaan)

m2 <- "
  # Loadings - Constrain 3rd loading for POS & 1st loading for NEG to 1
  POS =~ NA*lambda1*pos1 + lambda2*pos2 + 1*lambda3*pos3
  NEG =~  1*lambda4*neg1 + lambda5*neg2 + lambda6*neg3

  # Latent variances and covariance
  POS ~~ phi11*POS
  NEG ~~ phi22*NEG
  POS ~~ phi12*NEG

  # Latent means 
  POS ~ kappa1*1
  NEG ~ kappa2*1

  # Residual variances
  pos1 ~~ theta1*pos1
  pos2 ~~ theta2*pos2
  pos3 ~~ theta3*pos3
  neg1 ~~ theta4*neg1
  neg2 ~~ theta5*neg2
  neg3 ~~ theta6*neg3  

  # Intercepts - Constrain 3rd intercept for POS & 1st intercept for NEG to 0
  pos1 ~ tau1*1
  pos2 ~ tau2*1
  pos3 ~ 0*tau3*1
  neg1 ~ 0*tau4*1
  neg2 ~ tau5*1
  neg3 ~ tau6*1
"

m2_fit = sem(m2, sample.cov = mcov, sample.nobs = n, sample.mean = vmean)
OpenMx = coef(fit2)
lavaan = coef(m2_fit)

# Get coefs in same order
lavaan <- lavaan[match(names(OpenMx), names(lavaan))]      
t(rbind(OpenMx, lavaan))
########################



### Method 3
## See model diagram in TwoFactors3.svg
## Constrain sum of loadings to equal number of loadings
## Constrain sum of intercepts to equal 0

# Factor loadings
loadings1 <- mxPath(from = "POS", to = c("pos1", "pos2", "pos3"), arrows = 0.5,
   free = TRUE, values = 1,
   labels = c("lambda1", "lambda2", "lambda3"))

loadings2 <- mxPath(from = "NEG", to = c("neg1", "neg2", "neg3"), arrows = 0.5,
   free = TRUE, values = 1,
   labels = c("lambda4", "lambda5", "lambda6"))

# Factor variances and covariance
varFac <- mxPath(from = c("POS", "NEG"), arrows = 2, connect = "unique.pairs",
   free = TRUE, values = 1, labels = c("phi11", "phi12", "phi22"))

# Factor means
means <- mxPath(from = "one", to = c("POS", "NEG"), arrows = 1,
   free = TRUE, values = 1,
   labels = c("kappa1", "kappa2"))

# Residual variances
varRes <- mxPath(from = names, arrows = 2,
   free = TRUE, values = 1,
   labels = c("theta1", "theta2", "theta3", "theta4", "theta5", "theta6"))

# Intercepts
intercepts <- mxPath(from = "one", to = names, arrows = 1,
   free = TRUE, values = 1, 
   labels = c("tau1", "tau2", "tau3", "tau4", "tau5", "tau6"))
   
# Constraints
conLoadPOS <- mxConstraint(lambda1 + lambda2 + lambda3 == 3)
conLoadNEG <- mxConstraint(lambda4 + lambda5 + lambda6 == 3)
conInterPOS <- mxConstraint(tau1 + tau2 + tau3 == 0)
conInterNEG <- mxConstraint(tau4 + tau5 + tau6 == 0)


## Setup the model
model3 <- mxModel("Two Factor Model", type = "RAM",
   manifestVars = names, latentVars = c("POS", "NEG"),
   data, loadings1, loadings2, varFac, means, varRes, intercepts,
   conLoadPOS, conLoadNEG, conInterPOS, conInterNEG)


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
  NEG =~ NA*lambda4*neg1 + lambda5*neg2 + lambda6*neg3

  # Latent variances and covariance
  POS ~~ phi11*POS
  NEG ~~ phi22*NEG
  POS ~~ phi12*NEG

  # Latent means
  POS ~ kappa1*1
  NEG ~ kappa2*1

  # Residual variances
  pos1 ~~ theta1*pos1
  pos2 ~~ theta2*pos2
  pos3 ~~ theta3*pos3
  neg1 ~~ theta4*neg1
  neg2 ~~ theta5*neg2
  neg3 ~~ theta6*neg3

  # Intercepts 
  pos1 ~ tau1*1
  pos2 ~ tau2*1
  pos3 ~ tau3*1
  neg1 ~ tau4*1
  neg2 ~ tau5*1
  neg3 ~ tau6*1
  
  # Constraints
  lambda1 + lambda2 + lambda3 == 3
  lambda4 + lambda5 + lambda6 == 3

  tau1 + tau2 + tau3 == 0
  tau4 + tau5 + tau6 == 0
"

m3_fit = sem(m3, sample.cov = mcov, sample.nobs = n, sample.mean = vmean)
summary(m3_fit)
OpenMx = coef(fit3)
lavaan = coef(m3_fit)

# Get coefs in same order
lavaan <- lavaan[match(names(OpenMx), names(lavaan))]
t(rbind(OpenMx, lavaan))
########################
