

## Little, T., Slegers, D., & Card, N. (2006). A non-arbitrary method of
## identifying and scaling latent variables in SEM and MACS models.
## Structural Equation Modeling, 13(1), 59-72.


## Methods of Scaling and Identification

## Demonstrates three methods of scaling in two-factor, two-group model:
## 1. Reference-group method - Constrain both latent variables' variances and means;
## 2. Marker-variable method - Constrain one loading and that indicator's intercept in both factors;
## 3. Effect-scaling method - Constrain sums of loadings and intercepts for both factors.


## Little et al assume strong metric invariance:
## Corresponding loadings and intercepts constrained to equality across groups


## Compare with results given in Table 2 (pp. 64-65)
## Little et al show results for three versions of Method 2.
## Only the third is demonstrated here.


## Load package
library(OpenMx)


## Get the data from Appendix A
# 7th grade
cor7 <- c(
   1.00000,
   0.75854,  1.00000,
   0.76214,  0.78705,  1.00000,
   0.02766,  0.00973, -0.05762,  1.00000,
  -0.06112, -0.06105, -0.14060,  0.78501,  1.00000,
  -0.02222, -0.05180, -0.10250,  0.81616,  0.81076,  1.00000)

mean7 <- c(3.13552, 2.99061, 3.06945, 1.70069, 1.52705, 1.54483)
sd7   <- c(0.66770, 0.68506, 0.70672, 0.71418, 0.66320, 0.65276)
n7    <- 380

# 8th grade
cor8 <- c(
   1.00000,
   0.81366,  1.00000,
   0.84980,  0.83523,  1.00000,
  -0.18804, -0.15524, -0.21520,  1.00000,
  -0.28875, -0.24951, -0.33769,  0.78418,  1.00000,
  -0.29342, -0.21022, -0.30553,  0.79952,  0.83156,  1.00000)

mean8 <- c(3.07338, 2.84716, 2.97882, 1.71700, 1.57955, 1.55001)
sd8   <- c(0.70299, 0.71780, 0.76208, 0.65011, 0.60168, 0.61420)
n8    <- 379


## Get the variable names from Appendix A
names = c("pos1", "pos2", "pos3", "neg1", "neg2", "neg3")


# Get full correlation matrix for each Grade
mcor7 <- matrix( , 6, 6)                         # Empty matrix
mcor7[upper.tri(mcor7, diag = TRUE)] <- cor7     # Fill the upper triangle
mcor7 <- pmax(mcor7, t(mcor7), na.rm = TRUE)     # Fill the lower triangle

mcor8 <- matrix( , 6, 6)                          # Empty matrix
mcor8[upper.tri(mcor8, diag = TRUE)] <- cor8      # Fill the upper triangle
mcor8 <- pmax(mcor8, t(mcor8), na.rm = TRUE)      # Fill the lower triangle

# Get co/variance matrix for each grade
mcov7 <- outer(sd7, sd7) * mcor7
mcov8 <- outer(sd8, sd8) * mcor8

# Name the rows and columns
dimnames(mcov7) <- list(names, names)
dimnames(mcov8) <- list(names, names)
mcov7; mcov8

names(mean7) <- names   # OpenMx requires the means be named
names(mean8) <- names

# Get data into OpenMx format
data7 <- mxData(observed = mcov7, type = "cov", means = mean7, numObs = n7)
data8 <- mxData(observed = mcov8, type = "cov", means = mean8, numObs = n8)



### Reference-Group Method
## See model diagram in Scaling1.svg


## Constrain latent variances to 1
## Constrain latent means to 0
## These constraints apply to Grade 7 only.


## Collect the bits and pieces needed by OpenMx

# Factor loadings
loadings1 <- mxPath(from = "POS", to = c("pos1", "pos2", "pos3"), arrows = 1,
   free = TRUE, values = 0.5,
   labels = c("lambda1", "lambda2", "lambda3"))

loadings2 <- mxPath(from = "NEG", to = c("neg1", "neg2", "neg3"), arrows = 1,
   free = TRUE, values = 0.5,
   labels = c("lambda4", "lambda5", "lambda6"))

# Factor variances and covariance - constrain variances to 1 for Grade 7 only
varFac7 <- mxPath(from = c("POS", "NEG"), arrows = 2, connect = "unique.pairs",
   free = c(FALSE, TRUE, FALSE), values = 1,
   labels = c("phi7_11", "phi7_12", "phi7_22"))
 
varFac8 <- mxPath(from = c("POS", "NEG"), arrows = 2, connect = "unique.pairs",
   free = TRUE, values = 1,
   labels = c("phi8_11", "phi8_12", "phi8_22"))

# Factor means - constrain means to 0 for Grade 7 only
means7 <- mxPath(from = "one", to = c("POS", "NEG"), arrows = 1,
   free = FALSE, values = 0,
   labels = c("kappa7_1", "kappa7_2"))

means8 <- mxPath(from = "one", to = c("POS", "NEG"), arrows = 1,
   free = TRUE, values = 1,
   labels = c("kappa8_1", "kappa8_2"))

# Residual variances
varRes7 <- mxPath(from = names, arrows = 2,
   free = TRUE, values = 1,
   labels = c("theta7_1", "theta7_2", "theta7_3", "theta7_4", "theta7_5", "theta7_6"))

varRes8 <- mxPath(from = names, arrows = 2,
   free = TRUE, values = 1,
   labels = c("theta8_1", "theta8_2", "theta8_3", "theta8_4", "theta8_5", "theta8_6"))

# Intercepts
intercepts <- mxPath(from = "one", to = names, arrows = 1,
   free = TRUE, values = 1,
   labels = c("tau1", "tau2", "tau3", "tau4", "tau5", "tau6"))


## Setup models for each Grade
modGr7 <- mxModel("Grade7", type = "RAM",
   manifestVars = names, latentVars = c("POS", "NEG"),
   data7, loadings1, loadings2, varFac7, means7, varRes7, intercepts)

modGr8 <- mxModel("Grade8", type = "RAM",
   manifestVars = names, latentVars = c("POS", "NEG"),
   data8, loadings1, loadings2, varFac8, means8, varRes8, intercepts)


## Combine the two models
fun <- mxFitFunctionMultigroup(c("Grade7.fitfunction", "Grade8.fitfunction"))
model1 <- mxModel("Referemce Group Method", modGr7, modGr8, fun)


## Run the model and get summary
fit1 <- mxRun(model1)
summary1 <- summary(fit1, refModels = mxRefModels(fit1, run = TRUE))
summary1

# Compare with results for Method 1 in Table 2.



### Marker-Variable Method
## See model diagram in Scaling2.svg


## Constrain 3rd loading in POS to 1 in both groups
## Constrain 1st loading in NEG to 1 in both groups
## Constrain 3rd intercept in POS to 0 in both groups
## Constrain 1st intercept in NEG to 0 in both groups


## Collect the bits and pieces needed by OpenMx


# Factor loadings - 3rd loading for POS constrained to 1
#                 - 1st loading for NEG constrained to 1
loadings1 <- mxPath(from = "POS", to = c("pos1", "pos2", "pos3"), arrows = 1,
   free = c(TRUE, TRUE, FALSE), values = c(0.5, 0.5, 1),
   labels = c("lambda1", "lambda2", "lambda3"))

loadings2 <- mxPath(from = "NEG", to = c("neg1", "neg2", "neg3"), arrows = 1,
   free = c(FALSE, TRUE, TRUE), values = c(1, 0.5, 0.5),
   labels = c("lambda4", "lambda5", "lambda6"))

# Factor variances and covariance
varFac7 <- mxPath(from = c("POS", "NEG"), arrows = 2, connect = "unique.pairs",
   free = TRUE, values = 1, 
   labels = c("phi7_11", "phi7_12", "phi7_22"))

varFac8 <- mxPath(from = c("POS", "NEG"), arrows = 2, connect = "unique.pairs",
   free = TRUE, values = 1,
   labels = c("phi8_11", "phi8_12", "phi8_22"))

# Factor means
means7 <- mxPath(from = "one", to = c("POS", "NEG"), arrows = 1,
   free = TRUE, values = 1,
   labels = c("kappa7_1", "kappa7_2"))

means8 <- mxPath(from = "one", to = c("POS", "NEG"), arrows = 1,
   free = TRUE, values = 1,
   labels = c("kappa8_1", "kappa8_2"))

# Residual variances
varRes7 <- mxPath(from = names, arrows = 2,
   free = TRUE, values = 1,
   labels = c("theta7_1", "theta7_2", "theta7_3", "theta7_4", "theta7_5", "theta7_6"))

varRes8 <- mxPath(from = names, arrows = 2,
   free = TRUE, values = 1, 
   labels = c("theta8_1", "theta8_2", "theta8_3", "theta8_4", "theta8_5", "theta8_6"))

# Intercepts - 3rd intercept for POS & 1st intercept for NEG constrained to 0
intercepts <- mxPath(from = "one", to = names, arrows = 1,
   free = c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE), values = c(1, 1, 0, 0, 1, 1),
   labels = c("tau1", "tau2", "tau3", "tau4", "tau5", "tau6"))


## Setup models for each Grade
modGr7 <- mxModel("Grade7", type = "RAM",
   manifestVars = names, latentVars = c("POS", "NEG"),
   data7, loadings1, loadings2, varFac7, means7, varRes7, intercepts)

modGr8 <- mxModel("Grade8", type = "RAM",
   manifestVars = names, latentVars = c("POS", "NEG"),
   data8, loadings1, loadings2, varFac8, means8, varRes8, intercepts)


## Combine the two models
fun <- mxFitFunctionMultigroup(c("Grade7.fitfunction", "Grade8.fitfunction"))
model2 <- mxModel("Referemce Group Method", modGr7, modGr8, fun)


## Run the model and get summary
fit2 <- mxRun(model2)
summary2 <- summary(fit2, refModels = mxRefModels(fit2, run = TRUE))
summary2

# Compare with results for third Method 2 in Table 2.



### Effects-Scaling Method
## See model diagram in Scaling3.svg


## Constrain loadings to add to 3 in both factors for both Grades
## Constrain intercepts to add to 0 in both factors for both Grades


## Collect the bits and pieces needed by OpenMx

# Factor loadings
loadings1 <- mxPath(from = "POS", to = c("pos1", "pos2", "pos3"), arrows = 1,
   free = TRUE, values = 0.5,
   labels = c("lambda1", "lambda2", "lambda3"))

loadings2 <- mxPath(from = "NEG", to = c("neg1", "neg2", "neg3"), arrows = 1,
   free = TRUE, values = 0.5,
   labels = c("lambda4", "lambda5", "lambda6"))

# Factor variances and covariance
varFac7 <- mxPath(from = c("POS", "NEG"), arrows = 2, connect = "unique.pairs",
   free = TRUE, values = 1, 
   labels = c("phi7_11", "phi7_12", "phi7_22"))

varFac8 <- mxPath(from = c("POS", "NEG"), arrows = 2, connect = "unique.pairs",
   free = TRUE, values = 1,
   labels = c("phi8_11", "phi8_12", "phi8_22"))

# Factor means
means7 <- mxPath(from = "one", to = c("POS", "NEG"), arrows = 1,
   free = TRUE, values = 1,
   labels = c("kappa7_1", "kappa7_2"))

means8 <- mxPath(from = "one", to = c("POS", "NEG"), arrows = 1,
   free = TRUE, values = 1,
   labels = c("kappa8_1", "kappa8_2"))

# Residual variances
varRes7 <- mxPath(from = names, arrows = 2,
   free = TRUE, values = 1,
   labels = c("theta7_1", "theta7_2", "theta7_3", "theta7_4", "theta7_5", "theta7_6"))

varRes8 <- mxPath(from = names, arrows = 2,
   free = TRUE, values = 1,
   labels = c("theta8_1", "theta8_2", "theta8_3", "theta8_4", "theta8_5", "theta8_6"))

# Intercepts
intercepts <- mxPath(from = "one", to = names, arrows = 1,
   free = TRUE, values = 1,
   labels = c("tau1", "tau2", "tau3", "tau4", "tau5", "tau6"))


## Setup models for each Grade
modGr7 <- mxModel("Grade7", type = "RAM",
   manifestVars = names, latentVars = c("POS", "NEG"),
   data7, loadings1, loadings2, varFac7, means7, varRes7, intercepts)

modGr8 <- mxModel("Grade8", type = "RAM",
   manifestVars = names, latentVars = c("POS", "NEG"),
   data8, loadings1, loadings2, varFac8, means8, varRes8, intercepts)


## Constraints
conLoad1  <- mxConstraint(lambda1 + lambda2 + lambda3 == 3)
conLoad2  <- mxConstraint(lambda4 + lambda5 + lambda6 == 3)
conInter1 <- mxConstraint(tau1 + tau2 + tau3 == 0)
conInter2 <- mxConstraint(tau4 + tau5 + tau6 == 0)


## Combine the two models
fun <- mxFitFunctionMultigroup(c("Grade7.fitfunction", "Grade8.fitfunction"))
model3 <- mxModel("Referemce Group Method", modGr7, modGr8,
          conLoad1, conLoad2, conInter1, conInter2, fun)


## Run the model and get summary
fit3 <- mxRun(model3)
summary3 <- summary(fit3, refModels = mxRefModels(fit3, run = TRUE))
summary3

# Compare with results for Method 3 in Table 2.


## Get fit measures
#  Compare with fit measures given on page 66
models <- list(
   "Method 1" = summary1,
   "Method 2" = summary2,
   "Method 3" = summary3)

measures = c("Chi", "ChiDoF", "p", "CFI", "TLI", "RMSEA")

fit1 <- do.call(rbind, lapply(models, `[`, measures))
fit2 <- do.call(rbind, lapply(models, `[[`, "RMSEACI"))

cbind(fit1, fit2)


## Note: RMSEA given by OpenMx does not agree with the value given
#  by Little et al. OpenMx does not adjust RMSEA in multiple group models -
#  see the RMSEA section in ?mxSummary for brief explanation, and
# https://openmx.ssri.psu.edu/index.php/forums/opensem-forums/fit-functions/rmsea-multiple-group-analysis
#  for another discussion.
#  LISREL, lavaan, Mplus (and possibly other packages) do make the adjustment.

# Equation for RMSEA, as given by Steiger (1998, Eq. 25, p. 417):
#  RMSEA = sqrt((ChiSq/df - 1) / n)

#  Substitute values from above;
#  gives 'unadjusted RMSEA' as obtained by OpneMx.
sqrt((58.60023/24 - 1) / 759)

#  Steiger (1998, p. 417) states that 'adjusted RMSEA' can be obtained
#  by multiplying 'unadjusted RMSEA' by sqrt(g) (where g is the number
#  of groups).
sqrt(2) * sqrt((58.60023/24 - 1) / 759)

#  Gives 'adjusted RMSEA' as given by Little et al (ie, LISREL)
#  and by lavaan (see 'SEMs_with_lavaan').


# Steiger, J. (1998). A note on multiple sample extensions of the RMSEA 
# fit index. Structural Equation Modeling, 5(4), 411-419.
