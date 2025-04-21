

## Kurbanoglu, N. & Takunyaci, M. (2021). A structural equation modeling
## on relationship between self-efficacy, physics laboratory anxiety
## and attitudes. Journal of Family, Counseling and Education, 6(1), 47-56.


## Load package
library(OpenMx)


## Get the data from Table 1 (p. 50)
# Vectors of correlations, standard deviations, and means.
# Correlations entered row-by-row.
# (If entered column-by-column, use 'lower.tri' in line 32 below.)
vcor <- c(
   1,
   0.30,  1,
  -0.42, -0.32,  1)
vsd <- c(8.81, 7.95, 18.30)             # Standard deviations
vmean <- c(56.57, 40.39, 68.22)         # Means
n <- 513                                # Sample size


## Get the variable names (make sure order is the same as in Table 1)
names <- c("Att", "SE", "Anx")


## Get the (co)variance matrix
# First, get full correlation matrix
mcor = matrix( , 3, 3)                           # Empty matrix
mcor[upper.tri(mcor, diag = TRUE)] <- vcor       # Fill the upper triangle
mcor = pmax(mcor, t(mcor), na.rm = TRUE)         # Fill the lower triangle

# Get covariances
mcov <- outer(vsd, vsd) * mcor

# Name the rows and columns
dimnames(mcov) <- list(names, names); mcov

names(vmean) = names   # OpenMx requires the means be named


#### Collect the bits and pieces needed by OpenMx

## Get data into OpenMx format
dataCov <- mxData(observed = mcov, type = "cov", means = vmean, numObs = n)


## The model is shown in Fig 1 (p. 51); also see Kurbanoglu_2021.pdf in images folder


## Regressions
regPaths1 <- mxPath(from = c("SE", "Att"), to = "Anx",
   arrows = 1, values = 0.5, labels = c("cprime", "b"))

regPaths2 <- mxPath(from = "SE", to = "Att",
   arrows = 1, values = 0.5, labels = "a")


## Variances
varPaths <- mxPath(from = names, 
   arrows = 2, values = 1, labels = c("vAtt", "vSE", "vAnx"))


## Means and intercepts
means <- mxPath(from = "one", to = names,
   arrows = 1, values = 1, labels = c("iAtt", "mSE", "iAnx"))


## Indirect and total effects
indirect <- mxAlgebra(a * b, name = "indirect")
total <- mxAlgebra(a * b + cprime, name = "total")


## Setup the model with all the bits
medModel <- mxModel(model = "Mediation",
   type = "RAM",
   data = dataCov,
   manifestVars = names,
   varPaths,
   regPaths1, regPaths2,
   means, 
   indirect, total)	


## Run the model and get summary
fit <- mxRun(medModel)
summary(fit)


## Extract indirect and total effects (and their standard errors) from "fit" object
estimates <- mxEval(c(indirect, total), fit); estimates
SE <- sapply(c("indirect", "total"), function(x) mxSE(x, fit, silent = TRUE)); SE


## Get the standardised effects
mxStandardizeRAMpaths(fit)
estZ <- mxStandardizeRAMpaths(fit)[1:3, 8]
names(estZ) <- mxStandardizeRAMpaths(fit)[1:3, "label"]; estZ

# Calculate standardised effects for "indirect" and "total" by hand
estZ["indirect"] <- estZ["a"] * estZ["b"]
estZ["total"] <- estZ["indirect"] + estZ["cprime"]
estZ

# Compare with standardised estimates in Fig 1


## R squares - calculate by hand
# For variables implicated in regression,
# R square is given by matrix product of:
#   row matrix of correlations and
#   column matrix of standardised regression coefficients

# Get correlations (with names) and standardised regression coefficients
dimnames(mcor) = list(names, names); mcor
estZ

# R Squared for SE predicting Att:
# Correlation between Att and SE, and
# "a" standardised regression coefficient
RsqAtt = mcor["Att", "SE"] * estZ["a"]


# R Squared for SE and Att predicting Anx
# Correlations between Anx and SE, and Anx and Att
# "b" and "cprime" standardised regression coefficients
RsqAnx =  matrix(mcor["Anx", c("Att", "SE")], 1)  %*%
   matrix(estZ[c("b", "cprime")], 2)

# Combine them
Rsq = matrix(c(RsqAtt, RsqAnx), 
  dimnames = list(c("Att", "Anx"), "R Square"),
  2)
Rsq

# Compare with R squares given in Fig 1


## Likelihood-based CIs
ci <- mxCI(c("a", "b", "cprime", "indirect", "total"))

# Add to the model
medModel <- mxModel(medModel, ci)

# Run the model
fit <- mxRun(medModel, intervals = TRUE)
summary(fit)$CI
