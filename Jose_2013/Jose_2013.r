
## Jose, P. (2013). Doing statistical mediation and moderation. 
## New York, NY: Guilford Press.

## Chapter 3 describes a basic three-variable mediation analysis.

## Load packages
library(OpenMx)
library(haven)   # To read SPSS data files

# OpenMx is very verbose;
# it requires everything to be coded;
# means and intercepts must be included.

## Get the data from Guilford Press web site
url <- "http://www.guilford.com/add/jose/mediation_example.sav"
df <- data.frame(haven::read_spss(url))

str(df)
head(df)
summary(df)

## The model is shown in Fig 3.3 (p. 46);
## or see the `images` folder.

#### Collect the bits and pieces needed by OpenMx
## Get the names of the variables
manifest <- names(df)

## Get data into OpenMx format
dataRaw <- mxData(observed = df, type = "raw")

## Regressions
#  - Arrows from "ple" and "grat" to "shs".
#  - Arrows are single headed - arrow head at "shs".
#  - Starting values are 0.5.  When they are the same,
#    they need to be listed once only.
#  - Labels are "cprime" for arrow "ple" to "shs";
#    and "b" for arrow "grat" to "shs".
regPaths1 <- mxPath(from = c("ple", "grat"), to = "shs",
   arrows = 1, values = 0.5,
   labels = c("cprime", "b"))

#  - Arrow from "ple" to "grat".
#  - Arrow is single headed - arrow head at "grat".
#  - Starting value is 0.5.
#  - Label is "a".
regPaths2 <- mxPath(from = "ple", to = "grat",
   arrows = 1, values = 0.5,
   labels = "a")

## Variances
## Exogenous variables ("ple") have variances;
## Endogenous variables ("grat" and "shs") have residual or error variances.
## Variance for exogenous variable is not shown in the model diagram.
## The distinction does not matter to OpenMx, but I distinguish in the labels.

#  - Arrows are from manifest variables to manifest variables
#    (when "arrows = 2", the "to" argument can be omitted); 
#    ie, from "ple" to "ple"; from "grat" to "grat"; and from "shs" to "shs".
#  - Arrows are two headed.
#  - Starting values are 1.
#  - Labels are "vPLE" "eGRAT", and "eSHS".
varPaths <- mxPath(from = manifest,
   arrows = 2, values = 1,
   labels = c("vPLE", "eGRAT", "eSHS"))

## Means and intercepts
## Exogenous variables ("ple") have means;
## Endogenous variables ("grat" and "shs") have intercepts.
## Regress variables on a constant - in OpenMx, "one".
## Means and intercepts are not shown in the model diagram.
## The distinction does not matter to OpenMx, but I distinguish in the labels.

#  - Arrows from "one" to the manifest variables;
#    ie, from "one" to "ple"; from "one" to "grat"; from "one" to "shs".
#  - Arrows are single headed.
#  - Starting values are 1.
#  - labels are "mPLE", "iGRAT", "iSHS".
means <- mxPath(from = "one", to = manifest,
   arrows = 1, values = 1,
   labels = c("mPLE", "iGRAT", "iSHS"))

## Indirect and total effects
indirect <- mxAlgebra(a * b, name = "indirect")
total <- mxAlgebra(a * b + cprime, name = "total")

## Setup the model with all the bits
medModel <- mxModel(model = "Mediation",
   type = "RAM",
   data = dataRaw,
   manifestVars = manifest,
   varPaths,
   regPaths1, regPaths2,
   means, 
   indirect, total)

## Run the model and get summary
# Compare with regression outputs in Tables 3.4 and 3.5 (p. 52)
fit <- mxRun(medModel)
summary(fit)   # Note: summary() does not return indirect and total effects
coef(fit)      # Just the estimates

## Extract indirect and total effects (and their standard errors) from "fit" object
# Compare with unstandardised indirect and total effect in Fig 3.6 (p. 59).
estimates <- mxEval(c(indirect, total), fit); estimates
SE <- sapply(c("indirect", "total"), function(x) mxSE(x, fit, silent = TRUE)); SE

## To get the standardised effects, mxStandardizeRAMpaths(fit) will give standardised
## effects for free parameters, but not derived effects.

# Could use mxStandardizeRAMpaths(fit) to extract standardised "a", "b", and "cprime",
# then calculate standardised indirect effect (a * b),
# and standardised total effect (a * b + cprime) by hand.

# Compare with standardised indirect and total effect in Fig 3.6.
# Compare standardised regression coefficients given in Fig 3.6, or
# in Tables 3.4 and 3.5.
mxStandardizeRAMpaths(fit)
estZ <- mxStandardizeRAMpaths(fit)[1:3, 8]
names(estZ) <- mxStandardizeRAMpaths(fit)[1:3, "label"]; estZ

estZ["indirect"] <- estZ["a"] * estZ["b"]
estZ["total"] <- estZ["indirect"] + estZ["cprime"]
estZ

## To get bootstrap CIs
fitBoot <- mxBootstrap(fit, 2000)
# statistics <- summary(fitBoot, boot.quantile = c(0.025, 0.975),
#    boot.SummaryType = "bcbci"); statistics
# Note: No defined effects
ci <- mxBootstrapEval(c(a, b, cprime, indirect, total), fitBoot,
   bq = c(0.025, 0.975), method = "bcbci"); ci

## To get likelihood-based CIs
ci <- mxCI(c("a", "b", "cprime", "indirect", "total"))

# Add to the model
medModel <- mxModel(medModel, ci)

# Run the model
fit <- mxRun(medModel, intervals = TRUE)
summary(fit)$CI
