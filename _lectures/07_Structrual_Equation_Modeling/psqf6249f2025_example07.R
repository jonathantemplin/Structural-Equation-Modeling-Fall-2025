#Step 1: Find path of data file and copy path (if using windows)
#Step 2: In the console below, type readClipboard()
#Step 3: Copy and paste R's path to the line below in quotes

#AUTOMATING PACKAGES NEEDED FOR ANALYSES--------------------------------------------------------------------
needed_packages = c("lavaan","semPlot")

for (i in 1:length(needed_packages)){
  haspackage = require(needed_packages[i], character.only = TRUE)
  if (haspackage==FALSE){
    install.packages(needed_packages[i])
    library(needed_packages[i], character.only = TRUE)
  }
}

#FUNCTIONS FOR ANALYSES BELOW -------------------------------------------------------------------------------
data01 = read.csv(file = "gamblingdata.csv", na.strings="99")

#MODEL 01: Gambling GRI Single Factor Model ----------------------------------------------------------------
model01.syntax = "
  GAMBLING =~ GRI1 + GRI3 + GRI5
"

model01.fit = sem(
  model01.syntax, 
  data = data01, 
  estimator = "MLR", 
  mimic="Mplus"
)
summary(model01.fit, fit.measures=TRUE, rsquare=TRUE, standardized=TRUE)


#display normalized residual covariances
residuals(model01.fit, type="normalized")

#plot path diagram with standardized coefficients
semPaths(model01.fit,intercepts = FALSE, residuals = TRUE, style="mx", layout="tree", rotation=1, what = "std",
         optimizeLatRes=TRUE, whatLabels = "std", sizeLat = 5, sizeLat2=5, sizeMan=5, sizeMan2=5)

#model fitted values:
fitted.values(model01.fit)

#MODEL 02: Full Structural Equation Model ------------------------------------------------------------------
model02.syntax = "
  GAMBLING =~ GRI1 + GRI3 + GRI5
  GAMBLING ~ Student
"

model02.fit = sem(
  model02.syntax, 
  data=data01, 
  estimator = "MLR", 
  mimic="Mplus"
)
summary(model02.fit, fit.measures=TRUE, rsquare=TRUE, standardized=TRUE)

#display normalized residual covariances
residuals(model02.fit, type="normalized")

#display modification indices
modindices(model02.fit)


#MODEL Configural: Setting up invariance testing of GRI by student status ==========================

configural.syntax = "
#===================================================================================================
#Factor loadings all freely estimated in both groups with label for each group
  GAMBLING =~ c(SL1, NSL1)*GRI1 + c(SL3, NSL3)*GRI3 + c(SL5, NSL5)*GRI5
#===================================================================================================
#Item intercepts all freely estimated in both groups with label for each group
  GRI1 ~ c(SI1, NSI1)*1
  GRI3 ~ c(SI3, NSI3)*1
  GRI5 ~ c(SI5, NSI5)*1
#===================================================================================================
#Redidual variances all freely estimated with label for each group
  GRI1 ~~ c(SR1, NSR1)*GRI1
  GRI3 ~~ c(SR3, NSR3)*GRI3
  GRI5 ~~ c(SR5, NSR5)*GRI5
#===================================================================================================
#Factor variances all freely estimated in both groups with label for each group
  GAMBLING ~~ c(1, 1)*GAMBLING
#===================================================================================================
#Factor means all freely estimated in both groups with label for each group
  GAMBLING ~ c(0, 0)*1
#===================================================================================================  
"
modelConfigural.fit = lavaan(
  configural.syntax, 
  data=data01, 
  estimator = "MLR", 
  mimic="Mplus", 
  group = "Student"
)

summary(modelConfigural.fit, fit.measures=TRUE, rsquare=TRUE, standardized=TRUE)

#MODEL Configural: Testing equality of loadings of GRI by student status ===========================

metric.syntax = "
#===================================================================================================
#Factor loadings all freely estimated in both groups with label for each group
  GAMBLING =~ c(L1, L1)*GRI1 + c(L3, L3)*GRI3 + c(L5, L5)*GRI5
#===================================================================================================
#Item intercepts all freely estimated in both groups with label for each group
  GRI1 ~ c(SI1, NSI1)*1
  GRI3 ~ c(SI3, NSI3)*1
  GRI5 ~ c(SI5, NSI5)*1
#===================================================================================================
#Redidual variances all freely estimated with label for each group
  GRI1 ~~ c(SR1, NSR1)*GRI1
  GRI3 ~~ c(SR3, NSR3)*GRI3
  GRI5 ~~ c(SR5, NSR5)*GRI5
#===================================================================================================
#Factor variances all freely estimated in one group but not other for identification
  GAMBLING ~~ c(1, NA)*GAMBLING
#===================================================================================================
#Factor means all freely estimated in both groups with label for each group
  GAMBLING ~ c(0, 0)*1
#===================================================================================================  
"
modelMetric.fit = lavaan(
  metric.syntax, 
  data=data01, 
  estimator = "MLR", 
  mimic="Mplus", 
  group = "Student"
)

summary(modelMetric.fit, fit.measures=TRUE, rsquare=TRUE, standardized=TRUE)

# likelihood ratio test: all loadings tested at once:
anova(modelConfigural.fit, modelMetric.fit)

# looks like invariance--move to scalar

#MODEL Scalar: Testing equality of loadings of GRI by student status ===============================

scalar.syntax = "
#===================================================================================================
#Factor loadings all freely estimated in both groups with label for each group
  GAMBLING =~ c(L1, L1)*GRI1 + c(L3, L3)*GRI3 + c(L5, L5)*GRI5
#===================================================================================================
#Item intercepts all freely estimated in both groups with label for each group
  GRI1 ~ c(I1, I1)*1
  GRI3 ~ c(I3, I3)*1
  GRI5 ~ c(I5, I5)*1
#===================================================================================================
#Redidual variances all freely estimated with label for each group
  GRI1 ~~ c(SR1, NSR1)*GRI1
  GRI3 ~~ c(SR3, NSR3)*GRI3
  GRI5 ~~ c(SR5, NSR5)*GRI5
#===================================================================================================
#Factor variances all freely estimated in one group but not other for identification
  GAMBLING ~~ c(1, NA)*GAMBLING
#===================================================================================================
#Factor means all freely estimated in one group but not other for identification
  GAMBLING ~ c(0, NA)*1
#===================================================================================================  
"
modelScalar.fit = lavaan(
  scalar.syntax, 
  data=data01, 
  estimator = "MLR", 
  mimic="Mplus", 
  group = "Student"
)

summary(modelScalar.fit, fit.measures=TRUE, rsquare=TRUE, standardized=TRUE)

# likelihood ratio test: all loadings tested at once:
anova(modelMetric.fit, modelScalar.fit)

# modification indices inspection:
scalarMI1 = lavTestScore(modelScalar.fit)

#change label values
labelMap = data.frame(
  lhs = modelScalar.fit@ParTable$plabel,
  parameter = modelScalar.fit@ParTable$label
)

scalarMI1$uni = merge(x = scalarMI1$uni, y = labelMap, by = "lhs", all.x = TRUE)

# reorder by decreasing values of X2:
scalarMI1$uni = scalarMI1$uni[order(scalarMI1$uni$X2, decreasing = TRUE),]

#restrict to only means shown
scalarMI1$uni[grep(x = scalarMI1$uni$parameter, pattern = "I"),]

# free parameter for I3 first:

scalar.syntax1 = "
#===================================================================================================
#Factor loadings all freely estimated in both groups with label for each group
  GAMBLING =~ c(L1, L1)*GRI1 + c(L3, L3)*GRI3 + c(L5, L5)*GRI5
#===================================================================================================
#Item intercepts all freely estimated in both groups with label for each group
  GRI1 ~ c(I1, I1)*1
  GRI3 ~ c(SI3, NSI3)*1
  GRI5 ~ c(I5, I5)*1
#===================================================================================================
#Redidual variances all freely estimated with label for each group
  GRI1 ~~ c(SR1, NSR1)*GRI1
  GRI3 ~~ c(SR3, NSR3)*GRI3
  GRI5 ~~ c(SR5, NSR5)*GRI5
#===================================================================================================
#Factor variances all freely estimated in one group but not other for identification
  GAMBLING ~~ c(1, NA)*GAMBLING
#===================================================================================================
#Factor means all freely estimated in one group but not other for identification
  GAMBLING ~ c(0, NA)*1
#===================================================================================================  
"
modelScalar.fit1 = lavaan(
  scalar.syntax1, 
  data=data01, 
  estimator = "MLR", 
  mimic="Mplus", 
  group = "Student"
)

summary(modelScalar.fit1, fit.measures=TRUE, rsquare=TRUE, standardized=TRUE)

# likelihood ratio test: all loadings tested at once:
anova(modelMetric.fit, modelScalar.fit1)

# proceed to residual with GRI3 non-invariant for intercept and residual ============================


#MODEL Residual: Testing equality of residual variances of GRI by student status ====================

residual.syntax = "
#===================================================================================================
#Factor loadings all freely estimated in both groups with label for each group
  GAMBLING =~ c(L1, L1)*GRI1 + c(L3, L3)*GRI3 + c(L5, L5)*GRI5
#===================================================================================================
#Item intercepts all freely estimated in both groups with label for each group
  GRI1 ~ c(I1, I1)*1
  GRI3 ~ c(SI3, NSI3)*1
  GRI5 ~ c(I5, I5)*1
#===================================================================================================
#Redidual variances all freely estimated with label for each group
  GRI1 ~~ c(R1, R1)*GRI1
  GRI3 ~~ c(R3, R3)*GRI3
  GRI5 ~~ c(R5, R5)*GRI5
#===================================================================================================
#Factor variances all freely estimated in one group but not other for identification
  GAMBLING ~~ c(1, NA)*GAMBLING
#===================================================================================================
#Factor means all freely estimated in one group but not other for identification
  GAMBLING ~ c(0, NA)*1
#===================================================================================================  
"

modelResidual.fit = lavaan(
  residual.syntax, 
  data=data01, 
  estimator = "MLR", 
  mimic="Mplus", 
  group = "Student"
)

summary(modelResidual.fit, fit.measures=TRUE, rsquare=TRUE, standardized=TRUE)

# likelihood ratio test: all loadings tested at once:
anova(modelScalar.fit1, modelResidual.fit)

# looks like a problem...investigate MIs
residualMI1 = lavTestScore(modelResidual.fit)

#change label values
labelMap = data.frame(
  lhs = modelResidual.fit@ParTable$plabel,
  parameter = modelResidual.fit@ParTable$label
)

residualMI1$uni = merge(x = residualMI1$uni, y = labelMap, by = "lhs", all.x = TRUE)

# reorder by decreasing values of X2:
residualMI1$uni = residualMI1$uni[order(residualMI1$uni$X2, decreasing = TRUE),]

#restrict to only means shown
residualMI1$uni[grep(x = residualMI1$uni$parameter, pattern = "R"),]

# free parameter for R5 first:

residual.syntax1 = "
#===================================================================================================
#Factor loadings all freely estimated in both groups with label for each group
  GAMBLING =~ c(L1, L1)*GRI1 + c(L3, L3)*GRI3 + c(L5, L5)*GRI5
#===================================================================================================
#Item intercepts all freely estimated in both groups with label for each group
  GRI1 ~ c(I1, I1)*1
  GRI3 ~ c(SI3, NSI3)*1
  GRI5 ~ c(I5, I5)*1
#===================================================================================================
#Redidual variances all freely estimated with label for each group
  GRI1 ~~ c(R1, R1)*GRI1
  GRI3 ~~ c(R3, R3)*GRI3
  GRI5 ~~ c(SR5, NSR5)*GRI5
#===================================================================================================
#Factor variances all freely estimated in one group but not other for identification
  GAMBLING ~~ c(1, NA)*GAMBLING
#===================================================================================================
#Factor means all freely estimated in one group but not other for identification
  GAMBLING ~ c(0, NA)*1
#===================================================================================================  
"

modelResidual.fit1 = lavaan(
  residual.syntax1, 
  data=data01, 
  estimator = "MLR", 
  mimic="Mplus", 
  group = "Student"
)

summary(modelResidual.fit1, fit.measures=TRUE, rsquare=TRUE, standardized=TRUE)

# likelihood ratio test: all loadings tested at once:
anova(modelScalar.fit1, modelResidual.fit1)

# looks like a problem...investigate MIs
residualMI2 = lavTestScore(modelResidual.fit1)

#change label values
labelMap = data.frame(
  lhs = modelResidual.fit1@ParTable$plabel,
  parameter = modelResidual.fit1@ParTable$label
)

residualMI2$uni = merge(x = residualMI2$uni, y = labelMap, by = "lhs", all.x = TRUE)

# reorder by decreasing values of X2:
residualMI2$uni = residualMI2$uni[order(residualMI2$uni$X2, decreasing = TRUE),]

#restrict to only means shown
residualMI2$uni[grep(x = residualMI2$uni$parameter, pattern = "R"),]

# free parameter for R1:

residual.syntax2 = "
#===================================================================================================
#Factor loadings all freely estimated in both groups with label for each group
  GAMBLING =~ c(L1, L1)*GRI1 + c(L3, L3)*GRI3 + c(L5, L5)*GRI5
#===================================================================================================
#Item intercepts all freely estimated in both groups with label for each group
  GRI1 ~ c(I1, I1)*1
  GRI3 ~ c(SI3, NSI3)*1
  GRI5 ~ c(I5, I5)*1
#===================================================================================================
#Redidual variances all freely estimated with label for each group
  GRI1 ~~ c(SR1, NSR1)*GRI1
  GRI3 ~~ c(R3, R3)*GRI3
  GRI5 ~~ c(SR5, NSR5)*GRI5
#===================================================================================================
#Factor variances all freely estimated in one group but not other for identification
  GAMBLING ~~ c(1, NA)*GAMBLING
#===================================================================================================
#Factor means all freely estimated in one group but not other for identification
  GAMBLING ~ c(0, NA)*1
#===================================================================================================  
"

modelResidual.fit2 = lavaan(
  residual.syntax2, 
  data=data01, 
  estimator = "MLR", 
  mimic="Mplus", 
  group = "Student"
)

summary(modelResidual.fit2, fit.measures=TRUE, rsquare=TRUE, standardized=TRUE)

# likelihood ratio test: all loadings tested at once:
anova(modelScalar.fit1, modelResidual.fit2)

# ignore and continue with structural model


#MODEL Structural: Testing equality of mean/variance of GRI by student status ======================

structural.syntax = "
#===================================================================================================
#Factor loadings all freely estimated in both groups with label for each group
  GAMBLING =~ c(L1, L1)*GRI1 + c(L3, L3)*GRI3 + c(L5, L5)*GRI5
#===================================================================================================
#Item intercepts all freely estimated in both groups with label for each group
  GRI1 ~ c(I1, I1)*1
  GRI3 ~ c(SI3, NSI3)*1
  GRI5 ~ c(I5, I5)*1
#===================================================================================================
#Redidual variances all freely estimated with label for each group
  GRI1 ~~ c(SR1, NSR1)*GRI1
  GRI3 ~~ c(R3, R3)*GRI3
  GRI5 ~~ c(SR5, NSR5)*GRI5
#===================================================================================================
#Factor variances fixed to compare against previous model
  GAMBLING ~~ c(1, 1)*GAMBLING
#===================================================================================================
#Factor means fixed to compare against previous model
  GAMBLING ~ c(0,0)*1
#===================================================================================================  
"
modelStructural.fit = lavaan(
  structural.syntax, 
  data=data01, 
  estimator = "MLR", 
  mimic="Mplus", 
  group = "Student"
)

summary(modelStructural.fit, fit.measures=TRUE, rsquare=TRUE, standardized=TRUE)

# likelihood ratio test: all parameters tested at once:
anova(modelResidual.fit2, modelStructural.fit)

# conclusion: difference between students and non-students on GAMBLING factor ======================
summary(modelResidual.fit2, fit.measures=TRUE, rsquare=TRUE, standardized=TRUE)

# comparison models: 
#MODEL 03a: Structural Equation Model #2 --  prediction of GRI items by Student -------------------------------

model03a.syntax = "
  GRI5 ~ Student
  GRI3 ~ Student
  GAMBLING =~ GRI1 + GRI3 + GRI5
  GAMBLING ~ Student
"

model03a.fit = sem(model03a.syntax, data=data01, estimator = "MLR", mimic="Mplus", fixed.x=FALSE)
summary(model03a.fit, fit.measures=TRUE, rsquare=TRUE, standardized=TRUE)


#MODEL 03b: Structural Equation Model #2 -- NO prediction of GRI 3 by Student -------------------------------
model03b.syntax = "
  GRI3 ~ Student
  GAMBLING =~ GRI1 + GRI3 + GRI5
  GAMBLING ~ 0*Student
"

model03b.fit = sem(model03b.syntax, data=data01, estimator = "MLR", mimic="Mplus", fixed.x=FALSE)
summary(model03b.fit, fit.measures=TRUE, rsquare=TRUE, standardized=TRUE)


#MODEL 04a: Model 03b with Standardized Factor -------------------------------
model04a.syntax = "
  GRI3 ~ Student
  GAMBLING =~ GRI1 + GRI3 + GRI5
  GAMBLING ~ 0*Student
"

model04a.fit = sem(model04a.syntax, data=data01, estimator = "MLR", mimic="Mplus", fixed.x=FALSE, std.lv = TRUE)
summary(model04a.fit, fit.measures=TRUE, rsquare=TRUE, standardized=TRUE)


#MODEL 04b: Model 03a with Standardized Gambling Factor  -------------------------------
model04b.syntax = "
  GRI3 ~ Student
  GAMBLING =~ GRI1 + GRI3 + GRI5
  GAMBLING ~ Student
"

model04b.fit = sem(model04b.syntax, data=data01, estimator = "MLR", mimic="Mplus", fixed.x=FALSE, std.lv = TRUE)
summary(model04b.fit, fit.measures=TRUE, rsquare=TRUE, standardized=TRUE)

#MODEL 05: ANALYSIS WITH SUM SCORE INSTEAD OF FACTOR ----------------------------

#creating sum score for GAMBLING 3-item Survey
data02 = data01
data02$GRI135sum = data02$GRI1 + data02$GRI3 + data02$GRI5


model05.syntax = "
  GRI135sum ~ Student 
"
model05.fit = sem(model05.syntax, data=data02, estimator = "MLR", mimic="Mplus", fixed.x=FALSE)
summary(model05.fit, fit.measures=TRUE, rsquare=TRUE, standardized=TRUE)

#getting correct standardizedSolution value (nox as Student is a coded variable):
standardizedSolution(model05.fit, type="std.nox")[1,]

#getting similar standardizedSolution values for model 03 (best fit--but with predictor of GRI3)
standardizedSolution(model03a.fit, type="std.nox")[5,]

#getting similar standardizedSolution values for model 02 (terrible fit, but similar in in that no predictor of GRI3 is part of model)
standardizedSolution(model02.fit, type="std.nox")[4,]


#Model 06: Calculation of Alpha Reliablity with Tau-Equivalent CFA Model -----------------------------
model06.syntax = "
  GAMBLING =~ (loading)*GRI1 + (loading)*GRI3 + (loading)*GRI5
  GRI1 ~~ (U1)*GRI1
  GRI3 ~~ (U3)*GRI3
  GRI5 ~~ (U5)*GRI5

  GCalpha := (3*loading*loading)/( (3*loading*loading) + (U1 + U3 + U5))
"
model06.fit = sem(model06.syntax, data=data02, estimator = "MLR", mimic="Mplus", fixed.x=FALSE, std.lv = TRUE)
summary(model06.fit, fit.measures=TRUE, rsquare=TRUE, standardized=TRUE)

