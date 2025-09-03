# clear R global environment to ensure script runs contiguously from top to bottom
rm(list = ls())

#Set Working Directory to location of data file
#Step 1: Find path of data file and copy path (if using windows)
#Step 2: In the console below, type readClipboard()
#Step 3: Copy and paste R's path to the line below in quotes

#CHANGE BETWEEN QUOTES IN THIS LINE TO REFLECT DIRECTORY OF DATA:
# myDataLocation = "C:\\Dropbox\\!PRE 906\\Lectures\\05 Path Analysis\\R"

#SET WORKING DIRECTORY TO LOCATION OF DATA FILE
# setwd(myDataLocation)

#AUTOMATING PACKAGES NEEDED FOR ANALYSES--------------------------------------------------------------------------------
needed_packages = c("lavaan","semPlot")

for (i in 1:length(needed_packages)){
  haspackage = require(needed_packages[i], character.only = TRUE)
  if (haspackage==FALSE){
    install.packages(needed_packages[i])
    library(needed_packages[i], character.only = TRUE)
  }
}

if (!require("semPlot")){
  oldURL = "https://cran.r-project.org/src/contrib/Archive/semPlot/semPlot_1.1.4.tar.gz"
  install.packages(oldURL, repos=NULL, type="source")
  library(semPlot)
}

#FUNCTIONS USED WITHIN ANALYSIS:--------------------------------------------------------------------------------------
plot_bvn_surface = function(meanvec, covmat, type, xlab, ylab, zlab, main){
  require("fields")
  require("mnormt")
  #meanvec is 2 x 1, covmat is 2 x 2
  
  #creating values for x and y axes based on estimated values from model
  x = matrix(seq(meanvec[1,1]-4*sqrt(covmat[1,1]), meanvec[1,1]+4*sqrt(covmat[1,1]), .1*sqrt(covmat[1,1])), ncol=1)
  y = matrix(seq(meanvec[2,1]-4*sqrt(covmat[2,2]), meanvec[2,1]+4*sqrt(covmat[2,2]), .1*sqrt(covmat[2,2])), ncol=1)
  z = matrix(0,nrow = dim(x)[1],ncol = dim(y)[1])
  
  for (i in 1:dim(x)[1]){
    for (j in 1:dim(y)[1]){
      z[i,j] = dmnorm(c(x[i,1],y[j,1]),mean = t(meanvec), varcov = covmat, log = FALSE)
    }
  }
  
  grid.list=list(x = x, y = y)
  z1 = z[1:dim(x)[1]-1,1:dim(y)[1]-1]
  mygrid = make.surface.grid(grid.list)
  out = list(x = grid.list$x, y = grid.list$y, z = z)
  
  plot.surface(out,type=type, xlab=xlab, ylab=ylab, zlab=zlab, main=main)
}


#Model Examples---------------------------------------------------------------------------------------------------------

#READ IN JOB PERFORMANCE DATA SET: #note: data files with missing noted by . needs option na.strings
math_data = read.csv(file = "mathdata.csv", na.strings=".") 


#Model 0: Empty Model for Two Variables; Unstructured Covariance Matrix ------------------------------------------------

#note: syntax here is overly verbose in order to show all parts of the multivariate model
model00.syntax = "
#Variances:
  perf ~~ perf
  use  ~~ use

#Covariance:
  perf ~~ use

#Means: 
  perf ~ 1
  use  ~ 1
"

#empty multivariate model estimation:
model00.fit = sem(model00.syntax, data=math_data, mimic="MPLUS", fixed.x=TRUE, estimator = "MLR")

#display empty model output
summary(model00.fit, fit.measures=TRUE)

#plot what the model says the data should look like:
model00.estimates = fitted(model00.fit)

#correlation between PERF and USE:
est_cov = matrix(model00.estimates$cov, nrow=2, ncol=2)
corr_mat = solve(sqrt(diag(diag(est_cov))))%*%est_cov%*%solve(sqrt(diag(diag(est_cov))))

#PLOTTING MODEL-ESTIMATED DENSITY
meanvec = matrix(model00.estimates$mean)
covmat = matrix(model00.estimates$cov, nrow=2, ncol=2)
plot_bvn_surface(meanvec = meanvec, covmat = covmat, type = "p", xlab = "PERF",
                 ylab = "USE", zlab = "Relative Frequency", main = "Multivariate Regression Model Estimated Density")

plot_bvn_surface(meanvec = meanvec, covmat = covmat, type = "c", xlab = "PERF",
                 ylab = "USE", zlab = "Relative Frequency", main = "Multivariate Regression Model Estimated Density")

xdata = math_data[is.na(math_data$perf)==FALSE & is.na(math_data$use)==FALSE,]$perf
ydata = math_data[is.na(math_data$perf)==FALSE & is.na(math_data$use)==FALSE,]$use
points(xdata,ydata,pch=16)

#Model 01: Empty Model for Two Variables; Variance Components Covariance Matrix ----------------------------------------

#note: syntax here is overly verbose in order to show all parts of the multivariate model
model01.syntax = "
#Variances:
  perf ~~ (var)*perf
  use  ~~ (var)*use

#Covariance:
  perf ~~ 0*use

#Means: 
  perf ~ 1
  use  ~ 1
"

#empty multivariate model estimation:
model01.fit = sem(model01.syntax, data=math_data, mimic="MPLUS", fixed.x=TRUE, estimator = "MLR")

#display empty model output
summary(model01.fit, fit.measures=TRUE)


#Model 02: Empty Model for Two Variables; Independent Variables Covariance Matrix --------------------------------------

#note: syntax here is overly verbose in order to show all parts of the multivariate model
model02.syntax = "
#Variances:
  perf ~~ perf
  use  ~~ use

#Covariance:
  perf ~~ 0*use

#Means: 
  perf ~ 1
  use  ~ 1
"

#empty multivariate model estimation:
model02.fit = sem(model02.syntax, data=math_data, mimic="MPLUS", fixed.x=TRUE, estimator = "MLR")

#display empty model output
summary(model02.fit, fit.measures=TRUE)


#Model 03: Empty Model for Two Variables; Compound Symmetry Covariance Matrix ------------------------------------------

#note: syntax here is overly verbose in order to show all parts of the multivariate model
model03.syntax = "
#Variances:
  perf ~~ (var)*perf
  use  ~~ (var)*use

#Covariance:
  perf ~~ use

#Means: 
  perf ~ 1
  use  ~ 1
"

#empty multivariate model estimation:
model03.fit = sem(model03.syntax, data=math_data, mimic="MPLUS", fixed.x=TRUE, estimator = "MLR")

#display empty model output
summary(model03.fit, fit.measures=TRUE)

# likelihood ratio tests:
anova(model00.fit, model01.fit)
anova(model00.fit, model02.fit)
anova(model00.fit, model03.fit)

# examining estimates:
model02.estimates = fitted(model02.fit)
model02.estimates$cov

# examining residuals
model02.residuals = residuals(model02.fit, type="raw")
model02.residuals$cov

residuals(model02.fit, type="normalized")$cov

# modification indices
modindices(model02.fit)

#Model 05: Multivariate Model Predicting PERF and USE ---------------------------------------------------------------

model05.syntax = " 
#endogenous variable equations
  perf ~ hsl + cc
  use ~ hsl + cc

#endogenous variable intercepts
  perf ~ 1
  use ~ 1

#endogenous variable residual variances
  perf ~~ perf
  use ~~ use

#endogenous variable residual covariances
  perf~~ use

#exogeneous variables put into likelihood function:

#exogeneous means(intercepts)
hsl ~ 1
cc ~ 1

#exogeneous variances
hsl ~~ hsl
cc ~~ cc

#exogeneous covariances
hsl ~~ cc
"

#model estimation
model05.fit = sem(model05.syntax, data=math_data, estimator = "MLR", mimic="Mplus", fixed.x=FALSE)

#display model output
summary(model05.fit, fit.measures=TRUE, standardized=TRUE)
standardizedSolution(model05.fit, type = "std.all")

#r-squared for performance:
inspect(model05.fit, what="r2")

#Model 02: Full path model #1 ---------------------------------------------------------------


model06.syntax = " 
#endogenous variable equations
  perf ~ hsl + msc + mse
  use  ~ mse
  mse  ~ hsl + cc + gender 
  msc  ~ mse + cc + hsl
  cc   ~ hsl
  hsl  ~ gender

#endogenous variable intercepts
  perf ~ 1
  use  ~ 1
  mse  ~ 1
  msc  ~ 1
  cc   ~ 1
  hsl  ~ 1

#endogenous variable residual variances
  perf ~~ perf
  use  ~~ use
  mse  ~~ mse
  msc  ~~ msc
  cc   ~~ cc  
  hsl  ~~ hsl
  
#endogenous variable residual covariances
  #none specfied in the original model so these have zeros:
  perf ~~ 0*use + 0*mse + 0*msc + 0*cc + 0*hsl
  use  ~~ 0*mse + 0*msc + 0*cc + 0*hsl
  mse  ~~ 0*msc + 0*cc + 0*hsl
  msc  ~~ 0*cc + 0*hsl
  cc   ~~ 0*hsl

#exogeneous variables put into likelihood function:

  #means(intercepts)
  gender ~ 1

  #variances
  gender ~~ gender

"

#model estimation
model06.fit = sem(model06.syntax, data=math_data, estimator = "MLR", mimic="Mplus", fixed.x=FALSE)

#display model output
summary(model06.fit, fit.measures=TRUE, standardized=TRUE)

#display normalized residuals to inspect local model misfit
residuals(model06.fit, type="normalized")

#display modification indices
modindices(model06.fit)


#Model 03: Full path model  ---------------------------------------------------------------

model07.syntax = " 
#endogenous variable equations
  perf ~ hsl + msc + mse
  use  ~ mse
  mse  ~ hsl + cc + gender 
  msc  ~ mse + cc + hsl
  cc   ~ hsl
  hsl  ~ gender

#endogenous variable intercepts
  perf ~ 1
  use  ~ 1
  mse  ~ 1
  msc  ~ 1
  cc   ~ 1
  hsl  ~ 1

#endogenous variable residual variances
  perf ~~ perf
  use  ~~ use
  mse  ~~ mse
  msc  ~~ msc
  cc   ~~ cc  
  hsl  ~~ hsl
  
#endogenous variable residual covariances
  use ~~ msc
  perf ~~ 0*use + 0*mse + 0*msc + 0*cc + 0*hsl
  use  ~~ 0*mse + 0*cc + 0*hsl
  mse  ~~ 0*msc + 0*cc + 0*hsl
  msc  ~~ 0*cc + 0*hsl
  cc   ~~ 0*hsl

#exogeneous variables put into likelihood function:

  #means(intercepts)
  gender ~ 1

  #variances
  gender ~~ gender

"

#model estimation
model07.fit = sem(model07.syntax, data=math_data, estimator = "MLR", mimic="Mplus", fixed.x=FALSE)

#display model output
summary(model07.fit, fit.measures=TRUE, standardized=TRUE)

#display no.x standardized results for gender
standardizedSolution(model07.fit, type="std.nox")

#display r2
inspect(model07.fit, what="r2")


#Model 08: Adding indirect effects to full path model ---------------------------------------------------------------

model08.syntax = " 
#endogenous variable equations
  perf ~ hsl + msc + mse
  use  ~ mse
  mse  ~ b_hsl_mse*hsl + b_cc_mse*cc + gender 
  msc  ~ mse + cc + hsl
  cc   ~ b_hsl_cc*hsl
  hsl  ~ gender

#endogenous variable intercepts
  perf ~ 1
  use  ~ 1
  mse  ~ 1
  msc  ~ 1
  cc   ~ 1
  hsl  ~ 1

#endogenous variable residual variances
  perf ~~ perf
  use  ~~ use
  mse  ~~ mse
  msc  ~~ msc
  cc   ~~ cc  
  hsl  ~~ hsl
  
#endogenous variable residual covariances
  use ~~ msc
  perf ~~ 0*use + 0*mse + 0*msc + 0*cc + 0*hsl
  use  ~~ 0*mse + 0*cc + 0*hsl
  mse  ~~ 0*msc + 0*cc + 0*hsl
  msc  ~~ 0*cc + 0*hsl
  cc   ~~ 0*hsl

#exogeneous variables put into likelihood function:

  #means(intercepts)
  gender ~ 1

  #variances
  gender ~~ gender

#indirect effect of interest:
  ind_hsl_mse := b_hsl_cc*b_cc_mse

#total effect of interest:
  tot_hsl_mse := b_hsl_mse + (b_hsl_cc*b_cc_mse)
"

#as lavaan automatically bootstraps p-value, set the random seed for same results
set.seed(007)

#model estimation
model04.fit = sem(model04.syntax, data=math_data, estimator = "MLR", mimic="Mplus", fixed.x=FALSE)

#display model output
summary(model04.fit, fit.measures=TRUE, standardized=TRUE)

#plot path diagram with standardized coefficients
semPaths(model04.fit,intercepts = FALSE, residuals = TRUE, style="mx", layout="spring", rotation=1,
         optimizeLatRes = TRUE, whatLabels="std")
