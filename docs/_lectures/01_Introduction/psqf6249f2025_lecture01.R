#AUTOMATING PACKAGES NEEDED FOR ANALYSES--------------------------------------------------------------------
needed_packages = c("lavaan","semPlot")

for (i in 1:length(needed_packages)){
  haspackage = require(needed_packages[i], character.only = TRUE)
  if (haspackage==FALSE){
    install.packages(needed_packages[i])
    library(needed_packages[i], character.only = TRUE)
  }
}


#DATA SIMULATION FUNCTION: SIMPLIFIES SIMULATING DATA FOR EXAMPLE------------------------------------------

simdata_type1 = function(n, prop1, nitems_factor1, beta, factorvariance1, factorloading1, itemintercept1, uniquevar1){
  #data are generated with one-factor model and without missing data
  datamatrix = matrix(NA, nrow = n, ncol=nitems_factor1+2)
  
  datamatrix[,1] = rbinom(n, size=1, prob=prop1)
  datamatrix[,2] = -beta + beta*datamatrix[,1] + rnorm(n,mean=0, sd = sqrt(factorvariance1))
  
  for (i in 1:nitems_factor1){
    datamatrix[,i+2] = rnorm(n, mean=itemintercept1+factorloading1*datamatrix[,2], sd=sqrt(uniquevar1))
  }
  
  sumscore = apply(X=datamatrix[,3:(nitems_factor1+2)], MARGIN=1, FUN=sum)
  
  sum_score_reliability = ((nitems_factor1*factorloading1)^2)*factorvariance1
  sum_score_reliability = sum_score_reliability/(sum_score_reliability+(nitems_factor1*uniquevar1))
  
  simdata_type1 = list(data = data.frame(married = datamatrix[,1], happiness_sumscore = sumscore, datamatrix[,3:(nitems_factor1+2)]),
                       reliability = sum_score_reliability)
  
}

#SIMULATING DATA-----------------------------------------------------------------------------------
#SETTING RANDOM NUMBER SEED SO SIMUATIONS WILL HAVE SAME VALUES AS IN LECTURE
random_seed = 201590601 
set.seed(random_seed)

data01 = simdata_type1(n = 150, prop1 = .5, nitems_factor1 = 5, beta = .5, factorvariance1 = 1, 
                 factorloading1 = .25, itemintercept1 = 3, uniquevar1 = 1)
View(data01$data)
#report sum score reliability:
data01$data$id = rownames(data01$data)

#output data set to file (step is un-needed but is given to show how to export/import data:
write.csv(data01$data, file="sem25psqf6249_lecture01.csv", row.names=FALSE)

#importing data again
data02 = read.csv(file="sem25psqf6249_lecture01.csv")

#ANALYSIS MODEL #1: [W/O] SEM MARITAL STATUS PREDICTING A SUM SCORE FOR HAPPINESS------------------------------

#empty model syntax:
model00.syntax = "happiness_sumscore ~~ happiness_sumscore"

#empty model estimation
model00.fit = sem(model00.syntax, data=data02, mimic="MPLUS", fixed.x=TRUE)

#display empty model output
summary(model00.fit, fit.measures=TRUE)

#compare output to sample statistics:
#--mean is identical to intercept
mean(data02$happiness_sumscore)

#variance is not quite (lavaan is ML (divides by N)/var() function is unbiased (divides by N-1))
var(data02$happiness_sumscore)

#convert unbiased (N-1) variance to N version and get the same result:
(var(data02$happiness_sumscore)*(dim(data02)[1]-1))/(dim(data02)[1])

#now, add predictor of happiness...

#model 1 syntax
model01.syntax = "happiness_sumscore ~ married"

#model 1 estimation
model01.fit = sem(model01.syntax, data=data02, mimic="MPLUS", fixed.x=TRUE)

#display model 1 output
summary(model01.fit, fit.measures=TRUE)

#display standardized parameter output
standardizedSolution(model01.fit, type="std.nox")

#plot path diagram with standardized coefficients
semPaths(model01.fit,intercepts = FALSE, residuals = TRUE, style="mx", layout="tree", rotation=1,
         optimizeLatRes = TRUE, whatLabels="std")

#EXAMINING THE HAPPINESS ITEMS--------------------------------------------------------------------------------

#correlation matrix of 5 happiness items:
cor(data02[,3:7])

#covariance matrix of 5 happiness items:
cov(data02[,3:7])

#sum of the item variances:
itemvar = sum(diag(cov(data02[,3:7])))

#sum of the entire covariance matrix
totalvar = sum(cov(data02[,3:7]))

#ANALYSIS MODEL #2: [W/ SEM] MARITAL STATUS PREDICTING A SUM SCORE FOR HAPPINESS------------------------------

#model 2 syntax
model02.syntax = "happiness =~ X1 + X2 + X3 + X4 + X5
                  happiness ~ married"

#model 2 estimation
model02.fit = sem(model02.syntax, data=data02, mimic="MPLUS", fixed.x=TRUE)

#display model 2 output
summary(model02.fit, fit.measures=TRUE)

#display standardized parameter output
standardizedSolution(model02.fit)

#plot path diagram with standardized coefficients
semPaths(model02.fit,intercepts = FALSE, residuals = TRUE, style="mx", layout="tree", rotation=1,
         optimizeLatRes = TRUE, whatLabels="std")

