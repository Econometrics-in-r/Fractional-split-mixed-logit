rm(list = ls())
require(stats4)
require(maxLik)
require(randtoolbox)

# reading and storing data in a dataframe
dataset <- read.csv(file.choose(),header=T)
#Sample Size
N <- nrow(dataset) 

#Fractions of the dependent variable (for two categories; adjust as needed)
w1 <- dataset[,'fraction']#(fraction of speeding category)
w2 = 1-w1#(fraction of non-speeding category)

# Halton Draws 
preparedraws=function()
{
  d=1
  while(d<(length(normaldraws)+1))
  {
    draws1[,normaldraws[d]]<<- qnorm(draws1[,normaldraws[d]])
    d=d+1
  }
}

Ndraws=200      # set number of draws 
dimensions=2    # define number of random parameters in the model

# generate draws (using Halton)
draws1=as.matrix(halton(Ndraws*N,dimensions))

# assign names to individual sets of draws - need one entry per dimension
colnames(draws1)=c("HRbeta1","HRbeta2")
# define whether any draws should be transformed to Normals, which is also needed for e.g. lognormals (leave empty if not)
normaldraws=c("HRbeta1","HRbeta2")

# preparing draws for estimation - this may take a while
preparedraws()

# accounting for panel setting
#fixing parameters across grouped observations i.e. grouped random parameters
block = length(unique(dataset[,'ID']))
ngroup = length(unique(dataset[,'Group']))
for (i in 1:Ndraws){
  tempInd = ((i-1)*block*ngroup) + (1:block)
  for (ii in 2:ngroup){
    draws1[tempInd+(ii-1)*block,] = draws1[tempInd,]
  }
}

## data preparation
# separating the variables with fixed parameters 
dataF =  as.matrix(data.frame(1,log(dataset$Length)))
# separating the variables with random parameters 
dataR = as.matrix(data.frame(log(dataset$AADT),dataset$LWIDTH))

dataR2=NULL
for(i in 1:Ndraws){
  dataR2=rbind(dataR2,dataR)
}

draws1 = draws1[,1:dimensions]

# Likelihood function
LL <- function(params){  
  Fbeta <- params[1:2] # Fixed parameters in the mean Function (excluding the constant)
  MRbeta <- params[3:4]  # Mean of Random parameters in the mean function
  SDRbeta <- params[5:6]  # Std of Random parameters in the mean function
  
  # vector of indipendent variables with fixed parameters
  offset = rep.int(dataF%*%as.matrix(Fbeta,ncol=1),Ndraws)
  # simulating random parameters from their means and standard deviation
  beta = t( t(draws1)*SDRbeta + MRbeta )
  # constructing the utility function of the chosen alternative
  u <- offset+rowSums(dataR2*beta)
  
  # probability of the chosen alternative
  prob1 <- exp(u)/(1+exp(u))
  
  # simulated loglikelihood for fractional split model
  PR1 <- rowMeans(matrix(prob1, ncol = Ndraws))
  PR2 <- 1-PR1
  loglik <- sum(w1*log(PR1)+w2*log(PR2))
  
  return(loglik)
}

# initial values for optimization
init <- c(0,0.1,#fixed parameters
          0.1,0.1,#mean of random parameters
          0.1,0.1)#standard deviation of random parameters

# optimization (maximization of likelihood function)
fit1 <- maxLik(LL,start=init,method="BFGS")

summary(fit1)


