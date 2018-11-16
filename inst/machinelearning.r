library(twoplane)
# requires HiddenMarkov, functional, Rcpp, boot, expm, Matrix
library(mvtnorm)
# library(devtools); install_github("b-steve/nspp")
library(nspp)

library(caret)
library(caretEnsemble)

#  Simulate for training and testing of ML model.
#  Fixed parameters sigma, E1, D in testsims.r are replaced with loops here. 
seed=654321 # arbitrary fixed seed - for generation of exactly same set of data on separate occasions
set.seed(seed) # initialise random number sequence (for repeatability)

# Number of buckets.
for(bincount in c(160, 170, 180, 190, 200, 210, 220, 230, 240)) {
#  Set up data frame for ML input. Last column is the variable to predict. 
allsim=as.data.frame(matrix(0, ncol = bincount+2, nrow = 0))

nm2km=1.852 # multiplier to convert nautical miles to kilometres

E2=110  # expected length unavailable period (from Hiby&Lovell, 1998) in seconds
#E1=28   # expected length available period (from Hiby&Lovell, 1998) in seconds

Nsim=1000

#   In each simulation values are selected for E1, N and sigma with uniform probability in some range.
for(cursim in 1:Nsim) {
  E1=runif(1,18,36)
  
  E=c(E1,E2)
  Ec=sum(E) # mean dive cycle length in seconds
  p.up=E1/Ec # proportion of time up
  
  #N=500 # number animals in strip
  N=runif(1, 300, 750)
  
  L=350;w=0.3 # length and width of strip in km
  planeknots=100 # observer speed in knots
  planespd=planeknots*nm2km/(60^2) # observer speed in km/sec
  
  k=10 # time separation of observers in seconds
  
  # Set diffusion coefficient in 1D (f(x,t)=N(x; mu=0,var=sigma^2*k)), so that about 95%
  # of animals have moved a distance of no more than 2*sigma*sqrt(k) nm in time k
  #sigmaknots=2.5
  sigmaknots=runif(1, 1.5, 3.3)
  
  sigmarate=sigmaknots*nm2km/60/60 # std dev of movement distance per second (in m)
  sigma=sigmarate*sqrt(k) # std dev of movement distance after k seconds (in m)
  sigma*1000*2 # 95% within this many m of start point after k seconds
  sigma.mult=6 # multiplier to set bound for maximum distance moved
  
  p.see.up=c(1,1) # prob(see|up) for each observer
  
  # prob up in k seconds, given up now:
  p.k=calc.p.avail(k,E1,Ec,p=c(1,1),sigma,planespd);p.k
  round(p.k/p.up*100,1)-100 # % inflation in detection prob due to recent availability
  
  #halfw=(w/2)/planespd # strip half width in planespd units
  Dstrip=N/(L*w) # density in number per sq km
  Dstrip.t=Dstrip*(planespd^2) # density in planespd units
  D.line.t=Dstrip.t*w/planespd # density in planespd along LINE units (1-dimensional)
  Dbound=list(lower=-5*abs(log(D.line.t)),upper=5*abs(log(D.line.t))) # need bounds only if doing 1-dim estimation
  method="Nelder-Mead" # this is ignored if doing 1-dim estimation
  control.opt=list(trace=0,maxit=1000)
  
  plot.sample=FALSE
  plot.displacement=FALSE
  plot.cuts=FALSE
  checkdists=FALSE
  fromfile=FALSE
  segiotime=segtime=nspptime=rep(NA,Nsim)
  
  true=list(D=Dstrip,sigma=sigmarate,E1=E1) # parameters to use in simulation
  
  skip=c()
  
  sdat=sim.2plane(N,L,w,sigmarate,k,planespd,p.up,Ec,p=p.see.up,
                  sigma.mult=sigma.mult,movement=list(forward=TRUE,sideways=TRUE,correct=FALSE),
                  cdfnormvt=NULL)
  
  #  For a simulated data set, take set of all differences between pairs of observations (losing distinction between planes)
  obs=sort(c(sdat$s1, sdat$s2))
  diffs=sort(c(dist(obs)))
  
  binsize=(L/planespd)/bincount
  
  freqs=rep(0,bincount)
  for(i in 1:bincount) {
    freqs[i]=length(which( diffs>=(i-1)*binsize & diffs<i*binsize ))
  }
  
  # Put into dataframe for training.
  numobs=length(sdat$s1)+length(sdat$s2)
  
  freqsnorm=freqs
  #  freqsnorm=freqs/length(diffs)   normalise so buckets sum to 1. 
  
  allsim<-rbind(allsim, c(freqsnorm, numobs, sigma, E1, Dstrip))
}

colnames(allsim)=c(seq(1,bincount), "numobs", "sigma", "E1", "Dstrip")

inTrain <- createDataPartition(y = allsim$Dstrip, p=0.75, list=FALSE)

training <- allsim[ inTrain,]
testing  <- allsim[-inTrain,]

# Predict Dstrip

training$sigma=NULL
training$E1=NULL
testing$sigma=NULL
testing$E1=NULL

rfModel <- train(Dstrip ~ ., data = training, method = "parRF")
rfPredict <- predict(rfModel, newdata = testing)
cat(bincount, " rmse:", mean((rfPredict-testing$Dstrip)^2)^0.5, " mae:", mean(abs(rfPredict-testing$Dstrip)), " bias: ", (mean(rfPredict)/mean(testing$Dstrip))*100, "\n")
}

saveRDS(training, "training-ml-200buckets.rds")
saveRDS(testing, "testing-ml-200buckets.rds")

training <- readRDS("training-ml2.rds")
testing <- readRDS("testing-ml2.rds")

# Predict Dstrip

training$sigma=NULL
training$E1=NULL
testing$sigma=NULL
testing$E1=NULL

#########################################
#
#   A simple linear regression.

linearModel <- train(Dstrip ~ ., data = training, method = "lm")

lmPredict <- predict(linearModel, newdata = testing)

# Root mean square error.  0.7629109
mean((lmPredict-testing$Dstrip)^2)^0.5
# MAE   0.622036
mean(abs(lmPredict-testing$Dstrip))
(mean(lmPredict)/mean(testing$Dstrip))*100

##########################################
#
#   Fit a linear SVM

svmLinearModel <- train(Dstrip ~ ., data = training, method = "svmLinear")

svmLinearPredict <- predict(svmLinearModel, newdata = testing)

# Root mean square error.
mean((svmLinearPredict-testing$Dstrip)^2)^0.5
mean(abs(svmLinearPredict-testing$Dstrip))
(mean(svmLinearPredict)/mean(testing$Dstrip))*100



##########################################
#
#   Fit a SVM with polynomial kernel

svmPolyModel <- train(Dstrip ~ ., data = training, method = "svmPoly")

svmPolyPredict <- predict(svmPolyModel, newdata = testing)

# Root mean square error.
mean((svmPolyPredict-testing$Dstrip)^2)^0.5
# [1] 0.7261662
mean(abs(svmPolyPredict-testing$Dstrip))
# [1] 0.5772985
#  Bias
(mean(svmPolyPredict)/mean(testing$Dstrip))*100

###########################################
#
#   SVM Radial kernel

svmRadialModel <- train(Dstrip ~ ., data = training, method = "svmRadial")
svmRadialPredict <- predict(svmRadialModel, newdata = testing)
# Root mean square error.
mean((svmRadialPredict-testing$Dstrip)^2)^0.5
mean(abs(svmRadialPredict-testing$Dstrip))
#  Bias
(mean(svmRadialPredict)/mean(testing$Dstrip))*100


########
#   Random forest model

rfModel <- train(Dstrip ~ ., data = training, method = "parRF")

rfPredict <- predict(rfModel, newdata = testing)

# Root mean square error.
mean((rfPredict-testing$Dstrip)^2)^0.5
# 100 buckets [1] 0.6964453
mean(abs(rfPredict-testing$Dstrip))
# 100 buckets [1] 0.5576064
#  Bias
(mean(rfPredict)/mean(testing$Dstrip))*100

# 100  rmse: 0.7027615  mae: 0.5620901  bias:  100.9723 
# 150  rmse: 0.6804895  mae: 0.5570778  bias:  99.84452 
# 200  rmse: 0.6804317  mae: 0.5253935  bias:  100.0988 
# 250  rmse: 0.695125   mae: 0.5509367  bias:  99.49994 

# 100  rmse: 0.7325592  mae: 0.5947724  bias:  100.3819 
# 150  rmse: 0.7001101  mae: 0.5815442  bias:  99.73506 
# 200  rmse: 0.6458036  mae: 0.5238597  bias:  100.5008 
# 250  rmse: 0.7028836  mae: 0.5652548  bias:  100.428 

# 160  rmse: 0.7163458  mae: 0.5827304  bias:  100.7418 
# 170  rmse: 0.7441702  mae: 0.6054921  bias:  100.1989 
# 180  rmse: 0.6720945  mae: 0.5443028  bias:  101.0459 
# 190  rmse: 0.6807437  mae: 0.5359679  bias:  100.0277 
# 200  rmse: 0.6894268  mae: 0.5512045  bias:  99.76174 
# 210  rmse: 0.6824541  mae: 0.5487435  bias:  100.4503 
# 220  rmse: 0.67356    mae: 0.5447459  bias:  102.3329 
# 230  rmse: 0.7407297  mae: 0.5824032  bias:  100.1861 
# 240  rmse: 0.6869024  mae: 0.5481898  bias:  100.6676 


#############################################
#
#   Fit an ensemble model


#  Models that fit in a short time (less than about a minute):

mlmethods=c("enet", "foba", "icr", "lars", "lars2", "lasso", 
            "leapBackward", "leapForward", "leapSeq", "lm", 
            "nnls", "pcr", "ppr", "relaxo", "ridge", "rqlasso", 
            "rqnc", "superpc",
            
            
            "svmLinear", "svmPoly"
            )

#  Reached end of regression only list. 

#  Packages produced errors:  RWeka, neuralnet, rPython, rlm, rvmLinear, SBC, spikeslab, 

ensembleList <- caretList(Dstrip~., data=training, methodList=mlmethods)

ensembleStack <- caretStack(ensembleList, method="glm")

ensemblePredict <- predict(ensembleStack, newdata = testing)

# Root mean square error.
mean((ensemblePredict-testing$Dstrip)^2)^0.5

(mean(ensemblePredict)/mean(testing$Dstrip))*100

#  Notes from Kasim:  try gaussian or radial basis function kernel for SVM because 
#  the pseudolikelihood is convex -- suggested kernels should fit a convex function.


#  Upper bound on MAE/RMSE. Model is just the mean of training data. 
mean(abs(mean(training$Dstrip)-testing$Dstrip))
mean((mean(training$Dstrip)-testing$Dstrip)^2)^0.5
(mean(training$Dstrip)/mean(testing$Dstrip))*100

#  Baseline model -- learn function from number of observations to density

train2 <- data.frame(training$Dstrip, training$numobs)
test2 <- data.frame(testing$Dstrip, testing$numobs)
colnames(train2)<-c("Dstrip","numobs")
colnames(test2)<-c("Dstrip","numobs")

svmBaseModel <- train(Dstrip ~ numobs, data = train2, method = "svmPoly")
svmBasePredict <- predict(svmBaseModel, newdata = test2)
# Root mean square error and MAE
mean((svmBasePredict-test2$Dstrip)^2)^0.5
mean(abs(svmBasePredict-test2$Dstrip))
#  Bias
(mean(svmBasePredict)/mean(test2$Dstrip))*100

#   Baseline with random forest instead of polynomial SVM.
train2 <- data.frame(training$Dstrip, training$numobs)
test2 <- data.frame(testing$Dstrip, testing$numobs)
colnames(train2)<-c("Dstrip","numobs")
colnames(test2)<-c("Dstrip","numobs")

rfBaseModel <- train(Dstrip ~ numobs, data = train2, method = "parRF")
rfBasePredict <- predict(rfBaseModel, newdata = test2)
# Root mean square error and MAE
mean((rfBasePredict-test2$Dstrip)^2)^0.5
mean(abs(rfBasePredict-test2$Dstrip))
# [1] 0.6362949
#  Bias
(mean(rfBasePredict)/mean(test2$Dstrip))*100

