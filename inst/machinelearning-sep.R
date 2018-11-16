library(twoplane)
# requires HiddenMarkov, functional, Rcpp, boot, expm, Matrix
library(mvtnorm)
# library(devtools); install_github("b-steve/nspp")
library(nspp)

library(caret)
library(caretEnsemble)

# Set up parallel processing
library(doMC)
registerDoMC(cores=2)

#  Some functions to support ML

maeFunc=function(data, lev=NULL, model=NULL) { 
  mae=mean(abs(data$obs-data$pred))
  names(mae) <- "mae"
  mae
}

observationsToBins<-function(sdat,separate=FALSE, L, planespd, bincount) {
  if(!separate) {
    #  For a simulated data set, take set of all differences between pairs of observations (losing distinction between planes)
    obs=sort(c(sdat$s1, sdat$s2))
    diffs=sort(c(dist(obs)))
    
    binsize=(L/planespd)/bincount
    
    freqs=rep(0,bincount)
    for(i in 1:bincount) {
      freqs[i]=length(which( diffs>=(i-1)*binsize & diffs<i*binsize ))
    }
    return(freqs)
  }
  else {
    #  Separate non-recapture pairs from possible recapture pairs. 
    #  Non-recapture pairs first. 
    diffs=sort(c(dist(sdat$s1),dist(sdat$s2)))
    
    binsize=(L/planespd)/bincount
    
    freqsnonr=rep(0,bincount)
    for(i in 1:bincount) {
      freqsnonr[i]=length(which( diffs>=(i-1)*binsize & diffs<i*binsize ))
    }
    
    # Possible recap pairs.
    diffs=c()
    for(i in 1:length(sdat$s1)) {
      diffs=c(diffs, abs(sdat$s1[i]-sdat$s2))
    }
    diffs=sort(diffs)
    freqspossr=rep(0,bincount)
    for(i in 1:bincount) {
      freqspossr[i]=length(which( diffs>=(i-1)*binsize & diffs<i*binsize ))
    }
    return(c(freqsnonr, freqspossr))
  }
}


simulateTestInstance<-function(N, E1, sigmaknots, separate) {
  E2=110  # expected length unavailable period (from Hiby&Lovell, 1998) in seconds
  #E1=28   # expected length available period (from Hiby&Lovell, 1998) in seconds
  
  E=c(E1,E2)
  Ec=sum(E) # mean dive cycle length in seconds
  p.up=E1/Ec # proportion of time up
  
  L=350;w=0.3 # length and width of strip in km
  planeknots=100 # observer speed in knots
  planespd=planeknots*nm2km/(60^2) # observer speed in km/sec
  
  k=10 # time separation of observers in seconds
  
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
  
  sdat=sim.2plane(N,L,w,sigmarate,k,planespd,p.up,Ec,p=p.see.up,
                  sigma.mult=sigma.mult,movement=list(forward=TRUE,sideways=TRUE,correct=FALSE),
                  cdfnormvt=NULL)
  
  bins=observationsToBins(sdat, separate, L, planespd, bincount)
  numobs=length(sdat$s1)+length(sdat$s2)
  if(!separate) {
    kludge=bins[1]/(mean(bins[2:11]))
  }
  else {
    kludge=bins[bincount+1]/(mean(bins[(bincount+2):(bincount+11)]))
  }
  
  return(c(bins, numobs, kludge, sigma, E1, Dstrip))
}

evalModel<-function(m, Nsim, bincount, separate) {
  out=data.frame()
  #   27 points
  for(N in c(350,500,650)) {
    for(E1 in c(22,28,34)) {
      for(sigmaknots in c(1.8,2.5, 3.2)) {
        # Simulate test instances
        allsim=data.frame()
        for(i in 1:Nsim) {
          allsim<-rbind(allsim,simulateTestInstance(N, E1, sigmaknots, separate))
        }
        colnames(allsim)<-c(seq(bincount*(separate+1)), "numobs", "kludge", "sigma", "E1", "Dstrip")
        # Evaluate m
        
        rfPredict <- predict(m, newdata = allsim)
        cat("RF model, bincount: ", bincount, " rmse:", mean((rfPredict-allsim$Dstrip)^2)^0.5, " mae:", mean(abs(rfPredict-allsim$Dstrip)), " bias: ", (mean(rfPredict)/mean(allsim$Dstrip))*100-100, "\n")
        
        RMSE=mean((rfPredict-allsim$Dstrip)^2)^0.5
        MAE=mean(abs(rfPredict-allsim$Dstrip))
        bias=(mean(rfPredict)/mean(allsim$Dstrip))*100-100
        out<-rbind(out, c(N, E1, sigmaknots, RMSE, MAE, bias))
      }
    }
  }
  colnames(out) <- c("N","E1","sigmaknots","RMSE","MAE","Bias")
  return(out)
}



#  Simulate for training and testing of ML model.
#  Fixed parameters sigma, E1, D in testsims.r are replaced with loops here. 
seed=654321 # arbitrary fixed seed - for generation of exactly same set of data on separate occasions
set.seed(seed) # initialise random number sequence (for repeatability)

separate=TRUE  #  Separate possible recapture pairs and known non-recapture pairs.  

# Number of buckets.
for(bincount in c(160, 170, 180, 190, 200, 210, 220, 230, 240)) {
#  Set up data frame for ML input. Last column is the variable to predict. 
numcols=(bincount*(separate+1))+5
allsim=as.data.frame(matrix(0, ncol = numcols, nrow = 0))

nm2km=1.852 # multiplier to convert nautical miles to kilometres

Nsim=1000

#   In each simulation values are selected for E1, N and sigma with uniform probability in some range.
for(cursim in 1:Nsim) {
  E1=runif(1,18,36)
  
  #N=500 # number animals in strip
  N=runif(1, 300, 750)
  
  # Set diffusion coefficient in 1D (f(x,t)=N(x; mu=0,var=sigma^2*k)), so that about 95%
  # of animals have moved a distance of no more than 2*sigma*sqrt(k) nm in time k
  #sigmaknots=2.5
  sigmaknots=runif(1, 1.5, 3.3)
  
  allsim<-rbind(allsim, simulateTestInstance(N, E1, sigmaknots, separate))
}

colnames(allsim)=c(seq(bincount*(separate+1)), "numobs", "kludge", "sigma", "E1", "Dstrip")

inTrain <- createDataPartition(y = allsim$Dstrip, p=0.75, list=FALSE)

training <- allsim[ inTrain,]
testing  <- allsim[-inTrain,]

#saveRDS(training, "training-mlsep.rds")
#saveRDS(testing, "testing-mlsep.rds")

#training <- readRDS("training-mlsep.rds")
#testing <- readRDS("testing-mlsep.rds")

# Predict Dstrip

rfModel <- train(Dstrip ~ ., data = training[,c(seq(bincount*(separate+1)),"numobs","Dstrip")], method = "parRF")
rfPredict <- predict(rfModel, newdata = testing)
cat("RF model, bincount: ", bincount, " rmse:", mean((rfPredict-testing$Dstrip)^2)^0.5, " mae:", mean(abs(rfPredict-testing$Dstrip)), " bias: ", (mean(rfPredict)/mean(testing$Dstrip))*100, "\n")

svmPolyModel <- train(Dstrip ~ ., data = training, method = "svmPoly")
svmPolyPredict <- predict(svmPolyModel, newdata = testing)
cat("SVMPoly model, bincount: ", bincount, " rmse:", mean((svmPolyPredict-testing$Dstrip)^2)^0.5, " mae:", mean(abs(svmPolyPredict-testing$Dstrip)), " bias: ", (mean(svmPolyPredict)/mean(testing$Dstrip))*100, "\n")
}

#########################################
#
#   A simple linear regression.

linearModel <- train(Dstrip ~ ., data = training, method = "lm")

lmPredict <- predict(linearModel, newdata = testing)

# Root mean square error.  [1] 1.008075
mean((lmPredict-testing$Dstrip)^2)^0.5
# MAE  0.7870068
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
# [1] 0.6462105
mean(abs(svmPolyPredict-testing$Dstrip))
# [1] 0.5207277
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

rfModel <- train(Dstrip ~., data = training[,c(seq(bincount*(separate+1)),"numobs","Dstrip")], method = "parRF")
rfPredict <- predict(rfModel, newdata = testing)
# Root mean square error.
mean((rfPredict-testing$Dstrip)^2)^0.5
#  [1] 0.6404903
mean(abs(rfPredict-testing$Dstrip))
#  [1] 0.5252377
#  Bias
(mean(rfPredict)/mean(testing$Dstrip))*100


#  separate=FALSE:

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

# separate=TRUE:
# 160  rmse: 0.6969134  mae: 0.5582303  bias:  99.79475 
# 170  rmse: 0.6902852  mae: 0.5599847  bias:  99.12903 
# 180  rmse: 0.7249718  mae: 0.5789915  bias:  100.0074 
# 190  rmse: 0.7055481  mae: 0.5761494  bias:  101.7178 
# 200  rmse: 0.7153248  mae: 0.5484245  bias:  98.93955 
# 210  rmse: 0.7143464  mae: 0.5818792  bias:  99.88956 
# 220  rmse: 0.6823385  mae: 0.5676768  bias:  100.144 
# 230  rmse: 0.6893109  mae: 0.5565844  bias:  100.8707 
# 240  rmse: 0.6400187  mae: 0.5242953  bias:  101.0098

# separate=TRUE:
# RF model, bincount:  160  rmse: 0.6969134  mae: 0.5582303  bias:  99.79475 
# SVMPoly model, bincount:  160  rmse: 0.725116  mae: 0.577099  bias:  98.35131 
# RF model, bincount:  170  rmse: 0.7164821  mae: 0.5947566  bias:  100.6603 
# SVMPoly model, bincount:  170  rmse: 0.7208282  mae: 0.5940402  bias:  99.93745 
# RF model, bincount:  180  rmse: 0.6259714  mae: 0.5092271  bias:  100.8095 
# SVMPoly model, bincount:  180  rmse: 0.6412249  mae: 0.5143678  bias:  98.28491 
# RF model, bincount:  190  rmse: 0.745025  mae: 0.6030331  bias:  99.70847 
# SVMPoly model, bincount:  190  rmse: 0.7513083  mae: 0.6047218  bias:  98.68249 
# RF model, bincount:  200  rmse: 0.6975202  mae: 0.5745566  bias:  100.3148 
# SVMPoly model, bincount:  200  rmse: 0.7082026  mae: 0.5808597  bias:  99.00322 
# RF model, bincount:  210  rmse: 0.6970294  mae: 0.5562649  bias:  98.79167 
# SVMPoly model, bincount:  210  rmse: 0.7150884  mae: 0.5691768  bias:  97.332 
# RF model, bincount:  220  rmse: 0.7082407  mae: 0.5746368  bias:  99.88641 
# SVMPoly model, bincount:  220  rmse: 0.7127181  mae: 0.5693054  bias:  98.67426 
# RF model, bincount:  230  rmse: 0.7593422  mae: 0.6230728  bias:  100.1292 
# SVMPoly model, bincount:  230  rmse: 0.7853156  mae: 0.632459  bias:  98.36945 
# RF model, bincount:  240  rmse: 0.7007125  mae: 0.5722464  bias:  100.0042 
# SVMPoly model, bincount:  240  rmse: 0.7185276  mae: 0.5731973  bias:  98.28671 

#  Testing with evalModel, separate=FALSE, bincount=180
# N E1 sigmaknots      RMSE       MAE       Bias
# 1  350 22        1.8 0.1664809 0.1191564   3.522469
# 2  350 22        2.5 0.1704117 0.1210696   3.597216
# 3  350 22        3.2 0.1596300 0.1177950   3.497479
# 4  350 28        1.8 0.5405571 0.4678550  14.035651
# 5  350 28        2.5 0.5532880 0.4832332  14.496997
# 6  350 28        3.2 0.5519287 0.4844753  14.534260
# 7  350 34        1.8 1.0242271 0.9664344  28.993031
# 8  350 34        2.5 1.0281875 0.9789294  29.367881
# 9  350 34        3.2 1.0175731 0.9524530  28.573590
# 10 500 22        1.8 0.5568836 0.4752480  -8.480306
# 11 500 22        2.5 0.5502603 0.4693344  -8.572224
# 12 500 22        3.2 0.5526460 0.4808968  -8.893102
# 13 500 28        1.8 0.5248840 0.4604518   8.936237
# 14 500 28        2.5 0.5086834 0.4481383   8.605371
# 15 500 28        3.2 0.5490623 0.4876274   9.344023
# 16 500 34        1.8 0.9611041 0.9326067  19.584742
# 17 500 34        2.5 0.9570015 0.9302377  19.534991
# 18 500 34        3.2 0.9604658 0.9306464  19.543574
# 19 650 22        1.8 0.8005427 0.7583803 -12.250759
# 20 650 22        2.5 0.7936538 0.7443708 -12.024451
# 21 650 22        3.2 0.7989949 0.7464890 -12.058668
# 22 650 28        1.8 0.2617035 0.2158555  -2.828941
# 23 650 28        2.5 0.2758639 0.2298909  -3.007606
# 24 650 28        3.2 0.2754697 0.2223769  -2.926232
# 25 650 34        1.8 0.2865608 0.2540903   3.577148
# 26 650 34        2.5 0.2915225 0.2570682   3.910054
# 27 650 34        3.2 0.2901203 0.2580734   3.836047

# Same with separate=TRUE
# N E1 sigmaknots      RMSE       MAE       Bias
# 1  350 22        1.8 0.1716267 0.1152317   3.327328
# 2  350 22        2.5 0.1662742 0.1156742   3.361883
# 3  350 22        3.2 0.1558251 0.1027072   2.928293
# 4  350 28        1.8 0.5322795 0.4550937  13.652810
# 5  350 28        2.5 0.5450258 0.4708560  14.124292
# 6  350 28        3.2 0.5505745 0.4647503  13.936522
# 7  350 34        1.8 1.0091665 0.9422542  28.267625
# 8  350 34        2.5 1.0446372 0.9821951  29.465854
# 9  350 34        3.2 1.0247940 0.9703513  29.110539
# 10 500 22        1.8 0.5592049 0.4752118  -8.771711
# 11 500 22        2.5 0.5270185 0.4413740  -8.141071
# 12 500 22        3.2 0.5855653 0.5031882  -9.560878
# 13 500 28        1.8 0.4982140 0.4295861   8.159531
# 14 500 28        2.5 0.5527294 0.4892579   9.364825
# 15 500 28        3.2 0.5331038 0.4695264   9.076009
# 16 500 34        1.8 0.9651677 0.9352988  19.641276
# 17 500 34        2.5 0.9367318 0.9031983  18.967164
# 18 500 34        3.2 0.9563076 0.9271551  19.465377
# 19 650 22        1.8 0.8000131 0.7567852 -12.224992
# 20 650 22        2.5 0.8312455 0.7772122 -12.554967
# 21 650 22        3.2 0.8210832 0.7757801 -12.531833
# 22 650 28        1.8 0.2863042 0.2405860  -3.212277
# 23 650 28        2.5 0.3005720 0.2443985  -3.500773
# 24 650 28        3.2 0.3113075 0.2597785  -3.440604
# 25 650 34        1.8 0.2890498 0.2604604   3.975966
# 26 650 34        2.5 0.3040831 0.2702906   3.900969
# 27 650 34        3.2 0.2866691 0.2559061   3.681310



###########################################################

#  k-nearest neighbour model

knnModel <- train(Dstrip ~ ., data = training, method = "knn")
knnPredict <- predict(knnModel, newdata=testing)
mean((knnPredict-testing$Dstrip)^2)^0.5
# [1] 0.6996156
mean(abs(knnPredict-testing$Dstrip))
# [1] 0.5628021


#############################################
#
#   Fit an ensemble model


#  Models that fit in a short time (less than about a minute):

mlmethods=c("enet", "foba", "icr", "lars", "lars2", "lasso", 
            "leapBackward", "leapForward", "leapSeq", "lm", 
            "nnls", "pcr", "ppr", "relaxo", "ridge", "rqlasso", 
            "rqnc", "superpc", "parRF",
            
            
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


#  Upper bound on MAE/RMSE. Model is just the mean of training data. 
mean(abs(mean(training$Dstrip)-testing$Dstrip))
mean((mean(training$Dstrip)-testing$Dstrip)^2)^0.5
(mean(training$Dstrip)/mean(testing$Dstrip))*100


#  Baseline model -- learn function from number of observations to density
svmBaseModel <- train(Dstrip~numobs, data=training, method = "svmPoly")
svmBasePredict <- predict(svmBaseModel, newdata = testing)
# Root mean square error and MAE
mean((svmBasePredict-testing$Dstrip)^2)^0.5
# [1] 0.637783
mean(abs(svmBasePredict-testing$Dstrip))
# [1] 0.5149678
#  Bias
(mean(svmBasePredict)/mean(testing$Dstrip))*100

#   Baseline with random forest instead of polynomial SVM.

rfBaseModel <- train(Dstrip ~ numobs, data = training, method = "parRF")
rfBasePredict <- predict(rfBaseModel, newdata = testing)
# Root mean square error and MAE
mean((rfBasePredict-testing$Dstrip)^2)^0.5
# [1] 0.7097087
mean(abs(rfBasePredict-testing$Dstrip))
# [1] 0.5653511
#  Bias
(mean(rfBasePredict)/mean(testing$Dstrip))*100

