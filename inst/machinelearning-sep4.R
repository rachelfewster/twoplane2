library(twoplane)
# requires HiddenMarkov, functional, Rcpp, boot, expm, Matrix

library(caret)
library(caretEnsemble)


# Set up parallel processing
library(doMC)
registerDoMC(cores=2)

#  Some functions to support ML

#   Stupid quadratic algorithm.
splitToBins<-function(diffs, binstart, binend, bincount) {
  freqs=rep(0,bincount)
  binsize=(binend-binstart)/bincount
  for(i in 1:bincount) {
    #  Density (per plane second) not frequency now.
    freqs[i]=length(which( diffs>=(i-1)*binsize+binstart & diffs<i*binsize+binstart ))/binsize
  }
  return(freqs)
}

observationsToBins<-function(sdat, L, planespd, bincount, maxsigma) {
  #  Separate non-recapture pairs from possible recapture pairs. 
  #  Non-recapture pairs first. 
  
  #  maxsigma is in km. Convert to plane seconds.
  maxsigmas=maxsigma/planespd
  
  diffs=sort(c(dist(sdat$s1),dist(sdat$s2)))
  
  freqsnonr1=splitToBins(diffs, 0, 2*maxsigmas, bincount/2)
  freqsnonr2=splitToBins(diffs, 2*maxsigmas, L/planespd, bincount/2)
  freqsnonr=c(freqsnonr1, freqsnonr2)
  
  # Possible recap pairs.
  diffs=c()
  for(i in 1:length(sdat$s1)) {
    diffs=c(diffs, abs(sdat$s1[i]-sdat$s2))
  }
  diffs=sort(diffs)
  
  freqspossr1=splitToBins(diffs, 0, 2*maxsigmas, bincount/2)
  freqspossr2=splitToBins(diffs, 2*maxsigmas, L/planespd, bincount/2)
  
  freqspossr=c(freqspossr1, freqspossr2)
  
  #  Fit quadratic to freqspossr1 
  qmodel <- lm(freqspossr1 ~ poly(seq(bincount/2),2))
  smooth<-predict(qmodel, data.frame(x=seq(bincount/2)))
  #model <- smooth.spline(x=seq(bincount/2), y=freqspossr1)
  #smooth<-predict(model, data.frame(x=seq(bincount/2)))$y$x
  
  return(c(freqsnonr, freqspossr, smooth))
}


simulateTestInstance<-function(N, E1, sigmaknots, maxsigmaknots, bincount) {
  E2=110  # expected length unavailable period (from Hiby&Lovell, 1998) in seconds
  #E1=28   # expected length available period (from Hiby&Lovell, 1998) in seconds
  
  E=c(E1,E2)
  Ec=sum(E) # mean dive cycle length in seconds
  p.up=E1/Ec # proportion of time up
  
  L=350;w=0.3 # length and width of strip in km
  planeknots=100 # observer speed in knots
  planespd=planeknots*nm2km/(60^2) # observer speed in km/sec
  
  k=10 # time separation of observers in seconds
  
  sigmarate=sigmaknots*nm2km/60/60 # std dev of movement distance per second (in km)
  sigma=sigmarate*sqrt(k) # std dev of movement distance after k seconds (in km)
  sigma*1000*2 # 95% within this many m of start point after k seconds
  sigma.mult=6 # multiplier to set bound for maximum distance moved
  
  maxsigma=(maxsigmaknots*nm2km/60/60)*sqrt(k)
  
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
  
  bins=observationsToBins(sdat, L, planespd, bincount, maxsigma)
  numobs=length(sdat$s1)+length(sdat$s2)
  
  return(c(bins, numobs, sigma, E1, Dstrip))
}

evalModel<-function(m, Nsim, bincount, biasmodel=NULL) {
  out=data.frame()
  #   27 points
  for(N in c(350,500,650)) {
    for(E1 in c(22,28,34)) {
      for(sigmaknots in c(1.8,2.5, 3.2)) {
        # Simulate test instances
        allsim=data.frame()
        for(i in 1:Nsim) {
          allsim<-rbind(allsim,simulateTestInstance(N, E1, sigmaknots, 3.2, bincount))
        }
        colnames(allsim)<-c(seq(bincount*2.5), "numobs", "sigma", "E1", "Dstrip")
        # Evaluate m
        
        if(!is.null(biasmodel)) {
          rfPredict <- predict2(m, biasmodel, newdata=allsim)
        }
        else {
          rfPredict <- predict(m, newdata = allsim)
        }
        
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

learnBiasModel<-function(inModel, Nsimperpt, Nsim, bincount) {
  #  Learn a model of the bias of inModel
  allsim=data.frame()
  for(i in 1:Nsim) {
    # Pick a point in the parameter space.
    E1=rnorm(1,28,28/3)
    while(E1<=0) {
      E1=rnorm(1,28,28/3)
    }
    
    N=rnorm(1, 500, 500/3)
    while(N<=0) {
      N=rnorm(1, 500, 500/3)
    }
    
    sigmaknots=rnorm(1, 2.5, 2.5/3)
    while(sigmaknots<=0) {
      sigmaknots=rnorm(1, 2.5, 2.5/3)
    }
    
    registerDoMC(cores=8)
    ptsim<-foreach(j=1:Nsimperpt, .combine=rbind, .init=data.frame()) %dopar% simulateTestInstance(N, E1, sigmaknots, 3.2, bincount)
    
    colnames(ptsim)<-c(seq(bincount*2.5), "numobs", "sigma", "E1", "Dstrip")
    
    inModelPredict <- predict(inModel, newdata = ptsim)
    cat("RF model, bincount: ", bincount, " rmse:", mean((inModelPredict-ptsim$Dstrip)^2)^0.5, " mae:", mean(abs(inModelPredict-ptsim$Dstrip)), " bias: ", (mean(inModelPredict)/mean(ptsim$Dstrip))*100-100, "\n")
    ptsim<-cbind(ptsim, (mean(inModelPredict)/mean(ptsim$Dstrip))*100-100 )
    colnames(ptsim)<-c(seq(bincount*2.5), "numobs", "sigma", "E1", "Dstrip", "Bias")
    
    #  Bind ptsim with bias column onto the main data frame.
    allsim<-rbind(allsim, ptsim)
  }
  colnames(allsim)<-c(seq(bincount*2.5), "numobs", "sigma", "E1", "Dstrip", "Bias")
  
  # Smooth only model. 
  registerDoMC(cores=3)
  biasModel <- train(Bias ~ ., data = allsim[,c(seq(bincount*2+1,bincount*2.5),"numobs", "Bias")], method = "parRF")
  
  return(list(biasModel, allsim))
}

predict2<-function(baseModel, biasModel, newdata) {
  bmpred=predict(baseModel, newdata=newdata)
  biaspred=predict(biasModel, newdata=newdata)
  return(bmpred / ((biaspred+100)/100))
}


################################################################
#  Simulate for training and testing of ML model.

seed=654321 # arbitrary fixed seed - for generation of exactly same set of data on separate occasions
set.seed(seed) # initialise random number sequence (for repeatability)

# Number of buckets.
bincount=40

#  Set up data frame for ML input. Last column is the variable to predict. 
numcols=(bincount*2.5)+4
allsim=as.data.frame(matrix(0, ncol = numcols, nrow = 0))
colnames(allsim)=c(seq(bincount*2.5), "numobs", "sigma", "E1", "Dstrip")

nm2km=1.852 # multiplier to convert nautical miles to kilometres

Nsim=5000

#   In each simulation values are selected for E1, N and sigma from normal distribution
do1<-function() {
  E1=rnorm(1,28,28/3)
  while(E1<=0) {
    E1=rnorm(1,28,28/3)
  }
  
  N=rnorm(1, 500, 500/3)
  while(N<=0) {
    N=rnorm(1, 500, 500/3)
  }
  
  sigmaknots=rnorm(1, 2.5, 2.5/3)
  while(sigmaknots<=0) {
    sigmaknots=rnorm(1, 2.5, 2.5/3)
  }
  
  #  The rbind converts the vector returned by simulateTestInstance into a row of a dataframe.  Uses allsim from outer scope.
  blah=rbind(allsim, simulateTestInstance(N, E1, sigmaknots, 3.2, bincount))
  colnames(blah)=c(seq(bincount*2.5), "numobs", "sigma", "E1", "Dstrip")
  return(blah)
}
registerDoMC(cores=8)
allsim<-foreach(cursim=1:Nsim, .combine=rbind, .init=allsim) %dopar% do1() 

colnames(allsim)=c(seq(bincount*2.5), "numobs", "sigma", "E1", "Dstrip")

inTrain <- createDataPartition(y = allsim$Dstrip, p=0.75, list=FALSE)

training <- allsim[ inTrain,]
testing  <- allsim[-inTrain,]

saveRDS(training, "training-mlsep4.rds")
saveRDS(testing, "testing-mlsep4.rds")

##################################################################
##
##   Learn ML model

training <- readRDS("training-mlsep4.rds")
testing <- readRDS("testing-mlsep4.rds")

# Predict Dstrip

rfModel <- train(Dstrip ~ ., data = training[,c(seq(bincount*2.5),"numobs","Dstrip")], method = "parRF")

rfPredict <- predict(rfModel, newdata = testing)
cat("RF model, bincount: ", bincount, " rmse:", mean((rfPredict-testing$Dstrip)^2)^0.5, " mae:", mean(abs(rfPredict-testing$Dstrip)), " bias: ", (mean(rfPredict)/mean(testing$Dstrip))*100-100, "\n")

# Smooth only
rfModelSO <- train(Dstrip ~ ., data = training[,c(seq(bincount*2+1,bincount*2.5),"numobs","Dstrip")], method = "parRF")
rfPredictSO <- predict(rfModelSO, newdata = testing)
cat("RF model smooth only, bincount: ", bincount, " rmse:", mean((rfPredictSO-testing$Dstrip)^2)^0.5, " mae:", mean(abs(rfPredictSO-testing$Dstrip)), " bias: ", (mean(rfPredictSO)/mean(testing$Dstrip))*100-100, "\n")


svmPolyModel <- train(Dstrip ~ ., data = training[,c(seq(bincount*2+1,bincount*2.5),"numobs","Dstrip")], method = "svmPoly")
svmPolyPredict <- predict(svmPolyModel, newdata = testing)
cat("SVMPoly model smooth only, bincount: ", bincount, " rmse:", mean((svmPolyPredict-testing$Dstrip)^2)^0.5, " mae:", mean(abs(svmPolyPredict-testing$Dstrip)), " bias: ", (mean(svmPolyPredict)/mean(testing$Dstrip))*100-100, "\n")

saveRDS(svmPolyModel, "svmPolyModel-mlsep4.rds")

svmPolyModel <- readRDS("svmPolyModel-mlsep4.rds")


svmRBFModel <- train(Dstrip ~ ., data = training[,c(seq(bincount*2+1,bincount*2.5),"numobs","Dstrip")], method = "svmRadial")
svmRBFPredict <- predict(svmRBFModel, newdata = testing)
cat("SVMRBF model smooth only, bincount: ", bincount, " rmse:", mean((svmRBFPredict-testing$Dstrip)^2)^0.5, " mae:", mean(abs(svmRBFPredict-testing$Dstrip)), " bias: ", (mean(svmRBFPredict)/mean(testing$Dstrip))*100-100, "\n")

#  Fit bias model etc. 
svmPolyModelBias=learnBiasModel(svmPolyModel, 100,100,bincount)
svmPolyPredictB=predict2(svmPolyModel,svmPolyModelBias, testing)
cat("SVMPoly model smooth only, bincount: ", bincount, " rmse:", mean((svmPolyPredictB-testing$Dstrip)^2)^0.5, " mae:", mean(abs(svmPolyPredictB-testing$Dstrip)), " bias: ", (mean(svmPolyPredictB)/mean(testing$Dstrip))*100, "\n")


# With intervals E1=runif(1,18,36)    N=runif(1, 300, 750)    sigmaknots=runif(1, 1.5, 3.3)
# RF model, bincount:  40  rmse: 0.5765255  mae: 0.4497514  bias:  100.0243 
# RF model smooth only, bincount:  40  rmse: 0.5672845  mae: 0.4492101  bias:  100.231 

# With norm E1=rnorm(1,28,28/3), N=rnorm(1, 500, 500/3), sigmaknots=rnorm(1, 2.5, 2.5/3)
# RF model, bincount:  40  rmse: 0.705118  mae: 0.5443951  bias:  99.45249 
# RF model smooth only, bincount:  40  rmse: 0.6651206  mae: 0.5184902  bias:  99.42875 
# SVMPoly model smooth only, bincount:  40  rmse: 0.5959252  mae: 0.4542322  bias:  99.13197 


#  SVM poly on smooth only is better!  with the norm training data. 

#  Broader grid search with SVM poly, smooth only. 
svmParamGrid <-  expand.grid(degree = c(2, 3),
                        scale= c(0.1, 0.2, 0.3, 0.4, 0.5),
                        C = c(1.0)
                        )
registerDoMC(cores=8)
svmPolyModel2 <- train(Dstrip ~ ., data = training[,c(seq(bincount*2+1,bincount*2.5),"numobs","Dstrip")], method = "svmPoly", tuneGrid=svmParamGrid)
svmPolyPredict2 <- predict(svmPolyModel2, newdata = testing)
cat("SVMPoly model smooth only, bincount: ", bincount, " rmse:", mean((svmPolyPredict2-testing$Dstrip)^2)^0.5, " mae:", mean(abs(svmPolyPredict2-testing$Dstrip)), " bias: ", (mean(svmPolyPredict2)/mean(testing$Dstrip))*100, "\n")

tmp=learnBiasModel(svmPolyModel2, 100,100,bincount)
svmPolyModel2Bias<-tmp[[1]]
svmPolyPredictn <- predict2(svmPolyModel2, svmPolyModel2Bias, newdata = testing)
cat("SVMPoly model smooth only, bincount: ", bincount, " rmse:", mean((svmPolyPredictn-testing$Dstrip)^2)^0.5, " mae:", mean(abs(svmPolyPredictn-testing$Dstrip)), " bias: ", (mean(svmPolyPredictn)/mean(testing$Dstrip))*100, "\n")


#  Bincount 40, svmPolyModel, smooth curve added
# N E1 sigmaknots      RMSE       MAE      Bias
# 1  350 22        1.8 0.1873612 0.1446336  3.114052
# 2  350 22        2.5 0.1628363 0.1210350  1.980391
# 3  350 22        3.2 0.2006799 0.1333638  3.295395
# 4  350 28        1.8 0.3797013 0.3332474  9.749896
# 5  350 28        2.5 0.4278125 0.3556073 10.028240
# 6  350 28        3.2 0.5483576 0.4491334 12.963130
# 7  350 34        1.8 0.5937211 0.5399048 16.197145
# 8  350 34        2.5 0.6978416 0.6335385 19.006156
# 9  350 34        3.2 0.9181135 0.8422396 25.239753
# 10 500 22        1.8 0.5863783 0.4865168 -8.815438
# 11 500 22        2.5 0.5666954 0.4831494 -6.371459
# 12 500 22        3.2 0.5565069 0.4552279 -2.980423
# 13 500 28        1.8 0.4197177 0.3468611  2.964613
# 14 500 28        2.5 0.4614301 0.3686209  4.889639
# 15 500 28        3.2 0.6150679 0.5194305 10.060890
# 16 500 34        1.8 0.6651889 0.5862685 11.906241
# 17 500 34        2.5 0.7366150 0.6493750 13.443946
# 18 500 34        3.2 0.9495502 0.8910501 18.646460
# 19 650 22        1.8 0.6724773 0.5865128 -9.348095
# 20 650 22        2.5 0.6973602 0.5758094 -9.160471
# 21 650 22        3.2 0.5011685 0.4050958 -6.180489
# 22 650 28        1.8 0.3566253 0.2868719 -3.571246
# 23 650 28        2.5 0.2980606 0.2306436 -2.666495
# 24 650 28        3.2 0.2082375 0.1667609 -1.465912
# 25 650 34        1.8 0.3070148 0.2742343  3.398326
# 26 650 34        2.5 0.2996797 0.2641338  3.394030
# 27 650 34        3.2 0.3040230 0.2641794  3.876197



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

svmPolyModel <- train(Dstrip ~ ., data = training[,c(seq(bincount*2+1,bincount*2.5),"numobs","Dstrip")], method = "svmPoly")

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

rfModel <- train(Dstrip ~., data = training[,c(seq(bincount*2),"numobs","Dstrip")], method = "parRF")
rfPredict <- predict(rfModel, newdata = testing)
# Root mean square error.
mean((rfPredict-testing$Dstrip)^2)^0.5
#  [1] 0.6404903
mean(abs(rfPredict-testing$Dstrip))
#  [1] 0.5252377
#  Bias
(mean(rfPredict)/mean(testing$Dstrip))*100




###########################################################

#  k-nearest neighbour model

knnModel <- train(Dstrip ~ ., data = training[,c(seq(bincount*2+1,bincount*2.5),"numobs","Dstrip")], method = "knn")
knnPredict <- predict(knnModel, newdata=testing)
mean((knnPredict-testing$Dstrip)^2)^0.5
mean(abs(knnPredict-testing$Dstrip))

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

