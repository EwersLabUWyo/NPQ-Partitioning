##Function for partitioning the q values and outputting all important data in one table
library(segmented)
library(tidyverse)
library(lubridate)
library(dplyr)
library(plyr)
library(rjags)
library(coda)
###TEST####
#Here are some sample data that are able to run according to  the function below
datFLUO <- read.csv("Data/Fluorpennpq.csv", header = TRUE)
#subsetting the data to mitigate run time
dataFLUO <- datFLUO[1:35,]

datNPQlong <- read.csv("Data/NPQlong.csv", header = TRUE)
#subsetting the data to mitigate run time
datNPQlong <- datNPQlong[1:30,]
data <- datNPQlong

#Running the function on the data provided for either segmented or Bayes.
diditwork1 <- NPQ_Parts("Segmented", data = datFLUO)
diditwork <- NPQ_Parts("Segmented", data = datNPQlong)
diditworkbayse1 <- NPQ_Parts("Bayes", data = datFLUO)
diditworkbayse <- NPQ_Parts("Bayes", data = datNPQlong)
diditworkboth1 <- NPQ_Parts("Both", data = datFLUO)



##############Function for calculating the intercepts and energy partitions of Light Relaxation Curves########
NPQ_Parts <- function(method, data){
  #Making list of IDs to reference throughout this function.
  IDList <- as.list(unique(data$ID))
  #####Being segmented method of partitioning###########
  if(method == "Segmented"){
    print("Running Segmented Model")

    ##splitting each ID into is own data frame/list 
    mylist <- list()
    for (i in IDList){
      #reading each iteration into its own 'Data_i' list
      mylist[[paste("Data", i, sep="_")]]<- data %>% 
        filter(ID == i) 
      
    }
  
    ##Creating Intercepts for data file
    InterceptList <- list()
    #Making a vector of potential break numbers
    #This list can be altered depending on expectation
    breaks <- c(1,2,3,4,5)
    #forloop to calculate intercept values using lm and segmented function
    for(i in 1:length(mylist)){
      print(paste("SegmentedMeasurment", i, sep= " "))
      for(x in 1:length(breaks)){
        #making linear model for each unique curve of NPQ data
        model1<- lm(NPQ~Pulse,data = mylist[[i]])
        ####TryCatch Function
        #TryCatch works to identify `errors` and rerun with only one breaks
        tryCatch( expr = { 
          #placing linear model into the segmented model and running trycatch
          segmod <- segmented::segmented (model1, seg.Z= ~ Pulse,
                                npsi = breaks[[x]],
                                data = mylist[[i]])
          
          print(paste ("Data ran with", breaks[x], "break points.", sep = " "))
          
          #making a data frame that contains the intercept values that are used to calculate the partitions
          singleintercept<- data.frame(matrix(unlist(intercept(segmod)),
                                              nrow=length(intercept(segmod)), byrow=TRUE))
          #for loop to rename columns based upon the number of intercepts that are produced
          for (x in 1:(ncol(singleintercept))){
            colnames(singleintercept)[x] <- paste("SegmentedIntercept",x, sep = "")
          }  
          for(x in 1:ncol(singleintercept)) {
            singleintercept[x]<- (abs(singleintercept[x]))
          }
          singleintercept  <- round(singleintercept,5)
          InterceptList[[paste("SegmentedIntercepts", IDList[i], sep = "_")]] <- singleintercept
        }, #if trycatch runs an error, do following
        error = function(e){
          print(paste("Did not run with",breaks[[x]],"break points.", sep = " "))
        }, #if trycatch runs a warning, do the following
        warning = function(w){
          print(paste("Did not run with",breaks[[x]],"break points.", sep = " "))
        }, #finish the looping trycatch function
        finally = { 
          print("Executed.")}
        )
      }
    }
    
############add partitioning code here! with unique output ie `segmentedout`
    ####### and combine everything together...
    ###Partitioning of NPQ values##########
    Partitioning <- list()
    ###adding partition function  into code
    qvalseg <- function(input){
      Qvalues <- data.frame(matrix(NA, nrow = 1, ncol = length(input)*3))
      #premaking columns in singledat
      for (x in 1:(length(input))){
        colnames(Qvalues)[x] <- paste("SegmentedQabsolute",x, sep = "")
      }
      for (x in (length(input)+1):(length(input)*2)) {
        colnames(Qvalues)[x] <- paste("SegmentedQrelative",x-length(input), sep = "")
      }
      for(x in ((length(input)*2)+1):(length(input)*3)){
        colnames(Qvalues)[x] <- paste("SegmentedQpercentage",x-length(input)*2, sep = "")
      } 
      
      ##Code for filling in the Qvalues table.
      #calculating Q(E,T,I ...) absolute values
      for (i in 1:(ncol(input)-1)) {
        Qvalues[,i] <- abs(input[i]-input[i+1])
        Qvalues[,ncol(input)]<- input[ncol(input)]
        #calculating Q(E,T,I...) relative values
        Qvalues[,i+(length(input))] <- ((input[i]-input[i+1])/(1-input[i+1]))
        Qvalues[,ncol(input)*2]<- input[ncol(input)]
      }
      #calculating Q(E,T,I...) percentage
      for (i in 1:(ncol(input))) {
        Qvalues[,i+(length(input)*2)] <- Qvalues[i]/input[,1]
      }
      round(Qvalues,5) ## rounding to 5 decimal places for convenience
    }
    
    #For loop to partition the data into the calculated qvalues
    for (i in 1:length(InterceptList)){
      Partitioning[[paste("SegmentedPartition", IDList[i], sep = "_")]] <- qvalseg(InterceptList[[i]])
    }
    
      ######SegmentedOUTPUT######
    OutputSegmented <- data.frame(matrix(nrow= length(Partitioning),ncol= (length(Partitioning[[1]])+length(InterceptList[[1]]))))
    for (i in 1:length(Partitioning)){
      OutputSegmented[i,] <- data.frame(Partitioning[[i]], InterceptList[[i]])
      names(OutputSegmented)<- c(names(Partitioning[[1]]),names(InterceptList[[1]]))
    }
    #adding original IDs to partition list
    OutputSegmented$ID <- IDList
    
    OutputSegmented
    
    
  } else if (method == "Bayes")
    
  {
    ####Beginning of Bayes Model of partitioning#####
    print("Running Bayes Model")
    ####BayesNPQ function directly placed into the code#####
    BayesNPQ<-function(data, modelstr, nChains=nChains, adaptSteps=adaptSteps, burnInSteps=burnInSteps, 
                       parameters, nPerChain=nPerChain , thinSteps=thinSteps, DICsteps=DICsteps){
      ### Begin model
      modelrun <- jags.model(textConnection(modelstr),
                             data = data, n.chains=nChains , n.adapt=adaptSteps)
      ### Update model
      update(modelrun, burnInSteps)
      ### save samples
      mcmc_samples<- coda.samples(modelrun,
                                  variable.names=parameters,
                                  n.iter=nPerChain , thin=thinSteps )
      mcmcChain<- as.matrix(mcmc_samples)
      sigma <-1  / sqrt( mcmcChain[, "tau" ] )
      mcmcChain<-as.data.frame(cbind( mcmcChain, sigma ))
      ## median estimates - transposed to call by col names
      medsH<-as.data.frame(t(apply(mcmcChain,2,median)))
      ###gelman diagnostics
      gelman_stats<-gelman.diag(mcmc_samples, multivariate = FALSE)
      ### DIC samples
      DIC <- dic.samples(modelrun, DICsteps, "pD")
      ### list output
      outlist<-list(mcmcChain ,medsH,gelman_stats, DIC)
      return(outlist)
    }
    #source code for the Bayes 2 transition points
    source("Jons Work/Segment_Model_Bayes_2Trans.R")
    #defining parameters for 2 breaks
    parameters2 = c("yi", "Intercepta","Intercept", "slope","Trans", "tau", "QE","QT")### pars to be monitored
    ###### Set up for rjags #########
    adaptSteps = 1500             # Number of steps to "tune" the samplers.
    burnInSteps = 1500            # Number of steps to "burn-in" the samplers.
    nChains = 4                   # Number of chains to run.
    DICsteps= 20                  # Number of steps of sample DIC
    numSavedSteps= 50             # Total number of steps in chains to save.
    thinSteps=100                 # Number of steps to "thin" (1=keep every step).
    nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
    ####
    ##splitting each ID into is own data frame/list Data
    mylist2 <- list()
    for (i in IDList){
      print(i)
      mylist2[[paste( "Data",i, sep = "_")]]<- data %>% 
        filter(ID == i) %>% 
        dplyr::select(ID,NPQ,Pulse)  
      mylist3 <- list()
    for (i in 1:length(mylist2)){
      mylist3[[paste(i)]] <- list(N = length(mylist2[[i]][["NPQ"]]), NPQ = mylist2[[i]][["NPQ"]], Time = mylist2[[i]][["Pulse"]])
    }
    }
    #making a list to feed the Bayes NPQ function into 
    modellist <- list()
    for (i in 1:length(mylist2)){
      print(i)
   ######Running Bayesian Model####
     modellist[[paste(i)]] <- BayesNPQ(mylist3[[i]], Seg_model_2Trans, nChains ,adaptSteps, burnInSteps,parameters2,nPerChain, thinSteps, DICsteps)
    }
    InterceptList <- list()
    for(i in 1:length(modellist)){
      intlist <- list(Bayesintercep1 = modellist[[i]][[2]][[5]], Bayesintercept2 = modellist[[i]][[2]][[3]], Bayesintercept3 = modellist[[i]][[2]][[4]])
      InterceptList[[paste("BayesIntercepts",i,sep="_")]] <- as.data.frame(x = intlist, row.names = NULL) 
    }
###End of Bayes Code
    
#########ADD partitioning code here with unique name for objects.
    ###Partitioning of NPQ values##########
    Partitioning <- list()
    ###adding partition function  into code
    qvalbayes <- function(input){
      Qvalues <- data.frame(matrix(NA, nrow = 1, ncol = length(input)*3))
      #premaking columns in singledat
      for (x in 1:(length(input))){
        colnames(Qvalues)[x] <- paste("BayesQabsolute",x, sep = "")
      }
      for (x in (length(input)+1):(length(input)*2)) {
        colnames(Qvalues)[x] <- paste("BayesQrelative",x-length(input), sep = "")
      }
      for(x in ((length(input)*2)+1):(length(input)*3)){
        colnames(Qvalues)[x] <- paste("BayesQpercentage",x-length(input)*2, sep = "")
      } 
      
      ##Code for filling in the Qvalues table.
      #calculating Q(E,T,I ...) absolute values
      for (i in 1:(ncol(input)-1)) {
        Qvalues[,i] <- abs(input[i]-input[i+1])
        Qvalues[,ncol(input)]<- input[ncol(input)]
        #calculating Q(E,T,I...) relative values
        Qvalues[,i+(length(input))] <- ((input[i]-input[i+1])/(1-input[i+1]))
        Qvalues[,ncol(input)*2]<- input[ncol(input)]
      }
      #calculating Q(E,T,I...) percentage
      for (i in 1:(ncol(input))) {
        Qvalues[,i+(length(input)*2)] <- Qvalues[i]/input[,1]
      }
      round(Qvalues,5) ## rounding to 5 decimal places for convenience
    }
    
    #For loop to partition the data into the calculated qvalues
    for (i in 1:length(InterceptList)){
      Partitioning[[paste("BayesPartition", IDList[i], sep = "_")]] <- qvalbayes(InterceptList[[i]])
    }
    
    ##########OUTPUT#########
    OutputBayes <- data.frame(matrix(nrow= length(Partitioning),ncol= (length(Partitioning[[1]])+length(InterceptList[[1]]))))
    for (i in 1:length(Partitioning)){
      OutputBayes[i,] <- data.frame( Partitioning[[i]], InterceptList[[i]])
      names(OutputBayes)<- c(names(Partitioning[[1]]),names(InterceptList[[1]]))
    }
    #adding original IDs to partition list
    OutputBayes$ID <- IDList
    
    OutputBayes
  } 
  
 #########combining both Outputs######
else if (method == "Both"){
  print("Running Segmented Model")
  
  ##splitting each ID into is own data frame/list 
  mylist <- list()
  for (i in IDList){
    #reading each iteration into its own 'Data_i' list
    mylist[[paste("Data", i, sep="_")]]<- data %>% 
      filter(ID == i) 
    
  }
  
  ##Creating Intercepts for data file
  InterceptList <- list()
  #Making a vector of potential break numbers
  #This list can be altered depending on expectation
  breaks <- c(1,2,3,4,5)
  #forloop to calculate intercept values using lm and segmented function
  for(i in 1:length(mylist)){
    print(paste("SegmentedMeasurment", i, sep= " "))
    for(x in 1:length(breaks)){
      #making linear model for each unique curve of NPQ data
      model1<- lm(NPQ~Pulse,data = mylist[[i]])
      ####TryCatch Function
      #TryCatch works to identify `errors` and rerun with only one breaks
      tryCatch( expr = { 
        #placing linear model into the segmented model and running trycatch
        segmod <- segmented::segmented (model1, seg.Z= ~ Pulse,
                                        npsi = breaks[[x]],
                                        data = mylist[[i]])
        
        print(paste ("Data ran with", breaks[x], "break points.", sep = " "))
        
        #making a data frame that contains the intercept values that are used to calculate the partitions
        singleintercept<- data.frame(matrix(unlist(intercept(segmod)),
                                            nrow=length(intercept(segmod)), byrow=TRUE))
        #for loop to rename columns based upon the number of intercepts that are produced
        for (x in 1:(ncol(singleintercept))){
          colnames(singleintercept)[x] <- paste("SegmentedIntercept",x, sep = "")
        }  
        for(x in 1:ncol(singleintercept)) {
          singleintercept[x]<- (abs(singleintercept[x]))
        }
        singleintercept  <- round(singleintercept,5)
        InterceptList[[paste("SegmentedIntercepts", IDList[i], sep = "_")]] <- singleintercept
      }, #if trycatch runs an error, do following
      error = function(e){
        print(paste("Did not run with",breaks[[x]],"break points.", sep = " "))
      }, #if trycatch runs a warning, do the following
      warning = function(w){
        print(paste("Did not run with",breaks[[x]],"break points.", sep = " "))
      }, #finish the looping trycatch function
      finally = { 
        print("Executed.")}
      )
    }
  }
  
  ############add partitioning code here! with unique output ie `segmentedout`
  ####### and combine everything together...
  ###Partitioning of NPQ values##########
  Partitioning <- list()
  ###adding partition function  into code
  qvalseg <- function(input){
    Qvalues <- data.frame(matrix(NA, nrow = 1, ncol = length(input)*3))
    #premaking columns in singledat
    for (x in 1:(length(input))){
      colnames(Qvalues)[x] <- paste("SegmentedQabsolute",x, sep = "")
    }
    for (x in (length(input)+1):(length(input)*2)) {
      colnames(Qvalues)[x] <- paste("SegmentedQrelative",x-length(input), sep = "")
    }
    for(x in ((length(input)*2)+1):(length(input)*3)){
      colnames(Qvalues)[x] <- paste("SegmentedQpercentage",x-length(input)*2, sep = "")
    } 
    
    ##Code for filling in the Qvalues table.
    #calculating Q(E,T,I ...) absolute values
    for (i in 1:(ncol(input)-1)) {
      Qvalues[,i] <- abs(input[i]-input[i+1])
      Qvalues[,ncol(input)]<- input[ncol(input)]
      #calculating Q(E,T,I...) relative values
      Qvalues[,i+(length(input))] <- ((input[i]-input[i+1])/(1-input[i+1]))
      Qvalues[,ncol(input)*2]<- input[ncol(input)]
    }
    #calculating Q(E,T,I...) percentage
    for (i in 1:(ncol(input))) {
      Qvalues[,i+(length(input)*2)] <- Qvalues[i]/input[,1]
    }
    round(Qvalues,5) ## rounding to 5 decimal places for convenience
  }
  
  #For loop to partition the data into the calculated qvalues
  for (i in 1:length(InterceptList)){
    Partitioning[[paste("SegmentedPartition", IDList[i], sep = "_")]] <- qvalseg(InterceptList[[i]])
  }
  
     ######SegmentedOUTPUT######
  OutputSegmented <- data.frame(matrix(nrow= length(Partitioning),ncol= (length(Partitioning[[1]])+length(InterceptList[[1]]))))
  for (i in 1:length(Partitioning)){
    OutputSegmented[i,] <- data.frame(Partitioning[[i]], InterceptList[[i]])
    names(OutputSegmented)<- c(names(Partitioning[[1]]),names(InterceptList[[1]]))
  }
  #adding original IDs to partition list
  OutputSegmented$ID <- IDList
  
  OutputSegmented
  
  
  print("Running Bayes Model")
   ######BayesNPQ function directly placed into the code#####
  BayesNPQ<-function(data, modelstr, nChains=nChains, adaptSteps=adaptSteps, burnInSteps=burnInSteps, 
                     parameters, nPerChain=nPerChain , thinSteps=thinSteps, DICsteps=DICsteps){
    ### Begin model
    modelrun <- jags.model(textConnection(modelstr),
                           data = data, n.chains=nChains , n.adapt=adaptSteps)
    ### Update model
    update(modelrun, burnInSteps)
    ### save samples
    mcmc_samples<- coda.samples(modelrun,
                                variable.names=parameters,
                                n.iter=nPerChain , thin=thinSteps )
    mcmcChain<- as.matrix(mcmc_samples)
    sigma <-1  / sqrt( mcmcChain[, "tau" ] )
    mcmcChain<-as.data.frame(cbind( mcmcChain, sigma ))
    ## median estimates - transposed to call by col names
    medsH<-as.data.frame(t(apply(mcmcChain,2,median)))
    ###gelman diagnostics
    gelman_stats<-gelman.diag(mcmc_samples, multivariate = FALSE)
    ### DIC samples
    DIC <- dic.samples(modelrun, DICsteps, "pD")
    ### list output
    outlist<-list(mcmcChain ,medsH,gelman_stats, DIC)
    return(outlist)
  }
  #source code for the Bayes 2 transition points
  source("Jons Work/Segment_Model_Bayes_2Trans.R")
  #defining parameters for 2 breaks
  parameters2 = c("yi", "Intercepta","Intercept", "slope","Trans", "tau", "QE","QT")### pars to be monitored
  ###### Set up for rjags #########
  adaptSteps = 1500             # Number of steps to "tune" the samplers.
  burnInSteps = 1500            # Number of steps to "burn-in" the samplers.
  nChains = 4                   # Number of chains to run.
  DICsteps= 20                  # Number of steps of sample DIC
  numSavedSteps= 50             # Total number of steps in chains to save.
  thinSteps=100                 # Number of steps to "thin" (1=keep every step).
  nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.
  ####
  ##splitting each ID into is own data frame/list Data
  mylist2 <- list()
  for (i in IDList){
    print(i)
    mylist2[[paste( "Data",i, sep = "_")]]<- data %>% 
      filter(ID == i) %>% 
      dplyr::select(ID,NPQ,Pulse)  
    mylist3 <- list()
    for (i in 1:length(mylist2)){
      mylist3[[paste(i)]] <- list(N = length(mylist2[[i]][["NPQ"]]), NPQ = mylist2[[i]][["NPQ"]], Time = mylist2[[i]][["Pulse"]])
    }
  }
  #making a list to feed the Bayes NPQ function into 
  modellist <- list()
  for (i in 1:length(mylist2)){
    print(i)
    ######Running Bayesian Model####
    modellist[[paste(i)]] <- BayesNPQ(mylist3[[i]], Seg_model_2Trans, nChains ,adaptSteps, burnInSteps,parameters2,nPerChain, thinSteps, DICsteps)
  }
  InterceptList <- list()
  for(i in 1:length(modellist)){
    intlist <- list(Bayesintercep1 = modellist[[i]][[2]][[5]], Bayesintercept2 = modellist[[i]][[2]][[3]], Bayesintercept3 = modellist[[i]][[2]][[4]])
    InterceptList[[paste("BayesIntercepts",i,sep="_")]] <- as.data.frame(x = intlist, row.names = NULL) 
  }
  ###End of Bayes Code
  
    ######Partitioning of NPQ values##########
  Partitioning <- list()

  ###adding partition function  into code
  qvalbayes <- function(input){
    Qvalues <- data.frame(matrix(NA, nrow = 1, ncol = length(input)*3))
    #premaking columns in singledat
    for (x in 1:(length(input))){
      colnames(Qvalues)[x] <- paste("BayesQabsolute",x, sep = "")
    }
    for (x in (length(input)+1):(length(input)*2)) {
      colnames(Qvalues)[x] <- paste("BayesQrelative",x-length(input), sep = "")
    }
    for(x in ((length(input)*2)+1):(length(input)*3)){
      colnames(Qvalues)[x] <- paste("BayesQpercentage",x-length(input)*2, sep = "")
    } 
    
    ##Code for filling in the Qvalues table.
    #calculating Q(E,T,I ...) absolute values
    for (i in 1:(ncol(input)-1)) {
      Qvalues[,i] <- abs(input[i]-input[i+1])
      Qvalues[,ncol(input)]<- input[ncol(input)]
      #calculating Q(E,T,I...) relative values
      Qvalues[,i+(length(input))] <- ((input[i]-input[i+1])/(1-input[i+1]))
      Qvalues[,ncol(input)*2]<- input[ncol(input)]
    }
    #calculating Q(E,T,I...) percentage
    for (i in 1:(ncol(input))) {
      Qvalues[,i+(length(input)*2)] <- Qvalues[i]/input[,1]
    }
    round(Qvalues,5) ## rounding to 5 decimal places for convenience
  }
  #For loop to partition the data into the calculated qvalues
  for (i in 1:length(InterceptList)){
    Partitioning[[paste("BayesPartition", IDList[i], sep = "_")]] <- qvalbayes(InterceptList[[i]])
  }

    ######BayesOUTPUT##########
  OutputBayes <- data.frame(matrix(nrow= length(Partitioning),ncol= (length(Partitioning[[1]])+length(InterceptList[[1]]))))
  for (i in 1:length(Partitioning)){
    OutputBayes[i,] <- data.frame( Partitioning[[i]], InterceptList[[i]])
    names(OutputBayes)<- c(names(Partitioning[[1]]),names(InterceptList[[1]]))
  }
  #adding original IDs to partition list
  OutputBayes$ID <- IDList
  
  OutputBayes
  OutputAll <- cbind(OutputSegmented,OutputBayes)
}

    }
