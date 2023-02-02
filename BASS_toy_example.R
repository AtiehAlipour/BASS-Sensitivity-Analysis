# A toy example of the BASS sensitivity analysis
# Note: the Sobol analysis in this script directly uses the function in 
# "BASS" package instead of sensobol package

# Authors: Haochen Ye(hxy46@psu.edu), 
          # Atieh Alipour(atieg.alipour@dartmouth.edu), 
          # and Klaus Keller(Klaus.Keller@dartmouth.edu)
#############################################################################
# how to run:
# - save the file in a directory
# - open R-Studio and navigate to the folder with the file
# - make this directory the work directory
# - open the file
# - source the file

#############################################################################
            ####Install and load the required libraries####

# Install necessary Packages

list.of.packages <- c("BASS","lhs","caTools","Metrics")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Clear any existing variables and plots. 

rm(list = ls())
graphics.off()

# Loading the libraries

library(BASS)
library(lhs)
library(caTools)
library(Metrics)


#############################################################################
                      ####Define the toy example####

# A simple bivariate function example
f <- function (x1,x2){
  return (x1^2 + 2*x1*x2^2)
}

#############################################################################
#Fit and evaluate the BASS emulator for few sample points with defined seeds#

# Define variables to save the sensitivity indices

S_x1<- c()
S_x2<- c()
S_x1x2<- c()
T_x1<- c()
T_x2<- c()

# Choose number of samples
ns=20

# Calculate sensitivity indices for different seeds
for(i in 1:3){
  
set.seed(i) 
  
  
# Define inputs
a <-  randomLHS(ns, 2)
  
# Find the corresponding Outputs
b <- f(a[ ,1], a[ ,2])  


#Split data to train and test

ind = sample.split(Y = b, SplitRatio = 0.7)

# subsetting into Train data

train_input = a[ind,]
train_output = b[ind]

# subsetting into Test data

test_input = a[!ind,]
test_output = b[!ind]


#Fit the BASS emulator using the train data
#Users may set related parameters for the emulator structure and MCMC
# see more detailed information in bass() function documentation.
#The emulator consists of a series MCMC samples

bass_model <- bass(train_input, train_output)

# Test the emulator
pred=predict(bass_model,test_input)
RMSE_S=rmse(test_output, colMeans(pred))
NRMSE_S=RMSE_S/mean(test_output)
corr_S=cor(test_output, colMeans(pred),  method = "pearson")


#calculate Sobol' indices based on the emulator
#Sensitivity indices for each MCMC sample will be calculated.

S <- sobol(bass_model)

#Plot the sensitivity indices to check its convergence
#Currently, there isn't a quantitative method to check convergence easily,

#plot(c(1:1000),S$S[ ,1],type="l",xlab="MCMC iteration",ylab="First-order sensitivity of x1")
#plot(c(1:1000),S$S[ ,2],type="l",xlab="MCMC iteration",ylab="First-order sensitivity of x2")

#Users can approximate the mean of sensitivity indices (after burn-in period) as the true sensitivity.
#First-order sensitivity:
S_x1[i] <- mean(S$S[ ,1])
S_x2[i]<- mean(S$S[ ,2])

#Second-order sensitivity:
S_x1x2[i] <- mean(S$S[ ,3])

#Total-order sensitivity:
T_x1[i] <- mean(S$T[ ,1])
T_x2[i] <- mean(S$T[ ,2])

}

Small_sample <-  matrix(c(S_x1, S_x2, S_x1x2, T_x1, T_x2),nrow=3,ncol=5,byrow = FALSE)
colnames(Small_sample) <- c("X1-first-order-sensitivity", "X2-first-order-sensitivity", "second-order-sensitivity", "X1-total-effect", "X2-total-effect")
rownames(Small_sample) <- c("Seed = 1","Seed = 2","Seed = 3")

print(Small_sample)


#############################################################################
      #### Display the results and save them in a pdf file ####

barplot(Small_sample, beside = TRUE,cex.axis=0.5,legend= rownames(Small_sample),cex.names=0.5, main = "Sensitivity Indices for Small Sample Size",col=c("cornsilk4","red","orange"),ylim = c(0, 0.8))
pdf(file="Small_sample.pdf",10,6.17)  
barplot(Small_sample, beside = TRUE,cex.axis=0.5,legend= rownames(Small_sample),cex.names=0.5, main = "Sensitivity Indices for Small Sample Size",col=c("cornsilk4","red","orange"),ylim = c(0, 0.8))
dev.off() 


#############################################################################
          ######Repeat the process for larger sample size######

#Fit and evaluate the BASS emulator for large sample points with defined seeds#

# Define variables to save the sensitivity indices

S_x1<- c()
S_x2<- c()
S_x1x2<- c()
T_x1<- c()
T_x2<- c()

# Choose number of samples
nl=100

# Calculate sensitivity indices for different seeds
for(i in 1:3){
  
  set.seed(i) 
  
  
  # Define inputs
  a <-  randomLHS(nl, 2)
  
  # Find the corresponding Outputs
  b <- f(a[ ,1], a[ ,2])  
  
  
  #Split data to train and test
  
  ind = sample.split(Y = b, SplitRatio = 0.7)
  
  # subsetting into Train data
  
  train_input = a[ind,]
  train_output = b[ind]
  
  # subsetting into Test data
  
  test_input = a[!ind,]
  test_output = b[!ind]
  
  
  #Fit the BASS emulator using the train data
  #Users may set related parameters for the emulator structure and MCMC
  # see more detailed information in bass() function documentation.
  #The emulator consists of a series MCMC samples
  
  bass_model <- bass(train_input, train_output)
  
  # Test the emulator
  pred=predict(bass_model,test_input)
  RMSE_L=rmse(test_output, colMeans(pred))
  NRMSE_L=RMSE_L/mean(test_output)
  corr_L=cor(test_output, colMeans(pred),  method = "pearson")

  
  #calculate Sobol' indices based on the emulator
  #Sensitivity indices for each MCMC sample will be calculated.
  
  S <- sobol(bass_model)
  
  #Plot the sensitivity indices to check its convergence
  #Currently, there isn't a quantitative method to check convergence easily,
  
  #plot(c(1:1000),S$S[ ,1],type="l",xlab="MCMC iteration",ylab="First-order sensitivity of x1")
  #plot(c(1:1000),S$S[ ,2],type="l",xlab="MCMC iteration",ylab="First-order sensitivity of x2")
  
  #Users can approximate the mean of sensitivity indices (after burn-in period) as the true sensitivity.
  #First-order sensitivity:
  S_x1[i] <- mean(S$S[ ,1])
  S_x2[i]<- mean(S$S[ ,2])

  #Second-order sensitivity:
  S_x1x2[i] <- mean(S$S[ ,3])
  
  #Total-order sensitivity:
  T_x1[i] <- mean(S$T[ ,1])
  T_x2[i] <- mean(S$T[ ,2])
}

large_sample <- matrix(c(S_x1, S_x2, S_x1x2, T_x1, T_x2),nrow=3,ncol=5,byrow = FALSE)
colnames(large_sample) <- c("X1-first-order-sensitivity", "X2-first-order-sensitivity", "second-order-sensitivity", "X1-total-effect", "X2-total-effect")
rownames(large_sample) <- c("Seed = 1","Seed = 2","Seed = 3")
print(large_sample)


#############################################################################
#### Display the results and save them in a pdf file ###

barplot(large_sample, beside = TRUE,cex.axis=0.5,legend= rownames(large_sample),cex.names=0.5, main = "Sensitivity Indices for Large Sample Size",col=c("cornsilk4","red","orange"),ylim = c(0, 0.8))
pdf(file="Large_sample.pdf",10,6.17)  
barplot(large_sample, beside = TRUE,cex.axis=0.5,legend= rownames(large_sample),cex.names=0.5, main = "Sensitivity Indices for Large Sample Size",col=c("cornsilk4","red","orange"),ylim = c(0, 0.8))
dev.off() 


