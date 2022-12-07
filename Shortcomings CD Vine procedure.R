# Experiment to show that the score fails
# to recover the true DAG.

# Loading the required packages:
library("bnlearn")    # Gaussian and discrete BN
library("Rgraphviz")  # Plotting DAGs
library("gRain")
library("cubature")   # Numerical Integration
library("rlang")
library("rbmn")       # Reseaux bayesiens multinormaux
library("lattice")
library("tictoc")     # to measure the time
library("igraph")
library("penalized")  # Lasso and Ridge regression
library("corpcor")    # inverse correlation matrix
library("compare")    # comparing TP, FP and FN
library("GGMselect")  # Gaussian Graphical Models
library("VineCopula") # Vine Copula models
library("ggplot2")    # Plots high quality
library("network")    # Plotting networks of the Lasso
library("tibble")     # Data frame to Latex format
library("knitr")      # Data frame to Latex format
library("latex2exp")  # Latex for the plots
library("Rfast")      # Max lik estimator of the Mult.normal
library("MASS")       # simulations of the Mult.normal
library("GGally")     # Plotting pairs data using ggplot
library("network")    # undirected graph GGMselec
library("MVN")        # goodness of fit test MVN


# The objective is to show that the Score proposed by the VINE DAG fails
# to recover the true DAG.


X <- paste0("[X1][X2][X3|X1:X2]")
dag<-model2network(X)
dag
graphviz.plot(dag)

Y <- paste0("[X1][X3|X1][X2|X1:X3]")
dag2 <-model2network(Y)
dag2
graphviz.plot(dag2)


# Simulation
Matrix <- c(1, 3, 2, 
            0, 2, 3, 
            0, 0, 3)
Matrix <- matrix(Matrix, 3, 3)

shd(dag,dag2)

# Possible scenarios: combination of copula families and parameters:
family13.2 <- rep(c(3,4,5,6,3,4,5,3,4,3),each=4)
family23 <- rep(c(3,4,5,6,4,5,6,5,6,6),each=6)


par13.2 <- rep(c(0.25, 0.25, 0.75, 0.75),10)
par23 <- rep(c(0.25, 0.75, 0.25, 0.75),10)



DF <- data.frame(matrix(NA, nrow = length(par23), ncol = 4))
colnames(DF) = c("BIC BN1","BIC BN2", "IC BN1", "IC BN2")


for(j in 1:length(par23)){
  
  # families
  family <- c(0, family13.2[j], 0,  
              0, 0, family23[j], 
              0, 0, 0)
  family <- matrix(family, 3, 3)
  
  # parameters
  par <- c(0, BiCopTau2Par(family13.2[j],par13.2[j]), 0,
           0, 0,  BiCopTau2Par(family23[j],par23[j]),
           0, 0, 0)
  par <- matrix(par, 3, 3)
  
  # second parameters
  par2 <- matrix(5, 3, 3)
  
  # define the Vine
  RVM <- RVineMatrix(Matrix = Matrix, family = family,
                     par = par, par2 = par2,
                     names = c("X1", "X2", "X3"))
  
  N = 2500
  databbn <- RVineSim(N, RVM)
  
  
  # BIC BN1.
  RVMcop <- RVineCopSelect(databbn, familyset = c(0:6), Matrix, 
                           selectioncrit = "logLik",indeptest = TRUE,
                           level = 0.05, trunclevel = NA,
                           weights = NA, rotations = TRUE, se = FALSE,
                           presel = TRUE, method = "mle", cores = 1)
  
  BIC.BN1 <-  RVMcop$logLik-log(N)
  
  
  # BIC BN2.
  CMatrix1 <- c(1, 2, 3,
                0, 2, 3,
                0, 0, 3)
  CMatrix1 <- matrix(CMatrix1,3,3)
  RVMcop1 <- RVineCopSelect(databbn, familyset = c(0:6), CMatrix1, 
                            selectioncrit = "logLik",indeptest = FALSE,
                            level = 0.05, trunclevel = NA,
                            weights = NA, rotations = TRUE, se = FALSE,
                            presel = TRUE, method = "mle", cores = 1)
  
  BIC.BN2 <-  RVMcop1$logLik-log(N)*3/2
  
  
  
  # IC BN1
  CMatrix2 <- c(1, 2, 3,
                0, 2, 3,
                0, 0, 3)
  CMatrix2 <- matrix(CMatrix2,3,3)
  RVMcop2 <- RVineCopSelect(databbn, familyset = c(0:6), CMatrix2, 
                            selectioncrit = "logLik",indeptest = FALSE,
                            level = 0.05, trunclevel = NA,
                            weights = NA, rotations = TRUE, se = FALSE,
                            presel = TRUE, method = "mle", cores = 1)
  
  C12.approx <- BiCopSelect(databbn[,1], databbn[,2], family = c(0,1,3:6), 
                            selectioncrit = "logLik", indeptest = FALSE, 
                            level = 0.05, weights = NA, rotations = TRUE,
                            se = FALSE,presel = TRUE,method = "mle")
  
  IC.BN1 <-  (RVMcop2$logLik - C12.approx$logLik) - C12.approx$logLik/(2*log(N))
  
  
  # IC BN2
  CMatrix3 <- c(1, 3, 2, 
                0, 2, 3, 
                0, 0, 3)
  CMatrix3 <- matrix(CMatrix3, 3, 3)
  RVMcop3 <- RVineCopSelect(databbn, familyset = c(0:6), CMatrix3, 
                            selectioncrit = "logLik",indeptest = FALSE,
                            level = 0.05, trunclevel = NA,
                            weights = NA, rotations = TRUE, se = FALSE,
                            presel = TRUE, method = "mle", cores = 1)
  
  C13.approx <- BiCopSelect(databbn[,1], databbn[,3], family = c(0,1,3:6), 
                            selectioncrit = "logLik", indeptest = FALSE, 
                            level = 0.05, weights = NA, rotations = TRUE,
                            se = FALSE,presel = TRUE,method = "mle")
  
  IC.BN2 <-  (RVMcop3$logLik) - (C13.approx$logLik/(2*log(N))+sum(log(databbn[,1]))/log(N))
  
  
  # Storing the values
  DF[j,] = c(BIC.BN1, BIC.BN2, IC.BN1, IC.BN2)
  print(j)
  
}


sum(DF["BIC BN1"]>DF["BIC BN2"])
sum(DF["IC BN1"]>DF["IC BN2"])







