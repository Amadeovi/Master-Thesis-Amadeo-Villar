# Simulation study chapter 6: Adjacency Matrix
# Experiment 17/09 part 3
# Investigate adjacency matrices, for the  big network
# gaussian case.


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
library("pcalg")      # PC algorithm library
library("reshape2")   # Melting data frame
library("plyr")


# This case we will have a tree.
par(mfrow=c(1,1))
X <- paste0("[X1][X2|X1][X3|X2][X4|X3][X5|X2][X6|X5][X7|X6][X8|X6][X9|X5][X10|X9]")
dag <- model2network(X)
dag
graphviz.plot(dag)


# 64 different scenarios:
# families:
family12 <- rep(c(1,1,1,1,1,1,1,1),each=8)
family23 <- rep(c(1,1,1,1,1,1,1,1),each=8)
family34 <- rep(c(1,1,1,1,1,1,1,1),each=8)
family25 <- rep(c(1,1,1,1,1,1,1,1),each=8)
family56 <- rep(c(1,1,1,1,1,1,1,1),each=8)
family67 <- rep(c(1,1,1,1,1,1,1,1),each=8)
family68 <- rep(c(1,1,1,1,1,1,1,1),each=8)
family59 <- rep(c(1,1,1,1,1,1,1,1),each=8)
family910 <- rep(c(1,1,1,1,1,1,1,1),each=8)


# parameters:
par12 <- rep(c(0.25, 0.75, 0.25, 0.75, 0.25, 0.75, 0.25, 0.75),8)
par23 <- rep(c(0.25, 0.75, 0.25, 0.75, 0.25, 0.75, 0.75, 0.25),8)
par34 <- rep(c(0.25, 0.75, 0.25, 0.75, 0.75, 0.25, 0.25, 0.75),8)
par25 <- rep(c(0.25, 0.75, 0.25, 0.75, 0.75, 0.25, 0.75, 0.25),8)
par56 <- rep(c(0.25, 0.75, 0.75, 0.25, 0.25, 0.75, 0.25, 0.75),8)
par67 <- rep(c(0.25, 0.75, 0.75, 0.25, 0.25, 0.75, 0.75, 0.25),8)
par68 <- rep(c(0.25, 0.75, 0.75, 0.25, 0.75, 0.25, 0.25, 0.75),8)
par59 <- rep(c(0.25, 0.75, 0.75, 0.25, 0.75, 0.25, 0.75, 0.25),8)
par910 <- rep(c(0.25, 0.75, 0.75, 0.25, 0.75, 0.25, 0.25, 0.75),8)



Matrix <- c(8, 10, 4, 9, 7, 3, 1, 2, 5, 6,
            0, 10, 7, 4, 6, 3, 1, 2, 5, 9,
            0, 0, 4, 7, 6, 9, 5, 1, 2, 3,
            0, 0, 0, 7, 9, 3, 1, 2, 5, 6,
            0, 0, 0, 0, 9, 6, 3, 1, 2, 5,
            0, 0, 0, 0, 0, 6, 3, 1, 2, 5,
            0, 0, 0, 0, 0, 0, 3, 5, 1, 2,
            0, 0, 0, 0, 0, 0, 0, 5, 1, 2,
            0, 0, 0, 0, 0, 0, 0, 0, 1, 2,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 2)

Matrix <- matrix(Matrix, 10, 10) 

adjmatrix <- c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
               1, 0, 1, 0, 1, 0, 0, 0, 0, 0,
               0, 1, 0, 1, 0, 0, 0, 0, 0, 0,
               0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
               0, 1, 0, 0, 0, 1, 0, 0, 1, 0,
               0, 0, 0, 0, 1, 0, 1, 1, 0, 0,
               0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
               0, 0, 0, 0, 1, 0, 0, 0, 0, 1,
               0, 0, 0, 0, 0, 0, 0, 0, 1, 0)

adjmatrix <- matrix(adjmatrix, 10, 10)


A=list()
B=list()
C=list()

m <- length(par12) # Number of scenarios.
M <- 1 # repetition of seeds.
N <- 15000 # sample size.



for(j in 1:length(par12)){
  
  
  # define R-vine family matrix
  family <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, family68[j], 
              0, 0, 0, 0, 0, 0, 0, 0, 0, family910[j],
              0, 0, 0, 0, 0, 0, 0, 0, 0, family34[j],
              0, 0, 0, 0, 0, 0, 0, 0, 0, family67[j],
              0, 0, 0, 0, 0, 0, 0, 0, 0, family59[j],
              0, 0, 0, 0, 0, 0, 0, 0, 0, family56[j],
              0, 0, 0, 0, 0, 0, 0, 0, 0, family23[j],
              0, 0, 0, 0, 0, 0, 0, 0, 0, family25[j],
              0, 0, 0, 0, 0, 0, 0, 0, 0, family12[j],
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  
  family <- matrix(family, 10, 10)
  
  
  # define R-vine parameter matrix
  par <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, BiCopTau2Par(family68[j],par68[j]),
           0, 0, 0, 0, 0, 0, 0, 0, 0, BiCopTau2Par(family910[j],par910[j]),
           0, 0, 0, 0, 0, 0, 0, 0, 0, BiCopTau2Par(family34[j],par34[j]),
           0, 0, 0, 0, 0, 0, 0, 0, 0, BiCopTau2Par(family67[j],par67[j]),
           0, 0, 0, 0, 0, 0, 0, 0, 0, BiCopTau2Par(family59[j],par59[j]),
           0, 0, 0, 0, 0, 0, 0, 0, 0, BiCopTau2Par(family56[j],par56[j]),
           0, 0, 0, 0, 0, 0, 0, 0, 0, BiCopTau2Par(family23[j],par23[j]),
           0, 0, 0, 0, 0, 0, 0, 0, 0, BiCopTau2Par(family25[j],par25[j]),
           0, 0, 0, 0, 0, 0, 0, 0, 0, BiCopTau2Par(family12[j],par12[j]),
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  
  par <- matrix(par, 10, 10)
  
  # define second R-vine pair-copula parameter matrix
  par2 <- matrix(5, 10, 10)
  
  # define RVineMatrix object
  RVM <- RVineMatrix(Matrix = Matrix, family = family,
                     par = par, par2 = par2,
                     names = c("U1", "U2", "U3", "U4","U5",
                               "U6", "U7", "U8", "U9","U10"))
  
  
  for(i in 1:M){
    
    
    # simulate from the D-vine.
    set.seed(i+j*30)
    databbn2 <- RVineSim(N, RVM)
    
    
    # We then check what are the results that we obtain with Gaussian Data:
    # Transforming the margins, and computing logLik
    gdatabbn2 <- qnorm(data.matrix(databbn2),0,1)
    logLik <- sum(dmvnorm(gdatabbn2, colMeans(gdatabbn2), cov(gdatabbn2), log=TRUE))
    
    
    gdatabbn2 <- data.frame(gdatabbn2)
    colnames(gdatabbn2) = c("X1", "X2", "X3", "X4","X5",
                            "X6", "X7", "X8", "X9","X10")
    
    
    # Constrain-based algorithm: Hiton-PC:
    pc.start.time <- Sys.time()
    
    cpdag.pc <- si.hiton.pc(gdatabbn2, undirected = FALSE, alpha = 0.05, test="cor")
    pc.shd <- shd(cpdag(dag), cpdag.pc, wlbl = FALSE, debug = FALSE)
    pc.ham <- hamming(dag, cpdag.pc)
    
    pc.arcs <- as.numeric(gsub(" arcs:", "", head(grep("arcs", capture.output(cpdag.pc), value = TRUE), 1)))
    #graphviz.plot(cpdag.pc)
    pc.final.time <- Sys.time()
    pc.sim.time <- difftime(pc.final.time, pc.start.time, units="secs")
    
    
    
    
    # Score-based algorithms: Hill Climbing:
    hc.start.time <- Sys.time()
    
    cpdag.hc <- hc(gdatabbn2, score = "bic-g",  restart = 15)
    hc.shd <- shd(cpdag(dag),cpdag(cpdag.hc), wlbl = FALSE, debug = FALSE)
    hc.ham <- hamming(dag,cpdag.hc)
    
    hc.arcs <- as.numeric(gsub(" arcs:", "", head(grep("arcs", capture.output(cpdag.hc), value = TRUE), 1)))
    #graphviz.plot(cpdag.hc)
    hc.final.time <- Sys.time()
    hc.sim.time <- difftime(hc.final.time, hc.start.time, units="secs")
    
    
    
    
    # Graphical Lasso:
    glasso.start.time <- Sys.time()
    
    gdatabbn2.mat <- data.matrix(gdatabbn2)
    GRest <- selectFast(gdatabbn2.mat , family=c("LA"))
    
    #g <- network(GRest$LA$G)
    #plot(g, jitter=TRUE, usearrows = FALSE, label=1:p,displaylabels=TRUE)
    
    glasso.ham <- (100-sum(GRest$LA$G==adjmatrix))/2
    glasso.arcs <- sum(GRest$LA$G)/2
    
    glasso.final.time <- Sys.time()
    glasso.sim.time <- difftime(glasso.final.time, glasso.start.time, units="secs")
    
    
  }
  
  
  A[[j]]=GRest$LA$G
  B[[j]]=amat(skeleton(cpdag.pc))
  C[[j]]=amat(skeleton(cpdag.hc))
  
  print(j)
  
}


write.csv(A,"C:/Users/amade/OneDrive/Escritorio/Master Thesis Applied Maths/data files//A glasso big network gauss.csv", row.names = FALSE)
write.csv(B,"C:/Users/amade/OneDrive/Escritorio/Master Thesis Applied Maths/data files//B pc big network gauss.csv", row.names = FALSE)
write.csv(C,"C:/Users/amade/OneDrive/Escritorio/Master Thesis Applied Maths/data files//C hc big network gauss.csv", row.names = FALSE)


A.bignetworkgauss=A
B.bignetworkgauss=B
C.bignetworkgauss=C

