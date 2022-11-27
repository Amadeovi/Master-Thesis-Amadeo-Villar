# Simulation Study Chapter 6: big network Gaussian data.

# Experiment 11/10.
# Investigate how Gaussian methods perform for Gaussian Data.
# Bigger Network.


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


m <- length(par12) # Number of scenarios.
M <- 30 # repetition of seeds.
N <- 15000 # sample size.


DF <- data.frame(matrix(NA, nrow = length(par12), ncol = 30))
colnames(DF) <- c("C12", "C23", "C34", "C25", "C56", "C67", "C68", "C59","C910",
                  "$\tau_{1,2}$", "$\tau_{2,3}$", "$\tau_{3,4}$", "$\tau_{2,5}$",
                  "$\tau_{5,6}$","$\tau_{6,7}$","$\tau_{6,8}$", 
                  "$\tau_{5,9}$", "$\tau_{9,10}$",
                  "pc.shd", "pc.ham", "pc.arcs", "pc.sim.time",
                  "hc.shd", "hc.ham", "hc.arcs", "hc.sim.time",
                  "glasso.ham", "glasso.edges","glasso.sim.time","True.glasso")




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
  
  df <- data.frame(matrix(NA, nrow = M, ncol = 11))
  colnames(df) <- c("pc.shd", "pc.ham", "pc.arcs", "pc.sim.time",
                    "hc.shd", "hc.ham", "hc.arcs", "hc.sim.time",
                    "glasso.ham", "glasso.arcs", "glasso.sim.time")
  
  
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
    
    
    # Store the results.
    df[i,] = c(pc.shd, pc.ham, pc.arcs, pc.sim.time,
               hc.shd, hc.ham, hc.arcs, hc.sim.time,
               glasso.ham, glasso.arcs, glasso.sim.time)
    
    
  }
  
  DF[j,] = round(c(family12[j], family23[j], family34[j], family25[j], family56[j],
                   family67[j], family68[j], family59[j], family910[j],
                   par12[j], par23[j], par34[j], par25[j], par56[j],
                   par67[j], par68[j], par59[j], par910[j], 
                   sum(df["pc.shd"])/M, sum(df["pc.ham"])/M, sum(df["pc.arcs"])/M, sum(df["pc.sim.time"])/M,
                   sum(df["hc.shd"])/M, sum(df["hc.ham"])/M, sum(df["hc.arcs"])/M, sum(df["hc.sim.time"])/M,
                   sum(df["glasso.ham"])/M, sum(df["glasso.arcs"])/M,  sum(df["glasso.sim.time"])/M,
                   sum(df["glasso.ham"]==0)), 3)
  
  print(j)
  
  
}


DF[,1:9][DF[,1:9]==1]="N"
DF[,1:9][DF[,1:9]==2]="t"
DF[,1:9][DF[,1:9]==3]="C"
DF[,1:9][DF[,1:9]==4]="G"
DF[,1:9][DF[,1:9]==5]="F"
DF[,1:9][DF[,1:9]==6]="J"


write.csv(DF,"C:/Users/amade/OneDrive/Escritorio/Master Thesis Applied Maths/data files//bignetwork Gauss15000.csv", row.names = FALSE)






