# Simulation study chapter 7.
# Diamond shape BN.

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


# Simulation from a non chordal graph.
par(mfrow=c(1,1))
X <- paste0("[X1][X2|X1][X3|X1][X4|X2:X3]")
dag <- model2network(X)
dag
graphviz.plot(dag)


# 80 different scenarios:
# families:
family12 <- rep(c(2,3,4,5,6,2,4,2,4,6),each=8)
family13 <- rep(c(2,3,4,5,6,2,4,3,5,2),each=8)
family34 <- rep(c(2,3,4,5,6,3,5,4,6,3),each=8)
family24.3 <- rep(c(2,3,4,5,6,3,5,5,2,4),each=8)


# parameters:
par12 <- rep(c(0.25, 0.75, 0.25, 0.75, 0.25, 0.75, 0.25, 0.75),10)
par13 <- rep(c(0.25, 0.75, 0.25, 0.75, 0.25, 0.75, 0.75, 0.25),10)
par34 <- rep(c(0.25, 0.75, 0.75, 0.25, 0.25, 0.75, 0.25, 0.75),10)
par24.3 <- rep(c(0.25, 0.75, 0.75, 0.25, 0.75, 0.25, 0.75, 0.25),10)



m <- length(par12) # Number of scenarios.
M <- 1 # repetition of seeds.
N <- 2500 # sample size.

DF <- data.frame(matrix(NA, nrow = length(par12), ncol = 19))




for(i in 1:length(par12)){
  
  # SIMULATION.
  start.sim.time <- Sys.time()
  
  # We first simulate from the first 3 variables:
  Cmatrix <- c(3, 2, 1,
               0, 2, 1,
               0, 0, 1)
  Cmatrix <- matrix(Cmatrix, 3, 3)
  
  # with families:
  family <- c(0, 0, family13[i],  
              0, 0, family12[i],  
              0, 0, 0)
  family <- matrix(family, 3, 3)
  
  
  
  # and parameters:
  par <- c(0, 0, BiCopTau2Par(family13[i],par13[i]),
           0, 0, BiCopTau2Par(family12[i],par12[i]), 
           0, 0, 0)
  par <- matrix(par, 3, 3)
  par2 <- matrix(5, 3, 3)
  
  RVM <- RVineMatrix(Matrix = Cmatrix, family = family,
                     par = par, par2 = par2,
                     names = c("U1", "U2", "U3"))
  
  set.seed(i)
  databbn1 <- RVineSim(N, RVM)
  
  # Lets compute the last column:
  databbn1 <- cbind(databbn1,rep(NA,N))
  colnames(databbn1) <- c("U1","U2","U3","U4")
  
  NDIM <- 1
  NCOMP <- 1
  
  
  # Integrations for simulating:
  for(j in 1:N){
    
    integrand <- function(y1) {
      ff <- BiCopPDF(y1, databbn1[j,3], family=family13[i],
                     par= BiCopTau2Par(family13[i],par13[i]),
                     par2=5, obj = NULL, check.pars = TRUE)*
        BiCopHfunc1(y1, databbn1[j,2], family=family12[i], 
                    par=BiCopTau2Par(family12[i],par12[i]),
                    par2=5, obj = NULL, check.pars = TRUE)
      return(ff)
    } 
    
    a <- cuhre(f = integrand,lowerLimit = rep(0, NDIM),upperLimit = rep(1, NDIM),
               relTol = 1e-2, absTol= 1e-9,flags = list(verbose = 0, final = 0))$integral
    
    if(a>1){
      a=1
    }
    
    b <- BiCopHinv1(a,runif(1),family=family24.3[i], par=BiCopTau2Par(family24.3[i],par24.3[i]),
                    par2 = 5, obj = NULL, check.pars = TRUE)
    
    
    databbn1[j,4] <- BiCopHinv1(databbn1[j,3],b,family=family34[i], par=BiCopTau2Par(family34[i],par34[i]),
                                par2 = 5, obj = NULL, check.pars = TRUE)
  }
  
  
  final.sim.time <- Sys.time()
  Sim.time <- difftime(final.sim.time,start.sim.time,units="secs")
  
  
  
  # TRUE LOGLIKELIHOOD.
  # We compute the True logLikelihood with the True Order:
  start.TruelogLik <- Sys.time()
  
  set.seed(123)
  C12logLik <- sum(log(BiCopPDF(databbn1[,1],databbn1[,2],family=family12[i],
                                par=BiCopTau2Par(family12[i],par12[i]), par2=5)))
  C13logLik <- sum(log(BiCopPDF(databbn1[,1],databbn1[,3],family=family13[i],
                                par=BiCopTau2Par(family13[i],par13[i]), par2=5)))
  C34logLik <- sum(log(BiCopPDF(databbn1[,3],databbn1[,4],family=family34[i],
                                par=BiCopTau2Par(family34[i],par34[i]), par2=5)))
  
  u2.3.true <- c()
  
  for (j in 1:N){
    
    integrand <- function(y1) {
      ff <- BiCopPDF(y1, databbn1[j,3], family=family13[i], par=BiCopTau2Par(family13[i],par13[i]), par2=5)*
        BiCopHfunc1(y1, databbn1[j,2], family=family12[i], par=BiCopTau2Par(family12[i],par12[i]), par2=5,
                    check.pars = TRUE) 
      return(ff)
    } 
    
    u2.3.true[j] <- cuhre(f = integrand, lowerLimit = rep(0, NDIM),
                          upperLimit = rep(1, NDIM), relTol = 1e-2, 
                          absTol= 1e-9,flags = list(verbose = 0, final = 0))$integral
    
    if(u2.3.true[j]>1){
      u2.3.true[j] <- 1
    }
  } 
  
  u4.3.true <- BiCopHfunc1(databbn1[,3], databbn1[,4], family=family34[i], 
                           par=BiCopTau2Par(family34[i],par34[i]), par2=5, obj = NULL, check.pars = TRUE)
  
  C24.3logLik <- sum(log(BiCopPDF(u2.3.true,u4.3.true,family=family24.3[i],
                                  par=BiCopTau2Par(family24.3[i],par24.3[i]), par2=5, check.pars = TRUE)))
  
  TruelogLik <- C12logLik + C13logLik + C34logLik + C24.3logLik
  
  
  final.TruelogLik <- Sys.time()
  Comp.TruelogLik <- difftime(final.TruelogLik, start.TruelogLik, units="secs")
  
  
  
  # CD VINE DECOMPOSITION.
  # The following step is to use the CD-Vine decomposition score:
  start.CDVine <- Sys.time()
  
  C12.cd <- BiCopSelect(databbn1[,1], databbn1[,2], family = c(0,1,2,3:6), 
                        selectioncrit = "logLik", indeptest = TRUE, 
                        level = 0.05, weights = NA, rotations = TRUE,
                        se = FALSE,presel = TRUE,method = "mle")
  C13.cd <- BiCopSelect(databbn1[,1], databbn1[,3], family = c(0,1,2,3:6), 
                        selectioncrit = "logLik", indeptest = TRUE, 
                        level = 0.05, weights = NA, rotations = TRUE,
                        se = FALSE,presel = TRUE,method = "mle")
  C23.cd <- BiCopSelect(databbn1[,2], databbn1[,3], family = c(0,1,2,3:6), 
                        selectioncrit = "logLik", indeptest = TRUE, 
                        level = 0.05, weights = NA, rotations = TRUE,
                        se = FALSE,presel = TRUE,method = "mle")
  
  
  CMatrix.cd <- c(1, 2, 3,
                0, 2, 3,
                0, 0, 3)
  CMatrix.cd <- matrix(CMatrix.cd, 3, 3)
  RVMcop.cd <- RVineCopSelect(databbn1[,2:4], familyset = c(0,1,2,3:6), CMatrix.cd, 
                              selectioncrit = "logLik", indeptest = TRUE,
                              level = 0.05, trunclevel = NA,
                              weights = NA, rotations = TRUE, se = FALSE,
                              presel = TRUE, method = "mle", cores = 1)
  
  
  logLik.cd <-  (C12.cd$logLik + C13.cd$logLik +
                   RVMcop.cd$logLik - C23.cd$logLik)
  
  
  final.CDVine <- Sys.time()
  Comp.CDVine <- difftime(final.CDVine, start.CDVine, units="secs")
  
  
  
  
  
  # APPROXIMATED LOGLIKELIHOOD TRUE ORDER
  # Approximated loglikelihood dodging the integrals
  # and using the correct order:
  start.estcop.to <- Sys.time()
  C13.estcop.to <- BiCopSelect(databbn1[,1], databbn1[,3], family = c(0,1,2,3:6), 
                               selectioncrit = "logLik", indeptest = TRUE, 
                               level = 0.05, weights = NA, rotations = TRUE,
                               se = FALSE,presel = TRUE,method = "mle")
  C12.estcop.to <- BiCopSelect(databbn1[,1], databbn1[,2], family = c(0,1,2,3:6), 
                               selectioncrit = "logLik", indeptest = TRUE, 
                               level = 0.05, weights = NA, rotations = TRUE,
                               se = FALSE,presel = TRUE,method = "mle")
  C23.estcop.to <- BiCopSelect(databbn1[,2], databbn1[,3], family = c(0,1,2,3:6), 
                               selectioncrit = "logLik", indeptest = TRUE, 
                               level = 0.05, weights = NA, rotations = TRUE,
                               se = FALSE,presel = TRUE,method = "mle")
  
  CMatrix.estcop.to <- c(1, 3, 2,
                         0, 3, 2,
                         0, 0, 2)
  CMatrix.estcop.to <- matrix(CMatrix.estcop.to, 3, 3)
  RVMcop.estcop.to <- RVineCopSelect(databbn1[,2:4], familyset = c(0,1,2,3:6), 
                                     CMatrix.estcop.to, selectioncrit = "logLik",
                                     indeptest = TRUE, level = 0.05, trunclevel = NA,
                                     weights = NA, rotations = TRUE, se = FALSE,
                                     presel = TRUE, method = "mle", cores = 1)
  
  
  # checking the value of the empirical tau between 3 and 4:
  tau.34.to <- RVMcop.estcop.to$emptau[3,2]
  logLik.estcop.to <- (C13.estcop.to$logLik + C12.estcop.to$logLik +
                         RVMcop.estcop.to$logLik - C23.estcop.to$logLik)
  
  
  
  final.estcop.to <- Sys.time()
  Comp.estcop.to <- difftime(final.estcop.to, start.estcop.to, units="secs")
  
  
  
  # APPROXIMATED LOGLIKELIHOOD WRONG ORDER:
  # Approximated loglikelihood dodging the integrals
  # and using the wrong order:
  start.estcop.wo <- Sys.time()
  
  C13.estcop.wo <-  BiCopSelect(databbn1[,1], databbn1[,3], family = c(0,1,2,3:6), 
                                selectioncrit = "logLik", indeptest = TRUE, 
                                level = 0.05, weights = NA, rotations = TRUE,
                                se = FALSE,presel = TRUE,method = "mle")
  
  C12.estcop.wo <- BiCopSelect(databbn1[,1], databbn1[,2], family = c(0,1,2,3:6), 
                               selectioncrit = "logLik", indeptest = TRUE, 
                               level = 0.05, weights = NA, rotations = TRUE,
                               se = FALSE,presel = TRUE,method = "mle")
  
  C23.estcop.wo <- BiCopSelect(databbn1[,2], databbn1[,3], family = c(0,1,2,3:6), 
                               selectioncrit = "logLik", indeptest = TRUE, 
                               level = 0.05, weights = NA, rotations = TRUE,
                               se = FALSE,presel = TRUE,method = "mle")
  
  
  CMatrix.estcop.wo <- c(3, 2, 1,
                         0, 2, 1,
                         0, 0, 1)
  
  CMatrix.estcop.wo <- matrix(CMatrix.estcop.wo, 3, 3)
  RVMcop.estcop.wo <- RVineCopSelect(databbn1[,2:4], familyset = c(0,1,2,3:6), CMatrix.estcop.wo, 
                                     selectioncrit = "logLik", indeptest = TRUE,
                                     level = 0.05, trunclevel = NA,
                                     weights = NA, rotations = TRUE, se = FALSE,
                                     presel = TRUE, method = "mle", cores = 1)
  
  
  # checking the value of the empirical tau between 2 and 4
  # to check how good is the heuristic
  tau.24.wo  <- RVMcop.estcop.wo$emptau[3,1]
  logLik.estcop.wo <- (C13.estcop.wo$logLik + C12.estcop.wo$logLik +
                         RVMcop.estcop.wo$logLik - C23.estcop.wo$logLik)
  
  
  final.estcop.wo<- Sys.time()
  Comp.estcop.wo <- difftime(final.estcop.wo, start.estcop.wo, units="secs")
  
  
  
  # Storing all the previously computed values:
  DF[i,] = round(c(family12[i], family13[i], family34[i], family24.3[i],
                   par12[i], par13[i], par34[i], par24.3[i], Sim.time,
                   TruelogLik, Comp.TruelogLik,
                   logLik.cd, Comp.CDVine,
                   logLik.estcop.to, tau.34.to, Comp.estcop.to,
                   logLik.estcop.wo, tau.24.wo, Comp.estcop.wo),
                 3)
                  
                 
  print(i)
  
  
}

DF[,1:4][DF[,1:4]==1]="N"
DF[,1:4][DF[,1:4]==2]="t"
DF[,1:4][DF[,1:4]==3]="C"
DF[,1:4][DF[,1:4]==4]="G"
DF[,1:4][DF[,1:4]==5]="F"
DF[,1:4][DF[,1:4]==6]="J"


colnames(DF) <- c("C12", "C13", "C34",  "C24.3",
                  "$\tau_{1,2}$", "$\tau_{1,3}$", "$\tau_{3,4}$",
                  "\tau_{2,4|3}", "Sim.time", "TruelogLik", "Comp.TruelogLik",
                  "logLik.cd", "Comp.CDVine",
                  "logLik.estcop.to", "tau.34.to", "Comp.estcop.to",
                  "logLik.estcop.wo", "tau.24.wo", "Comp.estcop.wo")

# Storing the results
write.csv(DF,"C:/Users/amade/OneDrive/Escritorio/Master Thesis Applied Maths/data files//diamond shape simulation.csv", row.names = FALSE)

# We define new columns of the data frame:
DF["Rel.error.CD"]=abs(DF["TruelogLik"]-DF["logLik.cd"])/DF["TruelogLik"]*100
DF["Rel.error.estcop.to"]=abs(DF["TruelogLik"]-DF["logLik.estcop.to"])/DF["TruelogLik"]*100
DF["Rel.error.estcop.wo"]=abs(DF["TruelogLik"]-DF["logLik.estcop.wo"])/DF["TruelogLik"]*100
DF["Rel.error.estcop.po"]=abs(DF["TruelogLik"]-ifelse(DF$tau.34.to > DF$tau.24.wo, DF$logLik.estcop.to, DF$logLik.estcop.wo))/DF["TruelogLik"]*100


sum(DF$tau.34.to > DF$tau.24.wo)


colnames(DF) <- c("C12", "C13", "C34",  "C24.3",
                  "$\tau_{1,2}$", "$\tau_{1,3}$", "$\tau_{3,4}$",
                  "\tau_{2,4|3}", "Sim.time", "TruelogLik", "Comp.TruelogLik",
                  "logLik.cd", "Comp.CDVine",
                  "logLik.estcop.to", "tau.34.to", "Comp.estcop.to",
                  "logLik.estcop.wo", "tau.24.wo", "Comp.estcop.wo",
                  "Rel.error.CD","Rel.error.estcop.to","Rel.error.estcop.wo","Rel.error.estcop.po")

# Analysis of the results:
# Relative errors:
length(which(DF["Rel.error.CD"]<DF["Rel.error.estcop.to"]))
# Only 5 out of 80, which are the t-copula case.
length(which(DF["Rel.error.CD"]<DF["Rel.error.estcop.wo"])%%8)
# 45/80 scenarios.
length(which(DF["Rel.error.CD"]<DF["Rel.error.estcop.po"])%%8)
# 19/80 scenarios -> among them it seems that scenario 2, fails.
DF[seq(2,80,8),"tau.24.wo"]
# The reason behind it is that they always choose the wrong order for this case.
length(which(DF["Rel.error.estcop.wo"]<DF["Rel.error.estcop.to"]))



# Let's check spread of the distribution, and make a table out of it:
kable(round(rbind(quantile(DF$Rel.error.CD, probs = c(0.05,0.25,0.5,0.75,0.95)),
                  quantile(DF$Rel.error.estcop.to, probs = c(0.05,0.25,0.5,0.75,0.95)),
                  quantile(DF$Rel.error.estcop.wo, probs = c(0.05,0.25,0.5,0.75,0.95)),
                  quantile(DF$Rel.error.estcop.po, probs = c(0.05,0.25,0.5,0.75,0.95))),3),
      "latex")
  
  
colnames(DF)[20:23]=c("CD-Vine","True order",
                      "Wrong order","Parental order")

# Let's make density plots.
# First comparison between CD and the estimation of copulas true order
DF.relerror <- melt(DF[c("CD-Vine","True order")],id=NULL)
median.relerror <- ddply(DF.relerror, "variable", summarise, median=round(median(value),3))


ggplot(DF.relerror , aes(x=value, color=variable, fill=variable)) +
  geom_density(alpha=.2) + scale_x_continuous(name="Relative error (%)", limits=c(0, 10))+
  scale_fill_manual(values=c("#F8766D", "#00BFC4"))+
  scale_color_manual(values = c("#F8766D", "#00BFC4"))+
  geom_vline(data=median.relerror, aes(xintercept=median, color=variable),linetype="dashed")+
  ggtitle(paste("Relative error: CD-Vine vs Extra estimation \n copulas True order, N=",N))+
  theme(plot.title = element_text(hjust = 0.5))+
  annotate("text", color="red",x = 5, y = 0.5, label = paste ("Median ==", median.relerror$median[1]), parse = TRUE) +
  annotate("text", color="blue",x = 5, y = 0.43, label = paste ("Median ==",  median.relerror$median[2]), parse = TRUE) 


# Second comparison between CD and the estimation of copulas parental order
DF.relerror <- melt(DF[c("CD-Vine","Parental order")],id=NULL)
median.relerror <- ddply(DF.relerror, "variable", summarise, median=round(median(value),3))


ggplot(DF.relerror , aes(x=value, color=variable, fill=variable)) +
  geom_density(alpha=.2) + scale_x_continuous(name="Relative error (%)", limits=c(0, 10))+
  scale_fill_manual(values=c("#F8766D", "#00BA38"))+
  scale_color_manual(values = c("#F8766D", "#00BA38"))+
  geom_vline(data=median.relerror, aes(xintercept=median, color=variable),linetype="dashed")+
  ggtitle(paste("Relative error: CD-Vine vs Extra estimation \n copulas chosen Parental order, N=",N))+
  theme(plot.title = element_text(hjust = 0.5))+
  annotate("text", color="red",x = 7, y = 0.25, label = paste ("Median ==", median.relerror$median[1]), parse = TRUE) +
  annotate("text", color="#00BA38",x = 7, y = 0.22, label = paste ("Median ==",  median.relerror$median[2]), parse = TRUE) 



# Let's investigate computational times.
kable(round(rbind(quantile(DF$Comp.TruelogLik, probs = c(0.05,0.25,0.5,0.75,0.95)),
                  quantile(DF$Comp.CDVine, probs = c(0.05,0.25,0.5,0.75,0.95)),
                  quantile(DF$Comp.estcop.to, probs = c(0.05,0.25,0.5,0.75,0.95)),
                  quantile(DF$Comp.estcop.wo, probs = c(0.05,0.25,0.5,0.75,0.95))),3),
      "latex")





# First the CD-Vine and our approximation:
DF.times <- melt(DF[c("Comp.CDVine","Comp.estcop.to")],id=NULL)
median.times <- ddply(DF.times, "variable", summarise, median=round(median(value),3))

ggplot(DF.times , aes(x=value, color=variable, fill=variable)) +
  geom_density(alpha=.2) + scale_x_continuous(name="Time (sec.)", limits=c(4,16))+
  scale_fill_manual(values=c("#F8766D", "#00BFC4"))+
  scale_color_manual(values = c("#F8766D", "#00BFC4"))+
  geom_vline(data=median.times, aes(xintercept=median, color=variable),linetype="dashed")+
  ggtitle(paste("Computational Time: CD-Vine vs Extra estimation \n copulas True order, N=",N))+
  theme(plot.title = element_text(hjust = 0.5))+
  annotate("text", color="red",x = 12, y = 0.25, label = paste ("Median ==", median.times$median[1]), parse = TRUE) +
  annotate("text", color="blue",x = 12, y = 0.22, label = paste ("Median ==",  median.times$median[2]), parse = TRUE) 


# Second the true logLikelihood:
DF.times <- melt(DF[c("Comp.TruelogLik")],id=NULL)
median.times <- ddply(DF.times, "variable", summarise, median=round(median(value),3))

ggplot(DF.times , aes(x=value, color=variable, fill=variable)) +
  geom_density(alpha=.2) + scale_x_continuous(name="Time (sec.)", limits=c(50,500))+
  scale_fill_manual(values=c("#FF9900"))+
  scale_color_manual(values = c("#FF9900"))+
  geom_vline(data=median.times, aes(xintercept=median, color=variable),linetype="dashed")+
  ggtitle(paste("Computational Time True logLikelihood, N=",N))+
  theme(plot.title = element_text(hjust = 0.5))+
  annotate("text", color="#FF9900",x = 400, y = 0.004, label = paste ("Median ==", median.times$median), parse = TRUE) 
  





length(which(abs(DF["TruelogLik"]-DF["logLik.cd"])>abs(DF["TruelogLik"]-DF["logLik.estcop.to"])))
which(abs(DF["TruelogLik"]-DF["logLik.cd"])>abs(DF["TruelogLik"]-DF["logLik.estcop.wo"]))


# Let's see our proposed score:
which(abs(DF["TruelogLik"]-DF["logLik.cd"])<abs(DF["TruelogLik"]- ifelse(DF$tau.34.to > DF$tau.24.wo, DF$logLik.estcop.to, DF$logLik.estcop.wo)))
# We get that in 61 out of 80 scenarios performs better than the CD-vine approx.


abs(DF["TruelogLik"]-ifelse(DF$tau.34.to > DF$tau.24.wo, DF$logLik.estcop.to, DF$logLik.estcop.wo))/DF["TruelogLik"]*100
abs(DF["TruelogLik"]-DF["logLik.cd"])/DF["TruelogLik"]*100

