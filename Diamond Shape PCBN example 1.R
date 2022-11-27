# Diamond Shape PCBN: used in Chapter 4 PCBNs 
# and to introduction Chapter 6 and 7.

# Experiment 23/06
# We simulate from a non-chordal graph.
# We then analyse the data in depth.


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


# We construct the real BBN:
X <- paste0("[X1][X2|X1][X3|X1][X4|X2:X3]")
dag <- model2network(X)
dag
graphviz.plot(dag)

# We check whether the cond.ind holds
dsep(dag, x = "X2", y = "X3", z = "X1")
dsep(dag, x = "X1", y = "X4", z = c("X2","X3"))

# Only 1 scenario:
family12 <- 3; par12 <- c(BiCopTau2Par(3,0.5))
family13 <- 4; par13 <- c(BiCopTau2Par(4,0.75))
family34 <- 5; par34 <- c(BiCopTau2Par(5,0.5))
family24.3 <- 6; par24.3 <- c(BiCopTau2Par(6,0.75))


###############################
#Simulations:
start.time <- Sys.time()

Cmatrix <- c(3, 2, 1,
             0, 2, 1,
             0, 0, 1)

family <- c(0, 0, family13,  
            0, 0, family12,  
            0, 0, 0)
family <- matrix(family, 3, 3)

par <- c(0, 0, par13,
         0, 0, par12, 
         0, 0, 0)
par <- matrix(par, 3, 3)
par2 <- matrix(0, 3, 3)

RVM <- RVineMatrix(Matrix = Cmatrix, family = family,
                   par = par, par2 = par2,
                   names = c("U1", "U2", "U3"))

# simulate from the Bayesian Network:
set.seed(123)
n <- 2500
databbn1 <- RVineSim(n, RVM)


# Lets compute the last column:
databbn1 <- cbind(databbn1,rep(NA,n))
colnames(databbn1) <- c("U1","U2","U3","U4")

NDIM <- 1
NCOMP <- 1

for(j in 1:n){
  
  integrand <- function(y1) {
    ff <- BiCopPDF(y1, databbn1[j,3], family=family13,
                   par=par13, par2=0, obj = NULL, check.pars = TRUE)*
      BiCopHfunc1(y1, databbn1[j,2], family=family12, 
                  par=par12, par2=0, obj = NULL, check.pars = TRUE)
    return(ff)
  } 
  
  a <- cuhre(f = integrand,lowerLimit = rep(0, NDIM),upperLimit = rep(1, NDIM),
             relTol = 1e-2, absTol= 1e-9,flags = list(verbose = 0, final = 0))$integral
  
  if(a>1){
    a=1
  }
  
  b <- BiCopHinv1(a,runif(1),family=family24.3, par=par24.3,
                  par2 = 0, obj = NULL, check.pars = TRUE)
  
 
  databbn1[j,4] <- BiCopHinv1(databbn1[j,3],b,family=family34, par=par34,
                              par2 = 0, obj = NULL, check.pars = TRUE)
}

final.time <- Sys.time()
Sim.time <- difftime(final.time,start.time,units="secs")



########################################
# Analysis of the simulated data:
# First, visual analysis:

cdata <- as.copuladata(databbn1)
pairs(cdata)

data <- data.frame(databbn1)

# For plotting we define:
size.dot = 0.6
family.list <- list("0" = "Ind", "1" = "Gaussian", "3" ="Clayton",
          "4" = "Gumbel", "5" = "Frank", "6" = "Joe",
          "13" = "rotated Clayton 180", "14" = "rotated Gumbel 180",
          "16" = "rotated Joe 180", "23" = "rotated Clayton 90",
          "24" = "rotated Gumbel 90", "26" = "rotated Joe 90",
          "33" = "rotated Clayton 270", "34" = "rotated Gumbel 270",
          "36" = "rotated Joe 270")


# Plotting the non conditional copulas:
ggplot(data, aes(x = U1, y = U2))+geom_point(size=size.dot)+
   labs(title = TeX(paste("$C_{1,2}$",":",  family.list[as.character(family12)],",",
        '$\\tau$',"=",BiCopPar2Tau(family12,par12))),
        x = TeX("$U_{1}$"),y = TeX("$U_{2}$"))+theme(plot.title = element_text(hjust = 0.5))

ggplot(data, aes(x = U1, y = U3))+geom_point(size=size.dot)+
  labs(title = TeX(paste("$C_{1,3}$",":" , family.list[as.character(family13)],",",
       '$\\tau$',"=",BiCopPar2Tau(family13,par13))),
       x = TeX("$U_{1}$"),y = TeX("$U_{3}$"))+theme(plot.title = element_text(hjust = 0.5))

ggplot(data, aes(x = U3, y = U4))+geom_point(size=size.dot)+
  labs(title = TeX(paste("$C_{3,4}$",":",  family.list[as.character(family34)],",",
       '$\\tau$',"=",BiCopPar2Tau(family34,par34))),
       x = TeX("$U_{3}$"),y = TeX("$U_{4}$"))+theme(plot.title = element_text(hjust = 0.5))



# Second quantitative analysis: 
# We compute the True logLikelihood.
start.true <- Sys.time()
C12logLik <- sum(log(BiCopPDF(databbn1[,1],databbn1[,2],family=family12,
                              par=par12, par2=0)))
C13logLik <- sum(log(BiCopPDF(databbn1[,1],databbn1[,3],family=family13,
                              par=par13, par2=0)))
C34logLik <- sum(log(BiCopPDF(databbn1[,3],databbn1[,4],family=family34,
                              par=par34, par2=0)))

u2.3.true <- c()

for (j in 1:n){
  
  integrand <- function(y1) {
    ff <- BiCopPDF(y1,databbn1[j,3],family=family13,par=par13,par2=0)*
      BiCopHfunc1(y1, databbn1[j,2],family=family12,par=par12,par2=0,
                  check.pars = TRUE) 
    return(ff)
  } 
  
  u2.3.true[j] <- cuhre(f = integrand,lowerLimit = rep(0, NDIM),
                   upperLimit = rep(1, NDIM), relTol = 1e-2, 
                   absTol= 1e-9,flags = list(verbose = 0, final = 0))$integral
  
  if(u2.3.true[j]>1){
    u2.3.true[j] <- 1
  }
} 

u4.3.true <- BiCopHfunc1(databbn1[,3], databbn1[,4], family=family34, 
                       par=par34, par2 = 0,obj = NULL, check.pars = TRUE)

C24.3logLik <- sum(log(BiCopPDF(u2.3.true,u4.3.true,family=family24.3,
                                par=par24.3, par2=0,check.pars = TRUE)))

TruelogLik <- C12logLik + C13logLik + C34logLik + C24.3logLik

final.true <- Sys.time()
Comp.true <- difftime(final.true,start.true,units="secs")
  
  


# We also compute the Estimated True logLikelihood.
start.time1 <- Sys.time()
C12.approx <- BiCopSelect(databbn1[,1], databbn1[,2], family = c(0,1,3:6), 
                          selectioncrit = "logLik", indeptest = TRUE, 
                          level = 0.05, weights = NA, rotations = TRUE,
                          se = FALSE,presel = TRUE,method = "mle")
C13.approx <- BiCopSelect(databbn1[,1], databbn1[,3], family = c(0,1,3:6), 
                          selectioncrit = "logLik", indeptest = TRUE, 
                          level = 0.05, weights = NA, rotations = TRUE,
                          se = FALSE,presel = TRUE,method = "mle")
C34.approx <- BiCopSelect(databbn1[,3], databbn1[,4], family = c(0,1,3:6), 
                          selectioncrit = "logLik", indeptest = TRUE, 
                          level = 0.05, weights = NA, rotations = TRUE,
                          se = FALSE,presel = TRUE,method = "mle")

u2.3 <- c()

for (j in 1:n){
  
  integrand <- function(y1) {
    ff <- BiCopPDF(y1,databbn1[j,3],C13.approx)*
      BiCopHfunc1(y1, databbn1[j,2], family=C12.approx$family, 
                  par=C12.approx$par, par2 = C12.approx$par2,
                  obj = NULL, check.pars = TRUE)
    return(ff)
  } 
  
  u2.3[j] <- cuhre(f = integrand,lowerLimit = rep(0, NDIM),
                   upperLimit = rep(1, NDIM), relTol = 1e-2, 
                   absTol= 1e-9,flags = list(verbose = 0, final = 0))$integral
  
  if(u2.3[j]>1){
    u2.3[j]=1
  }
} 

u4.3 <- BiCopHfunc1(databbn1[,3], databbn1[,4], family=C34.approx$family, 
                    par=C34.approx$par, par2 = C34.approx$par2,
                    obj = NULL, check.pars = TRUE)

C24.3.approx <- BiCopSelect(u2.3, u4.3, family = c(0,1,3:6), 
                            selectioncrit = "logLik", indeptest = TRUE, 
                            level = 0.05, weights = NA, rotations = TRUE,
                            se = FALSE,presel = TRUE,method = "mle")

Estim.truelogLik <- (C13.approx$logLik + C12.approx$logLik +
                 C24.3.approx$logLik + C34.approx$logLik)
Estim.truelogLik24.3 <- C24.3.approx$logLik
final.time1 <- Sys.time()
Comp1.time <- difftime(final.time1,start.time1,units="secs")


# We also plot the conditional copula:
pseudo.obs <- data.frame(u2.3, u4.3)

ggplot(pseudo.obs, aes(x = u2.3, y = u4.3))+geom_point(size=size.dot)+
  labs(title = TeX(paste("$C_{2,4|3}$",":",  family.list[as.character(family24.3)],",",
                         '$\\tau$',"=",round(BiCopPar2Tau(family24.3,par24.3),3))),
       x = TeX("$U_{2|3}$"),y = TeX("$U_{4|3}$"))+
  theme(plot.title = element_text(hjust = 0.5))









###############################
# Approximated loglikelihood dodging the integrals
# and using the correct order:
start.time2 <- Sys.time()
C13.approx <- BiCopSelect(databbn1[,1], databbn1[,3], family = c(0,1,3:6), 
                          selectioncrit = "logLik", indeptest = TRUE, 
                          level = 0.05, weights = NA, rotations = TRUE,
                          se = FALSE,presel = TRUE,method = "mle")

C12.approx <- BiCopSelect(databbn1[,1], databbn1[,2], family = c(0,1,3:6), 
                          selectioncrit = "logLik", indeptest = TRUE, 
                          level = 0.05, weights = NA, rotations = TRUE,
                          se = FALSE,presel = TRUE,method = "mle")


CMatrix1 <- c(1, 3, 2,
              0, 3, 2,
              0, 0, 2)
CMatrix1 <- matrix(CMatrix1, 3, 3)
RVMcop1 <- RVineCopSelect(databbn1[,2:4], familyset = c(0,1,3:6), CMatrix1, 
                          selectioncrit = "logLik", indeptest = TRUE,
                          level = 0.05, trunclevel = NA,
                          weights = NA, rotations = TRUE, se = FALSE,
                          presel = TRUE, method = "mle", cores = 1)

# checking the value of the empirical tau between 3 and 4:
tau34 <- RVMcop1$emptau[3,2]
logLik34 <- RVMcop1$pair.logLik[3,2]

ApproxlogLik <- (C13.approx$logLik + C12.approx$logLik +
                   RVMcop1$pair.logLik[2,1] + RVMcop1$pair.logLik[3,2])

ApproxlogLik24.3 <- RVMcop1$pair.logLik[2,1]

final.time2 <- Sys.time()
Comp2.time <- difftime(final.time2,start.time2,units="secs")




# It is also interesting comparing the true margins 2-3
# to the simulated data using the estimated copula 2-3.

margins <- BiCopSim(n, RVMcop1$family[3,1], RVMcop1$par[3,1],
                    par2 = RVMcop1$par2[3,1],
                    obj = NULL, check.pars = TRUE)

margins <- data.frame(margins)


ggplot(data, aes(x = U2, y = U3))+geom_point(size=size.dot)+
  labs(title = TeX(paste("True margins","$U_{2}$","and","$U_{3}$",
       ",","$\\tau$","=",round(RVMcop1$emptau[3,1],3))),
       x = TeX("$U_{2}$"),y = TeX("$U_{3}$"))+
  theme(plot.title = element_text(hjust = 0.5))


ggplot(margins, aes(x = margins[,1], y = margins[,2]))+geom_point(size=size.dot)+
  labs(title = TeX(paste("Simulated data from","$\\hat{C}_{2,3}$",
       ":",family.list[as.character(RVMcop1$family[3,1])],",",
       "$\\tau$","=",round(BiCopPar2Tau(RVMcop1$family[3,1],RVMcop1$par[3,1]),3))), 
       x = TeX("$U_{2}$") , y = TeX("$U_{3}$"))+theme(plot.title = element_text(hjust = 0.5))
      
# We also compute a goodness of fit test for these data:
BiCopGofTest(databbn1[,2],databbn1[,3],family=RVMcop1$family[3,1],
             par=RVMcop1$par[3,1],par2=RVMcop1$par2[3,1],
             method = "kendall", B = 30) # it rejects the null hypoth.
   



# This will induce an error in the higher trees, let's check
# what happens with these errors:

high.margins <- BiCopSim(n, RVMcop1$family[2,1], RVMcop1$par[2,1],
                         par2 = RVMcop1$par2[2,1],
                         obj = NULL, check.pars = TRUE)

high.margins <- data.frame(high.margins)

ggplot(pseudo.obs, aes(x = u2.3, y = u4.3))+geom_point(size=size.dot)+
  labs(title = TeX(paste("$C_{2,4|3}$",":",  family.list[as.character(family24.3)],",",
                         '$\\tau$',"=",round(BiCopPar2Tau(family24.3,par24.3),3))),
       x = TeX("$U_{2|3}$"),y = TeX("$U_{4|3}$"))+
  theme(plot.title = element_text(hjust = 0.5))


ggplot(high.margins, aes(x = high.margins[,1], y = high.margins[,2]))+
  geom_point(size=size.dot)+labs(title = TeX(paste("Simulated data from","$\\hat{C}_{2,4|3}$",
    ":",family.list[as.character(RVMcop1$family[2,1])],",",
     "$\\tau$","=",round(RVMcop1$emptau[2,1],3))),
      x = TeX("$U_{2|3}$"),y = TeX("$U_{4|3}$"))+
  theme(plot.title = element_text(hjust = 0.5))

# We also compute a goodness of fit test for these data:
BiCopGofTest(u2.3,u4.3,family=RVMcop1$family[2,1],
             par=RVMcop1$par[2,1],par2=RVMcop1$par2[2,1],
             method = "kendall", B = 30) # it rejects the null hypoth.


# Let's also check whether the \hat{u}_{2|3}=\hat{C}_{2|3}(u_{2},u_{3})
# gives a uniform distribution:
hatu2.3 <- BiCopHfunc2(data$U2, data$U3, family = RVMcop1$family[3,1], 
                       par = RVMcop1$par[3,1], par2 = RVMcop1$par2[3,1],
                       obj = NULL, check.pars = TRUE)
            

compare2.3 <- data.frame(u2.3, hatu2.3)
compare2.3dens <- melt(compare2.3[c("u2.3","hatu2.3")],id=NULL)
  
# Comparison between both values:
# Checking uniformity:
ggplot(compare2.3dens, aes(x=value, color=variable, fill=variable)) +
  geom_density(alpha=0.2)+ ylim(0,1.2) + labs(title = TeX(
    paste("Density functions:","$u_{2|3}$"," and ",
          "$\\hat{u}_{2|3}=\\hat{C}_{2|3}(u_{2},u_{3})")))


# Checking difference in values:
ggplot(compare2.3, aes(x = u2.3, y=hatu2.3)) + geom_point(alpha=0.2) +
  labs(title = TeX(paste("Comparison:","$u_{2|3}$"," and ",
          "$\\hat{u}_{2|3}=\\hat{C}_{2|3}(u_{2},u_{3})")),
       x=TeX("$u_{2|3}$"), y=TeX("$\\hat{u}_{2|3}$"))




################################
# Approximation using the CD-Vine decomposition:
start.time4 <- Sys.time()
C13.approx <- BiCopSelect(databbn1[,1], databbn1[,3], family = c(0,1,3:6), 
                          selectioncrit = "logLik", indeptest = TRUE, 
                          level = 0.05, weights = NA, rotations = TRUE,
                          se = FALSE,presel = TRUE,method = "mle")

C12.approx <- BiCopSelect(databbn1[,1], databbn1[,2], family = c(0,1,3:6), 
                          selectioncrit = "logLik", indeptest = TRUE, 
                          level = 0.05, weights = NA, rotations = TRUE,
                          se = FALSE,presel = TRUE,method = "mle")



CMatrix2 <- c(1, 2, 3,
              0, 2, 3,
              0, 0, 3)
CMatrix2 <- matrix(CMatrix2, 3, 3)
RVMcop2 <- RVineCopSelect(databbn1[,2:4], familyset = c(0,1,2,3:6), CMatrix2, 
                          selectioncrit = "logLik", indeptest = TRUE,
                          level = 0.05, trunclevel = NA,
                          weights = NA, rotations = TRUE, se = FALSE,
                          presel = TRUE, method = "mle", cores = 1)

logLik.cd <-  (C13.approx$logLik + C12.approx$logLik +
                 RVMcop2$logLik - RVMcop1$pair.logLik[3,1])


final.time4 <- Sys.time()
Comp4.time <- difftime(final.time4,start.time4,units="secs")



#################################
# Approximated loglikelihood dodging the integrals
# and using the wrong order:
start.time3 <- Sys.time()

C13.approx <- BiCopSelect(databbn1[,1], databbn1[,3], family = c(0,1,3:6), 
                          selectioncrit = "logLik", indeptest = TRUE, 
                          level = 0.05, weights = NA, rotations = TRUE,
                          se = FALSE,presel = TRUE,method = "mle")

C12.approx <- BiCopSelect(databbn1[,1], databbn1[,2], family = c(0,1,3:6), 
                          selectioncrit = "logLik", indeptest = TRUE, 
                          level = 0.05, weights = NA, rotations = TRUE,
                          se = FALSE,presel = TRUE,method = "mle")

CMatrix3 <- c(3, 2, 1,
              0, 2, 1,
              0, 0, 1)
CMatrix3 <- matrix(CMatrix3, 3, 3)
RVMcop3 <- RVineCopSelect(databbn1[,2:4], familyset = c(0,1,3:6), CMatrix3, 
                          selectioncrit = "logLik", indeptest = TRUE,
                          level = 0.05, trunclevel = NA,
                          weights = NA, rotations = TRUE, se = FALSE,
                          presel = TRUE, method = "mle", cores = 1)

# checking the value of the empirical tau between 2 and 4
# to check how good is the heuristic
tau24 <- RVMcop3$emptau[3,1]
logLik24 <- RVMcop3$pair.logLik[3,1]


ApproxlogLik.wo <-  (C13.approx$logLik + C12.approx$logLik +
                      RVMcop3$logLik - RVMcop3$pair.logLik[3,2])

final.time3 <- Sys.time()
Comp3.time <- difftime(final.time3,start.time3,units="secs")



####################################
# Approximated using the Strongly chordal model 
# and then substracting the edge 2-3.
# sc: strongly chordal, asc: approx. using strongly chordal.
start.time5 <- Sys.time()
CMatrix4 <- c(1, 4, 3, 2,
              0, 2, 4, 3,
              0, 0, 3, 4,
              0, 0, 0, 4)
CMatrix4 <- matrix(CMatrix4, 4, 4)
RVMcop4 <- RVineCopSelect(databbn1, familyset = c(0,1,3:6), CMatrix4, 
                          selectioncrit = "logLik", indeptest = TRUE,
                          level = 0.005, trunclevel = NA,
                          weights = NA, rotations = TRUE, se = FALSE,
                          presel = TRUE, method = "mle", cores = 1)

logLik.asc <- RVMcop4$logLik - RVMcop4$pair.logLik[4,2]
logLik.sc <- RVMcop4$logLik

final.time5 <- Sys.time()
Comp5.time <- difftime(final.time5,start.time5,units="secs")





#############################
# Adding the copula 2,3|1 to the model.
start.time6 <- Sys.time()

RVMcop5 <- RVineCopSelect(databbn1[,1:3], familyset = c(0,1,3:6), Cmatrix, 
                          selectioncrit = "logLik",indeptest = TRUE,
                          level = 0.005, trunclevel = NA,
                          weights = NA, rotations = TRUE, se = FALSE,
                          presel = TRUE, method = "mle", cores = 1)

RVMcop1 <- RVineCopSelect(databbn1[,2:4], familyset = c(0,1,3:6), CMatrix1, 
                          selectioncrit = "logLik", indeptest = TRUE,
                          level = 0.005, trunclevel = NA,
                          weights = NA, rotations = TRUE, se = FALSE,
                          presel = TRUE, method = "mle", cores = 1)

approxloglik.23.1 <- (RVMcop5$logLik + RVMcop1$pair.logLik[2,1] + 
                        RVMcop1$pair.logLik[3,2])


final.time6 <- Sys.time()
Comp6.time <- difftime(final.time6,start.time6,units="secs")



##########################
# Adding the copula 1,4|2,3 to the model.
start.time7 <- Sys.time()

RVMcop4 <- RVineCopSelect(databbn1, familyset = c(0,1,3:6), CMatrix4, 
                          selectioncrit = "logLik", indeptest = TRUE,
                          level = 0.005, trunclevel = NA,
                          weights = NA, rotations = TRUE, se = FALSE,
                          presel = TRUE, method = "mle", cores = 1)

approxloglik.14.23 <- (RVMcop4$pair.logLik[2,1]+RVMcop4$pair.logLik[3,2]+
                         RVMcop4$pair.logLik[4,3]+C13.approx$logLik+
                         C12.approx$logLik)


final.time7 <- Sys.time()
Comp7.time <- difftime(final.time7, start.time7, units="secs")


BIC.to <- -2*ApproxlogLik + 4*log(n)
BIC.wo <- -2*ApproxlogLik.wo + 4*log(n)
BIC.sc <- -2*logLik.sc + 5*log(n)




###########################
# Lets now analyse the possible hidden conditional independences
# that might arise from simulations.
CMatrix5 <- c(2, 4, 1, 3,
              0, 4, 3, 1,
              0, 0, 3, 1,
              0, 0, 0, 1)
CMatrix5 <- matrix(CMatrix5, 4, 4)
RVMcop5 <- RVineCopSelect(databbn1, familyset = c(0,1,3:6), CMatrix5, 
                          selectioncrit = "logLik", indeptest = TRUE,
                          level = 0.005, trunclevel = NA,
                          weights = NA, rotations = TRUE, se = FALSE,
                          presel = TRUE, method = "mle", cores = 1)


cond.ind.matrix1 = which(RVMcop5$pair.logLik[lower.tri(RVMcop5$pair.logLik)]==0)
if(length(cond.ind.matrix1) == 0){cond.ind.matrix1 = 0}
cond.ind.matrix1 = cond.ind.matrix1[1]


CMatrix6 <- c(2, 1, 4, 3,
              0, 1, 3, 4,
              0, 0, 3, 4,
              0, 0, 0, 4)
CMatrix6 <- matrix(CMatrix6, 4, 4)
RVMcop6 <- RVineCopSelect(databbn1, familyset = c(0,1,3:6), CMatrix6, 
                          selectioncrit = "logLik", indeptest = TRUE,
                          level = 5*10^-5, trunclevel = NA,
                          weights = NA, rotations = TRUE, se = FALSE,
                          presel = TRUE, method = "mle", cores = 1)


cond.ind.matrix2 = which(RVMcop6$pair.logLik[lower.tri(RVMcop6$pair.logLik)]==0)
if(length(cond.ind.matrix2) == 0){cond.ind.matrix2 = 0}
cond.ind.matrix2 = cond.ind.matrix2[1]


CMatrix7 <- c(3, 4, 1, 2,
              0, 4, 2, 1,
              0, 0, 2, 1,
              0, 0, 0, 1)
CMatrix7 <- matrix(CMatrix7, 4, 4)
RVMcop7 <- RVineCopSelect(databbn1, familyset = c(0,1,3:6), CMatrix7, 
                          selectioncrit = "logLik", indeptest = TRUE,
                          level = 5*10^-5, trunclevel = NA,
                          weights = NA, rotations = TRUE, se = FALSE,
                          presel = TRUE, method = "mle", cores = 1)


cond.ind.matrix3 = which(RVMcop7$pair.logLik[lower.tri(RVMcop7$pair.logLik)]==0)
if(length(cond.ind.matrix3) == 0){cond.ind.matrix3 = 0}
cond.ind.matrix3 = cond.ind.matrix3[1]


CMatrix8 <- c(3, 1, 4, 2,
              0, 1, 2, 4,
              0, 0, 2, 4,
              0, 0, 0, 4)
CMatrix8 <- matrix(CMatrix8, 4, 4)
RVMcop8 <- RVineCopSelect(databbn1, familyset = c(0,1,3:6), CMatrix8, 
                          selectioncrit = "logLik", indeptest = TRUE,
                          level = 5*10^-5, trunclevel = NA,
                          weights = NA, rotations = TRUE, se = FALSE,
                          presel = TRUE, method = "mle", cores = 1)


cond.ind.matrix4 = which(RVMcop8$pair.logLik[lower.tri(RVMcop8$pair.logLik)]==0)
if(length(cond.ind.matrix4) == 0){cond.ind.matrix4 = 0}
cond.ind.matrix4 = cond.ind.matrix4[1]








#############################################
# Lastly we check what happen if we apply
# the learning algorithms for Gaussian Data:

# Transforming the margins:
gdatabbn1 <- qnorm(data.matrix(databbn1),0,1)
logLik.Gaussian <- sum(dmvnorm(gdatabbn1, colMeans(gdatabbn1), cov(gdatabbn1), log=TRUE))


gdatabbn1 <- data.frame(qnorm(data.matrix(databbn1),0,1))
colnames(gdatabbn1) = c("X1","X2","X3","X4")




# Fitting multivariate normal using Max.Lik.
X <- mvnorm.mle(data.matrix(gdatabbn1))
fitted_gdata <- mvrnorm(n = n, mu = X$mu, Sigma = X$sigma)
fitted_gdata <- data.frame(fitted_gdata)
X.fitted_gdata <- mvnorm.mle(data.matrix(fitted_gdata))


# Comparing margins of both models, we can check the
# big differences, and the possible errors that we can
# make assuming that our distribution is multivariate gaussian.

ggpairs(fitted_gdata, lower = list(
  continuous = wrap("points", alpha = 0.3, size=0.4)))

ggpairs(gdatabbn1,lower = list(
  continuous = wrap("points", alpha = 0.3, size=0.4)))

# We also include a goodness of fit test to asses
# the gaussianity assumptions.
mvn(gdatabbn1,subset = NULL,mvnTest = "mardia",covariance = TRUE,
    tol = 1e-25,alpha = 0.5,scale = FALSE,desc = TRUE,
    transform = "none", R = 1000,univariateTest = "AD",
    univariatePlot = "qq",multivariatePlot = "qq",
    multivariateOutlierMethod = "none",bc = FALSE,
    bcType = "rounded",showOutliers = FALSE,showNewData = FALSE)


mvn(fitted_gdata,subset = NULL,mvnTest = "mardia",covariance = TRUE,
    tol = 1e-25,alpha = 0.5,scale = FALSE,desc = TRUE,
    transform = "none", R = 1000,univariateTest = "AD",
    multivariateOutlierMethod = "none",bc = FALSE,
    bcType = "rounded",showOutliers = FALSE,showNewData = FALSE)



# Let's try a bunch of tests:
par(mfrow=c(1,2))
# Constraint-based algorithms: Hiton-PC.
cpdag1 <- si.hiton.pc(gdatabbn1, undirected = FALSE)
par(mfrow=c(1,1))
graphviz.plot(cpdag(cpdag1))


par(mfrow=c(1,2))
shd(cpdag1, cpdag(dag), wlbl = FALSE, debug = TRUE)
graphviz.compare(cpdag(dag),cpdag1)
hamming(cpdag1,dag)


# We check if the conditional independencies holds:
# We first use the Correlation matrix:
# checking 2ind3|1.
invcor  <- round(cor2pcor(cor(databbn1[, c("U1", "U2", "U3")])),3)
invcor[2,3] # the value of rho_{X2,X3|X1} 

invcor  <- round(cor2pcor(cor(databbn1[, c("U1", "U2", "U3","U4")])),3)
invcor[1,4] # the value of rho_{X1,X4|X2,X3} 

invcor  <- round(cor2pcor(cor(gdatabbn1[, c("X1", "X2", "X3")])),3)
invcor[2,3] # the value of rho_{X2,X3|X1} Gaussian Case

invcor  <- round(cor2pcor(cor(gdatabbn1[, c("X1", "X2", "X3","X4")])),3)
invcor[1,4] # the value of rho_{X1,X4|X2,X3} Gaussian Case


# We have that 3\indep4|1
# we then check that in our copula based model 
# it happened the same.

colnames(gdatabbn1) = c("X1","X2","X3","X4")
# Some of the conditional independence tests that we can use:
# Correlation:
gcind23.cor = ci.test("X2", "X3", c("X1"), test = "cor", data = gdatabbn1)
gcind14.cor = ci.test("X1", "X4", c("X2","X3"), test = "cor", data = gdatabbn1)
gcind34.cor = ci.test("X3", "X4", c("X1"), test = "cor", data = gdatabbn1)


CMatrix10 <- c(2, 4, 1, 3,
               0, 3, 4, 1,
               0, 0, 1, 4,
               0, 0, 0, 4)
CMatrix10 <- matrix(CMatrix10, 4, 4)
RVMcop10 <- RVineCopSelect(databbn1, familyset = c(0,1,3:6), CMatrix10, 
                           selectioncrit = "logLik", indeptest = FALSE,
                           level = 0.05, trunclevel = NA,
                           weights = NA, rotations = TRUE, se = FALSE,
                           presel = TRUE, method = "mle", cores = 1)

RVMcop10 <- RVineCopSelect(databbn1, familyset = c(0,1,3:6), CMatrix10, 
                           selectioncrit = "logLik", indeptest = TRUE,
                           level = 0.05, trunclevel = NA,
                           weights = NA, rotations = TRUE, se = FALSE,
                           presel = TRUE, method = "mle", cores = 1)



# Score-based algorithms: Hill Climbing:
set.seed(123)
cpdag2 = hc(gdatabbn1, score = "bic-g",  restart = 30)
score.result <- score(cpdag2, data = gdatabbn1, type = "bic-g")
truescore <- score(dag, data = gdatabbn1, type = "bic-g")
par(mfrow=c(1,1))
graphviz.plot(cpdag(cpdag2))
hamming(cpdag2,dag)

par(mfrow=c(1,2))
shd(cpdag(cpdag2), cpdag(dag), wlbl = FALSE, debug = TRUE)
graphviz.compare(cpdag(dag),cpdag(cpdag2))




# Graphical Lasso:
# we use both, CO1 and  Lasso families:
# The only problem is that we obtain a undirected graph.
gdatabbn1 <- data.matrix(gdatabbn1)
p=4
GRest <- selectFast(gdatabbn1, family=c("C01","LA"))

adjmatrix<-c(0, 1, 1, 0,
             1, 0, 1, 1,
             1, 1, 0, 1,
             0, 1, 1, 0)
adjmatrix <- matrix(adjmatrix, 4, 4)

# Comparing the adjacency matrices we get:
C01.result = all(GRest$C01$G == adjmatrix)
LA.result = all(GRest$LA$G == adjmatrix)

par(mfrow=c(1,1))
gV <- network(GRest$C01$G)
plot(gV,jitter=TRUE, usearrows = FALSE, label=1:p,displaylabels=TRUE)


g <- network(GRest$LA$G)
plot(g, jitter=TRUE, usearrows = FALSE, label=1:p,displaylabels=TRUE)

penalty(4,2500, dmax=min(3,n-3,p-1), K=2.5)



