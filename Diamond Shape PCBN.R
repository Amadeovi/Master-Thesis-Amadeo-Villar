# Diamond Shape PCBN: used in Chapter 4 PCBNs 
# and in Chapter 7 motivation our proposed method.

# Experiment 20/10:
# We simulate from a the diamond shape PCBN
# We study the different results that we get:


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
library("reshape2")

###############################
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
# We start with the simulations:
start.sim.time <- Sys.time()

# We first simulate from the first 3 variables:
Cmatrix <- c(3, 2, 1,
             0, 2, 1,
             0, 0, 1)
Cmatrix <- matrix(Cmatrix, 3, 3)

family <- c(0, 0, family13,  
            0, 0, family12,  
            0, 0, 0)
family <- matrix(family, 3, 3)

par <- c(0, 0, par13,
         0, 0, par12, 
         0, 0, 0)
par <- matrix(par, 3, 3)
par2 <- matrix(5, 3, 3)

RVM <- RVineMatrix(Matrix = Cmatrix, family = family,
                   par = par, par2 = par2,
                   names = c("U1", "U2", "U3"))

# simulate from the Bayesian Network:
set.seed(123)
N <- 2500
databbn1 <- RVineSim(N, RVM)


# Lets compute the last column:
databbn1 <- cbind(databbn1,rep(NA,N))
colnames(databbn1) <- c("U1","U2","U3","U4")

NDIM <- 1
NCOMP <- 1

for(j in 1:N){
  
  integrand <- function(y1) {
    ff <- BiCopPDF(y1, databbn1[j,3], family=family13,
                   par=par13, par2=5, obj = NULL, check.pars = TRUE)*
      BiCopHfunc1(y1, databbn1[j,2], family=family12, 
                  par=par12, par2=5, obj = NULL, check.pars = TRUE)
    return(ff)
  } 
  
  a <- cuhre(f = integrand,lowerLimit = rep(0, NDIM),upperLimit = rep(1, NDIM),
             relTol = 1e-2, absTol= 1e-9,flags = list(verbose = 0, final = 0))$integral
  
  if(a>1){
    a=1
  }
  
  b <- BiCopHinv1(a,runif(1),family=family24.3, par=par24.3,
                  par2 = 5, obj = NULL, check.pars = TRUE)
  
  
  databbn1[j,4] <- BiCopHinv1(databbn1[j,3],b,family=family34, par=par34,
                              par2 = 5, obj = NULL, check.pars = TRUE)
}

final.sim.time <- Sys.time()
Sim.time <- difftime(final.sim.time,start.sim.time,units="secs")



###############################
# Let's check whether we have simulated properly:
cdata <- as.copuladata(databbn1)
pairs(cdata)

data <- data.frame(databbn1)

# For plotting we define:
size.dot = 0.6
family.list <- list("0" = "Ind", "1" = "Gaussian","2"="t", "3" ="Clayton",
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




###############################
# Second quantitative analysis: 
# We compute the True logLikelihood with the True Order:
start.TruelogLik <- Sys.time()

set.seed(123)
C12logLik <- sum(log(BiCopPDF(databbn1[,1],databbn1[,2],family=family12,
                              par=par12, par2=5)))
C13logLik <- sum(log(BiCopPDF(databbn1[,1],databbn1[,3],family=family13,
                              par=par13, par2=5)))
C34logLik <- sum(log(BiCopPDF(databbn1[,3],databbn1[,4],family=family34,
                              par=par34, par2=5)))

u2.3.true <- c()

for (j in 1:N){
  
  integrand <- function(y1) {
    ff <- BiCopPDF(y1,databbn1[j,3],family=family13,par=par13,par2=5)*
      BiCopHfunc1(y1, databbn1[j,2],family=family12,par=par12,par2=5,
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


final.TruelogLik <- Sys.time()
Comp.TruelogLik <- difftime(final.TruelogLik, start.TruelogLik, units="secs")


# Moreover, we also plot the conditional copula, to make sure
# we are simulating from the correct distribution:
pseudo.obs <- data.frame(u2.3.true, u4.3.true)

ggplot(pseudo.obs, aes(x = u2.3.true, y = u4.3.true))+geom_point(size=size.dot)+
  labs(title = TeX(paste("$C_{2,4|3}$",":",  family.list[as.character(family24.3)],",",
                         '$\\tau$',"=",round(BiCopPar2Tau(family24.3,par24.3),3))),
       x = TeX("$U_{2|3}$"),y = TeX("$U_{4|3}$"))+
  theme(plot.title = element_text(hjust = 0.5))




###############################
# We also compute the Estimated logLikelihood with the true structure 
# and the true order.

start.estimTruelogLik <- Sys.time()
C12.estim <- BiCopSelect(databbn1[,1], databbn1[,2], family = c(0,1,2,3:6), 
                          selectioncrit = "logLik", indeptest = TRUE, 
                          level = 0.05, weights = NA, rotations = TRUE,
                          se = FALSE,presel = TRUE,method = "mle")
C13.estim <- BiCopSelect(databbn1[,1], databbn1[,3], family = c(0,1,2,3:6), 
                          selectioncrit = "logLik", indeptest = TRUE, 
                          level = 0.05, weights = NA, rotations = TRUE,
                          se = FALSE,presel = TRUE,method = "mle")
C34.estim <- BiCopSelect(databbn1[,3], databbn1[,4], family = c(0,1,2,3:6), 
                          selectioncrit = "logLik", indeptest = TRUE, 
                          level = 0.05, weights = NA, rotations = TRUE,
                          se = FALSE,presel = TRUE,method = "mle")

u2.3.estim <- c()

for (j in 1:N){
  
  integrand <- function(y1) {
    ff <- BiCopPDF(y1,databbn1[j,3],C13.estim)*
      BiCopHfunc1(y1, databbn1[j,2], family=C12.estim$family, 
                  par=C12.estim$par, par2 = C12.estim$par2,
                  obj = NULL, check.pars = TRUE)
    return(ff)
  } 
  
  u2.3.estim[j] <- cuhre(f = integrand,lowerLimit = rep(0, NDIM),
                   upperLimit = rep(1, NDIM), relTol = 1e-2, 
                   absTol= 1e-9,flags = list(verbose = 0, final = 0))$integral
  
  if(u2.3.estim[j]>1){
    u2.3.estim[j]=1
  }
} 

u4.3.estim <- BiCopHfunc1(databbn1[,3], databbn1[,4], family=C34.estim$family, 
                    par=C34.estim$par, par2 = C34.estim$par2,
                    obj = NULL, check.pars = TRUE)

C24.3.estim <- BiCopSelect(u2.3.estim, u4.3.estim, family = c(0,1,2,3:6), 
                            selectioncrit = "logLik", indeptest = TRUE, 
                            level = 0.05, weights = NA, rotations = TRUE,
                            se = FALSE,presel = TRUE,method = "mle")

Estim.truelogLik <- (C13.estim$logLik + C12.estim$logLik +
                       C24.3.estim$logLik + C34.estim$logLik)
Estim.truelogLik24.3 <- C24.3.estim$logLik


final.estimTruelogLik <- Sys.time()
Comp.estimTruelogLik <- difftime(final.estimTruelogLik, 
                                 start.estimTruelogLik, units="secs")

C12.estim
C13.estim
C34.estim
C24.3.estim



###############################
# Now let's estimate the logLik of network with the true structure
# but the wrong order of the parents:

start.estimtruelogLik.wo <- Sys.time()
C12.estim <- BiCopSelect(databbn1[,1], databbn1[,2], family = c(0,1,2,3:6), 
                         selectioncrit = "logLik", indeptest = TRUE, 
                         level = 0.05, weights = NA, rotations = TRUE,
                         se = FALSE,presel = TRUE,method = "mle")
C13.estim <- BiCopSelect(databbn1[,1], databbn1[,3], family = c(0,1,2,3:6), 
                         selectioncrit = "logLik", indeptest = TRUE, 
                         level = 0.05, weights = NA, rotations = TRUE,
                         se = FALSE,presel = TRUE,method = "mle")
C24.estim <- BiCopSelect(databbn1[,2], databbn1[,4], family = c(0,1,2,3:6), 
                         selectioncrit = "logLik", indeptest = TRUE, 
                         level = 0.05, weights = NA, rotations = TRUE,
                         se = FALSE,presel = TRUE,method = "mle")



u3.2.estim <- c()

for (j in 1:N){
  
  integrand <- function(y1) {
    ff <- BiCopPDF(y1,databbn1[j,2],C12.estim)*
           BiCopHfunc1(y1, databbn1[j,3], family=C13.estim$family, 
                       par=C13.estim$par, par2 = C13.estim$par2,
                       obj = NULL, check.pars = TRUE)
    return(ff)
  } 
  
  u3.2.estim[j] <- cuhre(f = integrand,lowerLimit = rep(0, NDIM),
                         upperLimit = rep(1, NDIM), relTol = 1e-2, 
                         absTol= 1e-9,flags = list(verbose = 0, final = 0))$integral

  
  if(u3.2.estim[j]>1){
    u3.2.estim[j]=1
  }
} 

u4.2.estim <- BiCopHfunc1(databbn1[,2], databbn1[,4], family=C24.estim$family, 
                          par=C24.estim$par, par2 = C24.estim$par2,
                          obj = NULL, check.pars = TRUE)


C34.2.estim <- BiCopSelect(u3.2.estim, u4.2.estim, family = c(0,1,2,3:6), 
                           selectioncrit = "logLik", indeptest = TRUE, 
                           level = 0.05, weights = NA, rotations = TRUE,
                           se = FALSE,presel = TRUE,method = "mle")


Estim.truelogLik.wo <- (C12.estim$logLik + C13.estim$logLik +
                          C24.estim$logLik + C34.2.estim$logLik)
Estim.truelogLik34.2 <- C34.2.estim$logLik


final.estimtruelogLik.wo <- Sys.time()
Comp.estimTruelogLik <- difftime(final.estimtruelogLik.wo, 
                                 start.estimtruelogLik.wo, units="secs")





###############################
# The following step is to use the CD-Vine decomposition score:

start.CDVine <- Sys.time()

C12.cd <- BiCopSelect(databbn1[,1], databbn1[,2], family = c(0,2,1,3:6), 
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


CMatrix2 <- c(1, 2, 3,
              0, 2, 3,
              0, 0, 3)
CMatrix2 <- matrix(CMatrix2, 3, 3)
RVMcop.cd <- RVineCopSelect(databbn1[,2:4], familyset = c(0,1,2,3:6), CMatrix2, 
                            selectioncrit = "logLik", indeptest = TRUE,
                            level = 0.05, trunclevel = NA,
                            weights = NA, rotations = TRUE, se = FALSE,
                            presel = TRUE, method = "mle", cores = 1)


logLik.cd <-  (C12.cd$logLik + C13.cd$logLik +
                 RVMcop.cd$logLik - C23.cd$logLik)
                 
                
final.CDVine <- Sys.time()
Comp.CDVine <- difftime(final.CDVine, start.CDVine, units="secs")




##############################
# We are going to study another network, given by:
Y <- model2network(paste0("[X2][X4|X2][X3|X2:X4][X1|X3]"))
# The logLik of this network is given by:
# The CD Vine score for this network is given by:
C13.Y <- BiCopSelect(databbn1[,1], databbn1[,3], family = c(0,2,1,3:6), 
                     selectioncrit = "logLik", indeptest = TRUE, 
                     level = 0.05, weights = NA, rotations = TRUE,
                     se = FALSE,presel = TRUE,method = "mle")

CMatrixY <- c(1, 3, 2,
              0, 3, 2,
              0, 0, 2)
CMatrixY <- matrix(CMatrixY, 3, 3)
RVMcop.CY <- RVineCopSelect(databbn1[,2:4], familyset = c(0,1,2,3:6), CMatrixY, 
                            selectioncrit = "logLik", indeptest = TRUE,
                            level = 0.05, trunclevel = NA,
                            weights = NA, rotations = TRUE, se = FALSE,
                            presel = TRUE, method = "mle", cores = 1)


logLik.cd.Y <- C13.Y$logLik +  RVMcop.CY$logLik 
# We obtain that using the CD-Vine score, this network 
# has higher score than the true one. Therefore the score in misspecified.
shd(dag,Y)
hamming(dag,Y)
# We can see how the score prioritize the false DAG over the true DAG,
# and the SHD between both is 4, which is pretty high for the
# number of edges that we are considering.






##############################
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
RVMcop.estcop.to <- RVineCopSelect(databbn1[,2:4], familyset = c(0,2,1,3:6), 
                                    CMatrix.estcop.to, selectioncrit = "logLik",
                                    indeptest = TRUE, level = 0.05, trunclevel = NA,
                                    weights = NA, rotations = TRUE, se = FALSE,
                                    presel = TRUE, method = "mle", cores = 1)


# checking the value of the empirical tau between 3 and 4:
tau.34.to <- RVMcop.estcop.to$emptau[3,2]
logLik.estcop.to <- (C13.estcop.to$logLik + C12.estcop.to$logLik +
                       RVMcop.estcop.to$logLik - C23.estcop.to$logLik)
logLik.estcop.to.24.3 <- RVMcop.estcop.to$pair.logLik[2,1]


final.estcop.to <- Sys.time()
Comp.estcop.to <- difftime(final.estcop.to, start.estcop.to, units="secs")




##############################
# It is also interesting comparing the true margins vs data generated by the
# estimated copula and how this affect the upper trees.
# to the simulated data using the estimated copula 2-3.
set.seed(1)
margins <- BiCopSim(N, RVMcop.estcop.to$family[3,1], RVMcop.estcop.to$par[3,1],
                    par2 = RVMcop.estcop.to$par2[3,1],
                    obj = NULL, check.pars = TRUE)

margins <- data.frame(margins)

# True margins:
ggplot(data, aes(x = U2, y = U3))+geom_point(size=size.dot)+
  labs(title = TeX(paste("True margins","$U_{2}$","and","$U_{3}$",
                         ",","$\\tau$","=",round(RVMcop.estcop.to$emptau[3,1],3))),
       x = TeX("$U_{2}$"),y = TeX("$U_{3}$"))+
  theme(plot.title = element_text(hjust = 0.5))

# Estimated margins:
ggplot(margins, aes(x = margins[,1], y = margins[,2]))+geom_point(size=size.dot)+
  labs(title = TeX(paste("Simulated data from","$\\hat{C}_{2,3}$",
                         ":",family.list[as.character(RVMcop.estcop.to$family[3,1])],",",
                         "$\\tau$","=",round(RVMcop.estcop.to$emptau[3,1],3))), 
       x = TeX("$U_{2}$") , y = TeX("$U_{3}$"))+theme(plot.title = element_text(hjust = 0.5))

# We also compute a goodness of fit test for these data:
BiCopGofTest(databbn1[,2],databbn1[,3],family=RVMcop.estcop.to$family[3,1],
             par=RVMcop.estcop.to$par[3,1],par2=RVMcop.estcop.to$par2[3,1],
             method = "kendall", B = 30) # it rejects the null hypoth.



# This will induce an error in the higher trees, let's check
# what happens with these errors:
set.seed(1)
high.margins <- BiCopSim(N, RVMcop.estcop.to$family[2,1], RVMcop.estcop.to$par[2,1],
                         par2 = RVMcop.estcop.to$par2[2,1],
                         obj = NULL, check.pars = TRUE)

high.margins <- data.frame(high.margins)

ggplot(pseudo.obs, aes(x = u2.3.true, y = u4.3.true))+geom_point(size=size.dot)+
  labs(title = TeX(paste("$C_{2,4|3}$",":",  family.list[as.character(family24.3)],",",
                         '$\\tau$',"=",round(BiCopPar2Tau(family24.3,par24.3),3))),
       x = TeX("$U_{2|3}$"),y = TeX("$U_{4|3}$"))+
  theme(plot.title = element_text(hjust = 0.5))


ggplot(high.margins, aes(x = high.margins[,1], y = high.margins[,2]))+
  geom_point(size=size.dot)+labs(title = TeX(paste("Simulated data from","$\\hat{C}_{2,4|3}$",
                                                   ":",family.list[as.character(RVMcop.estcop.to$family[2,1])],",",
                                                   "$\\tau$","=",round(RVMcop.estcop.to$emptau[2,1],3))),
                                 x = TeX("$U_{2|3}$"),y = TeX("$U_{4|3}$"))+
  theme(plot.title = element_text(hjust = 0.5))


# We also compute a goodness of fit test for these data:
BiCopGofTest(u2.3.true, u4.3.true, family=RVMcop.estcop.to$family[2,1],
             par=RVMcop.estcop.to$par[2,1], par2=RVMcop.estcop.to$par2[2,1],
             method = "kendall", B = 30)
# In this case it does not reject the null hypothesis.



# Let's also check whether the \hat{u}_{2|3}=\hat{C}_{2|3}(u_{2},u_{3})
# gives a uniform distribution:
hatu2.3 <- BiCopHfunc2(data$U2, data$U3, family = RVMcop.estcop.to$family[3,1], 
                       par = RVMcop.estcop.to$par[3,1], par2 = RVMcop.estcop.to$par2[3,1],
                       obj = NULL, check.pars = TRUE)


compare2.3 <- data.frame(u2.3.true, hatu2.3)
compare2.3dens <- melt(compare2.3[c("u2.3.true","hatu2.3")],id=NULL)

# Comparison between both values:
# Checking uniformity:
ggplot(compare2.3dens, aes(x=value, color=variable, fill=variable)) +
  geom_density(alpha=0.2)+ ylim(0,1.2) + labs(title = TeX(
    paste("Density functions:","$u_{2|3}$"," and ",
          "$\\hat{u}_{2|3}=\\hat{C}_{2|3}(u_{2},u_{3})")))


# Checking difference in values:
ggplot(compare2.3, aes(x = u2.3.true, y=hatu2.3)) + geom_point(alpha=0.2) +
  labs(title = TeX(paste("Comparison:","$u_{2|3}$"," and ",
                         "$\\hat{u}_{2|3}=\\hat{C}_{2|3}(u_{2},u_{3})")),
       x=TeX("$u_{2|3}$"), y=TeX("$\\hat{u}_{2|3}$"))


ks.test(hatu2.3, "punif")





##############################
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
RVMcop.estcop.wo <- RVineCopSelect(databbn1[,2:4], familyset = c(0,1,3:6), CMatrix.estcop.wo, 
                                   selectioncrit = "logLik", indeptest = TRUE,
                                   level = 0.05, trunclevel = NA,
                                   weights = NA, rotations = TRUE, se = FALSE,
                                   presel = TRUE, method = "mle", cores = 1)

# checking the value of the empirical tau between 2 and 4
# to check how good is the heuristic
tau.34.wo  <- RVMcop.estcop.wo$emptau[3,1]
logLik.estcop.wo <- (C13.estcop.wo$logLik + C12.estcop.wo$logLik +
                       RVMcop.estcop.wo$logLik - C23.estcop.wo$logLik)


final.estcop.wo<- Sys.time()
Comp.estcop.wo <- difftime(final.estcop.wo, start.estcop.wo, units="secs")






##############################
# Let's apply the different Gaussian Learning algorithms to 
# our PCBN data with normal margins.
gdatabbn1 <- qnorm(data.matrix(databbn1),0,1)
logLik.Gaussian <- sum(dmvnorm(gdatabbn1, colMeans(gdatabbn1), cov(gdatabbn1), log=TRUE))
colnames(gdatabbn1) = c("X1","X2","X3","X4")



# Fitting multivariate normal using Max.Lik.
X <- mvnorm.mle(data.matrix(gdatabbn1))
fitted_gdata <- mvrnorm(n = N, mu = X$mu, Sigma = X$sigma)
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



# Graphical Lasso.
set.seed(1)
gdatabbn1 <- data.matrix(gdatabbn1)
p=4
GRest <- selectFast(gdatabbn1, family=c("LA"),K=2.5)

adjmatrix<-c(0, 1, 1, 0,
             1, 0, 1, 1,
             1, 1, 0, 1,
             0, 1, 1, 0)
adjmatrix <- matrix(adjmatrix, 4, 4)

# Comparing the adjacency matrices we get:
GRest$LA$G
LA.result = all(GRest$LA$G == adjmatrix)


g <- network(GRest$LA$G)
plot(g, jitter=TRUE, usearrows = FALSE, label=1:p,displaylabels=TRUE)

penalty(4,2500, dmax=min(3,N-3,p-1), K=2.5)



gdatabbn1 <- data.frame(gdatabbn1)

# Constraint-based algorithms: Hiton-PC algorithm:
pc.start.time <- Sys.time()

cpdag.pc <- si.hiton.pc(gdatabbn1, undirected = FALSE)
pc.shd <- shd(cpdag(dag),cpdag(cpdag.pc), wlbl = FALSE, debug = FALSE)
pc.ham <- hamming(dag,cpdag.pc)
pc.arcs <- as.numeric(gsub(" arcs:", "", head(grep("arcs", capture.output(cpdag.pc), value = TRUE), 1)))
graphviz.plot(cpdag.pc)

pc.final.time <- Sys.time()
pc.sim.time <- difftime(pc.final.time, pc.start.time, units="secs")
graphviz.compare(cpdag(dag),cpdag.pc)


# For the output graph we get X_{3}\indep X_{4}|X_{1}, but
# we don't recover the true conditional independencies.
# we then check that in our copula based model 
# it happened the same.

colnames(gdatabbn1) = c("X1","X2","X3","X4")
# Some of the conditional independence tests that we can use:
# Correlation:
ci.test("X2", "X3", c("X1"), test = "cor", data = gdatabbn1)
ci.test("X1", "X4", c("X2","X3"), test = "cor", data = gdatabbn1)
ci.test("X3", "X4", c("X1"), test = "cor", data = gdatabbn1)


CMatrix10 <- c(2, 4, 1, 3,
               0, 3, 4, 1,
               0, 0, 1, 4,
               0, 0, 0, 4)
CMatrix10 <- matrix(CMatrix10, 4, 4)
# False independence test.
RVMcop10.indtest <- RVineCopSelect(databbn1, familyset = c(0,1,3:6), CMatrix10, 
                                   selectioncrit = "logLik", indeptest = FALSE,
                                   level = 0.05, trunclevel = NA,
                                   weights = NA, rotations = TRUE, se = FALSE,
                                   presel = TRUE, method = "mle", cores = 1) 
# True independence test.
RVMcop10 <- RVineCopSelect(databbn1, familyset = c(0,1,3:6), CMatrix10, 
                           selectioncrit = "logLik", indeptest = TRUE,
                           level = 0.05, trunclevel = NA,
                           weights = NA, rotations = TRUE, se = FALSE,
                           presel = TRUE, method = "mle", cores = 1)



# Scored-based algorithms: HC algorithm:

hc.start.time <- Sys.time()

cpdag.hc <- hc(gdatabbn1, score = "bic-g",  restart = 15)
hc.shd <- shd(cpdag(dag),cpdag(cpdag.hc), wlbl = FALSE, debug = FALSE)
hc.ham <- hamming(dag,cpdag.hc)
hc.arcs <- as.numeric(gsub(" arcs:", "", head(grep("arcs", capture.output(cpdag.hc), value = TRUE), 1)))
graphviz.plot(cpdag.hc)

hc.final.time <- Sys.time()
hc.sim.time <- difftime(hc.final.time, hc.start.time, units="secs")
graphviz.compare(cpdag(dag),cpdag.hc)



###############################
# Rehearsal of what out method would be:
C14.try <- BiCopSelect(databbn1[,1], databbn1[,4], family = c(0,1,2,3:6), 
                         selectioncrit = "logLik", indeptest = TRUE, 
                         level = 0.05, weights = NA, rotations = TRUE,
                         se = FALSE,presel = TRUE,method = "mle")
C24.try <- BiCopSelect(databbn1[,2], databbn1[,4], family = c(0,1,2,3:6), 
                         selectioncrit = "logLik", indeptest = TRUE, 
                         level = 0.05, weights = NA, rotations = TRUE,
                         se = FALSE,presel = TRUE,method = "mle")
C34.try<- BiCopSelect(databbn1[,3], databbn1[,4], family = c(0,1,2,3:6), 
                         selectioncrit = "logLik", indeptest = TRUE, 
                         level = 0.05, weights = NA, rotations = TRUE,
                         se = FALSE,presel = TRUE,method = "mle")



# Let's see how we would order the parents:
CMatrix.try <- c(4, 1, 3, 2,
                 0, 1, 2, 3,
                 0, 0, 2, 3,
                 0, 0, 0, 3)
CMatrix.try <- matrix(CMatrix.try, 4, 4)
RVMcop.try <- RVineCopSelect(databbn1, familyset = c(0,1,2,3:6), 
                                   CMatrix.try, selectioncrit = "logLik",
                                   indeptest = TRUE, level = 0.05, trunclevel = NA,
                                   weights = NA, rotations = TRUE, se = FALSE,
                                   presel = TRUE, method = "mle", cores = 1)



CMatrix.try.2 <- c(4, 3, 1, 2,
                   0, 3, 2, 1,
                   0, 0, 2, 1,
                   0, 0, 0, 1)
CMatrix.try.2 <- matrix(CMatrix.try.2, 4, 4)
RVMcop.try.2 <- RVineCopSelect(databbn1, familyset = c(0,1,2,3:6), 
                                   CMatrix.try.2, selectioncrit = "logLik",
                                   indeptest = TRUE, level = 0.05, trunclevel = 1,
                                   weights = NA, rotations = TRUE, se = FALSE,
                                   presel = TRUE, method = "mle", cores = 1)

# Result that C_{1,4|2} treated as independent, so order is: 2,3,4

# BIC computation:
# Initial Network:
(sum(RVMcop.try$pair.logLik[,1])+C12.estim$logLik+C13.estim$logLik)-5*log(N)/2
# Subtracting U1->U4
(sum(RVMcop.try$pair.logLik[,1])+C12.estim$logLik+C13.estim$logLik)-4*log(N)/2
# Substracting U1->U2
(sum(RVMcop.try$pair.logLik[,1])+C13.estim$logLik)-4*log(N)/2
# Substracting U1->U3
(sum(RVMcop.try$pair.logLik[,1])+C12.estim$logLik)-4*log(N)/2


# Adding the edge 2->3 for instance we have:
CMatrix.try.3 <- c(3, 2, 1,
                   0, 2, 1,
                   0, 0, 1)
CMatrix.try.3 <- matrix(CMatrix.try.3, 3, 3)
RVMcop.try.3 <- RVineCopSelect(databbn1[,1:3], familyset = c(0,1,2,3:6), 
                               CMatrix.try.3, selectioncrit = "logLik",
                               indeptest = TRUE, level = 0.05, trunclevel = NA,
                               weights = NA, rotations = TRUE, se = FALSE,
                               presel = TRUE, method = "mle", cores = 1)


CMatrix.try.4 <- c(3, 1, 2,
                   0, 1, 2,
                   0, 0, 2)
CMatrix.try.4 <- matrix(CMatrix.try.4, 3, 3)
RVMcop.try.4 <- RVineCopSelect(databbn1[,1:3], familyset = c(0,1,2,3:6), 
                               CMatrix.try.4, selectioncrit = "logLik",
                               indeptest = TRUE, level = 0.05, trunclevel = NA,
                               weights = NA, rotations = TRUE, se = FALSE,
                               presel = TRUE, method = "mle", cores = 1)


# Addign edge U2->U3:
(sum(RVMcop.try$pair.logLik[,1])+sum(RVMcop.try.3$pair.logLik[,1])+C12.estim$logLik)-6*log(N)/2
      


# What if we add the edge: 3->2
CMatrix.try.5 <- c(2, 1, 3,
                   0, 1, 3,
                   0, 0, 3)
CMatrix.try.5 <- matrix(CMatrix.try.5, 3, 3)
RVMcop.try.5 <- RVineCopSelect(databbn1[,1:3], familyset = c(0,1,2,3:6), 
                               CMatrix.try.5, selectioncrit = "logLik",
                               indeptest = TRUE, level = 0.05, trunclevel = NA,
                               weights = NA, rotations = TRUE, se = FALSE,
                               presel = TRUE, method = "mle", cores = 1)
      

CMatrix.try.6 <- c(2, 3, 1,
                   0, 3, 1,
                   0, 0, 1)
CMatrix.try.6 <- matrix(CMatrix.try.6, 3, 3)
RVMcop.try.6 <- RVineCopSelect(databbn1[,1:3], familyset = c(0,1,2,3:6), 
                               CMatrix.try.6, selectioncrit = "logLik",
                               indeptest = TRUE, level = 0.05, trunclevel = NA,
                               weights = NA, rotations = TRUE, se = FALSE,
                               presel = TRUE, method = "mle", cores = 1)
      
      
# Adding the edge U3->U2
(sum(RVMcop.try$pair.logLik[,1])+sum(RVMcop.try.6$pair.logLik[,1])+C13.estim$logLik)-6*log(N)/2



# Let's see what happens if we reverse 3->4
C13.try <- BiCopSelect(databbn1[,1], databbn1[,3], family = c(0,1,2,3:6), 
                       selectioncrit = "logLik", indeptest = TRUE, 
                       level = 0.05, weights = NA, rotations = TRUE,
                       se = FALSE,presel = TRUE,method = "mle")


CMatrix.try.7 <- c(2, 3, 1,
                   0, 3, 1,
                   0, 0, 1)
CMatrix.try.7 <- matrix(CMatrix.try.7, 3, 3)
RVMcop.try.7 <- RVineCopSelect(databbn1[,c(1,3,4)], familyset = c(0,1,2,3:6), 
                               CMatrix.try.7, selectioncrit = "logLik",
                               indeptest = TRUE, level = 0.05, trunclevel = NA,
                               weights = NA, rotations = TRUE, se = FALSE,
                               presel = TRUE, method = "mle", cores = 1)


C24.try <- BiCopSelect(databbn1[,2], databbn1[,4], family = c(0,1,2,3:6), 
                       selectioncrit = "logLik", indeptest = TRUE, 
                       level = 0.05, weights = NA, rotations = TRUE,
                       se = FALSE,presel = TRUE,method = "mle")

CMatrix.try.8 <- c(3, 1, 2,
                   0, 1, 2,
                   0, 0, 2)
CMatrix.try.8 <- matrix(CMatrix.try.8, 3, 3)
RVMcop.try.8 <- RVineCopSelect(databbn1[,c(1,2,4)], familyset = c(0,1,2,3:6), 
                               CMatrix.try.8, selectioncrit = "logLik",
                               indeptest = TRUE, level = 0.05, trunclevel = NA,
                               weights = NA, rotations = TRUE, se = FALSE,
                               presel = TRUE, method = "mle", cores = 1)

# Reversing U4->U3
(sum(RVMcop.try.8$pair.logLik[,1])+sum(RVMcop.try.7$pair.logLik[,1])+C12.estim$logLik)-5*log(N)/2


# Let's try to reverse U4->U2
C24.try <- BiCopSelect(databbn1[,2], databbn1[,4],family = c(0,1,2,3:6), 
                       selectioncrit = "logLik", indeptest = TRUE, 
                       level = 0.05, weights = NA, rotations = TRUE,
                       se = FALSE,presel = TRUE,method = "mle")


CMatrix.try.9 <- c(1, 2, 3,
                   0, 2, 3,
                   0, 0, 3)
CMatrix.try.9 <- matrix(CMatrix.try.9, 3, 3)
RVMcop.try.9 <- RVineCopSelect(databbn1[,c(1,2,4)], familyset = c(0,1,2,3:6), 
                               CMatrix.try.9, selectioncrit = "logLik",
                               indeptest = TRUE, level = 0.05, trunclevel = NA,
                               weights = NA, rotations = TRUE, se = FALSE,
                               presel = TRUE, method = "mle", cores = 1)

CMatrix.try.10 <- c(3, 2, 1,
                    0, 2, 1,
                    0, 0, 1)
CMatrix.try.10 <- matrix(CMatrix.try.10, 3, 3)
RVMcop.try.10 <- RVineCopSelect(databbn1[,c(1,3,4)], familyset = c(0,1,2,3:6), 
                               CMatrix.try.10, selectioncrit = "logLik",
                               indeptest = TRUE, level = 0.05, trunclevel = NA,
                               weights = NA, rotations = TRUE, se = FALSE,
                               presel = TRUE, method = "mle", cores = 1)


(RVMcop.try.10$logLik + RVMcop.try.9$pair.logLik[3,2] + RVMcop.try.9$pair.logLik[2,1])-5*log(N)/2
  
  
  

