# Some experiments using PCBNs, only using the Gaussian copula
# therefore recovering GBNs.

library("bnlearn")    # Gaussian and discrete BN
library("Rgraphviz")  # Plotting DAGs
library("gRain")
library("cubature")   # Numerical Integration
library("rlang")
library("rbmn")       # Réseaux bayésiens multinormaux
library("lattice")
library("igraph")
library("penalized")  # Lasso and Ridge regression
library("corpcor")    # inverse correlation matrix
library("compare")    # comparing TP, FP and FN
library("GGMselect")  # Gaussian Graphs Models Selection
library("VineCopula") # deal with copulas
library("Rfast")      # MLE of the gaussian distribution
library("MASS")

##############
# We construct BBN1-> not chordal graph:
X <- paste0("[X1][X2|X1][X3|X1][X4|X2:X3]")
dag<-model2network(X)
dag
graphviz.plot(dag)


# After showing the factorisation of this BBNs using copulas
# we proceed to simulate from it.
# In this case all our copulas are going to be gaussian
# such that the multivariate distribution is also gaussian.
C12 <- BiCop(family=1, par=0.5, par2 = 0, tau = NULL, check.pars = TRUE)
C13 <- BiCop(family=1, par=0.75, par2 = 0, tau = NULL, check.pars = TRUE)
C34 <- BiCop(family=1, par=-0.5, par2 = 0, tau = NULL, check.pars = TRUE)
C24.3 <- BiCop(family=1, par=-0.75, par2 = 0, tau = NULL, check.pars = TRUE)


# We are going to check 2 different ways of simulating from 
# this distribution:

# The first way is trough integration as we did
# when we simulated from BBN1:

NDIM <- 1
NCOMP <- 1

# We have to check whether there is a faster method
# to perform this integration or if it has an closed form solution
# Function for sampling:
samplebbn1 = function(n){
  
  databbn1 <- matrix(nrow=n,ncol=4)
  colnames(databbn1)=c("U1","U2","U3","U4")
  
  for (i in 1:n){
    
    databbn1[i,1] <- runif(1)
    
    databbn1[i,2] <- BiCopHinv1(databbn1[i,1],runif(1),family=C12$family, par=C12$par,
                                par2 = C12$par2, obj = NULL, check.pars = TRUE)
    
    databbn1[i,3] <- BiCopHinv1(databbn1[i,1],runif(1),family=C13$family, par=C13$par,
                                par2 = C13$par2, obj = NULL, check.pars = TRUE)
    
    
    integrand <- function(y1) {
      ff <- BiCopPDF(y1,databbn1[i,3],C13)*
        BiCopHfunc1(y1, databbn1[i,2], family=C12$family, 
                    par=C12$par, par2 = C12$par2,
                    obj = NULL, check.pars = TRUE)
      return(ff)
    } 
    
    a <- cuhre(f = integrand,lowerLimit = rep(0, NDIM),upperLimit = c(1),
               relTol = 1e-2, absTol= 1e-9,flags = list(verbose = 2, final = 0))$integral
    
    b <- BiCopHinv1(a,runif(1),family=C24.3$family, par=C24.3$par,
                    par2 = C24.3$par2, obj = NULL, check.pars = TRUE)
    
    databbn1[i,4] <- BiCopHinv1(databbn1[i,3],b,family=C34$family, par=C34$par,
                                par2 = C34$par2, obj = NULL, check.pars = TRUE)
    
  }
  return(databbn1)
}


# Performing the sampling with the number of samples that we want:
n = 300
databbn1 <- samplebbn1(n)

# We define the margins:
f1 <- list(mu = 3, sd = 1)
f2 <- list(mu = 5, sd = 2)
f3 <- list(mu = -1, sd = 1.5)
f4 <- list(mu = 0, sd = 1)

# We transform the margins:
gausscop.data <- matrix(nrow = n, ncol = 4)
colnames(gausscop.data) <- c("X1","X2","X3","X4")
for(i in 1:4){
  gausscop.data[,i] = qnorm(databbn1[,i],
                            eval(parse(text = paste("f",as.character(i),sep="")))$mu,
                            eval(parse(text = paste("f",as.character(i),sep="")))$sd)
}

# Checking the maximum likelihood estimators:
mvnorm.mle(gausscop.data)

# Now we are ready to use the bnlearn package:
# Markov blanket:
pdag1 <- iamb(data.frame(gausscop.data), test = "cor")
graphviz.plot(pdag1)
# PC algorithm:
pdag2 <- pc.stable(data.frame(gausscop.data), whitelist = NULL, blacklist = NULL, 
                   test = NULL, alpha = 0.05, B = NULL, max.sx = NULL,
                   debug = FALSE, undirected = FALSE)
graphviz.plot(pdag2)



# The other way of simulating is using the theorem:
# of equivalence between the correlation matrix and 
# the set of partial correlation of the vine.
r12 = C12$par
r13 = C13$par
r34 = C34$par
r24.3 = C24.3$par

r23 = r12*r13
r24 = r23*r34+r24.3*sqrt((1-r23^2)*(1-r34^2))
r14 = r13*r34+r12*r24.3*(1-r13^2)*sqrt(1-r34^2)/sqrt(1-r23^2)

R = matrix(data=c(1,r12,r13,r14,
                  r12,1,r23,r24,
                  r13,r23,1,r34,
                  r14,r24,r34,1),ncol=4)

cor2cov <- function(R, S) {
  sweep(sweep(R, 1, S, "*"), 2, S, "*")
}

Sigma <- cor2cov(R,c(f1$sd,f2$sd,f3$sd,f4$sd))
# It is pretty similar to the values that we obtained
# through maximum likelihood simulating in the first way.

gausscop.data2 <- mvrnorm(n=5000, mu = c(f1$mu,f2$mu,f3$mu,f4$mu),
                         Sigma = Sigma)

 
# Now we are ready to use the bnlearn package:
# Markov blanket:
pdag1 <- iamb(data.frame(gausscop.data2), test = "cor")
graphviz.plot(pdag1)
# PC algorithm:
pdag2 <- pc.stable(data.frame(gausscop.data2), whitelist = NULL, blacklist = NULL, 
                   test = NULL, alpha = 0.05, B = NULL, max.sx = NULL,
                   debug = FALSE, undirected = FALSE)
graphviz.plot(pdag2)





##############
# We construct BBN2-> chordal graph:
X <- paste0("[X1][X2|X1][X3|X1:X2][X4|X2:X3]")
dag<-model2network(X)
dag
graphviz.plot(dag)


Matrix <- c(1, 4, 3, 2,
            0, 2, 4, 3,
            0, 0, 3, 4,
            0, 0, 0, 4)
Matrix <- matrix(Matrix, 4, 4)
# define R-vine family matrix
family <- c(0, 0, 1, 1, 
            0, 0, 1, 1, 
            0, 0, 0, 1,
            0, 0, 0, 0)
family <- matrix(family, 4, 4)
# define R-vine parameter matrix
par <- c(0, 0, 0.5, 0.25,
         0, 0, 0.75, 0.5, 
         0, 0, 0, -0.5,
         0, 0, 0, 0)
par <- matrix(par, 4, 4)
# define second R-vine pair-copula parameter matrix
par2 <- matrix(0, 4, 4)

# define RVineMatrix object
RVM <- RVineMatrix(Matrix = Matrix, family = family,
                   par = par, par2 = par2,
                   names = c("U1", "U2", "U3", "U4"))
summary(RVM)
plot(RVM, tree = "ALL",type = 0)
contour(RVM, tree= "ALL")

set.seed(1)
n=15000
databbn2 <- RVineSim(n, RVM)

gausscop.data2 <- matrix(nrow = n, ncol = 4)
colnames(gausscop.data2) <- c("X1","X2","X3","X4")
for(i in 1:4){
  gausscop.data2[,i] = qnorm(databbn2[,i],
                            eval(parse(text = paste("f",as.character(i),sep="")))$mu,
                            eval(parse(text = paste("f",as.character(i),sep="")))$sd)
}

# Checking the maximum likelihood estimators:
mvnorm.mle(gausscop.data2)

# Now we are ready to use the bnlearn package:
# Markov blanket:
pdag1 <- iamb(data.frame(gausscop.data2), test = "cor")
graphviz.plot(pdag1)
# PC algorithm:
pdag2 <- pc.stable(data.frame(gausscop.data2), whitelist = NULL, blacklist = NULL, 
                   test = NULL, alpha = 0.05, B = NULL, max.sx = NULL,
                   debug = FALSE, undirected = FALSE)
graphviz.plot(pdag2)

pdag3 = hc(data.frame(gausscop.data2), score = "bic-g")
graphviz.plot(pdag3)
# The score based is the only one able to recover the v-structs 
# when set.seed(1).


# The other way of simulating is using the theorem:
# of equivalence between the correlation matrix and 
# the set of partial correlation of the vine.
r12 = par[4,1]
r23 = par[4,2]
r34 = par[4,3]
r13.2 = par[3,1]
r24.3 = par[3,2]

r13 = r12*r23+r13.2*sqrt((1-r12^2)*(1-r23^2))
r24 = r23*r34+r24.3*sqrt((1-r23^2)*(1-r34^2))
r14 = r13*r34+r24.3*(r12-r13*r23)*sqrt(1-r34^2)/sqrt(1-r23^2)


R = matrix(data=c(1,r12,r13,r14,
                  r12,1,r23,r24,
                  r13,r23,1,r34,
                  r14,r24,r34,1),ncol=4)

cor2cov <- function(R, S) {
  sweep(sweep(R, 1, S, "*"), 2, S, "*")
}

Sigma <- cor2cov(R,c(f1$sd,f2$sd,f3$sd,f4$sd))
# Result really close to the maximum likelihood estimator
# obtained for the simulated data using the 1st way.

set.seed(1)
gausscop.data3 <- mvrnorm(n=15000, mu = c(f1$mu,f2$mu,f3$mu,f4$mu),
                          Sigma = Sigma)


# Now we are ready to use the bnlearn package:
# Markov blanket:
pdag3 <- iamb(data.frame(gausscop.data3), test = "cor")
graphviz.plot(pdag3)
# PC algorithm:
pdag4 <- pc.stable(data.frame(gausscop.data3), whitelist = NULL, blacklist = NULL, 
                   test = NULL, alpha = 0.05, B = NULL, max.sx = NULL,
                   debug = FALSE, undirected = FALSE)
graphviz.plot(pdag4)

# Hill climbing algorithm:
pdag5 = hc(data.frame(gausscop.data3), score = "bic-g")
graphviz.plot(pdag5)
# The scored based is the only one able to recover the v-structures
# when setting the seed(1).



