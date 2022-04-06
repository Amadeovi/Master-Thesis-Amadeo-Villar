library("bnlearn")
library("Rgraphviz")
library("gRain")

# We construct BBN2.
X <- paste0("[X1][X2|X1][X3|X1:X2][X4|X2:X3]")
dag<-model2network(X)
dag
graphviz.plot(dag)

# After showing the factorisation of this BBNs using copulas
# we proceed to simulate from it.

# We first have to create the Vine object
# To do so, we have to specify the different matrices that define the Vine
# define the R-vine structure
Matrix <- c(1, 4, 3, 2,
            0, 2, 4, 3,
            0, 0, 3, 4,
            0, 0, 0, 4)
Matrix <- matrix(Matrix, 4, 4)
# define R-vine family matrix
family <- c(0, 0, 1, 3, 
            0, 0, 1, 3, 
            0, 0, 0, 4,
            0, 0, 0, 0)
family <- matrix(family, 4, 4)
# define R-vine parameter matrix
par <- c(0, 0, 0.5, 2,
         0, 0, 0.75, 3, 
         0, 0, 0, 1.5,
         0, 0, 0, 0)
par <- matrix(par, 4, 4)
# define second R-vine pair-copula parameter matrix
par2 <- matrix(0, 4, 4)

# define RVineMatrix object
RVM <- RVineMatrix(Matrix = Matrix, family = family,
                   par = par, par2 = par2,
                   names = c("X1", "X2", "X3", "X4"))
summary(RVM)
plot(RVM, tree = "ALL",type = 0)
contour(RVM, tree= "ALL")


# simulate a sample of size 500 from the R-vine copula model
set.seed(123)
databbn2 <- RVineSim(500, RVM)

# We have this nice graphs to get some information about the copulas.
cdata = as.copuladata(databbn2)
pairs(cdata)

# We can also compute the loglikelihood of this object:
loglik<-sum(log(RVinePDF(databbn2, RVM, verbose = TRUE)))
print(loglik)

# Moreover we can make different estimations depending on the 
# ammount of information that we have:

# If we know nothing about the vine model:
system.time(RVM1<-RVineStructureSelect(databbn2, familyset = c(0:4), type = 0,
                                       selectioncrit = "logLik", indeptest = FALSE, 
                                       level = 0.05, trunclevel = NA, progress = FALSE,
                                       weights = NA, treecrit = "tau", rotations = TRUE, 
                                       se = FALSE, presel = TRUE, method = "mle", cores = 1))
summary(RVM1)
print(RVM1$logLik)

# If we know the structure of the vine model:
# assuming that we know the structure, the results are far better:
system.time(RVM2<-RVineCopSelect(databbn2, familyset = c(1:4), Matrix,
                                 selectioncrit = "logLik",indeptest = FALSE, level = 0.05,
                                 trunclevel = NA,weights = NA, rotations = TRUE,
                                 se = FALSE,presel = TRUE,method = "mle", cores = 1))
summary(RVM2)
print(RVM2$logLik)




