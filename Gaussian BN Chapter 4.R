# Code chapter 4
# Section 4.1: Gaussian Bayesian Networks.

library("bnlearn")    # Gaussian and discrete BN
library("Rgraphviz")  # Plotting DAGs
library("gRain")
library("cubature")   # Numerical Integration
library("rlang")
library("rbmn")       # Reseaux bay?siens multinormaux
library("lattice")
library("igraph")
library("penalized")  # Lasso and Ridge regression
library("corpcor")    # inverse correlation matrix
library("compare")    # comparing TP, FP and FN
library("GGMselect")  # Gaussian Graphs Models Selection


# Graphical Representation:
par(mfrow=c(1,1))
X <- paste0("[X1][X2|X1][X3|X1][X4|X2:X3]")
dag.bnlearn <- model2network(X)
dag.bnlearn
graphviz.plot(dag.bnlearn)


# Probabilistic Representation:
X1.dist <- list(coef = c("(Intercept)" = 1), sd = 1)
X2.dist <- list(coef = c("(Intercept)" = 1, X1 = 1), sd = 1)
X3.dist <- list(coef = c("(Intercept)" = 1, X1 = 1), sd = 1)
X4.dist <- list(coef = c("(Intercept)" = 1, X2 = 1, X3 = 1), sd = 1)
dist.list <- list(X1 = X1.dist, X2 = X2.dist, X3 = X3.dist, X4 = X4.dist)


# Defining the Gaussian Bayesian Network:
gbn.bnlearn <- custom.fit(dag.bnlearn, dist = dist.list)

# We can check the local distributions using:
gbn.bnlearn$X1
gbn.bnlearn$X4

# Moreover we can check the global distribution
# using the rbmn package:
gbn.rbmn <- bnfit2nbn(gbn.bnlearn)
gema.rbmn <- nbn2gema(gbn.rbmn)
mn.rbmn <- gema2mn(gema.rbmn)
print8mn(mn.rbmn)



# We can also generate data from this BBN using the 
# following commands:
set.seed(1)
start.time <- Sys.time()
gaussian.data <- rbn(gbn.bnlearn, n = 2500)
end.time <- Sys.time()
end.time-start.time


# Knowing the structure, we can estimate parameters:
bbn.fitted <- bn.fit(dag.bnlearn, data = gaussian.data)
# The way of doing it is through maximum likelihood:
print(bbn.fitted)
print(bbn.fitted$X4)
print(gbn.bnlearn)


# We can check that using linear regression we obtain
# the same results.
bbn.fitted$X4 <- lm(X4 ~ X2 + X3, data = gaussian.data)

# We can also fit a lasso, ridge or elastic net
bbn.fitted$X4 <- penalized(X4 ~ X2 + X3, lambda1 = 1.6, lambda2 = 1.5,
                           data = gaussian.data)

bbn.fitted$X4

mvn(gaussian.data,subset = NULL,mvnTest = "mardia",covariance = TRUE,
    tol = 1e-25,alpha = 0.5,scale = FALSE,desc = TRUE,
    transform = "none", R = 1000,univariateTest = "AD",
    multivariateOutlierMethod = "none",bc = FALSE,
    bcType = "rounded",showOutliers = FALSE,
    showNewData = FALSE)$multivariateNormality



# The following section is Structure Learning in DAGs
###############
# Constraint based algorithms:
# Conditional indepence test based on the correlation matrix
# X2 indep X3|X1
cormat <- cor(gaussian.data[, c("X1", "X2", "X3")])
invcor <- cor2pcor(cormat)
dimnames(invcor) <- dimnames(cormat)
# rho_{X2,X3|X1} is then given by:
invcor[2,3]


# X1 indep X4|X2,X3
cormat <- cov(gaussian.data[, c("X1", "X2", "X3","X4")])
solve(cov(gaussian.data[, c("X1", "X2", "X3","X4")]))
solve(cov(gaussian.data[, c("X1", "X2", "X3")]))
invcor <- cor2pcor(cormat)
dimnames(invcor) <- dimnames(cormat)
# rho_{X1,X4|X2,X3} is then given by:
invcor[1,4]

# There are already tests implemented on the bnlearn package
# to deal with these issues:
ci.test("X2", "X3", "X1", test = "cor", data = gaussian.data)
ci.test("X1", "X2", "X3", test = "cor", data = gaussian.data)
ci.test("X1","X4",c("X2","X3"), test = "cor", data = gaussian.data)


ggpairs(gaussian.data,lower = list(
  continuous = wrap("points", alpha = 0.3, size=0.4)))


# PC algorithm:
set.seed(1)
cpdag1 <- pc.stable(gaussian.data, whitelist = NULL, blacklist = NULL, test = NULL,
          alpha = 0.05, B = NULL, max.sx = NULL, debug = TRUE, undirected = FALSE)

shd(cpdag(dag.bnlearn), cpdag(cpdag1), wlbl = FALSE, debug = TRUE)
par(mfrow=c(1,2))
graphviz.compare(cpdag(dag.bnlearn), cpdag(cpdag1))



# Hill climbing algorithm:
par(mfrow=c(1,1))
set.seed(1)
cpdag2 <- hc(gaussian.data , start = NULL, whitelist = NULL, blacklist = NULL, score = "bic-g",
             debug = TRUE, restart = 0, perturb = 1, max.iter = Inf, maxp = Inf, optimized = TRUE)
graphviz.plot(cpdag2)
debugging.output =
  capture.output(hc(gaussian.data , start = NULL, whitelist = NULL, blacklist = NULL, score = "bic-g",
                    debug = TRUE, restart = 0, perturb = 1, max.iter = Inf, maxp = Inf, optimized = TRUE))

head(grep("^\\* (best operation)", debugging.output, value = TRUE), 15)
head(grep("^\\* (current score)", debugging.output, value = TRUE), 15)

par(mfrow=c(2,3))
for(i in 2:7){
  graphviz.plot(model2network(gsub(
    " ","",debugging.output[which(grepl("nodes:", debugging.output))-1])[i]),
    main = sub(".*:", "",  grep("^\\* (best operation)", debugging.output, value = TRUE)[i-1]),
    sub = grep("^\\* (current score)", debugging.output, value = TRUE)[i])
}

par(mfrow=c(1,2))
set.seed(1)
cpdag3 <- hc(gaussian.data , start = NULL, whitelist = NULL, blacklist = NULL, score = "bic-g",
             debug = FALSE, restart = 20, perturb = 1, max.iter = Inf, maxp = Inf, optimized = TRUE)

par(mfrow=c(1,2))
graphviz.compare(cpdag(dag.bnlearn), cpdag(cpdag2),
                 main=c("Original graph","No random restarts"),
                 sub=c("",paste("score: ",as.character(score(cpdag2, data = gaussian.data, type = "bic-g")))))
                   
graphviz.compare(cpdag(dag.bnlearn), cpdag(cpdag3),
                 main=c("Original graph","20 random restarts"),
                 sub=c("",paste("score: ",as.character(score(cpdag3, data = gaussian.data, type = "bic-g")))))



# checking score equivalence
score(model2network("[X1][X2|X1][X3|X1][X4|X2:X3]"),
      data = gaussian.data, type = "bic-g")
score(model2network("[X2][X1|X2][X3|X1][X4|X2:X3]"),
      data = gaussian.data, type = "bic-g")
score(model2network("[X3][X1|X3][X2|X1][X4|X2:X3]"),
      data = gaussian.data, type = "bic-g")



# Learning with white list:
wl <- matrix(c("X1", "X1",
               "X2","X3"), ncol = 2, nrow = 2)
pdag2 <- iamb(gaussian.data200, test = "cor", whitelist = wl)
graphviz.plot(pdag2)
all.equal(dag.bnlearn, pdag2)



# Graphical Lasso:
# Still need to understand how does it work exactly.
g.data20k.mat <- data.matrix(gaussian.data20k)
GRest <- selectFast(g.data20k.mat, family=c("C01","LA"),K=2.5)
adjmatrix<-c(0, 1, 1, 0,
             1, 0, 1, 1,
             1, 1, 0, 1,
             0, 1, 1, 0)
adjmatrix <- matrix(adjmatrix, 4, 4)
p=30
n=20000
C01.result = all(GRest$C01$G == adjmatrix)
LA.result = all(GRest$LA$G == adjmatrix)

g <- network(GRest$LA$G)
plot(g, jitter=TRUE, usearrows = FALSE, label=1:p,displaylabels=TRUE)



GRest <- selectFast(g.data20k.mat, family=c("C01","LA"),K=2.5)
g <- network(GRest$LA$G)
plot(g, jitter=TRUE, usearrows = FALSE, label=1:p,displaylabels=TRUE)




##############3
# Another way using the glasso package:
# suspicious use....
library("glasso")
MyFamily <- NULL
for (j in 1:10){
  MyFamily[[j]] <-  abs(sign(glasso(cov(g.data20k.mat),rho=j/5)$wi))
  diag(MyFamily[[j]]) <- 0
}


GMF <- selectMyFam(g.data20k.mat,MyFamily)
gMyFam <- network(GMF$G)
plot(gMyFam, jitter=TRUE, usearrows = FALSE, label=1:p,displaylabels=TRUE)



p=7
n=30
# generate graph
eta=5
Gr <- simulateGraph(p,eta)
# generate data
X <- rmvnorm(n, mu=rep(0,p), sigma=Gr$C)
# generate a family of candidate graphs with glasso

library("glasso")
MyFamily <- NULL
for (j in 1:3){
  MyFamily[[j]] <- abs(sign(glasso(cov(X),rho=j/5)$wi))
  diag(MyFamily[[j]]) <- 0
}
# select a graph within MyFamily
GMF <- selectMyFam(X,MyFamily)
# plot the result
library(network)
par(mfrow=c(1,2))
gV <- network(Gr$G)
plot(gV,jitter=TRUE, usearrows = FALSE, label=1:p,displaylabels=TRUE)
gMyFam <- network(GMF$G)
plot(gMyFam, jitter=TRUE, usearrows = FALSE, label=1:p,displaylabels=TRUE)




