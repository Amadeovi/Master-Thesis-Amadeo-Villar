# Code chapter 2: Copula and Vine Copula Models.

family.list <- list("0" = "Independent", "1" = "Gaussian", "2"="t","3" ="Clayton",
                    "4" = "Gumbel", "5" = "Frank", "6" = "Joe")

for(i in c(1,2,3,4,5,6)){
  
  sim <- BiCopSim(1500, family=i, par=BiCopTau2Par(i,0.75), 
                  par2 = 5, obj = NULL, check.pars = TRUE)
  
  sim <- data.frame(sim)
  
  graph <- ggplot(sim, aes(x = sim[,1], y = sim[,2]))+geom_point(size=0.1)+
                  labs(title = TeX(paste("$C$",":",  family.list[as.character(i)],",",
                  '$\\tau$',"=",0.75)),
                  x = TeX("$U_{1}$"),y = TeX("$U_{2}$"))+theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(paste("copula",as.character(i),".png"), graph, path = "C:/Users/amade/OneDrive/Escritorio")
  
}



# Vine Copula in 4 Dimensions:
Matrix <- c(1, 4, 3, 2,
            0, 2, 4, 3,
            0, 0, 3, 4,
            0, 0, 0, 4)
Matrix <- matrix(Matrix, 4, 4)
# define R-vine family matrix
family <- c(0, 0, 5, 1, 
            0, 0, 6, 3, 
            0, 0, 0, 4,
            0, 0, 0, 0)
family <- matrix(family, 4, 4)
# define R-vine parameter matrix

# We choose some copulas high correlated and other low correlated
par <- c(0, 0, BiCopTau2Par(5,0.75), BiCopTau2Par(1,0.25),
         0, 0, BiCopTau2Par(6,0.75), BiCopTau2Par(3,0.25), 
         0, 0, 0, BiCopTau2Par(4,0.75),
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

set.seed(123)
n <- 2500

start.time1 <- Sys.time()
vinedata <- RVineSim(n, RVM)
end.time1 <- Sys.time()
end.time1-start.time1

# The logLik is given by:
logLik <- sum(log(RVinePDF(vinedata, RVM, verbose = TRUE)))


cdata <- as.copuladata(vinedata)
pairs(cdata)



# Estimation procedure:
start.time3 <- Sys.time()
RVMcop1 <- RVineCopSelect(vinedata, familyset = c(0,1,2,3:6), Matrix, 
                          selectioncrit = "BIC", indeptest = TRUE,
                          level = 0.005, trunclevel = NA,
                          weights = NA, rotations = TRUE, se = FALSE,
                          presel = TRUE, method = "mle", cores = 1)

end.time3 <- Sys.time()
end.time3-start.time3

summary(RVMcop1)
RVMcop1$logLik


# Selecting the structure:
start.time4 <- Sys.time()
RVMcop2 <- RVineStructureSelect(vinedata, familyset = c(0,1,3:6), type = 0, selectioncrit = "AIC",
                                indeptest = TRUE, level = 0.05, trunclevel = NA, progress = FALSE,
                                weights = NA, treecrit = "tau", rotations = TRUE, se = FALSE,
                                presel = TRUE, method = "mle", cores = 1)

end.time4 <- Sys.time()
end.time4-start.time4

summary(RVMcop2)

RVMcop2$logLik

RVineGofTest(vinedata, RVMcop2, method = "Breymann", statistic = "CvM",alpha = 2)


