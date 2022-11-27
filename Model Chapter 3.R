# Code chapter 3: Background on Bayesian Networks.

library("bnlearn")
library("Rgraphviz")
library("gRain")

# how to construct the BBN
X<-paste0("[X1][X2][X3|X1:X2][X4|X2][X5|X1:X3][X6|X3:X4:X5]",
          "[X7|X5][X8|X6:X7]")
dag<-model2network(X)
dag
graphviz.plot(dag)

graphviz.plot(cpdag(dag))

# Notice that if we include an arc from X8 to X5 we would have a cycle
graphviz.plot(dag)

# We now want to study d-separation and conditional independencies
hlight <- list(nodes = nodes(dag), arcs = arcs(dag), col = "grey",
               textCol = "grey")
pp = graphviz.plot(dag,highlight = hlight)

# We study the relation between X1 and X6, there are 3 possible trails:
#trail 1: 
edgeRenderInfo(pp)=list(col = c("X1~X5" = "black","X5~X6" = "black"),
                         lwd = c("X1~X5" = 3, "X5~X6" = 3))
nodeRenderInfo(pp)=list(col =c("X1" = "black", "X5" = "black",
                           "X6" = "black"),
                       textCol =c("X1" = "black", "X5" = "black",
                           "X6" = "black"),
                       fill = c("X5" = "grey"))
renderGraph(pp)

# trail 2:
edgeRenderInfo(pp)=list(col = c("X1~X3" = "black","X3~X6" = "black"),
                        lwd = c("X1~X3" = 3, "X3~X6" = 3))
nodeRenderInfo(pp)=list(col =c("X1" = "black", "X3" = "black",
                               "X6" = "black"),
                        textCol =c("X1" = "black", "X3" = "black",
                                   "X6" = "black"),
                        fill = c("X3" = "grey"))
renderGraph(pp)

# trail 3:
edgeRenderInfo(pp)=list(col = c("X1~X3" = "black","X2~X3" = "black",
                                "X2~X4" = "black","X4~X6" = "black" ),
                        lwd = c("X1~X3" = 3,"X2~X3" = 3,
                                "X2~X4" = 3,"X4~X6" = 3 ))
nodeRenderInfo(pp)=list(col =c("X1" = "black", "X3" = "black",
                        "X2" = "black","X4" = "black","X6" = "black"),
                        textCol =c("X1" = "black", "X3" = "black",
                        "X2" = "black","X4" = "black","X6" = "black"),
                        fill = c("X3" = "grey"))
renderGraph(pp)







# d-separation:
dsep(dag,"X1","X6",c("X3","X5"))
dsep(dag,"X1","X6",c("X3","X5","X2"))
dsep(dag,"X1","X6",c("X5","X2"))
dsep(dag,"X1","X6",c("X5","X4"))
# independence
cropnodes<-nodes(dag)
for(n1 in cropnodes){
  for(n2 in cropnodes){
    if(dsep(dag,n1,n2)){
      cat(n1,"and",n2,"are independent.\n")
    }
  }
}

# skeleton and v-structures
skel<-skeleton(dag)
graphviz.plot(skel)
v<-vstructs(dag)

moral<-moral(dag)
graphviz.plot(moral)

for(i in 1:dim(v)[1]){
  skel<-set.arc(skel,v[i,1],v[i,2])
  skel<-set.arc(skel,v[i,3],v[i,2])
}
hlight <- list(nodes = nodes(skel), arcs = arcs(skel), col = "grey",
               textCol = "grey")
ppskel = graphviz.plot(skel,highlight = hlight)


edgeRenderInfo(ppskel)=list(col=c("X1~X3"= "black","X2~X3"= "black",
                                  "X3~X6"= "black","X4~X6"= "black",
                                  "X4~X6"= "black","X5~X6"= "black",
                                  "X6~X8"= "black","X7~X8"= "black"),
                            lwd=c("X1~X3"= 3,"X2~X3"= 3,"X3~X6"= 3,
                                  "X4~X6"= 3,"X4~X6"= 3,"X5~X6"= 3,
                                  "X6~X8"= 3,"X7~X8"= 3))
renderGraph(ppskel)


# equivalence classes.
cp1<-cpdag(dag)
equivdag <- dag
equivdag = set.arc(equivdag, from = "X4", to = "X2")
cp2<-cpdag(equivdag)
all.equal(cp1,cp2)


#markov blanquets:
mb(dag,"X5")
parents(dag,"X5")
children(dag,"X5")



hlight <- list(nodes = setdiff(nodes(dag),mb(dag,"X5")),  col = "black",
               textCol = "grey")

ppmb = graphviz.plot(dag,highlight = hlight)

nodeRenderInfo(ppmb)=list(col =c("X5" = "black"),
                        textCol =c("X5" = "black"),
                        fill = c("X5" = "grey"))

renderGraph(ppmb)




