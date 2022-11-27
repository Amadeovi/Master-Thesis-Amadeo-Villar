# Simulation study chapter 6: Tables and Figur Small Network.

# This script can be used to build all the data frames
# and figures that we will show in simulation study
# Chapter 6 of the dissertation - > SMALL NETWORK.

# The name of the data frames will be like smallnetwork.Gauss500
# or bignetwork.nonGauss2500, to refer to which network,
# whether it is Gaussian or not and the sample size used.



# Building data frames for SHD and Hamming distance for pc, hc and glasso
# depending on the sample size used:
# This is more like a function that should be applied different times.
N = 15000
DF.Gauss <- bignetwork.Gauss15000
DF.nonGauss <- bignetwork.nonGauss15000

# Change the column names accordingly
colnames(DF.Gauss) <- c("C12", "C23", "C34", "C25", "C56", "C67", "C68", "C59","C910",
                        "$\tau_{1,2}$", "$\tau_{2,3}$", "$\tau_{3,4}$", "$\tau_{2,5}$",
                        "$\tau_{5,6}$","$\tau_{6,7}$","$\tau_{6,8}$", 
                        "$\tau_{5,9}$", "$\tau_{9,10}$",
                        "pc.shd", "pc.ham", "pc.arcs", "pc.sim.time",
                        "hc.shd", "hc.ham", "hc.arcs", "hc.sim.time",
                        "glasso.ham", "glasso.edges","glasso.sim.time","True.glasso")

colnames(DF.nonGauss) <- c("C12", "C23", "C34", "C25", "C56", "C67", "C68", "C59","C910",
                           "$\tau_{1,2}$", "$\tau_{2,3}$", "$\tau_{3,4}$", "$\tau_{2,5}$",
                           "$\tau_{5,6}$","$\tau_{6,7}$","$\tau_{6,8}$", 
                           "$\tau_{5,9}$", "$\tau_{9,10}$",
                           "pc.shd", "pc.ham", "pc.arcs", "pc.sim.time",
                           "hc.shd", "hc.ham", "hc.arcs", "hc.sim.time",
                           "glasso.ham", "glasso.edges","glasso.sim.time","True.glasso")

####### We first prepare the data:
########## Gaussian Case
# Preparing the data SHD:
data.gauss.shd <- melt(DF.Gauss[c("pc.shd","hc.shd")],id=NULL)
mu.gauss.shd <- ddply(data.gauss.shd, .(variable), summarise, mean=round(mean(value),3))  
median.gauss.shd <- ddply(data.gauss.shd, "variable", summarise, median=round(median(value),3))
gauss.shd <- merge(mu.gauss.shd,median.gauss.shd)

# Preparing the data Hamming:
data.gauss.ham <- melt(DF.Gauss[c("pc.ham","hc.ham","glasso.ham")],id=NULL)
mu.gauss.ham <- ddply(data.gauss.ham , "variable", summarise, mean=round(mean(value),3))
median.gauss.ham <- ddply(data.gauss.ham, "variable", summarise, median=round(median(value),3))
gauss.ham <- merge(mu.gauss.ham,median.gauss.ham)

# Preparing arcs:
data.gauss.arc <-  melt(DF.Gauss[c("pc.arcs","hc.arcs","glasso.edges")],id=NULL)
mu.gauss.arc <- ddply(data.gauss.arc , "variable", summarise, mean=round(mean(value),3))
median.gauss.arc <- ddply(data.gauss.arc, "variable", summarise, median=round(median(value),3))
gauss.arc <- merge(mu.gauss.arc,median.gauss.arc)




########## Non Gaussian Case
# Preparing the data SHD:
data.nongauss.shd <- melt(DF.nonGauss[c("pc.shd","hc.shd")],id=NULL)
mu.nongauss.shd <- ddply(data.nongauss.shd, .(variable), summarise, mean=round(mean(value),3))  
median.nongauss.shd <- ddply(data.nongauss.shd, "variable", summarise, median=round(median(value),3))
nongauss.shd <- merge(mu.nongauss.shd,median.nongauss.shd)

# Preparing the data Hamming:
data.nongauss.ham <- melt(DF.nonGauss[c("pc.ham","hc.ham","glasso.ham")],id=NULL)
mu.nongauss.ham <- ddply(data.nongauss.ham , "variable", summarise, mean=round(mean(value),3))
median.nongauss.ham <- ddply(data.nongauss.ham, "variable", summarise, median=round(median(value),3))
nongauss.ham <- merge(mu.nongauss.ham,median.nongauss.ham)

# Preparing arcs:
data.nongauss.arc <-  melt(DF.nonGauss[,c("pc.arcs","hc.arcs","glasso.edges")],id=NULL)
mu.nongauss.arc <- ddply(data.nongauss.arc , "variable", summarise, mean=round(mean(value),3))
median.nongauss.arc <- ddply(data.nongauss.arc, "variable", summarise, median=round(median(value),3))
gauss.arc <- merge(mu.nongauss.arc,median.nongauss.arc)





# Data Frame SHD and Hamming distance for glasso, constraint and scored showed below:
kable(cbind(rbind(gauss.shd,gauss.ham),
            rbind(nongauss.shd,nongauss.ham)[,2:3])[c(3,2,5,1,4),] ,"latex")  


###### Let's plot all these previous results:
# Gaussian case:
# Plotting SHD: pc and hc.
ggplot(data.gauss.shd , aes(x=value, color=variable, fill=variable)) +
  geom_density(alpha=.2) + scale_x_continuous(name="SHD", limits=c(0, 4))+
  scale_fill_manual(values=c("#F8766D", "#00BFC4"))+
  scale_color_manual(values = c("#F8766D", "#00BFC4"))+
  geom_vline(data=median.gauss.shd, aes(xintercept=median, color=variable),linetype="dashed")+
  ggtitle(paste("SHD: Constraint vs Scored, Gaussian data, N=",N))+
  annotate("text", color="red",x = 3, y = 0.55, label = paste ("Median ==", median.gauss.shd$median[1]), parse = TRUE) +
  annotate("text", color="blue",x = 3, y = 0.40, label = paste ("Median ==",  median.gauss.shd$median[2]), parse = TRUE) 


# Plotting Ham: pc and hc.
ggplot(data.gauss.ham, aes(x=value, color=variable, fill=variable)) +
  geom_density(alpha=.2) + scale_x_continuous(name="Hamming distance", limits=c(0, 0.5))+
  scale_fill_manual(values=c("#F8766D", "#00BFC4","#00BA38"))+
  scale_color_manual(values = c("#F8766D", "#00BFC4","#00BA38"))+
  geom_vline(data=median.gauss.ham, aes(xintercept=median, color=variable),linetype="dashed")+
  ggtitle(paste("Ham: glasso vs Constraint vs Scored, Gaussian data, N=",N))+
  annotate("text", color="red",x = 0.4, y = 15, label = paste ("Median ==", median.gauss.ham$median[1]), parse = TRUE) +
  annotate("text", color="blue",x = 0.4, y = 10, label = paste ("Median ==", median.gauss.ham$median[2]), parse = TRUE) +
  annotate("text", color="#00BA38",x = 0.4, y = 5, label = paste ("Median ==", median.gauss.ham$median[3]), parse = TRUE) 


# Plotting number of arcs:
ggplot(data.gauss.arc, aes(x=value, color=variable, fill=variable)) +
  geom_density(alpha=.2) + scale_x_continuous(name="edges", limits=c(8,10))+
  scale_fill_manual(values=c("#F8766D", "#00BFC4","#00BA38"))+
  scale_color_manual(values = c("#F8766D", "#00BFC4","#00BA38"))+
  geom_vline(data=median.gauss.arc, aes(xintercept=median, color=variable),linetype="dashed")+
  ggtitle(paste("|E| edges: glasso vs Constraint vs Scored, Gaussian data, N=",N))+
  annotate("text", color="red",x = 9.7, y = 12, label = paste ("Median ==", median.gauss.arc$median[1]), parse = TRUE) +
  annotate("text", color="blue",x = 9.7, y = 9, label = paste ("Median ==", median.gauss.arc$median[2]), parse = TRUE) +
  annotate("text", color="#00BA38",x = 9.7, y = 6, label = paste ("Median ==", median.gauss.arc$median[3]), parse = TRUE)




# Non-Gaussian case:
# Plotting SHD
ggplot(data.nongauss.shd , aes(x=value, color=variable, fill=variable)) +
  geom_density(alpha=.2) + scale_x_continuous(name="SHD", limits=c(0, 15))+
  scale_fill_manual(values=c("#F8766D", "#00BFC4"))+
  scale_color_manual(values = c("#F8766D", "#00BFC4"))+
  geom_vline(data=median.nongauss.shd, aes(xintercept=median, color=variable),linetype="dashed")+
  ggtitle(paste("SHD: Constraint vs Scored, non Gaussian data, N=",N))+
  annotate("text", color="red",x = 5, y = 0.50, label = paste ("Median ==", median.nongauss.shd$median[1]), parse = TRUE) +
  annotate("text", color="blue",x = 5, y = 0.46, label = paste ("Median ==",  median.nongauss.shd$median[2]), parse = TRUE) 


# Plotting Ham
ggplot(data.nongauss.ham, aes(x=value, color=variable, fill=variable)) +
  geom_density(alpha=.2) + scale_x_continuous(name="Hamming distance", limits=c(0, 10))+
  scale_fill_manual(values=c("#F8766D", "#00BFC4","#00BA38"))+
  scale_color_manual(values = c("#F8766D", "#00BFC4","#00BA38"))+
  geom_vline(data=median.nongauss.ham, aes(xintercept=median, color=variable),linetype="dashed")+
  ggtitle(paste("Ham.d, Constraint vs Scored, non-Gaussian data, N=",N))+
  annotate("text", color="red",x = 7, y = 0.43, label = paste ("Median ==", median.nongauss.ham$median[1]), parse = TRUE) +
  annotate("text", color="blue",x = 7, y = 0.40, label = paste ("Median ==", median.nongauss.ham$median[2]), parse = TRUE) +
  annotate("text", color="darkgreen",x = 7, y = 0.37, label = paste ("Median ==", median.nongauss.ham$median[3]), parse = TRUE)


# Number of arcs
ggplot(data.nongauss.arc, aes(x=value, color=variable, fill=variable)) +
  geom_density(alpha=.2) + scale_x_continuous(name="edges", limits=c(7,18))+
  scale_fill_manual(values=c("#F8766D", "#00BFC4","#00BA38"))+
  scale_color_manual(values = c("#F8766D", "#00BFC4","#00BA38"))+
  geom_vline(data=median.nongauss.arc, aes(xintercept=median, color=variable),linetype="dashed")+
  ggtitle(paste("|E| edges: glasso vs Constraint vs Scored, non-Gaussian data, N=",N))+
  annotate("text", color="red",x = 13, y = 0.65, label = paste ("Median ==", median.nongauss.arc$median[1]), parse = TRUE) +
  annotate("text", color="blue",x = 13, y = 0.6, label = paste ("Median ==", median.nongauss.arc$median[2]), parse = TRUE) +
  annotate("text", color="#00BA38",x = 13, y = 0.55, label = paste ("Median ==", median.nongauss.arc$median[3]), parse = TRUE) 





#### Building data frames 
#### Computing SHD and Ham distances for each of the considered cases:
#### SMALL NETWORK:
kable(cbind(ddply(melt(smallnetwork.Gauss500[c("glasso.ham","pc.shd","pc.ham","hc.shd","hc.ham")],id=NULL)
                  ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                  quantile.95=round(quantile(value,0.95),3)),
            ddply(melt(smallnetwork.nonGauss500[c("glasso.ham","pc.shd","pc.ham","hc.shd","hc.ham")],id=NULL)
                  ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                  quantile.95=round(quantile(value,0.95),3))[,c(2,3,4)]),"latex")



kable(cbind(ddply(melt(smallnetwork.Gauss2500[c("glasso.ham","pc.shd","pc.ham","hc.shd","hc.ham")],id=NULL)
                  ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                  quantile.95=round(quantile(value,0.95),3)),
            ddply(melt(smallnetwork.nonGauss2500[c("glasso.ham","pc.shd","pc.ham","hc.shd","hc.ham")],id=NULL)
                  ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                  quantile.95=round(quantile(value,0.95),3))[,c(2,3,4)]),"latex")


kable(cbind(ddply(melt(smallnetwork.Gauss7500[c("glasso.ham","pc.shd","pc.ham","hc.shd","hc.ham")],id=NULL)
                  ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                  quantile.95=round(quantile(value,0.95),3)),
            ddply(melt(smallnetwork.nonGauss7500[c("glasso.ham","pc.shd","pc.ham","hc.shd","hc.ham")],id=NULL)
                  ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                  quantile.95=round(quantile(value,0.95),3))[,c(2,3,4)]),"latex")


kable(cbind(ddply(melt(smallnetwork.Gauss15000[c("glasso.ham","pc.shd","pc.ham","hc.shd","hc.ham")],id=NULL)
                  ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                  quantile.95=round(quantile(value,0.95),3)),
            ddply(melt(smallnetwork.nonGauss15000[c("glasso.ham","pc.shd","pc.ham","hc.shd","hc.ham")],id=NULL)
                  ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                  quantile.95=round(quantile(value,0.95),3))[,c(2,3,4)]),"latex")



##### BIG NETWORK:

kable(cbind(ddply(melt(bignetwork.Gauss500[c("glasso.ham","pc.shd","pc.ham","hc.shd","hc.ham")],id=NULL)
                  ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                  quantile.95=round(quantile(value,0.95),3)),
            ddply(melt(bignetwork.nonGauss500[c("glasso.ham","pc.shd","pc.ham","hc.shd","hc.ham")],id=NULL)
                  ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                  quantile.95=round(quantile(value,0.95),3))[,c(2,3,4)]),"latex")


kable(cbind(ddply(melt(bignetwork.Gauss2500[c("glasso.ham","pc.shd","pc.ham","hc.shd","hc.ham")],id=NULL)
                  ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                  quantile.95=round(quantile(value,0.95),3)),
            ddply(melt(bignetwork.nonGauss2500[c("glasso.ham","pc.shd","pc.ham","hc.shd","hc.ham")],id=NULL)
                  ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                  quantile.95=round(quantile(value,0.95),3))[,c(2,3,4)]),"latex")


kable(cbind(ddply(melt(bignetwork.Gauss7500[c("glasso.ham","pc.shd","pc.ham","hc.shd","hc.ham")],id=NULL)
                  ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                  quantile.95=round(quantile(value,0.95),3)),
            ddply(melt(bignetwork.nonGauss7500[c("glasso.ham","pc.shd","pc.ham","hc.shd","hc.ham")],id=NULL)
                  ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                  quantile.95=round(quantile(value,0.95),3))[,c(2,3,4)]),"latex")


kable(cbind(ddply(melt(bignetwork.Gauss15000[c("glasso.ham","pc.shd","pc.ham","hc.shd","hc.ham")],id=NULL)
                  ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                  quantile.95=round(quantile(value,0.95),3)),
            ddply(melt(bignetwork.nonGauss15000[c("glasso.ham","pc.shd","pc.ham","hc.shd","hc.ham")],id=NULL)
                  ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                  quantile.95=round(quantile(value,0.95),3))[,c(2,3,4)]),"latex")









#### Building data frames for the Number of arcs.
# This is a script so no need of adjusting any parameters:
# Smaller network:
kable(rbind(cbind(ddply(melt(smallnetwork.Gauss500[c("glasso.arcs","pc.arcs","hc.arcs")],id=NULL)
                        ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                        quantile.95=round(quantile(value,0.95),3)),
                  ddply(melt(smallnetwork.nonGauss500[c("glasso.arcs","pc.arcs","hc.arcs")],id=NULL)
                        ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                        quantile.95=round(quantile(value,0.95),3))[,c(2,3,4)]),
            cbind(ddply(melt(smallnetwork.Gauss2500[c("glasso.arcs","pc.arcs","hc.arcs")],id=NULL)
                        ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                        quantile.95=round(quantile(value,0.95),3)),
                  ddply(melt(smallnetwork.nonGauss2500[c("glasso.arcs","pc.arcs","hc.arcs")],id=NULL)
                        ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                        quantile.95=round(quantile(value,0.95),3))[,c(2,3,4)]),
            cbind(ddply(melt(smallnetwork.Gauss7500[c("glasso.arcs","pc.arcs","hc.arcs")],id=NULL)
                        ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                        quantile.95=round(quantile(value,0.95),3)),
                  ddply(melt(smallnetwork.nonGauss7500[c("glasso.arcs","pc.arcs","hc.arcs")],id=NULL)
                        ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                        quantile.95=round(quantile(value,0.95),3))[,c(2,3,4)]),
            cbind(ddply(melt(smallnetwork.Gauss15000[c("glasso.arcs","pc.arcs","hc.arcs")],id=NULL)
                        ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                        quantile.95=round(quantile(value,0.95),3)),
                  ddply(melt(smallnetwork.nonGauss15000[c("glasso.arcs","pc.arcs","hc.arcs")],id=NULL)
                        ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                        quantile.95=round(quantile(value,0.95),3))[,c(2,3,4)])),"latex")



      
     


# Bigger Network:
kable(rbind(cbind(ddply(melt(bignetwork.Gauss500[c("glasso.arcs","pc.arcs","hc.arcs")],id=NULL)
                        ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                        quantile.95=round(quantile(value,0.95),3)),
                  ddply(melt(bignetwork.nonGauss500[c("glasso.arcs","pc.arcs","hc.arcs")],id=NULL)
                        ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                        quantile.95=round(quantile(value,0.95),3))[,c(2,3,4)]),
            cbind(ddply(melt(bignetwork.Gauss2500[c("glasso.arcs","pc.arcs","hc.arcs")],id=NULL)
                        ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                        quantile.95=round(quantile(value,0.95),3)),
                  ddply(melt(bignetwork.nonGauss2500[c("glasso.arcs","pc.arcs","hc.arcs")],id=NULL)
                        ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                        quantile.95=round(quantile(value,0.95),3))[,c(2,3,4)]),
            cbind(ddply(melt(bignetwork.Gauss7500[c("glasso.arcs","pc.arcs","hc.arcs")],id=NULL)
                        ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                        quantile.95=round(quantile(value,0.95),3)),
                  ddply(melt(bignetwork.nonGauss7500[c("glasso.arcs","pc.arcs","hc.arcs")],id=NULL)
                        ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                        quantile.95=round(quantile(value,0.95),3))[,c(2,3,4)]),
            cbind(ddply(melt(bignetwork.Gauss15000[c("glasso.arcs","pc.arcs","hc.arcs")],id=NULL)
                        ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                        quantile.95=round(quantile(value,0.95),3)),
                  ddply(melt(bignetwork.nonGauss15000[c("glasso.arcs","pc.arcs","hc.arcs")],id=NULL)
                        ,"variable", summarise, median=round(median(value),3), quantile.5=round(quantile(value,0.05),3),
                        quantile.95=round(quantile(value,0.95),3))[,c(2,3,4)]))[,c(2:7)],"latex")



#### Building data frames for the Running times.
# This is a script so no need of adjusting any parameters:
# Smaller network:
kable(cbind(rbind(cbind(ddply(melt(smallnetwork.Gauss500[,c("glasso.sim.time","pc.sim.time","hc.sim.time")]),
                              .(variable), summarise, mean=round(mean(value),3),  median=round(median(value),3)),
                        ddply(melt(smallnetwork.nonGauss500[,c("glasso.sim.time","pc.sim.time","hc.sim.time")]),
                              .(variable), summarise, mean=round(mean(value),3),  median=round(median(value),3))[,c(2:3)]),
                  cbind(ddply(melt(smallnetwork.Gauss2500[,c("glasso.sim.time","pc.sim.time","hc.sim.time")]),
                              .(variable), summarise, mean=round(mean(value),3),  median=round(median(value),3)),
                        ddply(melt(smallnetwork.nonGauss2500[,c("glasso.sim.time","pc.sim.time","hc.sim.time")]),
                              .(variable), summarise, mean=round(mean(value),3),  median=round(median(value),3))[,c(2:3)]),
                  cbind(ddply(melt(smallnetwork.Gauss7500[,c("glasso.sim.time","pc.sim.time","hc.sim.time")]),
                              .(variable), summarise, mean=round(mean(value),3),  median=round(median(value),3)),
                        ddply(melt(smallnetwork.nonGauss7500[,c("glasso.sim.time","pc.sim.time","hc.sim.time")]),
                              .(variable), summarise, mean=round(mean(value),3),  median=round(median(value),3))[,c(2:3)]),
                  cbind(ddply(melt(smallnetwork.Gauss15000[,c("glasso.sim.time","pc.sim.time","hc.sim.time")]),
                              .(variable), summarise, mean=round(mean(value),3),  median=round(median(value),3)),
                        ddply(melt(smallnetwork.nonGauss15000[,c("glasso.sim.time","pc.sim.time","hc.sim.time")]),
                              .(variable), summarise, mean=round(mean(value),3),  median=round(median(value),3))[,c(2:3)]))[,c(1,3,5)],
            rbind(cbind(ddply(melt(bignetwork.Gauss500[,c("glasso.sim.time","pc.sim.time","hc.sim.time")]),
                              .(variable), summarise, mean=round(mean(value),3),  median=round(median(value),3)),
                        ddply(melt(bignetwork.nonGauss500[,c("glasso.sim.time","pc.sim.time","hc.sim.time")]),
                              .(variable), summarise, mean=round(mean(value),3),  median=round(median(value),3))[,c(2:3)]),
                  cbind(ddply(melt(bignetwork.Gauss2500[,c("glasso.sim.time","pc.sim.time","hc.sim.time")]),
                              .(variable), summarise, mean=round(mean(value),3),  median=round(median(value),3)),
                        ddply(melt(bignetwork.nonGauss2500[,c("glasso.sim.time","pc.sim.time","hc.sim.time")]),
                              .(variable), summarise, mean=round(mean(value),3),  median=round(median(value),3))[,c(2:3)]),
                  cbind(ddply(melt(bignetwork.Gauss7500[,c("glasso.sim.time","pc.sim.time","hc.sim.time")]),
                              .(variable), summarise, mean=round(mean(value),3),  median=round(median(value),3)),
                        ddply(melt(bignetwork.nonGauss7500[,c("glasso.sim.time","pc.sim.time","hc.sim.time")]),
                              .(variable), summarise, mean=round(mean(value),3),  median=round(median(value),3))[,c(2:3)]),
                  cbind(ddply(melt(bignetwork.Gauss15000[,c("glasso.sim.time","pc.sim.time","hc.sim.time")]),
                              .(variable), summarise, mean=round(mean(value),3),  median=round(median(value),3)),
                        ddply(melt(bignetwork.nonGauss15000[,c("glasso.sim.time","pc.sim.time","hc.sim.time")]),
                              .(variable), summarise, mean=round(mean(value),3),  median=round(median(value),3))[,c(2:3)])
            )[,c(3,5)]),"latex")






# Lasso results:
ggplot(DF.Gauss, aes(x=True.glasso))+
  geom_histogram(color="blue", fill="lightblue",binwidth=1)+
  ggtitle(paste("Graphical Lasso, Gaussian data, N=",N)) + xlab("Nº recovered true graphs out of 30")+
  theme(plot.title = element_text(hjust = 0.5))

# Lasso results:
ggplot(DF.nonGauss, aes(x=True.glasso))+
  geom_histogram(color="blue", fill="lightblue",binwidth=1)+
  ggtitle(paste("Graphical Lasso, non Gaussian data, N=",N)) + xlab("Nº recovered true graphs out of 30")+
  theme(plot.title = element_text(hjust = 0.5))





################################
##### ANALYSIS OF THE RESULTS:

#### SMALL NETWORK:

# We start looking at the glasso Gaussian:
length(which(smallnetwork.Gauss15000["glasso.ham"]>0.5)%%6)
length(which(smallnetwork.Gauss7500["glasso.ham"]>0.5)%%6)
# We check that all the scenarios 0 and 1, fails
smallnetwork.Gauss15000[which(smallnetwork.Gauss15000["glasso.ham"]>0.5),"glasso.arcs"]
smallnetwork.Gauss7500[which(smallnetwork.Gauss7500["glasso.ham"]>0.5),"glasso.arcs"]
# We then check that they are underestimating-> majority time: 4 instead of 5, 
# which one is not recovered.
smallnetwork.Gauss15000[which(smallnetwork.Gauss15000["glasso.ham"]>0.5),"glasso.ham"]
# We then check the adjacency matrix, to do so we generate a list of matrices
# obtained using Experiment 14-10-> store it under the name:
# A glasso small network gauss.csv
sapply(A.smallnetworkgauss[which(smallnetwork.Gauss15000["glasso.ham"]>0.5)],'[',2,3)
sapply(A.smallnetworkgauss[which(smallnetwork.Gauss15000["glasso.ham"]<0.5)],'[',2,3)
sum(sapply(A.smallnetworkgauss,'[',1,4))
# So all the scenarios miss the edge (2,3)



# We then check glasso nonGaussian:
length(which(smallnetwork.nonGauss15000["glasso.ham"]>1)%%6)
smallnetwork.nonGauss15000[which(smallnetwork.nonGauss15000["glasso.ham"]>1),"glasso.arcs"]
smallnetwork.nonGauss7500[which(smallnetwork.nonGauss7500["glasso.ham"]>1),"glasso.arcs"]
# Both results: underestimation and overestimation, nothing really visible.
# Let's check the 1st parameters scenario.
smallnetwork.nonGauss15000[seq(1,60,6),"glasso.arcs"]
smallnetwork.nonGauss15000[seq(6,60,6),"glasso.arcs"]
# nothing really relevant, both overestimation and underestimation.
sum(sapply(A.smallnetworknongauss,'[',1,4))
sum(sapply(A.smallnetworknongauss,'[',2,3))
sum(sapply(A.smallnetworknongauss,'[',1,3))



# We then check the pc algorithm, which seems to fail for both 
# Gaussian and non-Gaussian.
# Guassian case:
length(which(smallnetwork.Gauss15000["pc.shd"]>2)%%6)
sum(which(smallnetwork.Gauss15000["pc.shd"]>1)%%6==1)
# We can check that it always fails 1, 3 and 4, in particular 3 and 4!!
# Let's study these cases
smallnetwork.Gauss15000[which(smallnetwork.Gauss15000["pc.shd"]>2),"pc.arcs"]
smallnetwork.Gauss15000[which(smallnetwork.Gauss15000["pc.shd"]>2),"pc.ham"]
smallnetwork.Gauss15000[which(smallnetwork.Gauss15000["pc.shd"]>1),"pc.arcs"]
smallnetwork.Gauss15000[which(smallnetwork.Gauss15000["pc.shd"]>1),"pc.ham"]
# We then check the adjacency matrix, to do so we generate a list of matrices
# obtained using Experiment 17-10-> store it under the name:
# B pc small network gauss.csv
# Scenarios 3 fail to spot the edge (3,4)
sapply(B.smallnetworkgauss[which(smallnetwork.Gauss15000["pc.shd"]>2)],'[',3,4)
# Whereas scenarios 4 fail to recover the edge (1,2)
sapply(B.smallnetworkgauss[which(smallnetwork.Gauss15000["pc.shd"]>2)],'[',1,2)
# We can also check that the majority of the graphs does not include (1,4)
sum(sapply(B.smallnetworkgauss,'[',1,4))



# Non-Guassian case:
length(which(smallnetwork.nonGauss15000["pc.shd"]>2)%%6)
length(which(smallnetwork.nonGauss15000["pc.ham"]>1)%%6)
# It seems that mainly scenarios 3,4  have problems:
smallnetwork.nonGauss15000[which(smallnetwork.nonGauss15000["pc.shd"]>2),"pc.arcs"]
smallnetwork.nonGauss15000[which(smallnetwork.nonGauss15000["pc.shd"]>2),"pc.ham"]
# It seems no pattern on overestimation or underestimation.
sum(sapply(B.smallnetworknongauss[which(smallnetwork.nonGauss15000["pc.shd"]>2)],'[',1,4))
# 8/13 put an arc between 1 and 4.
sum(sapply(B.smallnetworknongauss,'[',1,4))
sum(sapply(B.smallnetworknongauss,'[',2,3))
sum(sapply(B.smallnetworknongauss,'[',3,4))
sum(sapply(B.smallnetworknongauss,'[',1,2))
# 36/60 put an arc between 1 and 4.
# Not visible patterns...



# We then check the hc algorithm, which seems to give good results
# for the gaussian, but not that good for the non-gaussian.
# Guassian case:
length(which(smallnetwork.Gauss15000["hc.shd"]>0.5)%%6)
# It seems only to fail for scenarios 0 and 1
# But the number of arcs is the correct one:
smallnetwork.Gauss15000[which(smallnetwork.Gauss15000["hc.shd"]>0.5),"hc.arcs"]
smallnetwork.Gauss15000[which(smallnetwork.Gauss15000["hc.shd"]>0.5),"hc.shd"]
smallnetwork.Gauss15000[which(smallnetwork.Gauss15000["hc.shd"]>0.5),"hc.ham"]
# Adjacency matrix:
sum(sapply(C.smallnetworkgauss,'[',2,3))
sum(sapply(C.smallnetworkgauss,'[',1,4))
sum(sapply(C.smallnetworkgauss,'[',1,3))
sum(sapply(C.smallnetworkgauss,'[',3,4))

sapply(C.smallnetworkgauss[seq(1,60,6)],'[',1,4)
sapply(C.smallnetworkgauss[seq(6,60,6)],'[',2,3)
sapply(C.smallnetworkgauss[seq(6,60,6)],'[',1,4)
# Not really observable patterns in our data.
mean(smallnetwork.Gauss15000[seq(1,60,6),"hc.shd"])
mean(smallnetwork.Gauss15000[seq(2,60,6),"hc.shd"])
mean(smallnetwork.Gauss15000[seq(5,60,6),"hc.shd"])
mean(smallnetwork.Gauss15000[seq(6,60,6),"hc.shd"])



# Non-gaussian case:
length(which(smallnetwork.nonGauss15000["hc.shd"]>0.9)%%6)
length(which(smallnetwork.nonGauss15000["hc.shd"]>1)%%6)
length(which(smallnetwork.nonGauss15000["hc.ham"]>1)%%6)
# It mainly fails scenarios 1,3,4,0
smallnetwork.nonGauss15000[which(smallnetwork.nonGauss15000["hc.shd"]>1),"hc.arcs"]
smallnetwork.nonGauss15000[which(smallnetwork.nonGauss15000["hc.shd"]>1),"hc.shd"]
smallnetwork.nonGauss15000[which(smallnetwork.nonGauss15000["hc.shd"]>1),"hc.ham"]
# Lets check how many
sum(sapply(C.smallnetworknongauss,'[',1,4))
sum(sapply(C.smallnetworknongauss,'[',2,3))
sum(sapply(C.smallnetworknongauss,'[',1,2))
sum(sapply(C.smallnetworknongauss,'[',3,4))
# Not really visible patterns...
mean(smallnetwork.nonGauss15000[seq(1,60,6),"hc.shd"])
mean(smallnetwork.nonGauss15000[seq(2,60,6),"hc.shd"])
mean(smallnetwork.nonGauss15000[seq(3,60,6),"hc.shd"])
mean(smallnetwork.nonGauss15000[seq(4,60,6),"hc.shd"])
mean(smallnetwork.nonGauss15000[seq(5,60,6),"hc.shd"])
mean(smallnetwork.nonGauss15000[seq(6,60,6),"hc.shd"])
# But it seems that for cases 2 and 5 in general the SHD are lower
# these scenarios have the conditional copulas low correlates
# this might be a thing to study further...





################################
##### ANALYSIS OF THE RESULTS:

#### BIGGER NETWORK:


# We start looking at the glasso Gaussian:
length(which(bignetwork.Gauss15000["glasso.ham"]>0)%%8)
# It seems to work pretty good.

# Let's look at the non-gaussian case:
length(which(bignetwork.nonGauss15000["glasso.ham"]>10)%%8)


# it seems that for the t-copula and gumbell copula works fine.
# but for Clayton and Joe the results are pretty bad.
bignetwork.nonGauss15000[seq(1,8,1),"glasso.ham"]
bignetwork.nonGauss15000[seq(9,16,1),"glasso.ham"]
bignetwork.nonGauss15000[seq(17,24,1),"glasso.ham"]
bignetwork.nonGauss15000[seq(25,32,1),"glasso.ham"]
bignetwork.nonGauss15000[seq(33,40,1),"glasso.ham"]
# Also for the low correlated case the Hamming distance is way lower.
bignetwork.nonGauss15000[seq(2,64,8),"glasso.ham"]
bignetwork.nonGauss15000[seq(1,64,8),"glasso.ham"]
# Checking whether overestimation or underestimation, we can see
# that there is always overestimation, but they are able to recover
# the true edges always!!
bignetwork.nonGauss15000["glasso.arcs"]-bignetwork.nonGauss15000["glasso.ham"]



#  Let's see whether it overestimate or what happens:
# Overestimation:
bignetwork.nonGauss15000[seq(9,16,1),"glasso.arcs"]
plot(graph.adjacency(A.bignetworknongauss[[10]],mode = "undirected"))
plot(graph.adjacency(A.bignetworknongauss[[2]],mode = "undirected"))
plot(graph.adjacency(A.bignetworknongauss[[18]],mode = "undirected"))




# We then check the pc algorithm.
# Gaussian case, seems to work pretty good:
length(which(bignetwork.Gauss15000["pc.shd"]>0.9)%%8)
bignetwork.Gauss15000[which(bignetwork.Gauss15000["pc.shd"]>0.9),"pc.ham"]
bignetwork.Gauss15000[which(bignetwork.Gauss15000["pc.shd"]>0.9),"pc.arcs"]
bignetwork.Gauss15000["pc.arcs"]-bignetwork.Gauss15000["pc.ham"]
sum(bignetwork.Gauss15000["pc.arcs"]>9.9)
# no extra edges in a general sense! -> no patterns of added edges.

# We can see that they put one edge more than they should.
B.bignetworkgauss[[6]]-adjmatrix
B.bignetworkgauss[[7]]-adjmatrix
# It puts more arcs in the few cases where it fails. 
# But the performance is quite good, even better than in the smaller network.



# Non-Gaussian case:
# Lets analyse only the case N=15000, where it seems to get the worse results:
length(which(bignetwork.nonGauss15000["pc.shd"]>10)%%8)
# We can see that in this case depends on the family: best t, worse clayton and joe
bignetwork.nonGauss15000[seq(1,8,1),"pc.shd"]
bignetwork.nonGauss15000[seq(9,16,1),"pc.shd"]
bignetwork.nonGauss15000[seq(17,24,1),"pc.shd"]
bignetwork.nonGauss15000[seq(25,32,1),"pc.shd"]
bignetwork.nonGauss15000[seq(33,40,1),"pc.shd"]
# Lets see low correlated vs high correlated
bignetwork.nonGauss15000[seq(1,60,8),"pc.shd"]
bignetwork.nonGauss15000[seq(2,60,8),"pc.shd"]
# Lets see about overestimation or underestimation:
bignetwork.nonGauss15000[seq(9,16,1),"pc.arcs"]
bignetwork.nonGauss15000[seq(9,16,1),"pc.ham"]
bignetwork.nonGauss15000["pc.arcs"]-bignetwork.nonGauss15000["pc.ham"]
# It always overestimate, but it always recovers the true edges!!!
# We can see from the difference vector!

sum(bignetwork.nonGauss15000["glasso.arcs"]-bignetwork.nonGauss15000["pc.arcs"])


# Finally we check the hill climbing for gaussian and non-gaussian data
# Gaussian case:
length(which(bignetwork.Gauss15000["hc.shd"]>0.5)%%8)
length(which(bignetwork.Gauss15000["hc.ham"]>0.5)%%8)
bignetwork.Gauss15000["hc.arcs"]-bignetwork.Gauss15000["hc.ham"]
# It always recovers the true network:
plot(graph.adjacency(C.bignetworkgauss[[10]],mode = "undirected"))
plot(graph.adjacency(C.bignetworkgauss[[19]],mode = "undirected"))
# Works pretty fine!



# Let's study what happens with the non-gaussian case:
length(which(bignetwork.nonGauss15000["hc.shd"]>10)%%8)
# Same as before clayton and joe copula worksreally bad
bignetwork.nonGauss15000[seq(1,8,1),"hc.shd"]
bignetwork.nonGauss15000[seq(9,16,1),"hc.shd"]
bignetwork.nonGauss15000[seq(17,24,1),"hc.shd"]
bignetwork.nonGauss15000[seq(33,40,1),"hc.shd"]
# Also way worse for the high correlated than the low correlated one.
bignetwork.nonGauss15000[seq(1,60,8),"hc.shd"]
bignetwork.nonGauss15000[seq(2,60,8),"hc.shd"]

# We can see that the true edges are always recovered, but there
# are a lot of extra edges that are addedd.
bignetwork.nonGauss15000["hc.arcs"]-bignetwork.nonGauss15000["hc.ham"]

plot(graph.adjacency(C.bignetworknongauss[[2]],mode = "undirected"))
plot(graph.adjacency(C.bignetworknongauss[[10]],mode = "undirected"))
plot(graph.adjacency(C.bignetworknongauss[[11]],mode = "undirected"))









