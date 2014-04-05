library(ggplot2)
setwd("/Users/VGupta/Documents/Research/Julia Stuff/Demand")
d1 = read.csv("fit1_Demand.csv", header=TRUE)
dk1 = read.csv("fit1_Kernel.csv", header=TRUE)
d1 = data.frame(d1, Kernel=dk1$Fit, LB=dk1$LB, UB=dk1$UB, Sample=dk1$Sample)

ggplot(aes(x=Prices, y=MargRev), data=d1) + 
  geom_line(size=.5, linetype="dashed") + 
  geom_point(aes(y=Fit), color="red", size=3 ) +
  geom_line(aes(y=Fit), color="red" ) +
  geom_ribbon(aes(ymin=LB, ymax=UB), fill="grey", alpha=.5)+
  geom_line(aes(y=Kernel), linetype="dotted", size=.5, color="black") +
  geom_point(aes(y=Kernel), color="black", fill="black", shape=15, size=3) +  
  theme_bw(base_size=18) + 
  ylab("Marginal Revenue") + xlab("Firm 1 Prices")


d2 = read.csv("fit2_Demand.csv", header=TRUE)
dk2 = read.csv("fit2_Kernel.csv", header=TRUE)
d2 = data.frame(d2, Kernel=dk2$Fit, LB=dk2$LB, UB=dk2$UB, Sample=dk2$Sample)


ggplot(aes(x=Prices, y=MargRev), data=d2) + 
  geom_line(size=.5, linetype="dashed") + 
  geom_point(aes(y=Fit), color="red", size=3 ) +
  geom_line(aes(y=Fit), color="red" ) +
  geom_ribbon(aes(ymin=LB, ymax=UB), fill="grey", alpha=.5)+
  geom_line(aes(y=Kernel), linetype="dotted", size=.5, color="black") +
  geom_point(aes(y=Kernel), color="black", fill="black", shape=15, size=3) +  
  theme_bw(base_size=18) + 
  ylab("Marginal Revenue") + xlab("Firm 2 Prices")


#### Straightup out-of-sample residuals
resids = read.csv("DemandResids.csv")
names(resids)
ggplot(aes(x=x2), data=resids) + 
  geom_histogram(aes(y=..density..), fill="grey") + 
  xlab("Approximation Error") + ylab("") +
  geom_vline(xintercept = 1.0630507131067428, linetype="dotted") +
  theme_bw(base_size=18)

mean(resids$x2 < 1.06305)


##### The idealized problem #####
dat1.ideal = read.csv("fit1_idealized.csv")
dat2.ideal = read.csv("fit2_idealized.csv")
colnames(dat2.ideal)<- colnames(dat1.ideal)

ggplot(aes(x=Prices, y=True), data=dat2.ideal) + 
  geom_line(size=.5, linetype="dashed") + 
  geom_point(aes(y=Fit), color="black", size=3, shape=15) +
  geom_line(aes(y=Fit), color="black", size=.5) +
  ylab("Marginal Revenue") + 
  geom_line(aes(y=Sample), linetype="dotted", size=.5, color="black") +
  geom_point(aes(y=Sample), color="blue", shape=17, size=3) +  
  geom_ribbon(aes(ymin=LB, ymax=UB), fill="grey", alpha=.3) +
  coord_cartesian(ylim=c(-10, 6)) + 
  xlab("Firm 2 Price") +
  theme_bw(base_size=18)





