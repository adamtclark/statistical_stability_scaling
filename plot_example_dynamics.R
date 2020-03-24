#plots Figure 1 in the main text
#rm(list=ls())

#load functions and packages
source("functions.R")
require(mvtnorm)
require(nlme)
require(deSolve)
collst_attributes<-c("firebrick3", "darkgoldenrod3", "dodgerblue2", "mediumpurple2", "olivedrab4", "springgreen4")
#resistance; resilience; invariability; frequency; persistence; robustness

#parameters
K<-1 #carrying capacity
r<-1.5 #rate of recovery
d<-(0.5) #mean size of disturbance (mu in text)
d_sd<-sqrt(0.1) #SD of disturbances (sigma in text)
f<-1 #average time between disturbances (1/lambda in text)
tmax<-4.99 #maximum time for simulation
tm<-0.1 #sampling interval

set.seed(191002) #set random seed

pdf("figures/example_dynamics.pdf", width=6, height=8, colormodel = "cmyk", useDingbats = FALSE)
  #simulate simple dynamics
  tmsim<-symdyn(r, f, d, d_sd, sf=tm, tmax, stocht = FALSE, stochd = FALSE, oscillate_dist=TRUE)
  vr<-mean((tmsim[,"state"])^2) #variability (since E[x]=0)

  par(mar=c(3,4.2,1,4.8), mfrow=c(2,1), oma=c(1,0,0,0.5))
  plot(tmsim[,"time"], tmsim[,"state"], type="l", xlab="", ylab="", axes=FALSE, lwd=2, col=adjustcolor(1, 0.8))
  title(main = "a.", outer = FALSE, line = -0.5, adj=0.03, cex.main=1.5)
  axis(1, at=c(0:5)-0.1, c(0, rep("", 5))); axis(2, at=seq(-0.4, 0.4, length=5), c("", "", 0, "", ""))
  mtext(expression(paste(time, ", ", italic(t))), 1, line = 2.5, cex=1.5)
  mtext(expression(paste(italic(x), "(", italic(t), ")")), 2, line=2.5, cex=1.5)
  abline(h=0, lty=3, lwd=2, col="gray50")
  text(3.95, -0.07, expression(paste(equilibrium)), pos = 4, xpd=NA, cex=1.5, col="gray50")
  
  abline(h=c(sqrt(vr), -sqrt(vr)), lty=2, lwd=2, col=collst_attributes[3])
  arrows(tmsim[,"time"][(which(tmsim[,"disturbed"]==1))-1], tmsim[,"state"][(which(tmsim[,"disturbed"]==1))-1],
         tmsim[,"time"][(which(tmsim[,"disturbed"]==1))-1], tmsim[,"state"][(which(tmsim[,"disturbed"]==1))],
         col=collst_attributes[1], lwd=2, length = 0.1)
  arrows(1:3+0.6-0.1, -0.05, 1:3-0.1+0.03, -0.05, col=collst_attributes[4], lwd=2, length=0.1)
  arrows(1:3+0.4-0.1, -0.05, 1:3-0.1+1-0.03, -0.05, col=collst_attributes[4], lwd=2, length=0.1)
  
  for(i in 2:5) {
    ps<-which(tmsim[,"time"]<=(i-0.15) & tmsim[,"time"]>=(i-1))
    lines(tmsim[ps,"time"]+0.15, tmsim[ps,"state"], lwd=2, col=collst_attributes[2])
    arrows(tmsim[max(ps)-1,"time"]+0.15, tmsim[max(ps)-1,"state"],
           tmsim[max(ps),"time"]+0.15, tmsim[max(ps),"state"], lwd=2, col=collst_attributes[2], length=0.1)
  }
  
  text(1.17, 0.11, expression(paste(frequency)), pos = 1, cex=1.5, col=collst_attributes[4])
  text(0.9, -0.48, expression(paste(resistance^-1)), cex=1.5, pos=2, xpd=NA, col=collst_attributes[1])
  text(1.2, -0.42, expression(paste(resilience)), cex=1.5, pos=4, xpd=NA, col=collst_attributes[2])
  text(4.73, -0.29, expression(paste(invariability^-1)), cex=1.5, xpd=NA,
       col=collst_attributes[3])
  
  #example with non-static dynamics
  tmsim2<-symdyn(r, f, d, d_sd, sf=0.1, tmax=18.99, stocht = FALSE, stochd = FALSE, oscillate_dist=TRUE)
  x<-tmsim2[,"state"]
  tm<-tmsim2[,"time"]
  xstar<-rep(0, length(tm))
  xstar2<-rep(-0.85, length(tm)) #new mean value for x, after transition
  
  #transform simulated x to show transition
  xobs<-x
  transt<-7 #time step for transition
  xobs[tm<=transt]<-x[tm<=transt]+xstar[tm<=transt]
  xobs[tm>transt]<-x[tm>transt]+xstar2[tm>transt]
  xobs[tm>=transt & tm<=(transt+1)]<-xstar2[tm>=transt & tm<=(transt+1)]-0.1*exp(-seq(0, 5, length=11))
  xobs[tm>(2*transt+1)]<-(-1.5)
  
  rmps<-c(which(tm%in%c(transt, (2*transt+1))))
  rmps<-c(rmps, rmps+1)
  plot(tm[-rmps], xobs[-rmps], type="l", xlab="", ylab="", axes=F, lwd=2, ylim=c(-1.5, 0.5), col=adjustcolor(1, 0.8))
  title(main = "b.", outer = FALSE, line = -0.5, adj=0.03, cex.main=1.5)
  axis(1, at=seq(min(tm), max(tm), by=f), c(0, rep("", 18)))
  axis(2, at=seq(-1.5, 0.5, by=0.5), c(0, rep("", 4)))
  mtext(expression(paste(time, ", ", italic(t))), 1, line = 2.5, cex=1.5)
  mtext(expression(paste(italic(N), "(", italic(t), ")")), 2, line=2.5, cex=1.5)
  
  lines(tm, xstar, lty=3, lwd=2, col="gray50")
  lines(tm, xstar2, lty=3, lwd=2, col="gray50")
  lines(c(tm[1], tm[length(tm)]), c(-1.5, -1.5), lty=3, lwd=2, col="grey50")
  text(max(tm), xstar[length(xstar)], expression(paste(italic(N)[A],"*")), pos = 4, xpd=NA, cex=1.5, col="gray50")
  text(max(tm), xstar2[length(xstar2)], expression(paste(italic(N)[B],"*")), pos = 4, xpd=NA, cex=1.5, col="gray50")
  text(max(tm), -1.5, expression(paste(italic(N)[C],"*")), pos = 4, xpd=NA, cex=1.5, col="gray50")
  
  arrows(transt-0.1, xobs[which(tm==transt)-1], transt-0.1, xobs[which(tm==transt)+2], col=collst_attributes[6], lwd=2, length=0.1)
  arrows(transt*2+1-0.1, xobs[which(tm==(2*transt+1))-1], transt*2+1-0.1, -1.5, col=collst_attributes[5], lwd=2, length=0.1)
  
  text(transt+0.1, -0.75, expression(paste(robustness^-1)), pos = 2, cex=1.5, col=collst_attributes[6])
  text(2*transt+1, -1.4, expression(paste(persistence^-1)), cex=1.5, pos=2, xpd=NA, col=collst_attributes[5])
dev.off()

