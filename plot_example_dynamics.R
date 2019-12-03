error
rm(list=ls())
setwd("~/Dropbox/Projects/052_GFO_2019/src/simulation_study/")
source("functions.R")
source("~/Dropbox/Rfunctions/logit_funs.R")
require(mvtnorm)
require(nlme)
require(deSolve)
require(viridis)


#test diffusion approx
K<-1 #carrying capacity
r<-1.5 #rate of recovery
d<-(0.5) #mean size of disturbance
d_sd<-sqrt(0.1) #SD of disturbances
f<-1 #mean frequency of disturbance
sf<-0.01 #sampling frequency
tmax<-4.99 #maximum time

tm<-0.1
N<-16
M<-16

set.seed(191002)
tmsim<-symdyn(r, f, d, d_sd, sf=tm, tmax, stocht = FALSE, stochd = FALSE, oscillate_dist=TRUE)
vr<-mean((tmsim[,"state"])^2)

pdf("figures/example_dynamics.pdf", width=6, height=8, colormodel = "cmyk", useDingbats = FALSE)
par(mar=c(3,4.2,1,4.8), mfrow=c(2,1), oma=c(1,0,0,0.5))
plot(tmsim[,"time"], tmsim[,"state"], type="l", xlab="", ylab="", axes=FALSE, lwd=2)
title(main = "a.", outer = FALSE, line = -0.5, adj=0.03, cex.main=1.5)
axis(1, at=c(0:5)-0.1, c(0, rep("", 5))); axis(2, at=seq(-0.4, 0.4, length=5), c("", "", 0, "", ""))
mtext(expression(paste(time, ", ", italic(t))), 1, line = 2.5, cex=1.5)
mtext(expression(paste(italic(x), "(", italic(t), ")")), 2, line=2.5, cex=1.5)
abline(h=0, lty=3, lwd=2, col="gray50")
abline(h=c(sqrt(vr), -sqrt(vr)), lty=2, lwd=2, col="dodgerblue3")

text(3.95, -0.07, expression(paste(equilibrium)), pos = 4, xpd=NA, cex=1.5, col="gray50")

arrows(tmsim[,"time"][(which(tmsim[,"disturbed"]==1))-1], tmsim[,"state"][(which(tmsim[,"disturbed"]==1))-1],
       tmsim[,"time"][(which(tmsim[,"disturbed"]==1))-1], tmsim[,"state"][(which(tmsim[,"disturbed"]==1))],
       col=2, lwd=2, length = 0.1)

arrows(1:3+0.6-0.1, -0.05, 1:3-0.1+0.03, -0.05, col="darkgoldenrod3", lwd=2, length=0.1)
arrows(1:3+0.4-0.1, -0.05, 1:3-0.1+1-0.03, -0.05, col="darkgoldenrod3", lwd=2, length=0.1)

text(1.17, 0.11, expression(paste(frequency)), pos = 1, cex=1.5, col="darkgoldenrod3")

text(0.9, -0.48, expression(paste(resistance)), cex=1.5, pos=2, xpd=NA, col=2)

text(1.1, -0.42, expression(paste(resilience)), cex=1.5, pos=4, xpd=NA)

text(4.73, -0.29, expression(paste(variability)), cex=1.5, xpd=NA,
     col="dodgerblue3")
#dev.off()



tmsim2<-symdyn(r, f, d, d_sd, sf=0.1, tmax=19.99, stocht = FALSE, stochd = FALSE, oscillate_dist=TRUE)
x<-tmsim2[,"state"]
tm<-tmsim2[,"time"]
xstar<-sin(tm)*0.2
xstar2<-(1-sin(tm))*0.2-1.05

xobs<-x
xobs[tm<=11]<-x[tm<=11]+xstar[tm<=11]
xobs[tm>11]<-x[tm>11]+xstar2[tm>11]
xobs[tm>=11 & tm<=12]<-xstar2[tm>=11 & tm<=12]-0.1*exp(-seq(0, 5, length=11))

#pdf("figures/example_dynamics_complex.pdf", width=8, height=5, colormodel = "cmyk", useDingbats = FALSE)
par(mar=c(3,4.2,2,2.5))
plot(tm, xobs, type="l", xlab="", ylab="", axes=F, lwd=2, ylim=c(-1.5, 0.5))
title(main = "b.", outer = FALSE, line = -0.5, adj=0.03, cex.main=1.5)
axis(1, at=seq(0, 20, by=5), c(0, rep("", 4)))
axis(2, at=seq(-1.5, 0.5, by=0.5), c(0, rep("", 4)))
mtext(expression(paste(time, ", ", italic(t))), 1, line = 2.5, cex=1.5)
mtext(expression(paste(italic(N), "(", italic(t), ")")), 2, line=2.5, cex=1.5)

lines(tm, xstar, lty=3, lwd=2, col="gray50")
lines(tm, xstar2, lty=3, lwd=2, col="gray50")
text(19.95, xstar[length(xstar)], expression(paste(italic(N)[A],"*")), pos = 4, xpd=NA, cex=1.5, col="gray50")
text(19.95, xstar2[length(xstar2)], expression(paste(italic(N)[B],"*")), pos = 4, xpd=NA, cex=1.5, col="gray50")
dev.off()

