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
r<-1 #rate of recovery
d<-(0) #mean size of disturbance
d_sd<-sqrt(0.1) #SD of disturbances
f<-1 #mean frequency of disturbance
sf<-0.01 #sampling frequency
tmax<-120 #maximum time


##### Run simulations
tm<-0.1
N<-16
M<-16

set.seed(191002)
tmsim<-symdyn(r, f, d, d_sd, sf=tm, tmax)
tmsim<-tmsim[tmsim$time>=20,]
ecsim<-symdynN(r, amu=(-0.1), asd=0, f, d, d_sd, d_cov=(d_sd^2)/2, N=N, sf=tm, tmax)
ecsim<-ecsim[ecsim$time>=20,]
spsim<-symdynN(r, amu=0, asd=0, f, d, d_sd, d_cov=(d_sd^2)/2, N=M, sf=tm, tmax, Ifrac = 0.1, dffun = df_col)
spsim<-spsim[spsim$time>=20,]

#transform to N
tmsim$state[tmsim$state<(-K)]<-(-K)


##### Plot examples
y0<-(-120)
dy<-19
dx<-0.4
x0<-(-(dx*4)/2)
collst1<-adjustcolor(c(magma(max(c(N, M)))), alpha.f = 0.8)
collst2<-adjustcolor(c(viridis(max(c(N, M)))), alpha.f = 0.8)
col0<-adjustcolor(1, alpha.f = 0.8)
colmd<-4
ylm<-c(-1.3, 1.3)
yscl<-diff(range(ylm))/127

cifun<-function(r, x0, y0, yscl=1) {
  pisq<-seq(0, 2*pi, length=200)
  xsq<-cos(pisq)*r+x0
  ysq<-sin(pisq)*r*yscl+y0
  return(data.frame(x=xsq, y=ysq))
}

radsp<-dy*0.3
x0sp<-c(y0+(4*dy)*(c(0.48, 0.3, 0.5, 0.3, 0.44, 0.35, 0.54, 0.5, 0.7, 0.87, 0.73, 0.98, 0.72, 0.71, 0.9, 0.93)-0.09))
y0sp<-c(x0+(4*dx)*c(0.17, 0.22, 0.36, 0.42, 0.56, 0.75, 0.72, 0.95, 0.68, 0.75, 0.87, 0.57, 0.44, 0.22, 0.4, 0.23))


pdf("figures/trajectories.pdf", width=10, height=7.5, colormodel = "cmyk", useDingbats = FALSE)
  m<-t(matrix(nrow=3, data=1:9))
  layout(m)
  par(mar=c(2,12,2,2), oma=c(2,1,1.2,0))
  
  #singles
  plot(tmsim$time-20, tmsim$state, xlab="", ylab="", type="l", ylim=ylm, col=col0); abline(h=0, lty=3)
  title("a.", line=-1.2, adj=0.04, cex.main=1.4)
  mtext("Temporal Scale", side = 3, line=1, cex=1.2)
  mtext(expression(paste(italic(x))), side = 2, line=2.2)
  citmp<-cifun(r=2*dy, x0 = y0+dy*2, y0 = 0, yscl = yscl)
  lines(citmp, xpd=NA)
  citmp2<-cifun(r=2*dy*0.9, x0 = y0+dy*2, y0 = 0, yscl = yscl)
  nps<-round(seq(1, nrow(citmp2), length=13))[-13]
  text(x=citmp2$x[nps], y=citmp2$y[nps], c(3:1, 12:4), xpd=NA)
  arrows(y0+dy*2, 0, citmp2$x[nps], citmp2$y[nps], xpd=NA, length = 0.1, col="darkred", lwd=2, lend=4)
  
  matplot(ecsim$time-20, ecsim[,-c(1:2)], xlab="", ylab="", type="l", ylim=ylm, col=collst1, lty=1); abline(h=0, lty=3)
  mtext("Spatial Scale", side = 3, line=1, cex=1.2)
  mtext(expression(paste(italic(x))), side = 2, line=2.2)
  title("d.", line=-1.2, adj=0.04, cex.main=1.4)
  
  for(i in 1:4) {
    for(j in 1:4) {
      polygon(c(y0, y0+dy, y0+dy, y0)+dy*(i-1), c(x0, x0,x0+dx, x0+dx)+dx*(j-1), col=collst1[c(1,2,5,6,3,4,7,8,9,10,13,14,11,12,15,16)[i+4*(j-1)]], xpd=NA)
    }
  }
  
  matplot(spsim$time-20, spsim[,-c(1:2)], xlab="", ylab="", type="l", ylim=ylm, col=collst2, lty=1); abline(h=0, lty=3)
  mtext("Ecological Scale", side = 3, line=1, cex=1.2)
  mtext(expression(paste(italic(x))), side = 2, line=2.2)
  
  for(i in 1:N) {
    citmp3<-cifun(radsp, x0sp[i], y0sp[i], yscl = yscl)
    polygon(citmp3$x, citmp3$y, xpd=NA, col = collst2[i])
  }
  title("g.", line=-1.2, adj=0.04, cex.main=1.4)
  
  #halfs
  ps<-tmsim$time%in%seq(20, 120, by=1)
  plot(tmsim$time[ps]-20, tmsim$state[ps], xlab="", ylab="", type="l", ylim=ylm, col=col0); abline(h=0, lty=3)
  mtext(expression(paste(italic(x))), side = 2, line=2.2)
  citmp<-cifun(r=2*dy, x0 = y0+dy*2, y0 = 0, yscl = yscl)
  lines(citmp, xpd=NA)
  citmp2<-cifun(r=2*dy*0.9, x0 = y0+dy*2, y0 = 0, yscl = yscl)
  nps<-round(seq(1, nrow(citmp2), length=13))[-13]
  text(x=citmp2$x[nps], y=citmp2$y[nps], c(3:1, 12:4), xpd=NA)
  arrows(y0+dy*2, 0, citmp2$x[nps][c(1,4,7,10)], citmp2$y[nps][c(1,4,7,10)], xpd=NA, length = 0.1, col="darkred", lwd=2, lend=4)
  title("b.", line=-1.2, adj=0.04, cex.main=1.4)
  
  matplot(spsim$time-20, cbind(rowMeans(spsim[,1:4+2]), rowMeans(spsim[,5:8+2]), rowMeans(spsim[,9:12+2]), rowMeans(spsim[,13:16+2])),
          xlab="", ylab="", type="l", ylim=ylm, col=collst1[c(1,4,8,16)], lty=1); abline(h=0, lty=3)
  mtext(expression(paste(Sigma, italic(x), "/", italic(A))), side = 2, line=2.2)
  for(i in 1:2) {
    for(j in 1:2) {
      polygon(c(y0, y0+dy*2, y0+dy*2, y0)+dy*2*(i-1), c(x0, x0,x0+dx*2, x0+dx*2)+dx*2*(j-1), col=collst1[c(1,4,8,16)[i+2*(j-1)]], xpd=NA)
    }
  }
  title("e.", line=-1.2, adj=0.04, cex.main=1.4)
  
  matplot(ecsim$time-20, cbind(rowMeans(ecsim[,1:4+2]), rowMeans(ecsim[,5:8+2]), rowMeans(ecsim[,9:12+2]), rowMeans(ecsim[,13:16+2])),
          xlab="", ylab="", type="l", ylim=ylm, col=collst2[c(1,4,8,16)], lty=1); abline(h=0, lty=3)
  mtext(expression(paste(Sigma, italic(x), "/", italic(M))), side = 2, line=2.2)
  title("h.", line=-1.2, adj=0.04, cex.main=1.4)
  
  for(i in 1:N) {
    citmp3<-cifun(radsp, x0sp[i], y0sp[i], yscl = yscl)
    polygon(citmp3$x, citmp3$y, xpd=NA, col = adjustcolor(collst2[i], alpha.f = 0.2))
  }
  for(i in 1:4) {
    citmp3<-cifun(radsp*4, mean(x0sp[1:4+4*(i-1)]), mean(y0sp[1:4+4*(i-1)]), yscl = yscl)
    polygon(citmp3$x, citmp3$y, xpd=NA, col = collst2[i+4*(i-1)])
  }
  
  #fulls
  ps<-tmsim$time%in%seq(20, 120, by=5)
  plot(tmsim$time[ps]-20, tmsim$state[ps], xlab="", ylab="", type="l", ylim=ylm, col=col0); abline(h=0, lty=3)
  mtext(expression(paste(italic(x))), side = 2, line=2.2)
  mtext("time", side = 1, line=2.6)
  title("c.", line=-1.2, adj=0.04, cex.main=1.4)
  
  citmp<-cifun(r=2*dy, x0 = y0+dy*2, y0 = 0, yscl = yscl)
  lines(citmp, xpd=NA)
  citmp2<-cifun(r=2*dy*0.9, x0 = y0+dy*2, y0 = 0, yscl = yscl)
  nps<-round(seq(1, nrow(citmp2), length=13))[-13]
  text(x=citmp2$x[nps], y=citmp2$y[nps], c(3:1, 12:4), xpd=NA)
  arrows(y0+dy*2, 0, citmp2$x[nps][c(4,10)], citmp2$y[nps][c(4,10)], xpd=NA, length = 0.1, col="darkred", lwd=2, lend=4)
  
  matplot(spsim$time-20, cbind(rowMeans(spsim[,1:16+2])), xlab="", ylab="", type="l", ylim=ylm, col=collst1[colmd], lty=1); abline(h=0, lty=3)
  mtext(expression(paste(Sigma, italic(x), "/", italic(A))), side = 2, line=2.2)
  mtext("time", side = 1, line=2.6)
  title("f.", line=-1.2, adj=0.04, cex.main=1.4)
  polygon(c(y0, y0+dy*4, y0+dy*4, y0), c(x0, x0,x0+dx*4, x0+dx*4), col=collst1[colmd], xpd=NA)
  
  matplot(ecsim$time-20, cbind(rowMeans(ecsim[,1:16+2])), xlab="", ylab="", type="l", ylim=ylm, col=collst2[colmd], lty=1); abline(h=0, lty=3)
  mtext(expression(paste(Sigma, italic(x), "/", italic(M))), side = 2, line=2.2)
  mtext("time", side = 1, line=2.6)
  
  for(i in 1:N) {
    citmp3<-cifun(radsp, x0sp[i], y0sp[i], yscl = yscl)
    polygon(citmp3$x, citmp3$y, xpd=NA, col = adjustcolor(collst2[i], alpha.f = 0.2))
  }
  citmp3<-cifun(radsp*7.3, mean(x0sp)-dy*0.15, mean(y0sp), yscl = yscl)
  polygon(citmp3$x, citmp3$y, xpd=NA, col = adjustcolor(collst2[colmd], alpha.f = 0.8))
  title("i.", line=-1.2, adj=0.04, cex.main=1.4)
dev.off()


