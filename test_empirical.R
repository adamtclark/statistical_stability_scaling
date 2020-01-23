error
rm(list=ls())
setwd("~/Dropbox/Projects/052_GFO_2019/src/simulation_study/")

#d<-read.csv("~/Dropbox/Projects/old/025_T2_Succession/src/data/e054mega16_soildat.csv")
d<-read.csv("~/Dropbox/Projects/old/025_T2_Succession/src/data/e014_wide_covars_161017.csv")
pltlst<-paste(d$field, d$transect, d$plot)
upltlst<-sort(unique(pltlst))

spps<-5:179
vcvout<-matrix(nrow=length(upltlst), ncol=3)
for(i in 1:length(upltlst)) {
  ds<-d[pltlst==upltlst[i],]
  spds<-ds[,spps]
  spds<-spds[,colSums(spds)>0,drop=FALSE]
  
  if(ncol(spds)>20 & nrow(spds)>5) {
    cvmat<-cov(spds)
    muv<-mean(diag(cvmat),na.rm=T)
    mucv<-mean(cvmat[row(cvmat)!=col(cvmat)],na.rm=T)
    
    vcvout[i,1]<-muv
    vcvout[i,2]<-mucv
    vcvout[i,3]<-ncol(spds)
  }
}
crrlst<-vcvout[,2]/vcvout[,1]
hist(crrlst)

pmin<-which.min(crrlst)
pmax<-which.max(crrlst)
pmid<-which.min(abs(crrlst))


ordr<-order(crrlst)

#run...
niter<-1000
n<-1
varestlst<-NULL
vlst<-rep(0, 3)
cvlst<-rep(0, 3)
for(k in c(pmin, pmid, pmax)) {
  psuse<-k
  ds<-d[pltlst==upltlst[psuse],]
  dim(pmin)
  spds<-ds[,spps]
  spds<-spds[,colSums(spds)>0]
  dim(spds)
  
  cvmat<-cov(spds)
  muv<-mean(diag(cvmat),na.rm=T)
  mucv<-mean(cvmat[row(cvmat)!=col(cvmat)],na.rm=T)
  
  vlst[n]<-muv
  cvlst[n]<-mucv
  
  N<-ncol(spds)
  varest<-matrix(nrow=niter, ncol=N)
  for(i in 1:niter) {
    for(j in 1:N) {
      smp<-sample(1:N,j)
      varest[i,j]<-var(rowSums(spds[,smp,drop=FALSE]))
    }
  }
  varestlst[[n]]<-varest
  n<-n+1
}
  
pdf("figures/empirical_var_example.pdf", width=3, height=6, colormodel = "cmyk", useDingbats = FALSE)
par(mfrow=c(3,1), mar=c(2,4.5,2,2), oma=c(2,0,0,0))
for(i in 1:3) {
  varest<-varestlst[[i]]
  N<-ncol(varest)
  mns<-colMeans(varest)
  sds<-apply(varest, 2, sd)
  matplot(1:N, mns+sds*matrix(nrow=N, rep(c(-1, 0, 1), each=N)), type="l", lty=c(2,1,2), col=1, xlab="", ylab=expression(paste("var(", Sigma, italic(N), ")")))
  
  #text(1, max(mns+sds)-diff(range(mns+sds))*0.1, paste("rho =", round(crrlst[c(pmin, pmid, pmax)[i]],3)), pos = 4)
  if(i==1) {
    text(1, max(mns+sds)-diff(range(mns+sds))*0.1, expression(paste(bar(italic(rho[x])), " = -0.037")), pos = 4)
  } else if(i==2) {
    text(1, max(mns+sds)-diff(range(mns+sds))*0.1, expression(paste(bar(italic(rho[x])), " = 0.000")), pos = 4)
  } else if(i==3) {
    text(1, max(mns+sds)-diff(range(mns+sds))*0.1, expression(paste(bar(italic(rho[x])), " = 0.095")), pos = 4)
  }
  
  nsp<-1:N
  #mod<-lm(mns~nsp+I(nsp^2))
  nsp_lng<-seq(1, N, length=1000)
  #prd<-predict(mod, newdata=data.frame(nsp=nsp_lng))
  #lines(nsp_lng, prd, col=2)
  
  lines(nsp_lng, vlst[i]*nsp_lng+cvlst[i]*(nsp_lng^2-nsp_lng), col=2, lwd=2, lty=3)
}
mtext("ecological grain", side = 1, line=2.5)
dev.off()
