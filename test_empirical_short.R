#rm(list=ls())
collst_attributes<-c("firebrick3", "darkgoldenrod3", "dodgerblue2", "mediumpurple2", "olivedrab4", "springgreen4")
#resistance; resilience; invariability; frequency; persistence; robustness

#load data
d<-read.csv("old_field_empirical.csv")
pltlst<-paste(d$field, d$transect, d$plot) #list of plots
upltlst<-unique(pltlst) #unique plots

niter<-1000 #number of iterations
varestlst<-NULL
vlst<-rep(0, 3) #store variance
cvlst<-rep(0, 3) #store covariance
spps<-5:179 #columns with species abundance data
n<-1

set.seed(200324) #set random seed

if(sum(grep("emirical_out.rda", dir("datout/")))==0) {
  #run sampling routine
  for(k in 1:length(upltlst)) {
    ds<-d[pltlst==upltlst[k],]
    spds<-ds[,spps]
    spds<-spds[,colSums(spds)>0]
    
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
  save(list = c("varestlst", "vlst", "cvlst"), file = "datout/emirical_out.rda")
} else {
  load("datout/emirical_out.rda")
}

pdf("figures/empirical_var_example.pdf", width=3, height=6, colormodel = "cmyk", useDingbats = FALSE)
  par(mfrow=c(3,1), mar=c(2,4.5,2,2), oma=c(2,0,0,0))
  for(i in c(3,1,2)) {
    varest<-varestlst[[i]]
    N<-ncol(varest)
    mns<-colMeans((varest))
    sds<-apply((varest), 2, sd)
    
    #plot actual relationship
    matplot(1:N, (mns+sds*matrix(nrow=N, rep(c(-1, 0, 1), each=N))), lty=c(2,1,2), col=1, xlab="", ylab=expression(paste("var(", Sigma, italic(N), ")")), type="n")
    polygon(c(1:N, N:1), pmax((c(mns+sds, rev(mns-sds))), 0), border=NA, col=adjustcolor(collst_attributes[3], alpha.f = 0.5))
    
    #text(1, max(mns+sds)-diff(range(mns+sds))*0.1, paste("rho =", round(crrlst[c(pmin, pmid, pmax)[i]],3)), pos = 4)
    if(i==3) {
      text(1, max((mns+sds))-diff((range(mns+sds)))*0.1, expression(paste(bar(italic(rho[x])), " = -0.037")), pos = 4)
    } else if(i==1) {
      text(1, max((mns+sds))-diff((range(mns+sds)))*0.1, expression(paste(bar(italic(rho[x])), " = 0.000")), pos = 4)
    } else if(i==2) {
      text(1, max((mns+sds))-diff((range(mns+sds)))*0.1, expression(paste(bar(italic(rho[x])), " = 0.095")), pos = 4)
    }
    
    #analytial scaling relationship
    nsp<-1:N
    nsp_lng<-seq(1, N, length=1000)
    lines(nsp_lng, vlst[i]*nsp_lng+cvlst[i]*(nsp_lng^2-nsp_lng), col=1, lwd=1.5, lty=2)
  }
  mtext("ecological scale", side = 1, line=2.5)
dev.off()
