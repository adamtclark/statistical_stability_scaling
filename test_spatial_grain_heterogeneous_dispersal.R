#plots Figure 4 and Figure S1 in the main text

#load functions and packages
source("functions.R")
require(mvtnorm)
require(nlme)
require(deSolve)
collst_attributes<-c("firebrick3", "darkgoldenrod3", "dodgerblue2", "mediumpurple2", "olivedrab4", "springgreen4")
#resistance; resilience; invariability; frequency; persistence; robustness

set.seed(1202) #set random seed

K<-1 #carrying capacity
r<-1 #rate of recovery
d<-(0) #mean size of disturbance (mu in text)
d_sd<-sqrt(0.1) #SD of disturbances (sigma in text)
f<-1 #average time between disturbances (1/lambda in text)
sf<-0.1 #sampling interval
tmax<-120 #maximum time for simulation
d_cov<-d_cov0<-(d_sd)^2/2 #covariance in disturbance size among patches
rsd = 0.1 # standard deviation for r

##################################
#simulate varying spatial grain, asymmetrical dispersal
##################################
#set parameters
Alst<-1:30 #patch numbers to test
A<-c(2) #number of patches for dispersal test
Ilst<-seq(0, 2, by=0.1) #disperal rates to test
niter<-1000 #number of iterations

#only run if no saved simulation is available
if(sum(grep("estmat_sp_disp_hetero_210209.csv", dir("datout/")))==0) {
  set.seed(200209) #set random seed
  
  A<-c(2)
  Ilst<-seq(0, 2, by=0.2)
  d_cov<-0
  niter<-1000
  estmat_disp<-as.matrix(data.frame(iter=rep(niter, each=length(Ilst)), tsmp=rep(Ilst, niter),
                               n=NA, I=NA,
                               var=NA, f=NA, r_median=NA, d_sd_true=NA,
                               var_sp=NA, cor_sp=NA))
  n<-1
  for(i in 1:niter) {
    for(j in 1:length(Ilst)) {
      #simulate dynamics
      Iuse = Ilst[j]*c(0.5, 1.5)
      
      ruse = abs(rnorm(max(A), r, rsd))
      tmpout<-symdynN(r = ruse, amu=0, asd=0, f=f, d=d,
                    d_sd=d_sd, d_cov=d_cov, N=A,
                    sf=sf, tmax=tmax, stochd = TRUE, stocht = TRUE, as.matrix = TRUE, Ifrac=Iuse, dffun = df_col, fullout = TRUE)
      dtmp<-as.matrix(tmpout$datout)
      #remove first 20 steps for burn-in
      dtmp<-dtmp[dtmp[,"time"]>20,]
      
      estmat_disp[n,"I"]<-Ilst[j]
      estmat_disp[n,"n"]<-nrow(dtmp)
      
      #observation is summed abundance across the plots
      xsum<-rowSums(dtmp[,-c(1:2),drop=FALSE])
      
      #variance
      varest<-mean(xsum^2,na.rm=T)
      estmat_disp[n,"var"]<-varest
      #average time between disturbances
      estmat_disp[n,"f"]<-(max(dtmp[,"time"])-min(dtmp[,"time"]))/sum(dtmp[,"disturbed"][-1])
      
      #patch-level variance and between-patch covariance
      estmat_disp[n,"var_sp"]<-mean(apply(dtmp[,-c(1:2),drop=FALSE], 2, function(x) var(x)))
      cvtmp<-cov(dtmp[,-c(1:2),drop=FALSE])
      estmat_disp[n,"cor_sp"]<-mean(cvtmp[row(cvtmp)!=col(cvtmp)])
      
      # record true disturbance strength
      estmat_disp[n,"d_sd_true"]<-sd(rowSums(tmpout$dquant)) #achieved disturbance standard deviations
      
      dt<-mean(diff(dtmp[,"time"]))
      #estimate realized r
      disttm_forward<-which(c(dtmp[,"disturbed"][-1], 1)==0 & dtmp[,"disturbed"]==1)
      if(sum(disttm_forward)>0) {
        tmpr = suppressWarnings(log(xsum[disttm_forward+1]/xsum[disttm_forward])/dt)
        tmpr = tmpr[is.finite(tmpr)]
        estmat_disp[n,"r_median"]<-median(tmpr)
      }
      n<-n+1
    }
    if(i/10 == floor(i/10)) {
      cat(paste(round(i/niter,2), ";"))
    }
  }
  
  write.csv(estmat_disp, "datout/estmat_sp_disp_hetero_210209.csv", row.names=F)
} else {
  estmat_disp<-read.csv("datout/estmat_sp_disp_hetero_210209.csv")
}

pdf("figures/spatial_grain_dispersal_patch_level_hetro_210209.pdf", width=6, height=4, colormodel = "cmyk", useDingbats = FALSE)
  par(mfrow=c(2,2), mar=c(2,4,1,1), oma=c(2,1,0,0))
  
  #r
  ps<-is.finite(estmat_disp[,"r_median"])
  pltqt(estmat_disp[ps,"I"], -estmat_disp[ps,"r_median"], "", r, domod=FALSE, do_N = FALSE, plog = "", xlab = "dispersal rate", ylim = c(0.9, 1.12), qtp = c(-1, 1), jfac = 1, cluse = collst_attributes[2])
  mtext(expression(paste("meta. resilience")), 2, line=3.1)
  title("a.", line=-1, adj=0.04, cex.main=1.2)
  
  #d_sd
  sd_N<-sqrt((d_sd)^2*A) #recall, d_cov=0
  ps<-is.finite(estmat_disp[,"d_sd_true"])
  pltqt(estmat_disp[ps,"I"], estmat_disp[ps,"d_sd_true"], "", sd_N, domod=FALSE, do_N = FALSE, plog = "", xlab = "dispersal rate", ylim=c(0.38, 0.52), jfac = 1, cluse = collst_attributes[1])
  mtext(expression(paste("meta. resistance"^{-1})), 2, line=3.1)
  title("b.", line=-1, adj=0.04, cex.main=1.2)
  
  #variance
  var_N<-var_approx(r,f,d_sd)*A #recall, d_cov=0
  pltqt(estmat_disp[,"I"], estmat_disp[,"var"], "", truev = var_N, do_N = FALSE, domod=FALSE, plog = "", mlog="", xlab = "dispersal rate", ylim = c(0.07, 0.13), jfac = 1, cluse = collst_attributes[3])
  mtext(expression(paste("meta. invariance"^-1)), 2, line=3.1)
  title("c.", line=-1, adj=0.04, cex.main=1.2)
  
  
  pltqt(estmat_disp[,"I"], estmat_disp[,"var_sp"]*2, "", do_N = FALSE, domod=FALSE, plog = "", mlog="xy", xlab = "dispersal rate", modoffset = 1e-1, jfac = 1, ylim = c(-0.02, 0.14), cluse=collst_attributes[3])
  mtext(expression(paste("patch-level invariance"^-1)), 2, line=3.1)
  title("d.", line=-1, adj=0.02, cex.main=1.2)
  
  addqt(estmat_disp[,"I"], estmat_disp[,"cor_sp"]*(2^2-2), jfac = 1, cluse="darkblue", pltdens = NA, border=NA)
  abline(h=var_approx(r,f,d_sd), lwd=1, col="black", lty=3)
  
  text(1.5, 0.10, expression(paste(italic(A), paste("<var(", italic(x[i]), ")>"))), col=collst_attributes[3], cex=1.2)
  text(1.3, -0.03, expression(paste("(", italic(A)^2, "-", italic(A), ")", paste("<cov(", italic(x[i]), ",", italic(x[j]), ")>"))), col="darkblue", cex=1.2, pos = 3)
  
  mtext("dispersal rate", 1, line = 0.5, outer=TRUE)
dev.off()



##################################
#simulate varying spatial grain, dispersal in a non-closed system
##################################
#set parameters
Kuse = 1

#only run if no saved simulation is available
if(sum(grep("estmat_sp_disp_open_210209.csv", dir("datout/")))==0) {
  set.seed(200209) #set random seed
  
  A<-c(2)
  Ilst<-1
  d_cov<-0
  niter<-1000
  lossratelst = seq(0, 1, by =0.2) # rate for dispersal that leaves the simulated region
  estmat_disp<-as.matrix(data.frame(iter=rep(niter, each=length(lossratelst)), tsmp=rep(lossratelst, niter),
                                    n=NA, I=NA, lossrate = NA,
                                    meanx1 = NA, meanx2 = NA, meanxsum = NA,
                                    var=NA, distx = NA, f=NA, r_median=NA, d_sd_true=NA,
                                    var_sp=NA, cor_sp=NA))

  
  n<-1
  for(i in 1:niter) {
    for(j in 1:length(lossratelst)) {
      ruse = abs(rnorm(max(A), r, rsd))
      
      #simulate dynamics
      Iuse = Ilst
      lossrate = lossratelst[j]
      
      tmpout<-symdynN(r = ruse, amu=0, asd=0, f=f, d=d,
                      d_sd=d_sd, d_cov=d_cov, N=A,
                      sf=sf, tmax=tmax, stochd = TRUE, stocht = TRUE, as.matrix = TRUE, Ifrac=Iuse, Iloss = lossrate, dffun = df_col_loss, fullout = TRUE)
      dtmp<-as.matrix(tmpout$datout)
      #remove first 20 steps for burn-in
      dtmp<-dtmp[dtmp[,"time"]>20,]
      #matplot(dtmp[,1], dtmp[,-c(1:2)], type = "l")
      
      estmat_disp[n,"I"]<-Iuse
      estmat_disp[n,"lossrate"]<-lossrate
      estmat_disp[n,"n"]<-nrow(dtmp)
      
      #observation is summed abundance across the plots
      xsum<-rowSums(dtmp[,-c(1:2),drop=FALSE])
      
      #save abundances
      estmat_disp[n,"meanxsum"]<-mean(xsum)
      estmat_disp[n,"meanx1"]<-mean(dtmp[,3])
      estmat_disp[n,"meanx2"]<-mean(dtmp[,4])
 
      
      lossrate/(r-lossrate)
      estmat_disp[n,"meanx1"]
      estmat_disp[n,"meanx2"]
      
      #variance
      varest<-var(xsum,na.rm=T)
      estmat_disp[n,"var"]<-varest
      distx<-mean(xsum^2,na.rm=T)
      estmat_disp[n,"distx"]<-distx
      #average time between disturbances
      estmat_disp[n,"f"]<-(max(dtmp[,"time"])-min(dtmp[,"time"]))/sum(dtmp[,"disturbed"][-1])
      
      #patch-level variance and between-patch covariance
      estmat_disp[n,"var_sp"]<-mean(apply(dtmp[,-c(1:2),drop=FALSE], 2, function(x) var(x)))
      cvtmp<-cov(dtmp[,-c(1:2),drop=FALSE])
      estmat_disp[n,"cor_sp"]<-mean(cvtmp[row(cvtmp)!=col(cvtmp)])
      
      # record true disturbance strength
      estmat_disp[n,"d_sd_true"]<-sd(rowSums(tmpout$dquant)) #achieved disturbance standard deviations
      
      dt<-mean(diff(dtmp[,"time"]))
      #estimate realized r
      disttm_forward<-which(c(dtmp[,"disturbed"][-1], 1)==0 & dtmp[,"disturbed"]==1)
      if(sum(disttm_forward)>0) {
        tmpr = suppressWarnings(log(xsum[disttm_forward+1]/xsum[disttm_forward])/dt)
        tmpr = tmpr[is.finite(tmpr)]
        estmat_disp[n,"r_median"]<-median(tmpr)
      }
      n<-n+1
    }
    if(i/10 == floor(i/10)) {
      cat(paste(round(i/niter,2), ";"))
    }
  }
  
  write.csv(estmat_disp, "datout/estmat_sp_disp_open_210209.csv", row.names=F)
} else {
  estmat_disp<-read.csv("datout/estmat_sp_disp_open_210209.csv")
}

pdf("figures/spatial_grain_dispersal_patch_level_open_210209.pdf", width=6, height=4, colormodel = "cmyk", useDingbats = FALSE)
  par(mfrow=c(2,2), mar=c(2,4,1,1), oma=c(2,1,0,0))
  
  pltqt(estmat_disp[ps,"lossrate"], (estmat_disp[ps,"meanx1"]+estmat_disp[ps,"meanx2"])/2, "", r, domod=FALSE, do_N = FALSE, plog = "", xlab = "loss rate", ylim = c(-0.5, 0.05), qtp = c(-1, 1), jfac = 1, cluse = "darkgrey")
  mtext(expression(paste("patch-level x")), 2, line=3.1)
  title("a.", line=-1, adj=0.04, cex.main=1.2)
  xstar = r*Kuse/(r+lossratelst)-Kuse
  lines(lossratelst, xstar, lty = 2)
  
  
  #r
  ps<-is.finite(estmat_disp[,"r_median"])
  pltqt(estmat_disp[ps,"lossrate"], -estmat_disp[ps,"r_median"], "", r, domod=FALSE, do_N = FALSE, plog = "", xlab = "loss rate", ylim = c(-0.2, 1.12), qtp = c(-1, 1), jfac = 1, cluse = collst_attributes[2])
  mtext(expression(paste("meta. resilience")), 2, line=3.1)
  title("b.", line=-1, adj=0.04, cex.main=1.2)
  
  
  #variance
  var_N<-var_approx(r,f,d_sd)*A #recall, d_cov=0
  pltqt(estmat_disp[,"lossrate"], estmat_disp[,"var"], "", truev = var_N, do_N = FALSE, domod=FALSE, plog = "", mlog="", xlab = "dispersal rate", ylim = c(0.04, 0.12), jfac = 1, cluse = collst_attributes[3])
  mtext(expression(paste("meta. var(x)")), 2, line=3.1)
  title("c.", line=-1, adj=0.04, cex.main=1.2)
  
  
  pltqt(estmat_disp[,"lossrate"], estmat_disp[,"distx"], "", truev = var_N, do_N = FALSE, domod=FALSE, plog = "", mlog="", xlab = "dispersal rate", ylim = c(0.05, 1.2), jfac = 1, cluse = collst_attributes[3])
  mtext(expression(paste("meta. mean(x"^2, ")")), 2, line=3.1)
  title("d.", line=-1, adj=0.04, cex.main=1.2)
  mtext("loss rate", 1, line = 0.5, outer=TRUE)
dev.off()



