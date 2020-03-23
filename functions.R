xt <- function(t0, t1, B0, r) {
  dt<-(t1-t0)
  exp(-r*dt)*B0
}

#1 species
symdyn<-function(r, f, d, d_sd, sf, tmax, stochd=TRUE, stocht=TRUE, as.matrix=FALSE, oscillate_dist=FALSE) {
  # r is intrinsic growth rate (i.e. resilience)
  # f is the waiting time (or average waiting time) between disturbance events
  # d is the mean size of disturbances
  # d_sd is the standard deviation used to generate disturbances
  # sf is the waiting time between sampling events
  # tmax is the time series length to be simulated
  # stochd is a logical variable, indicating whether disturbance size should be stochastic - otherwise, all disturbances are of magnitude d
  # stocht is a logical variable, indicating whether waiting time between disturbance events should be stochastic - otherwise, waiting time is always f
  # as.matrix indicates whether results should be returned as matrix (potentially faster for some applications)
  # oscillate_dist is a logical variable indicating whether the sign of the disturbance should oscillate between positive and negative - ignored if stochd==TRUE
  
  st<-seq(0, tmax, by=sf)
  nobs<-length(st)
  
  datout<-matrix(nrow=length(st), ncol=3)
  colnames(datout)<-c("time", "state", "disturbed")
  datout[,"time"]<-st
  datout[1,"state"]<-0
  datout[,"disturbed"]<-0
  
  x<-0 #standardized abundance
  tm<-0 #time
  n<-2 #sample position
  md<-1 #number of disturbances
  
  while(n <= nobs) {
    if(n==2) {
      if(stocht) {
        tdist<-tm+rexp(1, 1/f) #time to next disturbance
      } else {
        tdist<-tm+f
      }
      tsamp<-st[n]
    }
    
    while((n <= nobs) & (tsamp<tdist)) {
      x<-xt(tm, tsamp, x, r)
      datout[n,"state"]<-x
      tm<-tsamp
      n<-n+1
      if(n <= length(st)) {
        tsamp<-st[n]
      }
    }
    
    while((n <= nobs) & (tdist<=tsamp)) {
      x<-xt(tm, tdist, x, r)
      tm<-tdist
      if(stochd) {
        rd<-rnorm(1, d, d_sd)
      } else {
        if(!oscillate_dist | (md%%2 == 0)) {
          rd<-d
        } else {
          rd<-(-d)
        }
        md<-md+1
      }
      x<-x+rd
      datout[n,"disturbed"]<-datout[n,"disturbed"]+1
      
      if(stocht) {
        tdist<-tm+rexp(1, 1/f) #time to next disturbance
      } else {
        tdist<-tm+f
      }
    }
  }
  
  if(!as.matrix) {
    data.frame(datout)
  } else {
    datout
  }
}

df<-function(time, state, pars) {
  list(c(t(as.matrix(state))%*%pars$A))
}

df_col<-function(time, state, pars) {
  list(diag(pars$A)*state-pmax(pars$Ifrac*(state+1),0)+pars$Ifrac*mean(pmax(state+1, 0)))
  #list(diag(pars$A)*state-pars$Ifrac*state+pars$Ifrac*mean(state))
}

xtN <- function(t0, t1, B0, odepars, dffun, nsteps=2) {
  out<-ode(y=B0, times=seq(t0, t1, length=nsteps), parms=odepars, func = dffun)
  out[-1,-1]
}


#n species or patches
symdynN<-function(r, amu, asd, f, d, d_sd, d_cov, N, sf, tmax, stochd=TRUE, stocht=TRUE, as.matrix=FALSE, amax=0, amin=-Inf, Ifrac=NULL, dffun=df, fullout=FALSE, xstart=NULL) {
  # r is intrinsic growth rate (i.e. resilience)
  # amu is the mean interaction strength
  # asd is the standard deviation used to generate interaction strengths
  # f is the waiting time (or average waiting time) between disturbance events
  # d is the mean size of disturbances
  # d_sd is the standard deviation used to generate disturbances
  # d_cov is the covariance for generating disturbances
  # N is the number of species
  # sf is the waiting time between sampling events
  # tmax is the time series length to be simulated
  # stochd is a logical variable, indicating whether disturbance size should be stochastic - otherwise, all disturbances are of magnitude d
  # stocht is a logical variable, indicating whether waiting time between disturbance events should be stochastic - otherwise, waiting time is always f
  # as.matrix indicates whether results should be returned as matrix (potentially faster for some applications)
  # oscillate_dist is a logical variable indicating whether the sign of the disturbance should oscillate between positive and negative - ignored if stochd==TRUE
  # amax is the maximum value for interaction coefficiets
  # amin is the minimum value for interaction coefficiets
  # Ifrac is the dispersal rate
  # dffun is the function handed to the ODE solver
  # fullout is a logical, determining whether the full output or just a summary is returned
  # xstart is an optional vector of starting abundances
  
  st<-seq(0, tmax, by=sf)
  nobs<-length(st)
  
  datout<-matrix(nrow=length(st), ncol=2+N)
  colnames(datout)<-c("time", "disturbed", paste("N", 1:N, sep="_"))
  datout[,"time"]<-st
  datout[,-1]<-0
  sppos<-1:N+2
  
  if(is.null(xstart)) {
    x<-rep(0, N) #standardized abundance
  } else {
    x<-xstart
    datout[1,-c(1:2)]<-x
  }
  tm<-0 #time
  n<-2 #sample position
  m<-1 #disturbance position
  
  #disturbance times
  if(stocht) {
    dtime<-cumsum(rexp(round(2*(tmax/f)), 1/f))
    mm<-3
    while(max(dtime)<tmax) {
      dtime<-cumsum(rexp(round(mm*(tmax/f)), 1/f))
      mm<-mm+1
    }
    dtime<-dtime[dtime<=tmax]
  } else {
    dtime<-seq(f, tmax, by=f)
  }
  
  #disturbance quantities
  if(stochd) {
    covmat<-diag(N)*d_sd^2
    covmat[row(covmat)!=col(covmat)]<-d_cov
    dquant<-rmvnorm(length(dtime), mean=rep(d, N), sigma = covmat)
  } else {
    dquant<-matrix(nrow=length(dtime), ncol=N, data=d)
  }
  
  #interaction matrix
  A<-(-diag(N)*r)
  ps<-which(row(A)!=col(A))
  A[ps]<-rnorm(N^2-N, amu, asd)
  while(any(A[ps]<amin) | any(A[ps]>amax)) {
    A[ps][A[ps]<amin]<-rnorm(sum(A[ps]<amin), amu, asd)
    A[ps][A[ps]>amax]<-rnorm(sum(A[ps]>amax), amu, asd)
  }
  
  odepars<-list(A=A, Ifrac=Ifrac)
  
  while(n <= nobs) {
    if(n==2) {
      tdist<-dtime[m] #time to next disturbance
      tsamp<-st[n]
    }
    
    while((n <= nobs) & (tsamp<tdist)) {
      x<-xtN(tm, tsamp, x, odepars, dffun)
      datout[n,sppos]<-x
      tm<-tsamp
      n<-n+1
      if(n <= length(st)) {
        tsamp<-st[n]
      }
    }
    
    while((n <= nobs) & (tdist<=tsamp)) {
      x<-xtN(tm, tdist, x, odepars, dffun)
      tm<-tdist
      
      x<-x+dquant[m,]
      datout[n,"disturbed"]<-datout[n,"disturbed"]+1
      
      m<-m+1
      if(m<=length(dtime)) {
        tdist<-dtime[m]
      } else {
        tdist<-tmax+1
      }
    }
  }
  
  if(fullout) {
    return(list(datout=datout, A=A, dquant=dquant))
  } else if(!as.matrix) {
    data.frame(datout)
  } else {
    datout
  }
}



#mean state immediately before next disturbance
mu_dist_fun<-function(r, K, f, d) {
  K*d*exp(-r*f)/(1-exp(-r*f))
}

#mean state
mu_state<-function(r, K, f, d) {
  (K*d)/(f*r)
}

#variance state
var_state<-function(r, K, f, d, d_sd) {
  ((K*d)/(f*r)+d_sd^2)/(2*r*f)
}

var_approx<-function(r, f, d_sd) {
  (d_sd^2)/(2*r*f)
}


#inverse functions
inv_var<-function(var, r=NULL, f=NULL, d_sd=NULL) {
  ## var = (d_sd^2)/(2*r*f)
  
  if(is.null(r)) {
    return(d_sd^2/(2*f*var))
  } else if(is.null(f)) {
    return(d_sd^2/(2*r*var))
  } else if(is.null(d_sd)) {
    return(sqrt(var*(2*r*f)))
  }
}


pltqt<-function(tmp, x, ylab="", truev=NULL, plog="y", mlog="", domod=TRUE, do_N=TRUE, plotqtl=c(0,1), modoffset=0, xlab="temporal grain", ylim=NULL, qtp=c(-1, 1), jfac=10, cluse="black", linecol="dodgerblue", ...) {
  if(var(x, na.rm=T)<1e-10) {
    #remove floating point error for plotting
    x<-round(x, 10)
  }
  
  #op<-par(no.readonly = TRUE)
  #par(mar=c(4,4,2,4))
  if(is.null(ylim)) {
    pylim<-quantile(x, plotqtl, na.rm=T)
  } else {
    pylim<-ylim
  }
  plot(jitter(tmp, factor = jfac), x, xlab=xlab, ylab=ylab, col=adjustcolor(cluse, alpha.f = 0.2), cex=0.3, log=plog, ylim=pylim, pch=16, axes=F, ...)
  axis(1); axis(2, las=2); box()
  
  qtl<-t(matrix(nrow=2, unlist(tapply(x, tmp, function(x) quantile(x, pnorm(c(qtp)),na.rm=T)))))
  tlst<-sort(unique(tmp))
  nobs<-tapply(x, tmp, function(x) sum(!is.na(x)))
  ps<-which(!is.na(rowSums(qtl)) & nobs>3)
  
  polygon(c(tlst[ps], rev(tlst[ps])), c(qtl[ps,1], rev(qtl[ps,2])), col=adjustcolor(cluse, alpha.f = 0.5), border = NA)
  if(!is.null(truev)) {
    if(length(truev)==1) {
      abline(h=truev, lty=2, col=linecol, lwd=2)
    } else {
      lines(truev, lty=2, col=linecol, lwd=2)
    }
  }
  
  if(domod) {
    if(sum(grep("y", mlog))>0) {
      xmod<-log(x)
    } else {
      xmod<-x
    }
    if(sum(grep("x", mlog))>0) {
      tmpmod<-log(tmp+modoffset)
      tlstmod<-log(tlst+modoffset)
    } else {
      tmpmod<-tmp+modoffset
      tlstmod<-tlst+modoffset
    }
    
    mod<-try(gls(xmod~tmpmod, weights=varPower(form=~fitted(.)), na.action = na.exclude), silent = TRUE)
    if(is.character(mod)) {
      wt_tmp<-1/tapply(xmod, tmpmod, var)
      wts<-wt_tmp[match(tmpmod, tlstmod)]
      mod<-lm(xmod~tmpmod, weights=wts)
    }
    
    prd<-predict(mod, newdata=data.frame(tmpmod=tlstmod))
    if(sum(grep("y", mlog))>0) {
      prd<-exp(prd)
    }
    lines(tlst, prd, col=2, lty=3, lwd=2)
  }
  
  if(do_N) {
    par(new=TRUE)
    nobs<-tapply(x, tmp, function(x) sum(!is.na(x))/length(x))
    plot(tlst, nobs, lty=2, lwd=2, col=1, type="l", axes="F", xlab="", ylab="")
    axis(4)
    mtext("frac. simulations", 4, line=2.2)
    abline(h=c(0, 1), lty=3, col=1)
  }
  #par(op)
}

addqt<-function(tmp, x, qtp=c(-1, 1), jfac=10, cluse="red") {
  if(var(x, na.rm=T)<1e-10) {
    #remove floating point error for plotting
    x<-round(x, 10)
  }
  
  points(jitter(tmp, factor = jfac), x, col=adjustcolor(cluse, alpha.f = 0.2), cex=0.3, pch=16)
  
  qtl<-t(matrix(nrow=2, unlist(tapply(x, tmp, function(x) quantile(x, pnorm(c(qtp)),na.rm=T)))))
  tlst<-sort(unique(tmp))
  nobs<-tapply(x, tmp, function(x) sum(!is.na(x)))
  ps<-which(!is.na(rowSums(qtl)) & nobs>3)
  
  polygon(c(tlst[ps], rev(tlst[ps])), c(qtl[ps,1], rev(qtl[ps,2])), col=adjustcolor(cluse, alpha.f = 0.5), border = NA)
}


xtfun_bound<-function(x0, r, d, dt, ndist) {
  r<-ilogit(r)*5
  
  if(length(dt)==1) {
    dt<-rep(dt, length(x0))
  }
  if(length(r)==1) {
    r<-rep(r, length(x0))
  }
  tstep<-dt/(ndist+1)
  stmp1<-x0*exp(-r*tstep)
  if(max(ndist)>0) {
    for(k in 1:max(ndist)) {
      ps<-which(ndist>=k)
      stmp1[ps]<-(stmp1[ps]+d)*exp(-r[ps]*tstep[ps])
    }
  }
  
  stmp1
}



xtfun<-function(x0, r, d, dt, ndist) {
  if(length(d)==1) {
    d<-rep(d, length(x0))
  }
  if(length(r)==1) {
    r<-rep(r, length(x0))
  }
  tstep<-dt/(ndist+1)
  stmp1<-x0*exp(-r*tstep)
  if(max(ndist)>0) {
    for(k in 1:max(ndist)) {
      ps<-which(ndist>=k)
      stmp1[ps]<-(stmp1[ps]+d[ps])*exp(-r[ps]*tstep[ps])
    }
  }
  
  stmp1
}

xt2fun<-function(x0, r, d, d_sd, dt, ndist) {
  if(length(dt)==1) {
    dt<-rep(dt, length(x0))
  }
  if(length(r)==1) {
    r<-rep(r, length(x0))
  }
  if(length(d_sd)==1) {
    d_sd<-rep(d_sd, length(x0))
  }
  
  tstep<-dt/(ndist+1)
  stmp1<-(x0*exp(-r*tstep))^2
  if(max(ndist)>0) {
    for(k in 1:max(ndist)) {
      ps<-which(ndist>=k)
      stmp1[ps]<-(stmp1[ps]+2*d*sqrt(stmp1[ps])+d^2+d_sd[ps]^2)*exp(-2*r[ps]*tstep[ps])
    }
  }
  
  stmp1
}


xt2fun_alt<-function(x0, r, d, d_sd, dt, ndist) {
  tdiffs<-seq(dt/(ndist+1), dt, length=(ndist+1))
  tdiffs<-tdiffs[-length(tdiffs)]
  (x0^2*exp(-2*r*dt))+
    d_sd^2*sum(exp(-2*r*(dt-tdiffs)))
}




#xtfunb<-function(x0, r, d, dt, ndist) {
#  stmp<-rep(0, length(x0))
#  if(length(dt)==1) {
#    dt<-rep(dt, length(x0))
#  }
#  if(length(r)==1) {
#    r<-rep(r, length(x0))
#  }
#  dttmp<-dt/(ndist+1)
#  for(k in 1:max(ndist)) {
#    xps<-ndist>=k
#    stmp[xps]<-stmp[xps]+exp(-k*r[xps]*dttmp[xps])
#  }
#  
#  x0*exp(-r*dt)+d*stmp
#}