logit<-function(x) {
  suppressWarnings(res<-(-log(1/x-1)))
  res[!is.finite(res)]<-NA
  res
}

ilogit<-function(x) {
  1/(1+exp(-x))
}

cloglog<-function(x) {
  log(-log(1-x))
}

icloglog<-function(x) {
  1-exp(-exp(x))
}
