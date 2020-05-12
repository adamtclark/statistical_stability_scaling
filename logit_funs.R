logit<-function(x, ...) qlogis(x, ...)

ilogit<-function(x, ...) plogis(x, ...)

cloglog<-function(x) {
  log(-log(1-x))
}

icloglog<-function(x) {
  1-exp(-exp(x))
}
