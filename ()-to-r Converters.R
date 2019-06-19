TtoR <-function(t, df) {
  return(sqrt((t^2)/(t^2+df)))
}
ZtoR <- function(z,N){
  return(sqrt((z^2)/(z^2+N)))
}
Fn1toR <- function(f, dfd){
  return(sqrt((f)/(f+dfd)))
}
Fn2toR <- function(f, dfd, dfn){
  return(sqrt((dfn*f)/((dfn*f)+dfd)))
}
X21toR <- function(X2, N){
  return(sqrt((X2^2)/N))
}
X22toR <- function(X2, N){
  return(sqrt((X2^2)/(X2^2+ N)))
}
DtoR <- function(d){
  return((d)/sqrt((4+d^2)))
}
TautoR <- function(tau){
  return(sin(.5*pi*tau))
}
UtoZ <-function(U, N1, N2){
  return((U-((N1*N2)/2))/(sqrt((N1*N2*(N1+N2+1))/12)))
}

X21toR(4.6, 76)
X22toR(.4, 62)
ZtoR(-1.64, 19)
DtoR(.358)
TautoR(.002)
Fn1toR(9.2, 19)
ZtoR(UtoZ(29,17,7), (17+7))
TtoR(-.4, 6)

ZtoR(-.31, 31)
