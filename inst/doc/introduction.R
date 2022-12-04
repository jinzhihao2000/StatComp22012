## ----eval=T-------------------------------------------------------------------
set.seed(1234)
n=500;p=300;r0=3
U=matrix(rnorm(n*r0),nrow=n)
V=matrix(rnorm(p*r0),nrow=p)
X=U%*%t(V)
Z=matrix(rnorm(n*p),nrow=n)
M=X+Z

## ----eval=T-------------------------------------------------------------------
library(StatComp22012)
rest=estimate(M,rmin=1,rmax =15,p1=0.05,p2=0.3,p3=0.3,p4=0.05,reptime = 20)
rest

## ----eval=FALSE---------------------------------------------------------------
#  gibbsR <- function(N,rho,mu1,mu2,sigma1,sigma2) {
#    X <- matrix(0, N, 2)
#    s1 = sqrt(1 - rho^2) * sigma1
#    s2 = sqrt(1 - rho^2) * sigma2
#    X[1, ] = c(mu1, mu2)
#    for (i in 2:N) {
#      x2 = X[i - 1, 2]
#      m1 = mu1 + rho * (x2 - mu2) * sigma1/sigma2
#      X[i, 1] = rnorm(1, m1, s1)
#      x1 = X[i, 1]
#      m2 = mu2 + rho * (x1 - mu1) * sigma2/sigma1
#      X[i, 2] = rnorm(1, m2, s2) }
#    X
#  }

## ----eval=FALSE---------------------------------------------------------------
#  NumericMatrix gibbsC(int N,int thin,double mu1, double mu2, double sigma1, double sigma2, double rho) {
#    NumericMatrix mat(N, 2);
#    double x=0,y=0;
#    double s1;
#    s1= sqrt(1 - rho*rho) * sigma1;
#    double s2;
#    s2= sqrt(1 - rho*rho) * sigma2;
#    for(int i = 0; i < N; i++) {
#      for(int j = 0; j < thin; j++) {
#        x = rnorm(1,mu1+rho*(y-mu2)*sigma1/sigma2,s1)[0];
#        y = rnorm(1, mu2+rho*(x-mu1)*sigma2/sigma1,s2)[0];
#      }
#      mat(i, 0) = x;
#      mat(i, 1) = y;
#    }
#    return(mat);
#  }

## ----eval=T-------------------------------------------------------------------
library(microbenchmark)
library(StatComp22012)
tm1 <- microbenchmark(
rnR = gibbsR(N=2000, mu1=0, mu2=0, sigma1=1, sigma2=1, rho=0.9),
rnC = gibbsC(2000,1,0,0,1,1,0.9))
knitr::kable(summary(tm1)[,c(1,3,5,6)])

