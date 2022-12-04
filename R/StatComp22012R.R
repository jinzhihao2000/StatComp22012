#' @title Benchmark R and Rcpp functions.
#' @name benchmarks
#' @description Use R package \code{microbenchmark} to compare the performance of C functions \code{gibbsR} and Cpp functions \code{gibbsC}.
#' @examples
#' \dontrun{
#' tm1 <- microbenchmark::microbenchmark(
#'   rnR = gibbsR(N=2000, mu1=0, mu2=0, sigma1=1, sigma2=1, rho=0.9),
#'   rnC = gibbsC(2000,1,0,0,1,1,0.9)
#' )
#' print(summary(tm1)[,c(1,3,5,6)])
#' }
#' @import microbenchmark
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm rbinom
#' @useDynLib StatComp22012
NULL

#' @title Search rank magnitude of the approximating matrix
#' @description Estimate a preferable matrix rank magnitude for fitting a low-rank matrix approximation to a matrix with missing values via the block diagonal missing mechanism.
#' @param M the matrix of estimation
#' @param rmin the start rank for searching
#' @param rmax the max rank for searching
#' @param p1 the observed probability of the upper left block
#' @param p2 the observed probability of the upper right block
#' @param p3 the observed probability of the bottom left block
#' @param p4 the observed probability of the bottom right block
#' @param reptime the number of times of replicated experiments 
#' @return the estimated rank of the matrix with missing data via the block diagonal missing mechanism
#' @examples
#' \dontrun{

#' n=500;p=300;r0=3
#' U=matrix(rnorm(n*r0),nrow=n)
#' V=matrix(rnorm(p*r0),nrow=p)
#' X=U%*%t(V)
#' Z=matrix(rnorm(n*p),nrow=n)
#' M=X+Z
#' rest=estimate(M,rmin=1,rmax =15,p1=0.05,p2=0.3,p3=0.3,p4=0.05,reptime = 20)
#' rest

#' }
#' @export
estimate=function(M,rmin,rmax,p1,p2,p3,p4,reptime){
  n=nrow(M);p=ncol(M);mu=log(log(p))*log(n)
  I_Omega=matrix(rep(0,n*p),nrow = n)  
  error=numeric()
  r.est=rep(0,reptime)
  for (j in 1:reptime) {
    set.seed(j*10)
    I1=matrix(rbinom(n*p/4,1,p1),nrow = n/2)
    I2=matrix(rbinom(n*p/4,1,p2),nrow = n/2)
    I3=matrix(rbinom(n*p/4,1,p3),nrow = n/2)
    I4=matrix(rbinom(n*p/4,1,p4),nrow = n/2)
    I_Omega=rbind(cbind(I1,I2),cbind(I3,I4))
    
    W1=svd(I_Omega)
    W=W1$d[1]*W1$u[,1]%*%t(W1$v[,1])
    Y=I_Omega*M
    
    for (r in rmin:rmax) {
      
      b=svd((W^{-1/2})*Y)
      if(r==1)
      {Y_debiased=(W^{-1/2})*(b$d[1]*b$u[,1:r]%*%t(b$v[,1:r]))
      error[r]=norm((W^{1/2})*(Y-Y_debiased)*I_Omega,type = "f")/norm(W^{1/2} * Y*I_Omega,"f")}
      else
      {Y_debiased=(W^{-1/2})*(b$u[,1:r]%*%diag(b$d[1:r])%*%t(b$v[,1:r]))
      error[r]=norm((W^{1/2})*(Y-Y_debiased)*I_Omega,type = "f")/norm(W^{1/2} * Y*I_Omega,"f")}
    }
    
    rank.seq=rmin:rmax
    res=vector()
    res=p*log(error) + rank.seq*mu
    r.est[j]=rank.seq[which.min(res)]
  }
  
  r.est1=unique(r.est)
  restimate=r.est1[which.max(tabulate(match(r.est, r.est1)))]
  return(restimate)
  
}

#' @title A Gibbs sampler to generate a bivariate normal chain
#' @description A Gibbs sampler to generate a bivariate normal chain
#' @param N the number of samples
#' @param mu1 the mean of the first dimension
#' @param mu2 the mean of the second dimension
#' @param sigma1 the standard deviation of the first dimension
#' @param sigma2 the standard deviation of the second dimension
#' @param rho the correlation of two dimensions
#' @return a random sample of size \code{N}
#' @examples
#' \dontrun{
#' rnR <- gibbsR(N=2000, mu1=0, mu2=0, sigma1=1, sigma2=1, rho=0.9)
#' par(mfrow=c(1,3));
#' plot(rnR[,1],type='l')
#' plot(rnR[,2],type='l')
#' plot(rnR)
#' }
#' @export
gibbsR <- function(N,rho,mu1,mu2,sigma1,sigma2) {
  X <- matrix(0, N, 2)
  s1 = sqrt(1 - rho^2) * sigma1
  s2 = sqrt(1 - rho^2) * sigma2
  X[1, ] = c(mu1, mu2)
  for (i in 2:N) {
    x2 = X[i - 1, 2]
    m1 = mu1 + rho * (x2 - mu2) * sigma1/sigma2
    X[i, 1] = rnorm(1, m1, s1)
    x1 = X[i, 1]
    m2 = mu2 + rho * (x1 - mu1) * sigma2/sigma1
    X[i, 2] = rnorm(1, m2, s2) }
  X
}
