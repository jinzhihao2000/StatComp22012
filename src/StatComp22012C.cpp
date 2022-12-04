#include <Rcpp.h>
using namespace Rcpp;

//' @title A Gibbs sampler to generate a bivariate normal chain
//' @description A Gibbs sampler to generate a bivariate normal chain
//' @param N the number of samples
//' @param thin the number of between-sample random numbers
//' @param mu1 the mean of the first dimension
//' @param mu2 the mean of the second dimension
//' @param sigma1 the standard deviation of the first dimension
//' @param sigma2 the standard deviation of the second dimension
//' @param rho the correlation of two dimensions
//' @return a random sample of size \code{n}
//' @examples
//' \dontrun{
//' rnC <- gibbsC(2000,1,0,0,1,1,0.9)
//' par(mfrow=c(1,3));
//' plot(rnC[,1],type='l')
//' plot(rnC[,2],type='l')
//' plot(rnC)
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix gibbsC(int N,int thin,double mu1, double mu2, double sigma1, double sigma2, double rho) {
  NumericMatrix mat(N, 2);
  double x=0,y=0;
  double s1;
  s1= sqrt(1 - rho*rho) * sigma1;
  double s2;
  s2= sqrt(1 - rho*rho) * sigma2;
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < thin; j++) {
      x = rnorm(1,mu1+rho*(y-mu2)*sigma1/sigma2,s1)[0];
      y = rnorm(1, mu2+rho*(x-mu1)*sigma2/sigma1,s2)[0];
    }
    mat(i, 0) = x;
    mat(i, 1) = y;
  }
  return(mat);
}
