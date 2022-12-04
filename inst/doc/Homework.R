## ----warning=FALSE------------------------------------------------------------
library(MASS)
mydata=data.frame(UScereal[,0],UScereal[,2],UScereal[,4],UScereal[,8],UScereal[,9])
colnames(mydata)=c("Calories","Fat","Sugars","Shelf")
variables=mydata[,1:3]
Shelf=factor(mydata$Shelf)
knitr::kable(head(mydata),align = "llcrr")

## -----------------------------------------------------------------------------
means=aggregate(variables,by=list(Shelf),FUN=mean)
knitr::kable(means,align = "lcrr")

## ----fig.align='center',fig.width=7,fig.height=7------------------------------
d=mahalanobis(variables,colMeans(variables),cov(variables))
qqplot(qchisq(ppoints(nrow(variables)),df=ncol(variables)),d,main="Q-Q plot")
abline(a=0,b=1)

## ----fig.align='center',fig.width=7,fig.height=7------------------------------
library(mvoutlier)
aq.plot(variables)

## -----------------------------------------------------------------------------
set.seed(1234)
n=500;a=2;b=2
u=runif(n)
x=b*((1-u)^(-1/a))
x[1:50]

## ---- fig.align='center',fig.height=7,fig.width=7-----------------------------
hist(x,pro=T,main = expression(paste(f(x)==8*x^(-3))))
y=seq(2,100,0.01)
lines(y,8*y^(-3))

## ----fig.align='center',fig.height=7,fig.width=7------------------------------
x=seq(0,1,0.001)
plot(x,dbeta(x,shape1=2,shape2=2),type = "l",xlab = "x",ylab = "f(x)",ylim = c(0,2),col="blue")
lines(x,dbeta(x,shape1=0.5,shape2=0.5),col="red")
legend("top",c("a<1,b<1","a>1,b>1"),lty=c(1,1),col = c("red","blue"),cex = 0.6)
plot(x,dbeta(x,shape1=0.5,shape2=2),type = "l",xlab = "x",ylab = "f(x)",col="red")
lines(x,dbeta(x,shape1=2,shape2=0.5),col="blue")
legend("top",c("a<1,b>=1","a>=1,b<1"),lty=c(1,1),col = c("red","blue"),cex = 0.6 )
plot(x,dbeta(x,shape1=1,shape2=3),type = "l",xlab = "x",ylab = "f(x)",col="red")
lines(x,dbeta(x,shape1=3,shape2=1),col="blue")
legend("top",c("a=1,b>1","a>1,b=1"),lty=c(1,1),col = c("red","blue") ,cex = 0.6)
plot(x,dbeta(x,shape1=1,shape2=1),type = "l",xlab = "x",ylab = "f(x)",col="red")
legend("top",c("a=1,b=1"),lty=1,col = "red",cex = 0.6 )

## -----------------------------------------------------------------------------
set.seed(1000)
f=function(n,a,b){
  k=0
  y=numeric(n)
  while(k<n)
  {
    u=runif(1)
    x=runif(1)
    j1=x^(a-1)
    j2={(1-x)^(b-1)}
    j=j1*j2
    z={(a-1)/(a+b-2)}^(a-1)*{(b-1)/(a+b-2)}^(b-1)
    if(u<(j/z)){
      k=k+1
      y[k]=x
    }
  }
  y
}

## -----------------------------------------------------------------------------
f(1000,3,2)[1:50]

## ----fig.align='center',fig.height=7,fig.width=7------------------------------
hist(f(1000,3,2),pro=T,main =expression(f(x)==12*x^2*(1-x)))
x=seq(0,1,0.01)
lines(x,12*x*x*(1-x))

## -----------------------------------------------------------------------------
set.seed(1234)
n=10^3;r=4;beta=2
lambda=rgamma(n,r,beta)
x=rexp(n,lambda)
x[1:50]

## -----------------------------------------------------------------------------
set.seed(1234)
n=10^3;r=4;beta=2
lambda=rgamma(n,r,beta)
x=rexp(n,lambda)
x[1:50]

## ----fig.align='center',fig.height=7,fig.width=7------------------------------
hist(x,pro=T,main =expression(f(x)==frac(64,(x+2)^5)))
y=seq(0,100,0.01)
lines(y,(r*(beta^r))/((beta+y)^(r+1)))

## -----------------------------------------------------------------------------
set.seed(1)
quick_sort<-function(x){
  num<-length(x)
  if(num==0||num==1){return(x)
  }else{
    a<-x[1]
    y<-x[-1]
    lower<-y[y<a]
    upper<-y[y>=a]
    return(c(quick_sort(lower),a,quick_sort(upper)))}
}
test1=sample(1:1e4)
head(test1)
q1=quick_sort(test1)
head(q1)
test2=sample(1:2e4)
head(test2)
q2=quick_sort(test2)
head(q2)
test3=sample(1:4e4)
head(test3)
q3=quick_sort(test3)
head(q3)
test4=sample(1:6e4)
head(test4)
q4=quick_sort(test4)
head(q4)
test5=sample(1:8e4)
head(test5)
q5=quick_sort(test5)
head(q5)
q=list(test1,q1,test2,q2,test3,q3,test4,q4,test5,q5)
data=do.call(cbind, lapply(lapply(q, unlist), `length<-`, max(lengths(q))))

## -----------------------------------------------------------------------------
t=matrix(rep(0,8*100),nrow = 100)
for (i in c(1,2,4,6,8)) {
  t[,i]=replicate(100, {
    test=sample(1:(i*10^4))
    system.time(quick_sort(test))[1]
  })
}
an=colMeans(t)[colMeans(t)>0]
an

## ----fig.align='center',fig.height=7,fig.width=7------------------------------
an=colMeans(t)[colMeans(t)>0]
x=c(1e4*log(1e4),2e4*log(2e4),4e4*log(4e4),6e4*log(6e4),8e4*log(8e4))
y=an
plot(x,y,main = "scatter plot and regression line",xlab="tn",ylab="an")
data=data.frame(x,y)
fit=lm(y~x, data=data)
summary(fit)
abline(fit)

## -----------------------------------------------------------------------------
set.seed(1)
m=1e4
thetahat=replicate(500,{
  x=runif(m)
  mean(exp(x))
})
mean1=mean(thetahat)
mean1

## -----------------------------------------------------------------------------
set.seed(1)
m=1e4
thetahat2=replicate(500,{
  u=runif(m/2)
  v=1-u
  mean((exp(u)+exp(v))/2)
})
mean2=mean(thetahat2)
mean2

## -----------------------------------------------------------------------------
var1=var(thetahat)
var2=var(thetahat2)
(var1-var2)/var1

## ----fig.align='center',fig.width=7,fig.height=7------------------------------
x=seq(1,10,0.001)
y=(2*pi)^{-1/2}*x*x*exp(-x*x/2)
plot(x,y,main = expression(g(x)==frac(x^2,sqrt(2*pi))*exp(-x^2/2)),type = "l")

## ----fig.align='center',fig.width=7,fig.height=7------------------------------
x=seq(1,10,0.001)
y=(2*pi)^{-1/2}*x*x*exp(-x*x/2)
plot(x,y,main = expression(g(x)==frac(x^2,sqrt(2*pi))*exp(-x^2/2)),type = "l", ylim = c(0, 1))
lines(x, 2 * dnorm(x,1,1), lty = 2)
lines(x, dgamma(x - 1, sqrt(2), 1),lty=3)
legend("topright", inset = 0.02, legend = c("g(x)", "f1", "f2"), lty = 1:3)

## ----fig.align='center',fig.width=7,fig.height=7------------------------------
x=seq(1,10,0.001)
y=(2*pi)^{-1/2}*x*x*exp(-x*x/2)
plot(x, y/(dgamma(x - 1, sqrt(2), 1)), type = "l", col="red",ylim = c(0,1),ylab = "g/f")
lines(x, y/(2 * dnorm(x,1,1)), col="blue")
legend("topright",c("f1", "f2"),col=c("blue","red"),lty = c(1,1))

## -----------------------------------------------------------------------------
set.seed(1234)
m =10000
f_1 = replicate(1000, expr = {
x = sqrt(rchisq(m, 1)) + 1
f = 2 * dnorm(x,1,1)
g = x^2 * exp(-x^2/2)/sqrt(2 * pi)
mean(g/f)})
f_2 = replicate(1000, expr = {
 x = rgamma(m, sqrt(2), 1) + 1
 f = dgamma(x - 1, sqrt(2), 1)
 g = x^2 * exp(-x^2/2)/sqrt(2 * pi)
 mean(g/f)
 })
c(mean(f_1), mean(f_2))
c(var(f_1), var(f_2))
g <- function(x) x^2 * exp(-x^2/2)/sqrt(2 * pi)
integrate(g, lower = 1, upper = Inf)

## -----------------------------------------------------------------------------
set.seed(1234)
M = 10000
k = 5
m = M/k
thetahat = rep(0,k)
variance =rep(0,k)
g = function(x){exp(-x)/(1 + x^2)} 
f = function(x){(k/(1 - exp(-1))) * exp(-x)} 
for (j in 0:(k-1)) {
 u = runif(m, j/k, (j+1)/k)#inverse transform method
 x = -log(1 - (1 - exp(-1)) * u)
 fg = g(x)/f(x)
 thetahat[j+1] = mean(fg)
 variance[j+1] = var(fg)
 }
sum(thetahat)
mean(variance)
sqrt(mean(variance))
#True value
integrate(g,0,1)

## -----------------------------------------------------------------------------
set.seed(1234)
m = 10000
g = function(x) {
exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
}
u =runif(m) #inverse transform method
x = - log(1 - u * (1 - exp(-1)))
fg = g(x) / (exp(-x) / (1 - exp(-1)))
theta.hat = mean(fg)
theta.hat
se = sd(fg)
se

## -----------------------------------------------------------------------------
set.seed(1234)
M = 10000
k = 5
m = M/k
N=50
est <- matrix(0, N, 1)
for (i in 1:N) {
thetahat = rep(0,k)
g = function(x){exp(-x)/(1 + x^2)} 
f = function(x){(k/(1 - exp(-1))) * exp(-x)} 
for (j in 0:(k-1)) {
 u = runif(m, j/k, (j+1)/k)#inverse transform method
 x = -log(1 - (1 - exp(-1)) * u)
 fg = g(x)/f(x)
 thetahat[j+1] = mean(fg)
 }
est[i,1]=sum(thetahat)
}
round(apply(est,2,mean),4)
round(apply(est,2,sd),5)
#True value
integrate(g,0,1)

## -----------------------------------------------------------------------------
set.seed(1234)
M = 10000
k = 5
m = M/k
N=50
est <- matrix(0, N, 1)
for (i in 1:N) {
thetahat = rep(0,k)
g = function(x){exp(-x)/(1 + x^2)} 
f = function(x){(1/(1 - exp(-1))) * exp(-x)} 
for (j in 0:(k-1)) {
 u = runif(m, j/k, (j+1)/k)#inverse transform method
 x = -log(1 - (1 - exp(-1)) * u)
 fg = g(x)/f(x)
 thetahat[j+1] = mean(fg)
 }
est[i,1]=mean(thetahat)
}
round(apply(est,2,mean),4)
round(apply(est,2,sd),5)
#True value
integrate(g,0,1)

## -----------------------------------------------------------------------------
#clear up memory
rm(list=ls())
#set seed
set.seed(1)
#data generation
data=function(n)
  {
   x = rlnorm(n)
   y = log(x)
   return(y)
  }
#data analysis
CI=function(m){
 replicate(m, expr = {
 n1=30
 y=data(n1)
 ybar = mean(y)
 se = sd(y)/sqrt(n1)
 ybar + se * qt(c(0.025, 0.975),n1-1)
 })}
#result reporting
result=function(m){
  ci=CI(m)
  LCL= ci[1, ]
  UCL= ci[2, ]
  print(sum(LCL < 0 & UCL > 0))
  print(mean(LCL < 0 & UCL > 0))
}
result(10000)

## -----------------------------------------------------------------------------
#clear up memory
rm(list=ls())
#set seed
set.seed(1234)
#data generation
data1=function(n)
{ 
  sigma1 = 1
  x = rnorm(n, 0, sigma1)
  return(x)
}
data2=function(n)
{ 
  sigma2 = 1.5
  y = rnorm(n, 0, sigma2)
  return(y)
}
#data analysis
power1 = function(m){replicate(m, expr={
 count5test = function(x, y) {
 X = x - mean(x)
 Y = y - mean(y)
 outx = sum(X > max(Y)) + sum(X < min(Y))
 outy = sum(Y > max(X)) + sum(Y < min(X))
 return(as.integer(max(c(outx, outy)) > 5))
 }
 x=data1(20)
 y=data2(20)
count5test(x,y)
})}
#result reporting
result=function(m)
{
  mean(power1(m))
}
result(10000)

## -----------------------------------------------------------------------------
#clear up memory
rm(list=ls())
#set seed
set.seed(1234)
#data generation
data1=function(n)
{ 
  sigma1 = 1
  x = rnorm(n, 0, sigma1)
  return(x)
}
data2=function(n)
{ 
  sigma2 = 1.5
  y = rnorm(n, 0, sigma2)
  return(y)
}
#data analysis
power2 = function(m){replicate(m, expr={
 x=data1(20)
 y=data2(20)
 Fp = var.test(x, y)$p.value
 Ftest = as.integer(Fp <= 0.055)
})}
#result reporting
result=function(m)
{
  mean(power2(m))
}
result(10000)

## -----------------------------------------------------------------------------
#clear up memory
rm(list=ls())
#set seed
set.seed(1234)
#data generation
data=function(n){
  sigma1 = 1
  sigma2 = 1.5
  x = rnorm(n, 0, sigma1)
  y = rnorm(n, 0, sigma2)
  data=c(x,y)
  data
}
#data analysis
power1 = function(m,n)
  {
  replicate(m, expr={
  count5test = function(x, y)
  {
  X = x - mean(x)
  Y = y - mean(y)
  outx = sum(X > max(Y)) + sum(X < min(Y))
  outy = sum(Y > max(X)) + sum(Y < min(X))
  return(as.integer(max(c(outx, outy)) > 5))
  }
  x = data(n)[1:n]
  y = data(n)[(n+1):(2*n)]
  count5test(x, y)
  })
  }
power2=function(m,n)
  {replicate(m, expr={
  x = data(n)[1:n]
  y = data(n)[(n+1):(2*n)]
  Fp = var.test(x, y)$p.value
  Ftest = as.integer(Fp <= 0.055)
  })
}
#result reporting
result=function(m,n)
{ 
  c(mean(power1(m,n)),mean(power2(m,n)))
}
rbind(result(10000,30),result(10000,100),result(10000,500))

## -----------------------------------------------------------------------------
#clear up memory
rm(list=ls())
#set seed
set.seed(1234)
#data generation
data=function(){
  library(boot)
  data=aircondit
  data
}
#data analysis
obj=function(x){
  library(boot)
  rate = function(x, i) 1/mean(as.matrix(x[i, ]))
  obj1=boot(data(), statistic = rate, R = x)
  output=list(obj1$t0,mean(obj1$t),sd(obj1$t))
  return(output)
}
#result reporting
result=function(x){
  a=obj(1e4)
  round(c(original=a[[1]],bias=a[[2]]-a[[1]],
            se=a[[3]]),x)
}
result(4)

## -----------------------------------------------------------------------------
#clear up memory
rm(list=ls())
#set seed
set.seed(1)
#data generation
data=function(){
  library(boot)
  data=aircondit
  data
}
#data analysis
obj=function(x){
  means = function(x, i) mean(as.matrix(x[i, ]))
  obj1=boot(data(), statistic = means, R = x)
  obj2=boot.ci(obj1, type = c("norm","basic","perc","bca"))
  return(obj2)
}

#result reporting
result=function(x){
 obj(x)
}
result(2000)

## ----fig.align='center',fig.width=7,fig.height=7------------------------------
#clear up memory
rm(list=ls())
#set seed
set.seed(1)
#data generation
data=function(){
  library(boot)
  data=aircondit
  data
}
#data analysis
obj=function(x){
  means = function(x, i) mean(as.matrix(x[i, ]))
  obj1=boot(data(), statistic = means, R = x)
  return(obj1)
}

#result reporting
result=function(x){
 hist(obj(x)$t, prob = TRUE, main = "")
 points(obj(x)$t0, 0, cex = 1, pch = 16)
}
result(2000)

## -----------------------------------------------------------------------------
#clear up memory
rm(list=ls())
#set seed
set.seed(1234)
#data generation
data1=function(n){
  rnorm(n, 0, 1)
}
data2=function(m){
  matrix(0, m, 2)
}
nor.norm =data2(10000)
nor.basic=data2(10000)
nor.perc =data2(10000)
#data analysis
cp=function(m){
  library(boot)
  for (i in 1:m) {
  data.nor =data1(30)
  means = function(x,i) mean(x[i])
  nor.means <- boot(data.nor, statistic = means, R=1000)
  nor <- boot.ci(nor.means, type=c("norm","basic","perc"))
  nor.norm[i,] <- nor$norm[2:3]
  nor.basic[i,] <- nor$basic[4:5]
  nor.perc[i,] <- nor$percent[4:5]
  }
  mu = 0
  norm1 = mean(nor.norm[,1] <= mu & nor.norm[,2] >= mu)
  basic1 = mean(nor.basic[,1] <= mu & nor.basic[,2] >= mu)
  perc1 = mean(nor.perc[,1] <= mu & nor.perc[,2] >= mu)
  return(list(norm1,basic1,perc1))
}

#result reporting
result=function(m){cp(m)}
a=result(10000)
Distribution = c("N(0,1)","N(0,1)","N(0,1)")
Type = c("basic", "norm", "perc")
P.coverage = c(a[[1]],a[[2]],a[[3]])
data.frame(Distribution,Type,P.coverage)


## -----------------------------------------------------------------------------
#clear up memory
rm(list=ls())
#set seed
set.seed(1234)
#data generation
data1=function(n){
  rnorm(n, 0, 1)
}
data2=function(m){
  matrix(0, m, 2)
}
nor.norm =data2(10000)
nor.basic=data2(10000)
nor.perc =data2(10000)
#data analysis
cp=function(m){
  library(boot)
  for (i in 1:m) {
  data.nor =data1(30)
  means = function(x,i) mean(x[i])
  nor.means <- boot(data.nor, statistic = means, R=1000)
  nor <- boot.ci(nor.means, type=c("norm","basic","perc"))
  nor.norm[i,] <- nor$norm[2:3]
  nor.basic[i,] <- nor$basic[4:5]
  nor.perc[i,] <- nor$percent[4:5]
  }
  mu = 0
#Calculate the probability of the left side of the normal distribution
    norm1.left <- mean(nor.norm[,1] >= mu )
    basic1.left <- mean(nor.basic[,1] >= mu )
    perc1.left <- mean(nor.perc[,1] >= mu )
#Calculate the right side probability of a normal distribution
    norm1.right <- mean(nor.norm[,2] <= mu )
    basic1.right <- mean(nor.basic[,2] <= mu )
    perc1.right <- mean(nor.perc[,2] <= mu)
  return(list(norm1.left,  basic1.left, perc1.left,norm1.right,basic1.right, perc1.right))
}
#result reporting
result=function(m){cp(m)}
a=result(10000)

Distribution = c("N(0,1)")
Type = c("basic", "norm", "perc")
Left = c(a[[1]],a[[2]],a[[3]])
Right = c(a[[4]],a[[5]],a[[6]])
data.frame(Distribution, Type, Left, Right)

## -----------------------------------------------------------------------------
#clear up memory
rm(list=ls())
#set seed
set.seed(1234)
#data generation
data=function(){
  library(bootstrap)
  x = as.matrix(scor)
  x
}
#data analysis
jack=function(n){
  theta.jack = numeric(n)
  lambda = eigen(cov(data()))$values
  theta.hat = max(lambda/sum(lambda))
  for (i in 1:n) {
    y = data()[-i, ]
    s = cov(y)
    lambda = eigen(s)$values
    theta.jack[i] = max(lambda/sum(lambda))
 }
  bias.jack = (n - 1) * (mean(theta.jack) - theta.hat)
  se.jack = sqrt((n - 1)*mean((theta.jack - mean(theta.jack))^2))
  return(list(est = theta.hat, bias = bias.jack, se = se.jack))
}
#result reporting
result=function(n)
{
  jack(n)
}
result(nrow(data()))

## -----------------------------------------------------------------------------
#clear up memory
rm(list=ls())
#set seed
set.seed(1234)
#data generation
data=function(){
  library(DAAG)
  data=ironslag
  data
}
#data analysis
pe=function(n){

N = choose(n, 2)
e1 = e2 = e3 = e4 = e5 = numeric(N)
ij = 1
for (i in 1:(n - 1)) for (j in (i + 1):n) {
 k = c(i, j)
 y = data()$magnetic[-k]
 x = data()$chemical[-k]
 J1 = lm(y ~ x)
 yhat1 = J1$coef[1] + J1$coef[2] * data()$chemical[k]
 e1[ij] = sum((data()$magnetic[k] - yhat1)^2)
 J2 = lm(y ~ x + I(x^2))
 yhat2 = J2$coef[1] + J2$coef[2] * data()$chemical[k] + J2$coef[3] * data()$chemical[k]^2
 e2[ij] = sum((data()$magnetic[k] - yhat2)^2)
 J3 = lm(log(y) ~ x)
 logyhat3 = J3$coef[1] + J3$coef[2] * data()$chemical[k]
 yhat3 = exp(logyhat3)
 e3[ij] = sum((data()$magnetic[k] - yhat3)^2)
 J4 = lm(log(y) ~ log(x))
 logyhat4 = J4$coef[1] + J4$coef[2] * log(data()$chemical[k])
 yhat4 = exp(logyhat4)
 e4[ij] = sum((data()$magnetic[k] - yhat4)^2)
 c2 = x^2
 c3 = x^3
 J5 = lm(y ~ x + c2 + c3)
 yhat5 = J5$coef[1] + J5$coef[2] * data()$chemical[k] + J5$coef[3] * data()$chemical[k]^2 + J5$coef[4] * data()$chemical[k]^3
 e5[ij] = sum((data()$magnetic[k] - yhat5)^2)
 ij = ij + 1
}
return(list(sum(e1),sum(e2),sum(e3),sum(e4),sum(e5)))
}
#result reporting
result=function(n)
{ 
  N = choose(n, 2)
  c(pe(n)[[1]],pe(n)[[2]],pe(n)[[3]],pe(n)[[4]],pe(n)[[5]])/N
}
result(length(data()$magnetic))

## -----------------------------------------------------------------------------
#clear up memory
rm(list=ls())
#set seed
set.seed(0)
#data generation
data=function(n){
   set.seed(0)
   library(MASS)
   mu = c(0, 0)
   Sigma = matrix(c(1, 0.5, 0.5, 1), 2, 2)
   x = mvrnorm(n, mu, Sigma)
   x
}
#data analysis
spearman.permutaion = function(x, y) {
  stest = cor.test(x, y, method = "spearman")
  n = length(x)
  R = 499
  rs = replicate(R, expr = {
  k = sample(1:n)
  cor.test(x, y[k], method = "spearman")$estimate})
  rs1 = c(stest$estimate, rs)
  pval = mean(as.integer(stest$estimate <=
  rs1))
  return(list(rho.s = stest$estimate, pvalue = pval))
}
#result reporting
result1=function(n)
{
  cor.test(data(n)[, 1], data(n)[, 2], method = "spearman")
}
result2=function(n){
  spearman.permutaion(data(n)[, 1], data(n)[, 2])
}
result1(40)
result2(40)

## ----fig.align='center',fig.width=7,fig.height=7------------------------------
rm(list = ls())
set.seed(12)
rw.Laplace <- function(N, x0, sigma) {
x = numeric(N)
x[1] = x0
u = runif(N)
k = 0
for (i in 2:N) {
  xt = x[i - 1]
  y = rnorm(1, xt, sigma)
  if (u[i] <= exp(abs(xt) - abs(y)))
  x[i] = y
  else {
  x[i] = x[i - 1]
  k = k + 1
  }
 }
 return(list(x = x, k = k))
}
N = 10000
sigma = c(0.5, 1, 1.5, 4)
x0 = 0.5
rw1 = rw.Laplace(N, x0, sigma[1])
rw2 = rw.Laplace(N, x0, sigma[2])
rw3 = rw.Laplace(N, x0, sigma[3])
rw4 = rw.Laplace(N, x0, sigma[4])
#trace plot
plot(rw1$x, type = "l",xlab=bquote(sigma == .(round(sigma[1],3))))
plot(rw2$x, type = "l",xlab=bquote(sigma == .(round(sigma[2],3))))
plot(rw3$x, type = "l",xlab=bquote(sigma == .(round(sigma[3],3))))
plot(rw4$x, type = "l",xlab=bquote(sigma == .(round(sigma[4],3))))

#histogram

p = ppoints(100)
y = qexp(p, 1)
z = c(-rev(y), y)
fx = 0.5 * exp(-abs(z))
hist(rw1$x, breaks = "Scott", freq = FALSE, ylim = c(0,0.5),main="")
lines(z, fx,col="red")
hist(rw2$x, breaks = "Scott", freq = FALSE, ylim = c(0,0.5),main="")
lines(z, fx,col="red")
hist(rw3$x, breaks = "Scott", freq = FALSE, ylim = c(0,0.5),main="")
lines(z, fx,col="red")
hist(rw4$x, breaks = "Scott", freq = FALSE, ylim = c(0,0.5),main="")
lines(z, fx,col="red")

#qqplot

Q1 = quantile(rw1$x, p)
qqplot(z, Q1,xlab="Theoretical Quantiles", ylab="Sample Quantiles")
abline(0, 1,col="red")
Q2 = quantile(rw2$x, p)
qqplot(z, Q2,xlab="Theoretical Quantiles", ylab="Sample Quantiles")
abline(0, 1,col="red")
Q3 = quantile(rw3$x, p)
qqplot(z, Q3,xlab="Theoretical Quantiles", ylab="Sample Quantiles")
abline(0, 1,col="red")
Q4 = quantile(rw4$x, p)
qqplot(z, Q4,xlab="Theoretical Quantiles", ylab="Sample Quantiles")
abline(0, 1,col="red")


## -----------------------------------------------------------------------------
reject=c(rw1$k, rw2$k, rw3$k, rw4$k)/N
accept=1-reject
cat("acceptance rates",accept)

## ----fig.align='center',fig.width=7,fig.height=7------------------------------
X=cbind(rw1$x, rw2$x, rw3$x, rw4$x)
#compute R-hat statistics
Gelman.Rubin = function(psi) {
# psi[i,j] is the statistic psi(X[i,1:j])
# for chain in i-th row of X
psi = as.matrix(psi)
n = ncol(psi)
k = nrow(psi)
psi.means = rowMeans(psi) #row means
B = n * var(psi.means) #between variance est.
psi.w = apply(psi, 1, "var") #within variances
W = mean(psi.w) #within est.
v.hat = W*(n-1)/n + (B/n) #upper variance est.
r.hat = v.hat / W #G-R statistic
return(r.hat)
}

#compute diagnostic statistics
psi <- apply(X, 1, cumsum)
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))

#plot psi for the four chains
for (i in 1:4)
plot(psi[i, 1:N], type="l",
xlab=i, ylab=bquote(psi))


#plot the sequence of R-hat statistics
rhat <- rep(0, N)
for (j in 1:N)
rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[1:N], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
rm(list=ls())
set.seed(1234)
N = 5000
X = matrix(0, N, 2)
rho = 0.9;mu1 = 0;mu2 = 0;sigma1 = 1;sigma2 = 1
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
head(X)

## ----fig.align='center',fig.height=7,fig.width=7------------------------------
burn=500
b=burn+1
x=X[b:N, ]
x1=x[, 1]
x2=x[, 2]
plot(x1,x2)
abline(h = 0, v = 0)

## -----------------------------------------------------------------------------
lm1=lm(x2 ~ x1)
summary(lm1)

## ----fig.align='center',fig.height=7,fig.width=7------------------------------
qqnorm(lm1$res)
qqline(lm1$res)

## ----fig.align='center',fig.height=7,fig.width=7------------------------------
plot(lm1$fit, lm1$res)
abline(h = 0)

## -----------------------------------------------------------------------------
rm(list=ls())
set.seed(1234)
N = 5000
#chain1
X1 = matrix(0, N, 2)
rho = 0.9;mu1 = 0;mu2 = 0;sigma1 = 1;sigma2 = 1
s1 = sqrt(1 - rho^2) * sigma1
s2 = sqrt(1 - rho^2) * sigma2
X1[1, ] = c(mu1, mu2)
for (i in 2:N) {
x2 = X1[i - 1, 2]
m1 = mu1 + rho * (x2 - mu2) * sigma1/sigma2
X1[i, 1] = rnorm(1, m1, s1)
x1 = X1[i, 1]
m2 = mu2 + rho * (x1 - mu1) * sigma2/sigma1
X1[i, 2] = rnorm(1, m2, s2) }
#chain2
X2 = matrix(0, N, 2)
X2[1, ] =c(-1,1)
for (i in 2:N) {
x2 = X2[i - 1, 2]
m1 = mu1 + rho * (x2 - mu2) * sigma1/sigma2
X2[i, 1] = rnorm(1, m1, s1)
x1 = X2[i, 1]
m2 = mu2 + rho * (x1 - mu1) * sigma2/sigma1
X2[i, 2] = rnorm(1, m2, s2) }
#chain3
X3 = matrix(0, N, 2)
X3[1, ] =c(2,3)
for (i in 2:N) {
x2 = X3[i - 1, 2]
m1 = mu1 + rho * (x2 - mu2) * sigma1/sigma2
X3[i, 1] = rnorm(1, m1, s1)
x1 = X3[i, 1]
m2 = mu2 + rho * (x1 - mu1) * sigma2/sigma1
X3[i, 2] = rnorm(1, m2, s2) }
#chain4
X4 = matrix(0, N, 2)
X4[1, ] =c(-3,-2)
for (i in 2:N) {
x2 = X4[i - 1, 2]
m1 = mu1 + rho * (x2 - mu2) * sigma1/sigma2
X4[i, 1] = rnorm(1, m1, s1)
x1 = X4[i, 1]
m2 = mu2 + rho * (x1 - mu1) * sigma2/sigma1
X4[i, 2] = rnorm(1, m2, s2) }

## -----------------------------------------------------------------------------
XX=cbind(X1[,1], X2[,1], X3[,1], X4[,1])
YY=cbind(X1[,2], X2[,2], X3[,2], X4[,2])

## ----fig.align='center',fig.height=7,fig.width=7------------------------------
Gelman.Rubin = function(psi) {
psi = as.matrix(psi)
n = ncol(psi)
k = nrow(psi)
psi.means = rowMeans(psi) #row means
B = n * var(psi.means) #between variance est.
psi.w = apply(psi, 1, "var") #within variances
W = mean(psi.w) #within est.
v.hat = W*(n-1)/n + (B/n) #upper variance est.
r.hat = v.hat / W #G-R statistic
return(r.hat)
}

#compute diagnostic statistics
psi <- apply(XX, 1, cumsum)
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))

#plot psi for the four chains
for (i in 1:4)
plot(psi[i, 1:N], type="l",
xlab=i, ylab=bquote(psi))

#plot the sequence of R-hat statistics
rhat <- rep(0, N)
for (j in 1:N)
rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[1:N], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)


## ----fig.align='center',fig.height=7,fig.width=7------------------------------
#compute diagnostic statistics
psi <- apply(YY, 1, cumsum)
for (i in 1:nrow(psi))
psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))

#plot psi for the four chains
for (i in 1:4)
plot(psi[i, 1:N], type="l",
xlab=i, ylab=bquote(psi))

#plot the sequence of R-hat statistics
rhat <- rep(0, N)
for (j in 1:N)
rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[1:N], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
rm(list=ls())
set.seed(99999)
#generate sample X,M,Y
N=10000;alpha=0;beta=0;gamma=1
X=rnorm(N,0,1);em=rnorm(N,0,1);ey=rnorm(N,0,1)
M=alpha*X+em
Y=beta*M+gamma*X+ey
#original test
library(bda)
#mediation.test(M,X,Y)
t01=mediation.test(M,X,Y)[1,1]
p1.origin=mediation.test(M,X,Y)[2,1]
p1.origin
#permutation test
m=1000
XX=matrix(0,nrow = N,ncol = m)
YY=matrix(0,nrow = N,ncol = m)
t=numeric(0)
for (i in 1:m) {
  XX[,i]=  X[sample(1:N,N)]
  YY[,i]=  Y[sample(1:N,N)]
  t[i]=mediation.test(M,XX[,i],YY[,i])[1,1]
}
p1.permutation=sum(abs(t)>abs(t01))/m
p1.permutation

## -----------------------------------------------------------------------------
rm(list=ls())
set.seed(99999)
#generate sample X,M,Y
N=10000;alpha=0;beta=1;gamma=1
X=rnorm(N,0,1);em=rnorm(N,0,1);ey=rnorm(N,0,1)
M=alpha*X+em
Y=beta*M+gamma*X+ey
#original test
library(bda)
#mediation.test(M,X,Y)
t02=mediation.test(M,X,Y)[1,1]
p2.origin=mediation.test(M,X,Y)[2,1]
p2.origin
#permutation test
m=1000
XX=matrix(0,nrow = N,ncol = m)
t=numeric(0)
for (i in 1:m) {
  XX[,i]=  X[sample(1:N,N)]
  t[i]=mediation.test(M,XX[,i],Y)[1,1]
}
p2.permutation=sum(abs(t)>abs(t02))/m
p2.permutation

## -----------------------------------------------------------------------------
rm(list=ls())
set.seed(99999)
#generate sample X,M,Y
N=10000;alpha=1;beta=0;gamma=1
X=rnorm(N,0,1);em=rnorm(N,0,1);ey=rnorm(N,0,1)
M=alpha*X+em
Y=beta*M+gamma*X+ey
#original test
library(bda)
#mediation.test(M,X,Y)
t03=mediation.test(M,X,Y)[1,1]
p3.origin=mediation.test(M,X,Y)[2,1]
p3.origin
#permutation test
m=1000
YY=matrix(0,nrow = N,ncol = m)
t=numeric(0)
for (i in 1:m) {
  YY[,i]=  X[sample(1:N,N)]
  t[i]=mediation.test(M,X,YY[,i])[1,1]
}
p3.permutation=sum(abs(t)>abs(t03))/m
p3.permutation

## -----------------------------------------------------------------------------
rm(list = ls())
f <- function(N,b1,b2,b3,f0){
  g <- function(alpha){
  tmp <- exp(-alpha-b1*x1-b2*x2-b3*x3)
  p<-1/(1+tmp)
  mean(p) - f0
  }
  solution <- uniroot(g,c(-20,0))
  alpha <- solution$root
  alpha
}

## -----------------------------------------------------------------------------
set.seed(1234)
N=10^6;b1=0;b2=1;b3=-1;
f0=c(0.1,0.01,0.001,0.0001)
x1 <- rpois(N,1);x2=rexp(N,1); x3 <- sample(0:1,N,replace=TRUE)
alpha=c(f(N,b1,b2,b3,f0[1]),f(N,b1,b2,b3,f0[2]),f(N,b1,b2,b3,f0[3]),f(N,b1,b2,b3,f0[4]))
data=data.frame(f0=f0,alpha=alpha)
data

## ----fig.align='center',fig.width=7,fig.height=7------------------------------
plot(x=-log10(f0),alpha)


## -----------------------------------------------------------------------------
rm(list = ls())
u=c(11,8,27,13,16,0,23,10,24,2)
v=c(12,9,28,14,17,1,24,11,25,3)
f=function(lambda){
  I=numeric(10)
  J=numeric(10)
  for (i in 1:10) {
    I[i]=v[i]*exp((-1)*(lambda*v[i]))-u[i]*exp((-1)*(lambda*u[i]))
    J[i]=exp((-1)*(lambda*u[i]))-exp((-1)*(lambda*v[i]))
  }
  sum(I/J)
}
solution=uniroot(f,c(0,10))
solution

## -----------------------------------------------------------------------------
rm(list = ls())
u=c(11,8,27,13,16,0,23,10,24,2)
v=c(12,9,28,14,17,1,24,11,25,3)
f=function(lambda){
  I=numeric(10)
  J=numeric(10)
  for (i in 1:10) {
    I[i]=u[i]*exp((-1)*(lambda*u[i]))-v[i]*exp((-1)*(lambda*v[i]))
    J[i]=exp((-1)*(lambda*u[i]))-exp((-1)*(lambda*v[i]))
  }
  sum(I/J)
}
EM=function(lambda2,n,max.it=10000,eps=1e-5){

  i=1
  lambda1=1

  
  while( abs(lambda1 - lambda2) >= eps){
    lambda1 = lambda2
    lambda2 = n/((n/lambda1)+f(lambda1))
    if(i == max.it) break
    i = i + 1    
  }
  return(lambda2)
}

initial=c(0.5,0.1,0.01)
estimate=c(EM(lambda2=0.5,n=10,max.it=10000,eps=1e-5),
EM(lambda2=0.1,n=10,max.it=10000,eps=1e-5),
EM(lambda2=0.01,n=10,max.it=10000,eps=1e-5))
data.frame(initial,estimate)

## -----------------------------------------------------------------------------
rm(list = ls())
a=c(1,2,3)
dim(a)

## -----------------------------------------------------------------------------
rm(list = ls())
a=matrix(1:9,nrow=3)
is.matrix(a)
is.array(a)

## -----------------------------------------------------------------------------
rm(list = ls())
a=1:3
b=c("1","2","3")
c=c(1.233,2.455,88.888)
d=c(TRUE,TRUE,FALSE)
e=c(FALSE,FALSE,TRUE)
data1=data.frame(a,b,c,d)
as.matrix(data1)
data2=data.frame(d,e)
as.matrix(data2)
data3=data.frame(a,d,e)
as.matrix(data3)

## -----------------------------------------------------------------------------
data("iris")
iris[FALSE,]
iris[ , FALSE]
iris[FALSE, FALSE]

## -----------------------------------------------------------------------------
scale01 <- function(x) {
rng <- range(x, na.rm = TRUE)
(x - rng[1]) / (rng[2] - rng[1])
}

## -----------------------------------------------------------------------------
scale01 <- function(x) {
rng <- range(x, na.rm = TRUE)
(x - rng[1]) / (rng[2] - rng[1])
}
library(MASS)
head(UScereal)
head(data.frame(lapply(UScereal, function(x) if (is.numeric(x)) scale01(x) else x)))

## -----------------------------------------------------------------------------
head(quakes)
vapply(quakes, sd,numeric(1))

## -----------------------------------------------------------------------------
head(warpbreaks)
vapply(warpbreaks[vapply(warpbreaks, is.numeric, logical(1))],sd, numeric(1))

## ----eval=FALSE---------------------------------------------------------------
#  #include <Rcpp.h>
#  using namespace Rcpp;
#  // [[Rcpp::export]]
#  NumericMatrix gibbsC(int N,int thin) {
#    NumericMatrix mat(N, 2);
#    double rho = 0.9, mu1 = 0,mu2=0,sigma1=1,sigma2=1,x=0,y=0;
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

## -----------------------------------------------------------------------------
rm(list=ls())
set.seed(1234)
gibbsR <- function(N) {
  X <- matrix(0, N, 2)
  rho = 0.9;mu1 = 0;mu2 = 0;sigma1 = 1;sigma2 = 1
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
gibbsR=gibbsR(2000)

## ----eval=FALSE---------------------------------------------------------------
#  library(Rcpp)
#  cppFunction('NumericMatrix gibbsC(int N,int thin) {
#    NumericMatrix mat(N, 2);
#    double rho = 0.9, mu1 = 0,mu2=0,sigma1=1,sigma2=1,x=0,y=0;
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
#  }')
#  gibbsC=gibbsC(2000,1)

## ----fig.align='center',eval=FALSE--------------------------------------------
#  par(mfrow=c(1,3))
#  qqplot(gibbsC,gibbsR)
#  qqplot(gibbsC[,1],gibbsR[,1])
#  qqplot(gibbsC[,2],gibbsR[,2])

## ----eval=FALSE---------------------------------------------------------------
#  library(microbenchmark)
#  library(Rcpp)
#  dir_cpp <- 'D:/研究生/研一上/统计计算/A-22012-2022-11-18/Rcpp/'
#  source(paste0(dir_cpp,'gibbsR.R'))
#  sourceCpp(paste0(dir_cpp,'gibbsC.cpp'))
#  ts = microbenchmark(gibbsR=gibbsR(2000),
#                           gibbsC=gibbsC(2000,1))
#  summary(ts)[,c(1,3,5,6)]

