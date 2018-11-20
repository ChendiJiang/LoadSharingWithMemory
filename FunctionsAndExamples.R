#***********************************************************#
# Preliminary. Important functions
# ***********************************************************#
# Part I: Data generation from LSMFM
# input:
# n: sample size
# gamma: the load-share parameter gamma=c(1, gamma_2,\dots, gamma_J)
# theta: the parameter in the initial hazard rate function r_1(u)=r_1(u|theta)
# iniF.gr: a function to generate random variable from the distribution
# with r_1(u|theta)
# iniF.cdf: the cumulative distribution function corresponding to r_1(u|theta)
# iniF.qf: the quantile function corresponding to r_1(u|theta)
# output: A n times J matrix. For i=1,\dots, n, the $i$th row of the output is
# c(S_{i1}, \dots, S_{ij}).
Data.generation=function(n,gamma,theta,iniF.gr,iniF.cdf,iniF.qf)
{
J=length(gamma)
res=matrix(-99,n,J)
latent=matrix(iniF.gr(n*J,theta),n,J)
res[,1]=apply(latent,1,min)
# in terms of surivaval probability exp(-cum.hazard)
cum.hazard=1-iniF.cdf(res[,1],theta)
for(j in 2:J)
{
kappa=iniF.qf(1-(cum.hazard)^(1/gamma[j]),theta)
kappaU=rep(kappa,J-j+1)
sa=(1-iniF.cdf(kappaU,theta))^gamma[j]
Fa=1-sa
qa=runif(n*(J-j+1))*sa+Fa
latentU=iniF.qf(1-(1-qa)^(1/gamma[j]),theta)-kappaU
if(j<J)
{
latentU=apply(matrix(latentU,n,J-j+1),1,min)
}
res[,j]=res[,j-1]+latentU
cum.hazard=(1-iniF.cdf(kappa+latentU,theta))^gamma[j]
}
return(res)
}


#***********************************************************#
# Part II: The sequential search for the initial values
# input:
# Data: the observed sample (output from Data.generation)
# iniF.hr: the hazard rate function r_1(u|theta)
# iniF.cdf: the cumulative distribution function corresponding to r_1(u|theta)
# start.theta: starting point for the search for the intial value of theta.
# LB: lower bound of parameter theta, required by R function nlminb
# output:
# ini.theta: initial value of theta
# ini.gamma: initial value of gamma
INI.values=function(Data,iniF.hr,iniF.cdf,start.theta,LB)
{
J=ncol(Data)
INI.theta=function(tem.theta)
{
sum(-J*log(1-iniF.cdf(Data[,1],tem.theta))-
log(iniF.hr(Data[,1],tem.theta)))
}
ini.theta=nlminb(start.theta,INI.theta,lower=LB)$par
cum.hazard=1-iniF.cdf(Data[,1],ini.theta)
ini.gamma=rep(1,J)
for(j in 2:J)
{
Tlag=Data[,j]-Data[,j-1] # this is about from kappa to kappa+Tlag
INI.gamma=function(tem.gamma)
{
kappa=iniF.qf(1-(cum.hazard)^(1/tem.gamma),ini.theta)
sum(-log(tem.gamma)-log(iniF.hr(kappa+Tlag,ini.theta))+
(J-j+1)*tem.gamma*
(log(1-iniF.cdf(kappa,ini.theta))-
log(1-iniF.cdf(kappa+Tlag,ini.theta))))
}
ini.gamma[j]=optimize(INI.gamma,interval=c(0.001,2*J))$minimum
kappa=iniF.qf(1-(cum.hazard)^(1/ini.gamma[j]),ini.theta)
cum.hazard=(1-iniF.cdf(kappa+Tlag,ini.theta))^ini.gamma[j]
}
return(list(ini.theta=ini.theta,ini.gamma=ini.gamma))
}



#***********************************************************#
# Part III: The log-likelihood function l_n(beta)
# input:
# tem.beta: the value of beta
# Data: the observed sample (output from Data.generation)
# iniF.hr: the hazard rate function r_1(u|theta)
# iniF.cdf: the cumulative distribution function corresponding to r_1(u|theta)
# Memory: type of models. Could be "Memoryless" (for the MLSM); or
# "Partial" (for the LSMRM); or "Full" (for the LSMFM)
# Output: the value of -l_n(tem.beta).
MLE=function(tem.beta, Data,iniF.hr, iniF.cdf ,Memory)
{
J=ncol(Data)
p=length(tem.beta)
n=nrow(Data)
tem.theta=tem.beta[1:(p-J+1)]
tem.gamma=c(1,tem.beta[(p-J+2):p])
res=sum(-J*log(1-iniF.cdf(Data[,1],tem.theta))-
log(iniF.hr(Data[,1],tem.theta)))
if(Memory=="Full"){
cum.hazard=1-iniF.cdf(Data[,1],tem.theta)}
for(j in 2:J){
Tlag=Data[,j]-Data[,j-1] # this is about from kappa to kappa+Tlag
if(Memory=="Full"){
kappa=iniF.qf(1-(cum.hazard)^(1/tem.gamma[j]),tem.theta)}
if(Memory=="Partial"){
kappa=Data[,j-1]}
if(Memory=="Memoryless"){
kappa=rep(0,n)}
res=res+sum(-log(tem.gamma[j])-log(iniF.hr(kappa+Tlag,tem.theta))+
(J-j+1)*tem.gamma[j]*
(log(1-iniF.cdf(kappa,tem.theta))-
log(1-iniF.cdf(kappa+Tlag,tem.theta))))
if(Memory=="Full")
{
cum.hazard=(1-iniF.cdf(kappa+Tlag,tem.theta))^tem.gamma[j]
}
}
return(res)
}


##############################################################
# An illustrative example:
# Set the sample size
n=200
# Set the true load-sharing system
# number of components
J=3
# true value of load-share parameters
gamma=rep(1,J)
# the true initial hazard rate r_1(u|theta)
# we consider the gamma distribution in this example.
theta=c(2,10)
#Generate random variables from initial hazard
iniF.gr=function(J,theta)
{
return(rgamma(J,shape=theta[1],rate=theta[2]))
}
# initial hazard rate function
iniF.hr=function(t, theta)
{
return(
dgamma(t,shape=theta[1],rate=theta[2])/
(1-pgamma(t,shape=theta[1],rate=theta[2])))
}
# cdf of the intial hazard
iniF.cdf=function(t, theta)
{
return(pgamma(t,shape=theta[1],rate=theta[2]))
}
# qunatile of the intial hazard
iniF.qf=function(t, theta)
{
return(qgamma(t,shape=theta[1],rate=theta[2]))
}


#################################################################
# Data generation
set.seed(100)
Data=Data.generation(n,gamma,theta,iniF.gr,iniF.cdf,iniF.qf)
head(Data)
# > head(Data)
# [,1] [,2] [,3]
# [1,] 0.03849823 0.11170043 0.1306754
# [2,] 0.01864550 0.14703415 0.2543772
# [3,] 0.03647510 0.09423964 0.1572879
# [4,] 0.16466805 0.22399873 0.2680226
# [5,] 0.07858919 0.10243634 0.1680540
# [6,] 0.02508926 0.02563251 0.4824652
###############################################################
# Find the inital value:
# Since it is a gamma distribution.
# We set the lower bound of theta to be c(0.000001,0.000001)
require(numDeriv)
LB=rep(1e-6,length(theta))
ini.values=INI.values(Data,iniF.hr,iniF.cdf,rep(1,length(theta)),LB)
ini.values
# > ini.values
# $ini.theta
# [1] 1.926317 9.528780
#
# $ini.gamma
# [1] 1.000000 1.040698 1.097800
ini.theta=ini.values$ini.theta
ini.gamma=ini.values$ini.gamma
ini.beta=c(ini.theta,ini.gamma[-1])
ini.beta
# > ini.beta
# [1] 1.926317 9.528780 1.040698 1.097800


###########################################################
# Computing MLE under LSMFM
mle.full=nlminb(ini.beta,MLE,Data=Data,iniF.hr=iniF.hr,iniF.cdf=iniF.cdf,
Memory="Full", lower=LB)
# MLE under LSMFM
mle.LSMFM=mle.full$par
# estimate standard error of MLE under LSMFM
mle.LSMFM.se=sqrt(diag(solve(hessian(MLE,mle.LSMFM,
Data=Data,iniF.hr=iniF.hr,iniF.cdf=iniF.cdf, Memory="Full"))))
###########################################################
# Computing MLE under LSMRM
mle.partial=nlminb(ini.beta,MLE,Data=Data,iniF.hr=iniF.hr,iniF.cdf=iniF.cdf,
Memory="Partial", lower=LB)
# MLE under LSMRM
mle.LSMRM=mle.partial$par
# estimate standard error of MLE under LSMRM
mle.LSMRM.se=sqrt(diag(solve(hessian(MLE,mle.LSMRM,
Data=Data,iniF.hr=iniF.hr,iniF.cdf=iniF.cdf, Memory="Partial"))))
###########################################################
# Computing MLE under MLSM
mle.memoryless=nlminb(ini.beta,MLE,Data=Data,iniF.hr=iniF.hr,iniF.cdf=iniF.cdf,
Memory="Memoryless", lower=LB)
# MLE under MLSM
mle.MLSM=mle.memoryless$par
# estimate standard error of MLE under MLSM
mle.MLSM.se=sqrt(diag(solve(hessian(MLE,mle.MLSM,Data=Data,
iniF.hr=iniF.hr,iniF.cdf=iniF.cdf, Memory="Memoryless"))))
###########################################################
# Indicator: whether the search for MLE convergned.
# If 0 means a successful convergence
convergence=c(mle.full$convergence,
mle.partial$convergence,
mle.memoryless$convergence)
convergence
# > convergence
# [1] 0 0 0
###########################################################
# Summary of results:
Truth=c(theta,gamma[-1])
res=rbind(Truth,ini.beta,mle.LSMFM,mle.LSMFM.se,
mle.LSMRM,mle.LSMRM.se,mle.MLSM,mle.MLSM.se)

###########################################################
# chenck numerical optimazation
# If this is 0 0 0, then all the searches for mle under LSMFM, LSMRM,
# MLSM are successfull, respectively.
print(convergence)
# > print(convergence)
# [1] 0 0 0
###########################################################
# Final results
# Truth: true value of beta
# ini.beta: initial value found by the sequential algorithm
# mle.LSMFM: mle under LSMFM
# mle.LSMFM.se: estiamted standard error for the mle under LSMFM
# mle.LSMRM: mle under LSMRM
# mle.LSMRM.se: estiamted standard error for the mle under LSMRM
# mle.MLSM: mle under MLSM
# mle.MLSM.se: estiamted standard error for the mle under MLSM
print(res)
# > print(res)
# [,1] [,2] [,3] [,4]
# Truth 2.00000000 10.0000000 1.0000000 1.0000000
# ini.beta 1.92631749 9.5287803 1.0406981 1.0977997
# mle.LSMFM 1.87297712 9.1729677 1.0663132 1.1308257
# mle.LSMFM.se 0.13620426 1.0503823 0.1352474 0.1478414
# mle.LSMRM 1.87697555 9.2016456 1.0550601 1.1170797
# mle.LSMRM.se 0.13429387 1.0357154 0.1163636 0.1336029
# mle.MLSM 1.34740647 5.6980729 1.6828107 1.9582801
# mle.MLSM.se 0.06022699 0.5300502 0.1680700 0.1971233