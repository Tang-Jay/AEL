#########################################
#   function L=ELratio(z)
#  find the root lambda of sum(u/(1+lambda*u))=0
# n=numel(z)
# b1=-1./(max(z)-0.0000001);b2=-1./(min(z)+0.0000001)
# dif=1
# tol=1e-4
# while (dif>tol)
#     lambda0=(b2+b1)/2
#        gg=(1/n)*sum(z./(1+lambda0*z+0.00000001))
#        if (gg>0) 
#             b1=lambda0
#         else b2=lambda0
#         end
#         %dif=abs(gg);
#        dif=b2-b1
#  end
#  L=2*sum(log(1+lambda0*z))
#####################################################
# R codes for profiling the EL confidence interval #
#####################################################
lambda<-function(x,mu)
{
  L<--1/max(x-mu)
  R<--1/min(x-mu)
  dif<-1
  tol<-1e-08
  while(dif>tol){
    M<-(L+R)/2
    glam<-sum((x-mu)/(1+M*(x-mu)))
    if(glam>0) L<-M
    if(glam<0) R<-M
    dif<-abs(glam)
  }
  return(M)
}
#=============================================
# Repeated simulation runs start from here!!
#=============================================
nsim<-1000
a<-0.95
cut<-qchisq(a,1)
mu<-0
f<-0
L<-0
U<-0
m=1
for(m in 1:nsim){
  x<-rnorm(100)
  #------------
  tol<-1e-08
  t1<-mean(x)
  t2<-max(x) 
  dif<-t2-t1
  while(dif>tol){
    tau<-(t1+t2)/2
    M<-lambda(x,tau)
    elratio<-2*sum(log(1+M*(x-tau)))
    if(elratio>cut) t2<-tau
    if(elratio<=cut) t1<-tau
    dif<-t2-t1
  }
  UB<-(t1+t2)/2
  #------------
  t1<-mean(x)
  t2<-min(x) 
  dif<-t1-t2
  while(dif>tol){
    tau<-(t1+t2)/2
    M<-lambda(x,tau)
    elratio<-2* sum(log(1+M*(x-tau)))
    if(elratio>cut) t2<-tau
    if(elratio<=cut) t1<-tau
    dif<-t1-t2
  }
  LB<-(t1+t2)/2
  el<-2*sum(log(1+lambda(x,mu)*(x-mu)))
  if(el<=cut) f=f+1
  el
  cut
  L=LB+1
  U=UB+1
  m=m+1
}
fL=L/nsim
fU=U/nsim
ff=f/nsim
length=fU-fL
fL
fU
ff
length

