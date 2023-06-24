#===============================================
#             SARAR H0 C1 C2 C3
#===============================================

# source('Main.R')
# install.packages("raster",type='binary')
# install.packages("spData",type='binary')
# install.packages("terra",type='binary')
# library('spData')

rm(list = ls()) 	
source('GlambdaChen.R')
library('sp')
library('terra')
library('spdep')
nsim = 5
beta = 3.5
rou1 = 0.85
rou2 = 0.15
tol=1e-6
a=0.95
k=1
cut = qchisq(a,k+3)
size = c(3,4,5,6,7)
size = c(10,13,15)

ff_elH0 = c()
ff_aelH0 = c()
ff_elC1 = c()
ff_aelC1 = c()
ff_elC2 = c()
ff_aelC2 = c()
ff_elC3 = c()
ff_aelC3 = c()
for (m in size){
  # SARAR模型
  n = m*m
  i = 1:n;Xn = i/(n+1)
  Wnb = cell2nb(m,m,type='queen')
  Ws = nb2listw(Wnb)
  Wn = listw2mat(Ws)
  
  Mn = Wn
  In = diag(n)
  An = In - rou1*Wn
  Bn = In - rou2*Mn
  Ani = solve(An)
  Bni = solve(Bn)
  # 估计方程准备
  Gn = Bn%*%Wn%*%Ani%*%Bni
  Gnn = 1/2*(Gn + t(Gn))
  Hn = Mn%*%Bni
  Hnn = 1/2*(Hn + t(Hn))
  # 简化符号
  b = t(Xn)%*%t(Bn)
  g = diag(Gnn)
  h = diag(Hnn)
  s = Bn%*%Wn%*%Ani%*%Xn%*%beta
  if(1>log(n)/2) an=1 else an=log(n)/2
  f<-function(Matrix,Vector){
    irow = nrow(Matrix)
    v =c(0)
    for(i in 2:irow){
      v[i] =Matrix[i,][1:(i-1)]%*%Vector[1:i-1]
    }
    return(v)
  }
  score<-function(sigmC2){
    zz = matrix(NA,nrow=n,ncol=k+3)
    zz[,k] =  b*e
    zz[,k+1] = g*(e^2-sigmC2) + 2*e*f(Gnn,e) + s*e
    zz[,k+2] = h*(e^2-sigmC2) + 2*e*f(Hnn,e)
    zz[,k+3] = e*e - rep(sigmC2, n)
    return(zz)
  }
  # 启动模拟
  f1H0 = 0
  f2H0 = 0
  f1C1 = 0
  f2C1 = 0
  f1C2 = 0
  f2C2 = 0
  f1C3 = 0
  f2C3 = 0
  for(m in 1:nsim){
    En = rnorm(n);
    e = En
    
    # H0估计方程赋值
    z=score(sigmC2=1)
    az=rbind(z,-an*colMeans(z))
    # 计算EL_H0值
    lam = lambdaChen(z)
    el = 2*sum( log(1+t(lam)%*%t(z) ) )
    if(el>cut) f1H0=f1H0+1
    # 计算AEL_H0值
    alam=lambdaChen(az)
    ael=2*sum( log(1+t(alam)%*%t(az) ) )
    if(ael>cut) f2H0=f2H0+1
    
    # C1估计方程赋值
    z=score(sigmC2=1.001)
    az=rbind(z,-an*colMeans(z))
    # 计算EL_C1值	
    lam = lambdaChen(z)
    el = 2*sum( log(1+t(lam)%*%t(z) ) )
    if(el>cut) f1C1=f1C1+1
    # 计算AEL_C1值
    alam=lambdaChen(az)
    ael=2*sum( log(1+t(alam)%*%t(az) ) )  		
    if(ael>cut) f2C1=f2C1+1
    
    # C2估计方程赋值
    z=score(sigmC2=1.1)
    az=rbind(z,-an*colMeans(z))
    # 计算EL值
    lam = lambdaChen(z)
    el = 2*sum( log(1+t(lam)%*%t(z) ) )
    if(el>cut) f1C2=f1C2+1
    # 计算AEL值
    alam=lambdaChen(az)
    ael=2*sum( log(1+t(alam)%*%t(az) ) )
    if(ael>cut) f2C2=f2C2+1
    
    # C3估计方程赋值
    z=score(sigmC2=2)
    az=rbind(z,-an*colMeans(z))
    # 计算EL值
    lam = lambdaChen(z)
    el = 2*sum( log(1+t(lam)%*%t(z) ) )
    if(el>cut) f1C3=f1C3+1
    # 计算AEL值
    alam=lambdaChen(az)
    ael=2*sum( log(1+t(alam)%*%t(az) ) )
    if(ael>cut) f2C3=f2C3+1
    
  }
  cat('样本个数为',n,'完成模拟',m,'次','\n')
  ff_elH0=append(ff_elH0,f1H0/nsim)
  ff_aelH0=append(ff_aelH0,f2H0/nsim)
  ff_elC1=append(ff_elC1,f1C1/nsim)
  ff_aelC1=append(ff_aelC1,f2C1/nsim)
  ff_elC2=append(ff_elC2,f1C2/nsim)
  ff_aelC2=append(ff_aelC2,f2C2/nsim)
  ff_elC3=append(ff_elC3,f1C3/nsim)
  ff_aelC3=append(ff_aelC3,f2C3/nsim)
  Sys.sleep(30)
}
# ff = matrix(NA,nrow=length(size),ncol=8)
# ff[,1]=ff_elH0
# ff[,2]=ff_aelH0
# ff[,3]=ff_elC1
# ff[,4]=ff_aelC1
# ff[,5]=ff_elC2
# ff[,6]=ff_aelC2
# ff[,7]=ff_elC3
# ff[,8]=ff_aelC3
# rownames(ff) <- size^2
# colnames(ff) <- c("E_0","A_0","E_1","A_1","E_2","A_2","E_3","A_3")
# cat('--------------------------------------','\n')
# print(ff)

ff = matrix(NA,nrow=length(size),ncol=8)
ff[,1]=ff_elH0
ff[,2]=ff_elC1
ff[,3]=ff_elC2
ff[,4]=ff_elC3
ff[,5]=ff_aelH0
ff[,6]=ff_aelC1
ff[,7]=ff_aelC2
ff[,8]=ff_aelC3
rownames(ff) <- size^2
colnames(ff) <- c("E_0","E_1","E_2","E_3","A_0","A_1","A_2","A_3")
cat('--------------------------------------','\n')
print(ff)

# ff = matrix(NA,nrow=length(size),ncol=4)
# ff[,1]=ff_elH0
# ff[,2]=ff_elC1
# ff[,3]=ff_elC2
# ff[,4]=ff_elC3
# rownames(ff) <- size^2
# colnames(ff) <- c("E_0","E_1","E_2","E_3")
# cat('--------------------------------------','\n')
# print(ff)
# 
# ff = matrix(NA,nrow=length(size),ncol=4)
# ff[,1]=ff_aelH0
# ff[,2]=ff_aelC1
# ff[,3]=ff_aelC2
# ff[,4]=ff_aelC3
# rownames(ff) <- size^2
# colnames(ff) <- c("A_0","A_1","A_2","A_3")
# cat('--------------------------------------','\n')
# print(ff)
