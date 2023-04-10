#===============================================
#             SARAR H0 A1 A2
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
nsim = 1000
beta = 3.5
rou1 = 0.85
rou2 = 0.15

tol=1e-6
a=0.95
k=1
cut = qchisq(a,k+3)
size = c(3,4,5,6,7,10,13,20)
size = c(20)

ff_elH0 = c()
ff_aelH0 = c()
ff_elA1 = c()
ff_aelA1 = c()
ff_elA2 = c()
ff_aelA2 = c()
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
  f<-function(Matrix,Vector){
    irow = nrow(Matrix)
    v =c(0)
    for(i in 2:irow){
      v[i] =Matrix[i,][1:(i-1)]%*%Vector[1:i-1]
    }
    return(v)
  }
  if(1>log(n)/2) an=1 else an=log(n)/2
  # 启动模拟
  f1H0 = 0
  f2H0 = 0
  f1A1 = 0
  f2A1 = 0
  f1A2 = 0
  f2A2 = 0
  for(m in 1:nsim){
    En = rnorm(n);
    e = En
    # # H0估计方程赋值
    # sigma2=1
    # z = matrix(NA,nrow=n,ncol=k+3)
    # z[,k] =  b*e
    # z[,k+1] = g*(e^2-sigma2) + 2*e*f(Gnn,e) + s*e
    # z[,k+2] = h*(e^2-sigma2) + 2*e*f(Hnn,e)
    # z[,k+3] = e*e - rep(sigma2, n)
    # az=rbind(z,-an*colMeans(z))	
    # # 计算EL_H0值	
    # lam = lambdaChen(z)
    # el = 2*sum( log(1+t(lam)%*%t(z) ) )
    # if(el>cut) f1H0=f1H0+1
    # # 计算AEL_H0值
    # alam=lambdaChen(az)
    # ael=2*sum( log(1+t(alam)%*%t(az) ) )  		
    # if(ael>cut) f2H0=f2H0+1
    
    # A1估计方程赋值
    # sigma2=1+n^(-1/2)
    # sigma2=1+n^(-3)
    sigma2=1
    z = matrix(NA,nrow=n,ncol=k+3)
    z[,k] =  b*e
    z[,k+1] = g*(e^2-sigma2) + 2*e*f(Gnn,e) + s*e
    z[,k+2] = h*(e^2-sigma2) + 2*e*f(Hnn,e)
    z[,k+3] = e*e - rep(sigma2, n)
    az=rbind(z,-an*colMeans(z))	
    # 计算EL_A1值	
    lam = lambdaChen(z)
    el = 2*sum( log(1+t(lam)%*%t(z) ) )
    if(el>cut) f1A1=f1A1+1
    # 计算AEL_A1值
    alam=lambdaChen(az)
    ael=2*sum( log(1+t(alam)%*%t(az) ) )  		
    if(ael>cut) f2A1=f2A1+1
    
    # # A2估计方程赋值
    # sigma2=2
    # z = matrix(NA,nrow=n,ncol=k+3)
    # z[,k] =  b*e
    # z[,k+1] = g*(e^2-sigma2) + 2*e*f(Gnn,e) + s*e
    # z[,k+2] = h*(e^2-sigma2) + 2*e*f(Hnn,e)
    # z[,k+3] = e*e - rep(sigma2, n)
    # az=rbind(z,-an*colMeans(z))
    # # 计算EL值	
    # lam = lambdaChen(z)
    # el = 2*sum( log(1+t(lam)%*%t(z) ) )
    # if(el>cut) f1A2=f1A2+1
    # # 计算AEL值
    # alam=lambdaChen(az)
    # ael=2*sum( log(1+t(alam)%*%t(az) ) )  		
    # if(ael>cut) f2A2=f2A2+1
    
  }
  cat('样本个数为',n,'完成模拟',m,'次','\n')
  ff_elH0=append(ff_elH0,f1H0/nsim)
  ff_aelH0=append(ff_aelH0,f2H0/nsim)
  ff_elA1=append(ff_elA1,f1A1/nsim)
  ff_aelA1=append(ff_aelA1,f2A1/nsim)
  ff_elA2=append(ff_elA2,f1A2/nsim)
  ff_aelA2=append(ff_aelA2,f2A2/nsim)
  Sys.sleep(20)
}
ff = matrix(NA,nrow=length(size),ncol=6)
ff[,1]=ff_elH0
ff[,2]=ff_aelH0
ff[,3]=ff_elA1
ff[,4]=ff_aelA1
ff[,5]=ff_elA2
ff[,6]=ff_aelA2
rownames(ff) <- size^2
colnames(ff) <- c("EL_H0","AEL_H0","EL_A1","AEL_A1","EL_A2","AEL_A2")
cat('--------------------------------------','\n')
print(ff)


