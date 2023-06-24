#==================================================#
#                     SARAR model EL
#==================================================#
# Remark: 
# 2023/06/24
# 一些包需要提前下载，见如下命令：
# install.packages("raster",type='binary')
# install.packages("spData",type='binary')
# install.packages("terra",type='binary')
# library('spData')
# 该程序为 Tables 1-7 的核心代码 
# Tables 1-7 皆为该段代码的不同变体
# 了解该段代码即掌握模拟的核心技术
rm(list = ls())
#===================================================#
#                    求解lambda数值解
#===================================================#
lambdaChen<-function(u){
  n=dim(u)[1]
  df=dim(u)[2]
  M=rep(0,df)
  k=0
  gama=1
  tol=1e-11
  dif=1
  R<-function(lam){R0=sum(log(1+t(lam)%*%t(u)));return(R0)}
  R1=rep(0,df)
  R2=R1%*%t(R1)
  
  while(dif>tol && k<=300){
    # 计算R1、R2
    aa=1+t(M)%*%t(u)
    for(i in 1:df){
      R1[i]=sum(t(u[,i])/aa)
      for(j in 1:df){
        R2[i,j]=-sum(u[,i]*u[,j]/aa^2)
      }
    }
    delta=-solve(R2)%*%R1
    dif=c(sqrt(t(delta)%*%delta))
    sigma=gama*delta
    while(min(1+t(M+sigma)%*%t(u))<=0){
      gama=gama/2
      sigma=gama*delta
    }
    # print(k)
    M=M+sigma
    gama=1/sqrt(k+1)
    k=k+1
  }
  # cat(k,'\n')
  return(M)
} 
#===================================================#
#                    加载必要的包
#===================================================#
library('sp')
library('terra')
library('spdep')
#===================================================#
#                    修改模拟参数
#===================================================#
nsim = 100
beta = 3.5
rou1 = 0.85
rou2 = 0.15
tol=1e-6
a=0.95
k=1
cut = qchisq(a,k+3)
size = c(7,10,13,40)
#===================================================#
#                    模拟正式运行
#===================================================#
ff_el = c()
for (m in size){
  # SARAR模型
  n = m*m
  # i = 1:n;Xn = i/(n+1)
  Xn = rnorm(n)^2
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
  # 启动模拟
  f1 = 0
  for(m in 1:nsim){
    # cat('样本个数为',n,'正在模拟第 ',m,'次','\n')
    En = rnorm(n);sigma2=1
    # En = rt(n,5);sigma2=5/3
    # En = rchisq(n,4)-4;sigma2=8
    e = En
    
    z = matrix(NA,nrow=n,ncol=k+3)
    z[,k] =  b*e
    z[,k+1] = g*(e^2-sigma2) + 2*e*f(Gnn,e) + s*e
    z[,k+2] = h*(e^2-sigma2) + 2*e*f(Hnn,e)
    z[,k+3] = e*e - rep(sigma2, n)
    
    # 计算EL值
    lam = lambdaChen(z)
    el = 2*sum( log(1+t(lam)%*%t(z) ) )
    if(el<=cut) f1=f1+1
  }
  cat('样本个数为',n,'覆盖率',f1/nsim,'\n')
}

