#===============================================
#                    Table 5                   #
#===============================================
# Remark:
# 这是用来统计两种方法运行时间的程序

# install.packages("raster",type='binary')
# install.packages("spData",type='binary')
# install.packages("terra",type='binary')
# library('spData')

rm(list = ls())
library('sp')
library('sf')
library('spData')
library('terra')
library('spdep')

#===============================================
#                 计算lambda程序
#===============================================
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
#                      修改模拟参数
#===================================================#
nsim = 1000
beta = 3.5
rou1 = 0.85
rou2 = 0.15
tol = 1e-6
a = 0.95
k = length(beta)
cut = qchisq(a,k+3)
size = c(3,7,10,13,16)
#===================================================#
#                 EL、AEL计时程序启动
#===================================================#
zero = 0
azero = 0
tt1=0
tt2=0
ff_el = c()
ff_time = c()
ff_ael = c()
ff_atime = c()
T1<-lubridate::now()
for (m in size){
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
  Gn = Bn%*%Wn%*%Ani%*%Bni
  Gnn = 1/2*(Gn + t(Gn))
  Hn = Mn%*%Bni
  Hnn = 1/2*(Hn + t(Hn))
  b = t(Xn)%*%t(Bn)
  g = diag(Gnn)
  h = diag(Hnn)
  s = Bn%*%Wn%*%Ani%*%Xn%*%beta
  f<-function(Matrix,Vector){
    irow = nrow(Matrix)
    v = c(0)
    for(i in 2:irow){
      v[i] = Matrix[i,][1:(i-1)]%*%Vector[1:i-1]
    }
    return(v)
  }
  if(1>log(n)/2) an=1 else an=log(n)/2
  # 启动模拟
  f1 = 0
  f2 = 0
  for(m in 1:nsim){
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
    t1<-lubridate::now()
    lam = lambdaChen(z)
    el = 2*sum( log(1+t(lam)%*%t(z) ) )
    t2<-lubridate::now()
    if(el<=cut) f1=f1+1
    tt1=tt1+t2-t1
    aa=1+t(lam)%*%t(z)
    glam=rowSums(t(z)/t(matrix(rep(aa,2),n,k+3)))
    if(max(abs(glam))>tol) zero=zero+1
    
    # 计算AEL值
    az=rbind(z,-an*colMeans(z))
    at1<-lubridate::now()
    alam=lambdaChen(az)
    ael=2*sum( log(1+t(alam)%*%t(az) ) )
    at2<-lubridate::now()
    if(ael<=cut) f2=f2+1
    tt2=tt2+at2-at1
    aa=1+t(alam)%*%t(az)
    aglam=rowSums(t(az)/t(matrix(rep(aa,2),n+1,k+3)))
    if(max(abs(aglam))>tol) azero=azero+1
  }
  # cat('样本个数为',n,'完成模拟',m,'次',tol,zero,azero,'\n')
  ff_el=append(ff_el,f1/nsim)
  ff_time=append(ff_time,tt1)
  ff_ael=append(ff_ael,f2/nsim)
  ff_atime=append(ff_atime,tt2)
  tt1=0
  tt2=0
  zero=0
  azero=0
}
T2<-lubridate::now()
TT<-T2-T1
#===================================================#
#                   模拟结果可视化
#===================================================#
ff = matrix(NA,nrow=length(size),ncol=4)
ff[,1]=ff_el
ff[,2]=ff_ael
ff[,3]=ff_time
ff[,4]=ff_atime
rownames(ff) <- size^2
colnames(ff) <- c("EL","AEL","times","atimes")

print(TT)
cat(sum(ff_time),"+",sum(ff_atime),'\n')
print(ff)


