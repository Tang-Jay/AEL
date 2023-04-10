#===============================================
#             SARAR model W897
#===============================================

rm(list = ls()) 	
source('GlambdaChen.R')
library('sp')
library('terra')
library('rgdal')
library('spdep')

nsim = 100
tol = 1e-4
a = 0.95
beta = 3.5
k = 1
cut = qchisq(a,k+3)
m = 897
n = m*m
i = 1:n;Xn = i/(n+1)
chi.poly <- rgdal::readOGR('foreclosures.shp')
list.queen=poly2nb(chi.poly,queen=TRUE)
Ws=nb2listw(list.queen,style="W",zero.policy=TRUE)
Wn = listw2mat(Ws) 
Mn = Wn
In = diag(n)
if(1>log(n)/2) an=1 else an=log(n)/2
#====================N===========================
ff_el = c()
ff_ael = c()
rou_1 = c()
rou_2 = c()
# 启动模拟
cat('N','\n')
for(rou1 in c(-0.85, 0.85 )){
  for(rou2 in c(-0.15, 0.15)){
    cat(rou1,rou2,'\n')
    f1 = 0
    f2 = 0
    zero = 0
    azero = 0
    for(m in 1:nsim){
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
      En = rnorm(n);sigma2=1
      # df=5;En = rt(n,df);sigma2=df/(df-2)
      # df=4;En = rchisq(n,df)-df;sigma2=2*df
      e = En

      # 模拟Yi(程序运行不需要Yi值)
      # Yn = Ani%*%Xn%*%beta + Ani%*%Bni%*%En
      # 估计方程赋值
      z = matrix(NA,nrow=n,ncol=k+3)
      z[,k  ] = b*e
      z[,k+1] = g*(e^2-sigma2) + 2*e*f(Gnn,e) + s*e
      z[,k+2] = h*(e^2-sigma2) + 2*e*f(Hnn,e)
      z[,k+3] = e*e - rep(sigma2, n)

      # 计算EL值
      lam = lambdaChen(z)
      el = 2*sum( log(1+t(lam)%*%t(z)) )
      if(el<=cut) f1=f1+1

      # 计算AEL值
      az=rbind(z,-an*colMeans(z))
      alam=lambdaChen(az)
      ael=2*sum( log(1+t(alam)%*%t(az) ) )
      if(ael<=cut) f2=f2+1
    }
    rou_1 = append(rou_1,rou1)
    rou_2 = append(rou_2,rou2)
    ff_el = append(ff_el,f1/nsim)
    ff_ael = append(ff_ael,f2/nsim)
  }
}
ff = matrix(NA,nrow=4,ncol=2)
ff[,1]=ff_el
ff[,2]=ff_ael
rownames(ff) <- c('[-0.85, -0.15]','[-0.85,  0.15]','[0.85,  -0.15]','[0.85,   0.15]')
colnames(ff) <- c("EL","AEL")
print(ff)










