#===============================================
#             SARAR H1 
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
nsim = 100
beta = 3.5
rou1 = 0.85
rou2 = 0.15

tol=1e-7
a=0.95
k=1
cut = qchisq(a,k+3)
size = c(3,4,5,7,10,13)

zero = 0
azero = 0
ff_elH1 = c()
ff_aelH1 = c()
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
  f1H1 = 0
  f2H1 = 0
  for(m in 1:nsim){
    # cat('样本个数为',n,'正在模拟第 ',m,'次','\n')
    En = rnorm(n)
    sigma2=2
    e = En
    # 模拟Yi(程序运行不需要Yi值)
    # Yn = Ani%*%Xn%*%beta + Ani%*%Bni%*%En
    # 估计方程赋值
    z = matrix(NA,nrow=n,ncol=k+3)
    z[,k] =  b*e
    z[,k+1] = g*(e^2-sigma2) + 2*e*f(Gnn,e) + s*e
    z[,k+2] = h*(e^2-sigma2) + 2*e*f(Hnn,e)
    z[,k+3] = e*e - rep(sigma2, n)
    az=rbind(z,-an*colMeans(z))
    
    # 计算EL值	
    lam = lambdaChen(z)
    el = 2*sum( log(1+t(lam)%*%t(z) ) )
    if(el>cut) f1H1=f1H1+1
    
    # 计算AEL值
    alam=lambdaChen(az)
    ael=2*sum( log(1+t(alam)%*%t(az) ) )  		
    if(ael>cut) f2H1=f2H1+1
  }
  cat('样本个数为',n,'完成模拟',m,'次',zero,azero,'\n')
  ff_elH1=append(ff_elH1,f1H1/nsim)
  ff_aelH1=append(ff_aelH1,f2H1/nsim)
  zero = 0
  azero = 0
}
ff = matrix(NA,nrow=length(size),ncol=2)
ff[,1]=ff_elH1
ff[,2]=ff_aelH1
rownames(ff) <- size^2
colnames(ff) <- c("EL","AEL")
cat('------------H1------------------------','\n')
print(ff)

# cat('ff_elH1 ', ff_elH1 ,'\n')
# cat('ff_aelH1 ', ff_aelH1 ,'\n')

# 可视化
# nums=size^2
# data=t(data.frame( EL = ff_elH1, AEL = ff_aelH1))
# colnames(data)=nums
# par(font = 2, lwd = 2)
# barplot(data, col = c("black", "red"), beside = T, ylim = c(0, 1.3), font = 2, legend.text=c('el','ael'))
# axis(side = 2, lwd = 3, font = 2)
# abline(h=0.95)
# title(main = list("SARAR model", font = 2))
