#===============================================
#             SARAR model t(df)
#===============================================

SARAR.model.t<-function(nsim,rou1,rou2,df,size){
  
  beta=3.5
  tol=1e-4
  a=0.95
  k=length(beta)
  cut = qchisq(a,k+3)
  
  zero = 0
  azero = 0
  ff_el = c()
  ff_ael = c()
  for (m in size){
    # SARAR模型
    n = m*m
    i = 1:n;Xn = i/(n+1)
    # Xn = rnorm(n)^2
    Wnb = cell2nb(m,m,type='queen')
    Ws = nb2listw(Wnb)
    Wn = listw2mat(Ws)
    # Wn = matrix(NA,nrow=n,ncol=n)
    # for(i in 1:n){
    #   for(j in 1:n){
    #     if(i==j) Wn[i,j]=0 else Wn[i,j]=1/(n-1)
    #   }
    # }
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
    f1 = 0
    f2 = 0
    for(i in 1:nsim){
      # cat('样本个数为',n,'正在模拟第 ',i,'次','\n')
      En = rt(n,df);sigma2=df/(df-2)
      e = En
      # 模拟Yi(程序运行不需要Yi值)
      # Yn = Ani%*%Xn%*%beta + Ani%*%Bni%*%En
      # 估计方程赋值
      z = matrix(NA,nrow=n,ncol=k+3)
      z[,k] =  b*e
      z[,k+1] = g*(e^2-sigma2) + 2*e*f(Gnn,e) + s*e
      z[,k+2] = h*(e^2-sigma2) + 2*e*f(Hnn,e)
      z[,k+3] = e*e - rep(sigma2, n)
      
      # 计算EL值	
      lam = lambdaChen(z)
      el = 2*sum( log(1+t(lam)%*%t(z) ) )
      if(el<=cut) f1=f1+1
      # aa=1+t(lam)%*%t(z)
      # glam=rowSums(t(z)/t(matrix(rep(aa,2),n,k+3)))
      # if(max(abs(glam))>tol) zero=zero+1
      
      
      # 计算AEL值
      az=rbind(z,-an*colMeans(z))	 		
      alam=lambdaChen(az)
      ael=2*sum( log(1+t(alam)%*%t(az) ) )  		
      if(ael<=cut) f2=f2+1
      # aa=1+t(alam)%*%t(az)
      # glam=rowSums(t(az)/t(matrix(rep(aa,2),n+1,k+3)))
      # if(max(abs(glam))>tol) azero=azero+1
    }
    # cat('样本个数为',n,'完成模拟',i,'次',zero,azero,'\n')
    # cat(paste0('t(',df,') ',n),' ',f1/nsim,f2/nsim,'\n')
    ff_el=append(ff_el,f1/nsim)
    ff_ael=append(ff_ael,f2/nsim)
    zero = 0
    azero = 0
    # Sys.sleep(10)
  }
  ff=matrix(NA,nrow=length(size),ncol=2)
  ff[,1]=ff_el
  ff[,2]=ff_ael
  rownames(ff) <- size^2
  colnames(ff) <- c("EL","AEL")
  print(ff)
}