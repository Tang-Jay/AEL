#===============================================
#           Chen lm
#===============================================
rm(list = ls()) 
source('GlambdaChen.R')

nsim=2000
a=0.95
df=2
cut=qchisq(a,df)
tol=1e-6
zero=0
azero=0
ff_el = c()
ff_ael = c()
size=c(10,100)
for(n in size){
  f=0
  f2=0
  if(1>log(n)/2) an=1 else an=log(n)/2
  
  for(m in 1:nsim){
    i = 1:n;x0 = i/(n+1)
    x1 = c(x0,x0^2)               						
    x = matrix(x1,nrow=n, ncol=2)  			
    e = runif(n,-0.5,0.5)          						
    e = t(t(e))
    e1=matrix(c(e,e),nrow=n, ncol=2) 		
    
    z=x*e1                       		
    lam=lambdaChen(z)                						
    aa=1+t(lam)%*%t(z)
    glam=rowSums(t(z)/t(matrix(rep(aa,2),n,df)))
    if(max(abs(glam))>tol)zero=zero+1
    el=sum(2*log(1+x%*%lam*e))
    if(el<=cut) f=f+1
    
    az=rbind(z,-an*colMeans(z))
    alam=lambdaChen(az)
    aa=1+t(alam)%*%t(az)
    glam=rowSums(t(az)/t(matrix(rep(aa,2),n+1,df)))
    ael=2*sum(log(1+t(alam)%*%t(az)))
    if(max(abs(glam))>tol)azero=azero+1
    if(ael<=cut) f2=f2+1
  }
  cat(n,'个样本',m,'次模拟','\n')
  ff_el=append(ff_el,f/nsim)
  ff_ael=append(ff_ael,f2/nsim)
}

cat('模拟',nsim,'次Glambda>',tol,'的次数',zero,azero,'\n')
cat('ff_el ', ff_el ,'\n')
cat('ff_ael ', ff_ael ,'\n')

# nums=size
# data=t(data.frame( EL = ff_el, AEL = ff_ael))
# colnames(data)=nums
# par(font = 2)
# barplot(data,col = c("black", "red"), beside = T, ylim = c(0, 1.2), font = 2)
# legend('top', legend=c('el','ael'), fill = c("black", "red"),bty='n',cex=1,horiz=T)
# axis(side = 2, lwd = 3, font = 2)
# abline(h=0.95)
# title(main = list("Linear model", font = 2))
