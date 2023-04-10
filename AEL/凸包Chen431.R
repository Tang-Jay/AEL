#===============================================
#                 Plot_ConvexHull
#===============================================
Plot_ConVexHull<-function(u,au){
  library(geometry)
  g1 <- u[,1]
  g2 <- u[,2]
  
  ag1 <- au[,1]
  ag2 <- au[,2]
  
  ConVexHull<-convhulln(u,"FA")
  aConVexHull<-convhulln(au,"FA")
  plot_ConvexHull<-function(xcoord, ycoord, lcolor){
    hpts <- chull(x = xcoord, y = ycoord)
    hpts <- c(hpts, hpts[1])
    lines(xcoord[hpts], ycoord[hpts], col = lcolor)
  } 

  axrange <- range(c(ag1))
  ayrange <- range(c(ag2))
  
  par(mfrow=c(1,2))
  plot(g1, g2, type="p", pch=1, col="black", xlim=c(axrange), ylim=c(ayrange),main='EL',xlab='g1',ylab='g2')
  points(0,0,pch=20)
  plot_ConvexHull(xcoord = g1, ycoord = g2, lcolor = "black")
  
  plot(ag1, ag2, type="p", pch=1, col="black", xlim=c(axrange), ylim=c(ayrange),main='AEL',xlab='g1',ylab='g2')
  points(0,0,pch=20)
  points(au[n+1,][1],au[n+1,][2],col='4')
  plot_ConvexHull(xcoord = ag1, ycoord = ag2, lcolor = "black")
}
#===============================================
#                 Main
#===============================================
n=50 # 样本个数
a=0.95
df=2
cut=qchisq(a,df)
theta=c(4,4)
mu=rep(1,n)%*%t(theta)
if(1>log(n)/2) an=1 else an=log(n)/2

x0<-rnorm(n,2)
y0<-rnorm(n,2)
x1<-c(x0,y0)
x=matrix(x1,nrow=n, ncol=2)

u=x-mu # EL估计方程
au=rbind(u,-an*colMeans(u)) # AEL估计方程

Plot_ConVexHull(u,au)