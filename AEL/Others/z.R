# 调试 
# Y(n=5,beta=3.5,rou1=0.85,rou2=1.5)
n=5;beta=3.5;rou1=0.85;rou2=1.5
sigma2=1
k=length(beta)
i = 1:n
Xn = i/(n+1)
En = rnorm(n)
In = diag(n)
Wn = matrix(NA,nrow=n,ncol=n)
for(i in 1:n){
  for(j in 1:n){
    if(i==j) Wn[i,j]=0 else Wn[i,j]=1/(n-1)
  }
}
Mn = Wn
An = In - rou1*Wn
Bn = In - rou2*Mn
# det(An)==0 || det(Bn)==0)
Ani = solve(An)
Bni = solve(Bn)

# 模拟Yi
Yn = Ani%*%Xn%*%beta + Ani%*%Bni%*%En
# 估计方程准备
Gn = Bn%*%Wn%*%Ani%*%Bni
Gnn = 1/2*(Gn + t(Gn))
Hn = Mn%*%Bni
Hnn = 1/2*(Hn + t(Hn))
# 简化符号
b=t(Xn)%*%t(Bn)
e=En
g=diag(Gnn)
h=diag(Hnn)
s=Bn%*%Wn%*%Ani%*%Xn%*%beta
f<-function(Matrix,Vector){
  irow = nrow(Matrix)
  v =c(0)
  for(i in 2:irow){
    v[i] =Matrix[i,][1:(i-1)]%*%Vector[1:i-1]
  }
  return(v)
}
# 估计方程赋值
z = matrix(NA,nrow=n,ncol=k+3)
z[,k] =  b*e
z[,k+1] = g*(e^2-sigma2) + 2*e*f(Gnn,e) + s*e
z[,k+2] = h*(e^2-sigma2) + 2*e*f(Hnn,e)
z[,k+3] = t(e)*e - t(rep(sigma2, n))
cat('-------------------------------','\n')
print(z)

# 估计方程 唐洁-错版（EE第二行和第三行写错的估计方程）
EE = matrix(NA,nrow=k+3,ncol=n)
EE[k,] =  b*e
EE[k+1,] = t(s)*En + t(En)%*%Gnn*En - sigma2*t(diag(Gnn))
EE[k+2,] = t(En)%*%Hnn*En - sigma2*t(diag(Hnn))
EE[k+3,] = e*e - rep(sigma2, n)
cat('-------------------------------','\n')
print(t(EE))


# # 得分函数
# tr<-function(X){X=sum(diag(X));return(X)} 	
# g1 = matrix(NA,nrow=k+3,ncol=1)
# g1[k,] =  t(Xn)%*%t(Bn)%*%En
# g1[k+1,] = t(Bn%*%Wn%*%Ani%*%Xn%*%beta)%*%En + t(En)%*%Gnn%*%En - sigma2*tr(Gnn)
# g1[k+2,] = t(En)%*%Hnn%*%En - sigma2*tr(Hnn)
# g1[k+3,] = t(En)%*%En - n*sigma2
# # 检测
# if( all( rowSums(g) == c(g1))  ) print('True')

# # 估计方程 师姐-原版向量版 写法一
# z1 = c()
# for(i in 1:n){
# z1[i] = b[i]%*%e[i]
# }
# z2 = c()
# z2[1] = Gnn[1,1]*e[1]^2-Gnn[1,1]*sigma2+s[1]*e[1]
# for(i in 2:n){
# z2[i] = Gnn[i,i]*(e[i]^2-sigma2)+2*e[i]*(Gnn[i,][1:(i-1)]%*%e[1:i-1])+s[i]*e[i]
# }
# z3 = c()
# z3[1] = Hnn[1,1]*e[1]^2-Hnn[1,1]*sigma2
# for(i in 2:n){
# z3[i] = Hnn[i,i]*(e[i]^2-sigma2)+2*e[i]*(Hn[i,][1:(i-1)]%*%e[1:i-1])
# }
# z4 = c()
# for(i in 1:n){
# z4[i] = e[i]^2 - sigma2
# }
# cat('-------------------------------','\n')
# print(z1)
# print(z2)
# print(z3)
# print(z4)


# # 估计方程 师姐-改写矩阵版 写法二
# z1=z
# z = matrix(NA,nrow=n,ncol=k+3)
# for(i in 1:n){
# z[i,1] = b[i]%*%e[i]
# }
# z[1,2] = Gnn[1,1]*e[1]^2-Gnn[1,1]*sigma2+s[1]*e[1]
# for(i in 2:n){
# z[i,2] = Gnn[i,i]*(e[i]^2-sigma2)+2*e[i]*(Gnn[i,][1:(i-1)]%*%e[1:i-1])+s[i]*e[i]
# }
# z[1,3] = Hnn[1,1]*e[1]^2-Hnn[1,1]*sigma2
# for(i in 2:n){
# z[i,3] = Hnn[i,i]*(e[i]^2-sigma2)+2*e[i]*(Hn[i,][1:(i-1)]%*%e[1:i-1])
# }
# for(i in 1:n){
# z[i,4] = e[i]^2 - sigma2
# }
# # 检验
# z2=z
# print(z1==z2)

# # 估计方程 唐洁-错版（EE第二行和第三行写错的估计方程）
# EE = matrix(NA,nrow=k+3,ncol=n)
# #colnames(g) = paste('g',1:n,sep='')
# #rownames(g) = paste(1:(k+3),sep='')
# EE[k,] =  t(Xn)%*%t(Bn)*En
# EE[k+1,] = t(Bn%*%Wn%*%Ani%*%Xn%*%beta)*En + t(En)%*%Gnn*En - sigma2*t(diag(Gnn))
# EE[k+2,] = t(En)%*%Hnn*En - sigma2*t(diag(Hnn))
# EE[k+3,] = t(En)*En - t(rep(sigma2, n))

# # 估计方程 唐洁-正确版
# z = matrix(NA,nrow=n,ncol=k+3)
# z[,k] =  t(Xn)%*%t(Bn)*En
# z[,k+1] = g*(e^2-sigma2) + 2*En*f(Gnn,e) + s*e
# z[,k+2] = h*(e^2-sigma2) + 2*En*f(Hnn,e)
# z[,k+3] = t(En)*En - t(rep(sigma2, n))