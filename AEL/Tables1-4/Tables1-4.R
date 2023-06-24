rm(list = ls())
source('GlambdaChen.R')
source('N.R')
source('t.R')
source('Chi.R')
library('sp')
library('sf')
library('spData')
library('terra')
library('spdep')

nsim = 5000
size = c(3,5,7)
#===============================================
#                    Table 1
#===============================================
cat('----------------------------------','\n')
for(rou1 in c(-0.85, 0.85 )){
  for(rou2 in c(-0.15, 0.15)){
    cat('标准正态',rou1,rou2,'\n')
    SARAR.model.N(nsim,rou1,rou2,1,size)
    cat('\n')
  }
}

#===============================================
#                    Table 2
#===============================================
cat('----------------------------------','\n')
for(rou1 in c(-0.85, 0.85 )){
  for(rou2 in c(-0.15, 0.15)){
    cat('薄尾正态',rou1,rou2,'\n')
    SARAR.model.N(nsim,rou1,rou2,0.75,size)
    cat('\n')
  }
}

#===============================================
#                    Table 3
#===============================================
cat('----------------------------------','\n')
for(rou1 in c(-0.85, 0.85)){
  for(rou2 in c(-0.15, 0.15)){
    cat('厚尾t5',rou1,rou2,'\n')
    SARAR.model.t(nsim,rou1,rou2,5,size)
    cat('\n')
  }
}

#===============================================
#                    Table 4
#===============================================
cat('----------------------------------','\n')
for(rou1 in c(-0.85, 0.85)){
  for(rou2 in c(-0.15, 0.15)){
    cat('Chi4',rou1,rou2,'\n')
    SARAR.model.chi(nsim,rou1,rou2,4,size)
    cat('\n')
  }
}



