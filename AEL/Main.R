# rm(list = ls()) 	
source('GlambdaChen.R')
source('SARAR model N.R')
source('SARAR model t.R')
source('SARAR model chi.R')
library('sp')
library('sf')
library('spData')
library('terra')
library('spdep')

# SARAR.model.N(1,0.85, 0.15, 1)

# cat('----------------------------------','\n')
# for(rou1 in c(-0.85, 0.85 )){
#   for(rou2 in c(-0.15, 0.15)){
#     cat('标准正态',rou1,rou2,'\n')
#     SARAR.model.N(nsim=5000,rou1,rou2,1)
#     cat('\n')
#   }
# }
# 
# cat('----------------------------------','\n')
# for(rou1 in c(-0.85, 0.85 )){
#   for(rou2 in c(-0.15, 0.15)){
#     cat('薄尾正态',rou1,rou2,'\n')
#     SARAR.model.N(nsim=5000,rou1,rou2,0.75)
#     cat('\n')
#   }
# }
# 
# cat('----------------------------------','\n')
# for(rou1 in c(-0.85, 0.85)){
#   for(rou2 in c(-0.15, 0.15)){
#     cat('厚尾t5',rou1,rou2,'\n')
#     SARAR.model.t(nsim=5000,rou1,rou2,5)
#     cat('\n')
#   }
# }

cat('----------------------------------','\n')
for(rou1 in c(-0.85, 0.85)){
  for(rou2 in c(-0.15, 0.15)){
    cat('chi4',rou1,rou2,'\n')
    SARAR.model.chi(nsim=5000,rou1,rou2,4)
    cat('\n')
  }
}



