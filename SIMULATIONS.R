rm(list=ls())
set.seed(111)
setwd("E:/git")

load("SEQ15")#15,50,100,200
longueur=15

# # set.seed(111) for longueur=15
# # set.seed(222) for longueur=50
# # set.seed(333) for longueur=100
# # set.seed(444) for longueur=200
# 
 iter11=1000
 LL1=0
 fct_calls1=0
 parameters<-array(0,c(1,6))
 AetPi<-array(0,c(3,2))
 ptm <- proc.time()
 datahist=array(0,c(2,longueur+1))
 
#     
#     
#     comp1<-function(prev_obs){
#         mu=1+0.2*prev_obs
#         new_obs=mu+rnorm(1)*0.5
#         return(new_obs)
#     }
#     comp2<-function(prev_obs){
#         mu=3+0.6*prev_obs
#         new_obs=mu+rnorm(1)*2
#         return(new_obs)
#     }
#     
#     
#     generating_process<-function(prev_obs,prev_state){
#         if(prev_state==1){
#             if(runif(1)<0.75){
#                 new_state<-1
#                 new_obs<-comp1(prev_obs)
#             } else{ new_state<-2
#             new_obs<-comp2(prev_obs) }
#         }else if(prev_state==2){
#             if(runif(1)<0.6){
#                 new_state<-2
#                 new_obs<-comp2(prev_obs)
#             } else{ new_state<-1
#             new_obs<-comp1(prev_obs) }
#         }
#         return(c(new_state,new_obs))
#     }
#     
#     
#     SEQ=array(0,c(2,longueur))
#     
# for(i in 1:iter11){
#     
#     seq=array(0,c(2,longueur))   
#     
#     # First value
#     if(runif(1)<0.75){
#         seq[1,1]<-1
#         seq[2,1]=comp1(0)
#     } else{ seq[1,1]<-2
#     seq[2,1]=comp2(0) }
#     
#     for(m in 2:longueur){
#         seq[1:2,m]=generating_process(seq[2,m-1],seq[1,m-1])
#     }
#     
#     seq=round(seq,3)
#     
#     SEQ=rbind(SEQ,seq)
#     
#     }
# 
# SEQ=SEQ[3:nrow(SEQ),]
# 
# save(SEQ,ascii=FALSE,file="SEQ200")


for(i in 1:iter11){
    
   
    
    Data=SEQ[2*i,,drop=FALSE]
    
    
    
    Param_hist<-array(0,c(1,7))
    Vit_hist<-array(0,c(1,1034))
    LL_hist=0
    total=2
    # pb <- tkProgressBar(title = "ça prends beaucoup trop de temps non?", min = 0,
    #                     max = total, width = 500)
    source("OPTIMISATION1_universal_NEW_Adiag_covminmax.R")
    
    envir<-new.env()
    
    individuals<-nrow(Data)
    
    q=0
    p=1
    k=2
    OHC=1
    Adiag=0 # 1 if A is diagonal, 0 otherwise
    if(Adiag==1){OHC=1} # si A diagonale -> OHC=1
    
    P<-array(p,dim=c(k,1))   
    Q<-array(q,dim=c(k,1))
    
    
    
    X_cov=array(0)
    Phi_cov_constr=array(0,c(k,nrow(X_cov)))   # constraints on Phi_cov   Phi_cov_constr=array(1,c(k,nrow(X_cov)))
    phi_cov_places=which(Phi_cov_constr!=0)
    
    possibleError <- tryCatch(
        OPTIMISATION1_universal_NEW_Adiag_covminmax(Data,X_cov,Phi_cov,phi_cov_places,k,OHC, P,Q,
                                                    procedure=c("SEMP_v2"),it=500,envir,Adiag,covmin=c(-6,-8,-2,-6,-5),covmax=c(5,5,15,6,8),picov),
        error=function(e) e
    )
    
    if(inherits(possibleError, "error")) next
    
    LL<-envir$LL
    Param<-envir$param
    fct_calls<-envir$fct_calls
    A<-envir$A
    Pi<-envir$Pi
    #proc.time() - ptm
    
    Param
    LL
    A
    Pi
    fct_calls
    
    Phi=array(Param[(k+1):(k+k*(p+1))],dim=c(k,p+1))
    Theta=array(Param[1:k],c(k,1))
    
    LL1=rbind(LL1,LL)
    fct_calls1=rbind(fct_calls1,fct_calls)
    parameters<-rbind(parameters,Param)
    AetPi<-rbind(AetPi,A,Pi)
    datahist=rbind(datahist,c(i,Data))
}# 
proc.time() - ptm
# 
Results <- setClass(
    # Name of the class
    "Results",
    # Define the slots
    slots = c(
        parameters = "matrix",
        AetPi = "matrix",
        LL1 = "matrix",
        fct_calls1 = "matrix",
        datahist = "matrix"
    ),
    
    # Default values for the slots
    prototype=list(
        parameters = matrix(0),
        AetPi = matrix(0),
        LL1 = matrix(0),
        fct_calls1 = matrix(0),
        datahist = matrix(0)
    )
)

setwd("D:/Article Optimization")
output<-Results(parameters =parameters,AetPi =AetPi,LL1 = LL1,fct_calls1 =fct_calls1,datahist =datahist)
#save(output,ascii=FALSE,file="DERNIER_SEMP_v2_500it_t100")
