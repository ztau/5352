gamma1 <- function(n, r, k, OHC, kOHC, SAlpha, SAlog, SBeta, SBlog, LLAlpha, Epsilon){
    
    OHC=max(1,OHC)
    Gamma=array(0,c(kOHC,n))
    if(OHC>1){
    for (t in (r+1):(r+OHC-1)){
        for (i in 1:k^(t-r)){
            Gamma[i,t]=Epsilon[[1]][i,,t]%*%array(1,c(k,1))
        }
    }
    }
    
    for (t in max((r+OHC),1):(n-1)){
        for (i in 1:k^OHC){
            Gamma[i,t]=Epsilon[[1]][i,,t]%*%array(1,c(k,1))
        }
    }
    
    for (i in 1:kOHC){
        if (SAlpha[i,n]!=0){
            Gamma[i,n]=exp(log(SAlpha[i,n])+ SAlog[1:n]%*%array(1,c(n,1)) + SBlog[n] - LLAlpha)
        }
    }
   return=list(Gamma) 
}