SEM_v2<- function(fct, VECTPAR,lengthPhi,lengthTheta,k, target=NA, it, envir, lower= VECT_limits_min ,upper= VECT_limits_max){
    
    # SEM
    

    
    amplitude=upper-lower # amplitude of the acceptable values of the parameters (between the lower and upper bounds)
    

    
    # initialisation of some objects
    LLAlpha_vis_last=fct(VECTPAR)
    LL_last=LLAlpha_vis_last
    VECTPAR_elite=cbind(t(VECTPAR),0)
    ii=3     # iterations
    jj=0     # nbr of calls
    lengthVECTPAR=length(VECTPAR)
    changes<-array(0.2, c(1,lengthVECTPAR))
    VECThist=array(0,c(5000,lengthVECTPAR))
    LLAlphahist_up<- array(0,c(1,lengthVECTPAR))
    LLAlphahist_down<- array(0,c(1,lengthVECTPAR))
    Lvisible=array(0,c(1,it))
    Lvisible[ii]=-1
    Lvisible[ii-1]=-0.5
    nbr_fcalls=0
    direction=array("up",c(1,length(VECTPAR)))
    
    
    
    min_param_change=amplitude/100  # the minimum and the maximum of the step are calculated as fractions of the amplitude between the limits (parameter space) 
    max_param_change=amplitude/10
    
    
    while(nbr_fcalls<it){
        LLAlpha_vis_last<-LL_last
        
        
        VECTPAR_up=VECTPAR
        VECTPAR_down=VECTPAR
        
        for (m in 1:(lengthVECTPAR)){
            
            

            if(VECTPAR[m]<0.1 & m<=lengthTheta){
                VECTPAR[m]=lower[m]                    # introduce the "hard" limits that cannot be exceeded (in this case just the variance parameters cannot be negative)
            }                                          # all limits that are not specified here are considered as "floating" (not compulsary)
            
            
            VECTPAR_up[m]=VECTPAR[m]+max(changes[m],min_param_change[m])             # increase the parameters
            VECTPAR_down[m]=VECTPAR[m]-max(changes[m],min_param_change[m])           # diminish the parameters
            
            LL_up=fct(VECTPAR_up)     # logL values obtaned if we increase the parameters (for each parameter separately)
            jj=jj+1
            if(LL_up>LL_last){
                LL_last=LL_up
                VECTPAR[m]=VECTPAR_up[m]
                if(direction[m]=="up"){changes[m]=min(changes[m]*2,max_param_change)}else{changes[m]=max(changes[m]/2,min_param_change)} # if the last direction las the same increase the step, if not diminish it because we may have gone too far.
                direction[m]="up"
            } else{
                LL_down=fct(VECTPAR_down)       # logL values obtaned if we reduce the parameters (for each parameter separately)
                jj=jj+1
                if(LL_down>LL_last){
                    LL_last=LL_down
                    VECTPAR[m]=VECTPAR_down[m]
                    changes[m]=min(changes[m]*2,max_param_change) 
                    if(direction[m]!="up"){changes[m]=min(changes[m]*2,max_param_change)}else{changes[m]=max(changes[m]/2,min_param_change)} # if the last direction was the same increase the step, if not diminish it because we may have gone too far.
                    direction[m]="down"
                } else{
                    changes[m]<-max(changes[m]/3,min_param_change)
                }
            }
            
        }
        
        
        
        Lvisible[ii]=LL_last

        
        
        if(Lvisible[ii]==Lvisible[ii-1]){VECTPAR_elite=rbind(VECTPAR_elite,cbind(t(VECTPAR),ii)) # save the last VECTPAR (AND the nbr of its corresponding function call) before the "jitter" to escape the local optimum
        VECTPAR=VECTPAR+rnorm(length(VECTPAR))*0.3*VECTPAR
        VECTPAR[1:lengthTheta]=abs(VECTPAR[1:lengthTheta])
        LL_last=fct(VECTPAR) 
        jj=jj+1 }
        
        ii=ii+1
        nbr_fcalls=jj
        if ((ii/10)%%1==0){    # save every 10 parameter proposition (just for information how did we explore the parameter space)
            VECThist[ii/10,]=VECTPAR
        }

        
        
    }
    envir$LL<-max(Lvisible[Lvisible<(-0.5)])
    envir$param<-if(length(as.array(VECTPAR_elite[which(VECTPAR_elite[,length(VECTPAR)+1]==  max(which(Lvisible==max(Lvisible[Lvisible<(-0.5)]))) ), 1:length(VECTPAR) ] )  )==0){as.array(VECTPAR)} else{as.array(VECTPAR_elite[which(VECTPAR_elite[,length(VECTPAR)+1]==  max(which(Lvisible==max(Lvisible[Lvisible<(-0.5)]))) ), 1:length(VECTPAR) ] ) }
    envir$iter<-ii-3
    envir$fct_calls<- nbr_fcalls
    
}