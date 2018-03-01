#'Maximization procedure "SEM_v2"
#'
#'@param fct function to be maximized (logLikelihood)
#'@param VECTPAR vector of visiible parameters Phi and Theta
#'@return logL, parameters, iterations and number of function calls


SEM_v2<- function(fct, VECTPAR,lengthPhi,lengthTheta,k, target=NA, it, envir, lower= VECT_limits_min ,upper= VECT_limits_max,mutheta,sigtheta,muphi,sigphi,initialVECTPAR=NULL){

    # SEM fct=LL_10


    amplitude=upper-lower # amplitude of the acceptable values of the parameters (between the lower and upper bounds)

    min_param_change=amplitude/300  # the minimum and the maximum of the step are calculated as fractions of the amplitude between the limits (parameter space)
    max_param_change=amplitude/20

    # initialisation of some objects
    LLAlpha_vis_last=fct(VECTPAR)
    LL_last=LLAlpha_vis_last
    VECTPAR_elite=cbind(t(VECTPAR),0)
    ii=3     # iterations
    jj=0     # nbr of calls
    lengthVECTPAR=length(VECTPAR)
    changes<-array(min_param_change, c(1,length(min_param_change)))
    VECThist=array(0,c(5000,lengthVECTPAR))
    LLAlphahist_up<- array(0,c(1,lengthVECTPAR))
    LLAlphahist_down<- array(0,c(1,lengthVECTPAR))
    Lvisible=array(0,c(1,it))
    Lvisible[ii]=-1
    Lvisible[ii-1]=-0.5
    nbr_fcalls=0
    direction=array("up",c(1,length(VECTPAR)))



    while(nbr_fcalls<it){
        LLAlpha_vis_last<-LL_last




        for (m in 1:length(VECTPAR)){

            VECTPAR_up=VECTPAR
            VECTPAR_down=VECTPAR

            if(VECTPAR[m]<0.1 & m<=lengthTheta){
                VECTPAR[m]=lower[m]                    # introduce the "hard" limits that cannot be exceeded (in this case just the variance parameters cannot be negative)
            }                                          # all limits that are not specified here are considered as "floating" (not compulsary)


            VECTPAR_up[m]=VECTPAR[m]+max(changes[m],min_param_change[m])             # increase the parameters
            VECTPAR_down[m]=VECTPAR[m]-max(changes[m],min_param_change[m])           # diminish the parameters

            LL_up=fct(VECTPAR_up)     # logL values obtaned if we increase the parameters (for each parameter separately)
            jj=jj+1
            if(LL_up>LL_last){
                LL_last<<-LL_up
                VECTPAR[m]=VECTPAR_up[m]
                if(direction[length(direction)]=="up"){changes[m]=min(changes[m]*2,max_param_change[m])}else{changes[m]=max(changes[m]/3,min_param_change[m])} # if the last direction las the same increase the step, if not diminish it because we may have gone too far.
                direction[length(direction)]="up"
            } else{
                LL_down=fct(VECTPAR_down)       # logL values obtaned if we reduce the parameters (for each parameter separately)
                jj=jj+1
                if(LL_down>LL_last){
                    LL_last<<-LL_down
                    VECTPAR[m]=VECTPAR_down[m]
                    changes[m]=min(changes[m]*2,max_param_change[m])
                    if(direction[length(direction)]!="up"){changes[m]=min(changes[m]*2,max_param_change[m])}else{changes[m]=max(changes[m]/3,min_param_change[m])} # if the last direction was the same increase the step, if not diminish it because we may have gone too far.
                    direction[length(direction)]="down"
                } else{
                    changes[m]<-max(changes[m]/3,min_param_change[m])
                }
            }

        }

        # print(VECTPAR)
        # print(LL_last)
        # print(nbr_fcalls)

        Lvisible[ii]=LL_last



        if(Lvisible[ii]==Lvisible[ii-2]){VECTPAR_elite=rbind(VECTPAR_elite,cbind(t(VECTPAR),ii)) # save the last VECTPAR (AND the nbr of its corresponding function call) before the "jitter" to escape the local optimum
        if(any(is.na(lower))){lower[which(is.na(lower))]=upper[which(is.na(lower))]*0.5}
        if(any(is.na(upper))){upper[which(is.na(upper))]=abs(lower[which(is.na(upper))])*0.5}
        #VECTPAR=c(Theta,Phi,Phicov)
        if(!is.null(initialVECTPAR)){VECTPAR=initialVECTPAR+rnorm(length(initialVECTPAR),0,0.02)}  # if we have a best initial solution, explore around it
        else{VECTPAR =c(matrix(abs(rnorm(lengthTheta,mutheta,sigtheta)), c(1,lengthTheta)),         matrix(cbind(matrix(rnorm(k,muphi,sigphi),1,k),matrix(runif(lengthPhi-k,-0.2,0.3),1,lengthPhi-k) ),c(1,lengthPhi)),     array(runif(length(VECTPAR)-lengthPhi-lengthTheta,min=-0.5,max=1), c(1,(length(VECTPAR)-lengthPhi-lengthTheta))))}
        #VECTPAR=runif(length(VECTPAR),lower,upper)#    VECTPAR+rnorm(length(VECTPAR))*0.1*VECTPAR#
        VECTPAR[1:lengthTheta]=abs(VECTPAR[1:lengthTheta])
        VECTPAR[which(VECTPAR>upper)]=upper[which(VECTPAR>upper)]
        #VECTPAR[which(VECTPAR<lower)]=lower[which(VECTPAR<lower)]
        LL_last=fct(VECTPAR)
        jj=jj+1 }

        ii=ii+1
        nbr_fcalls=jj
        if ((ii/10)%%1==0){    # save every 10 parameter proposition (just for information how did we explore the parameter space)
            VECThist[ii/10,]=VECTPAR
        }



    }
    #browser()
    envir$LL<-max(Lvisible[Lvisible<(-0.5)])
    if(length(as.array(VECTPAR_elite[which(VECTPAR_elite[,ncol(VECTPAR_elite)]==  max(which(Lvisible==max(Lvisible[Lvisible<(-0.5)]))) ), 1:(ncol(VECTPAR_elite)-1) ] )  )==0){envir$param<-as.array(VECTPAR)}else{envir$param<-as.array(VECTPAR_elite[which(VECTPAR_elite[,ncol(VECTPAR_elite)]==  max(which(Lvisible==max(Lvisible[Lvisible<(-0.5)]))) ), 1:(ncol(VECTPAR_elite)-1) ] ) }
    envir$iter<-ii-3
    envir$fct_calls<- nbr_fcalls

}
