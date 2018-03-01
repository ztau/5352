#'Maximization procedure "S"
#'
#'@param fct function to be maximized (logLikelihood)
#'@param VECTPAR vector of visiible parameters Phi and Theta
#'@return logL, parameters, iterations and number of function calls

la_nouvelle_approche_max_AVEC_CHANGEMENTS_PRECEDENTS_maxiterations <- function(fct, VECTPAR, target, maxdelta, mindelta, it, envir...){

    # S

    # maxdelta=0.3 - max percentage change (ex: 0,25)
    # mindelta=0.005 - min percentage change (ex: 0,01)
    # fct=LL_10
    # target=0

    LLAlpha_vis_last=fct(VECTPAR)
    LLAlpha_vis=LLAlpha_vis_last
    VECTPAR_elite=cbind(t(VECTPAR),0)
    ii=3
    lengthVECTPAR=length(VECTPAR)
    changes<-array(0.2, c(1,lengthVECTPAR))
    VECThist=array(0,c(5000,lengthVECTPAR))
    LLAlphahist_up<- array(0,c(1,lengthVECTPAR))
    LLAlphahist_down<- array(0,c(1,lengthVECTPAR))
    Lvisible=array(0,c(1,it))
    Lvisible[ii]=-1
    Lvisible[ii-1]=-0.5
    nbr_fcalls=0

    while(nbr_fcalls<it){
        LLAlpha_vis_last<-LLAlpha_vis


        VECTPAR_up=VECTPAR
        VECTPAR_down=VECTPAR

        for (m in 1:lengthVECTPAR){


            if(abs(VECTPAR[m])<0.05 & m>lengthTheta){
                VECTPAR[m]=-VECTPAR[m]                    # if parameter too small AND IS NOT a VARIANCE PARAMETER, try NEGATIVE value
            }




            VECTPAR_up[m]=VECTPAR[m]+max(VECTPAR[m]*changes[m],0.005)                                                       # positive change


            LLAlphahist_up[m]=fct(VECTPAR_up)     # vecteur de L pour chaque param modifi? vers le haut




            VECTPAR_down[m]=VECTPAR[m]-max(VECTPAR[m]*changes[m],0.005)                                                       # negative change


            LLAlphahist_down[m]=fct(VECTPAR_down)   # vecteur de L pour chaque param modifi? vers le bas

            if(LLAlphahist_up[m]>LLAlphahist_down[m]){VECTPAR[m]=VECTPAR_up[m]}
            else if(LLAlphahist_up[m]==LLAlphahist_down[m]){VECTPAR[m]=VECTPAR[m]+runif(1,min=-0.5*abs(VECTPAR[m]),max=0.5*abs(VECTPAR[m]))} ######################################################################## new 13.03.2015
            else{VECTPAR[m]=VECTPAR_down[m]}

        }

        negative_deltaL_place<-which((LLAlphahist_up-LLAlpha_vis_last)<0 & (LLAlphahist_down-LLAlpha_vis_last)<0)     # when the L doesn't grow up by changing 1 param
        changes[negative_deltaL_place]<-max(changes[negative_deltaL_place]/3,mindelta)                      # diminish its change step

        positive_deltaL_place_up<-which((LLAlphahist_up-LLAlpha_vis_last)>0)                                             #
        positive_deltaL_place_down<-which((LLAlphahist_down-LLAlpha_vis_last)>0)                                         #

        #         # we modify only the parameters that raise L (once upwards and downwards)                       #
        #         VECTPAR[positive_deltaL_place_up]=VECTPAR[positive_deltaL_place_up]#*(1+changes[positive_deltaL_place_up])
        #         VECTPAR[positive_deltaL_place_down]=VECTPAR[positive_deltaL_place_down]#*(1-changes[positive_deltaL_place_down])
        # we increase the  changes in this case
        changes[positive_deltaL_place_up]=min(changes[positive_deltaL_place_up]*2,maxdelta)                              #
        changes[positive_deltaL_place_down]=min(changes[positive_deltaL_place_down]*2,maxdelta)                          #




        LLAlpha_vis=fct(VECTPAR)
        #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



        Lvisible[ii]=LLAlpha_vis
        print(LLAlpha_vis)

        if(Lvisible[ii]>max(Lvisible[Lvisible<(-0.5)])){VECTPAR_elite=rbind(VECTPAR_elite,cbind(t(VECTPAR),ii)) # if the last logL is the best so far, save the last VECTPAR !!
        #VECTPAR=VECTPAR+rnorm(length(VECTPAR))*0.2*VECTPAR
        #LLAlpha_vis=fct(VECTPAR)
        }

        ii=ii+1
        nbr_fcalls=ii*(lengthVECTPAR*2+1)
        if ((ii/10)%%1==0){    # si l'iteration est divisable par 10 sans d?cimales (chaque 10e iteration on note VECTPAR)
            VECThist[ii/10,]=VECTPAR
        }
        # print(changes)

    }
    #return=list(max(Lvisible[Lvisible<(-0.5)]), VECTPAR_elite[nrow(VECTPAR_elite),1:(ncol(VECTPAR_elite)-1)], ii-3, nbr_fcalls)  # car le 2e est tjs =-0.5
    envir$LL<-max(Lvisible[Lvisible<(-0.5)])
    envir$param<-VECTPAR_elite[nrow(VECTPAR_elite),1:(ncol(VECTPAR_elite)-1)]
    envir$iter<-ii-3
    envir$fct_calls<- nbr_fcalls

}
