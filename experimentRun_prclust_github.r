library("fossil")
library("prclust")

getParamtersFrmPrclust <- function(data){
    data=t(data)
    lambda1=1 #fixed
    lambda2=NULL
    tau=NULL
    stab_max=-1
    for(tau_tmp in seq(1,10,2)){
        tau_tmp=tau_tmp/10
        
        for(lambda2_tmp in c(0.01, 0.05, 0.1, 0.2,0.5,0.7,1,1.5,2)){
            stab_tmp=NULL
            stab_tmp=stability(data,lambda1,lambda2_tmp,tau_tmp,n.times=50)
            
            if(is.na(stab_tmp)){
                next
            }
            
            cat("stab=",stab_tmp," ")
            if (stab_tmp>stab_max){
                lambda2=lambda2_tmp
                tau=tau_tmp
            }
        }
    }
    res=list()
    res=list("lambda1"=1,"lambda2"=lambda2,"tau"=tau)
    return(res)
}

resFrmPrclust <- function(data,lambda1,lambda2,tau,lables){
    data=t(data)
    res1=NULL
    res1=PRclust(data,lambda1,lambda2,tau)
    
    res2=clusterStat(res1$group,lables)
    
    randIdx=NULL;AdjustedRandIdx=NULL;
    randIdx=round(res2$Rand,digits = 2)
    AdjustedRandIdx=round(res2$AdjustedRand,digits = 2)
    
    res=NULL
    res=list("randIdx"=randIdx,"AdjustedRandIdx"=AdjustedRandIdx)
    return(res)
}

### for prclust end
main <- function(simulationData,savePath){
    
    finalRes4backup_prclust=list()
    #finalRes=list()
    
    for (i1 in 1:length(simulationData)){
        testType=NULL
        testType=names(simulationData[i1]) # simu_k,simu_n,simu_sigma,....
        
        k=2
        
        allTestTypeRes=list()
        for(i2 in 1:length(simulationData[[i1]])){ # test each parameter in one test type
            parameter_value=NULL
            parameter_value=names(simulationData[[i1]][i2]) # sig=0.1
            
            if(testType=="simu_k"){
                k=NULL
                k=strtoi(strsplit(parameter_value,'=')[[1]][2])
                print(paste("The new k:",k,sep = " "))
            }
            
            oneTestTypeRes=list()
            
            
            print(paste("start test type:",parameter_value,'method:prclust',sep = ' '))
            tmpRes=c()
            
            data=NULL
            data=simulationData[[i1]][[i2]][[1]]$mat
            parameters_prclust=NULL
            parameters_prclust=getParamtersFrmPrclust(data)
            for(i3 in 1:length(simulationData[[i1]][[i2]])){
                print(paste("The data Id:",i3))
                data=NULL;xy=NULL;labels_true=NULL;beta_true=NULL;
                data=simulationData[[i1]][[i2]][[i3]]$mat
                #data=data.frame(data)
                xy=simulationData[[i1]][[i2]][[i3]]$xy
                labels_true=simulationData[[i1]][[i2]][[i3]]$cl_new
                beta_true=simulationData[[i1]][[i2]][[i3]]$beta
                
                # get the results
                res_prclust=NULL
                res_prclust=resFrmPrclust(data,parameters_prclust$lambda1,parameters_prclust$lambda2,parameters_prclust$tau,labels_true)
                
                # get the outputs rand inex and ARandIndex
                randIdx=NULL;AdjustedRandIdx=NULL;
                randIdx=res_prclust$randIdx
                AdjustedRandIdx=res_prclust$AdjustedRandIdx
                
                outlierAccuracy="NaN";beta_euDist="NaN";
                tmpRes=append(tmpRes,c(randIdx,AdjustedRandIdx,outlierAccuracy,beta_euDist))
            }
            
            # convert list to matrix
            tmpRes=matrix(tmpRes,ncol=4,byrow = TRUE)
            colnames(tmpRes)=c("RI","ARI","outlierAcc","beta_euDist")
            oneTestTypeRes["res_prclust"]=list(tmpRes)
            print(paste("End test type:",parameter_value,'| method:prclust',sep = ' '))
            cat("\n")# new line
            
            # save the results of 4 methos to allTestTypeRes
            allTestTypeRes[parameter_value]=list(oneTestTypeRes)
        }
        # save the results of all types to finalRes4backup
        finalRes4backup_prclust[testType]=list(allTestTypeRes)
        save(finalRes4backup_prclust,file = paste(savePath,"finalRes4backup_prclust.rdata",sep = ''))
    }
    return(TRUE)
}
savePath <- "prclust/"
load("simulated_data_list.RData")
main(simu_listAndList,savePath)
