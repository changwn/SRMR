library("fossil")
library("ClustGeo")
library(robustbase)

getParamtersFrmClustGeo <- function(data_euDist,xy_euDist,k){
    range.alpha <- seq(0,1,0.1)
    cr=NULL
    cr = choicealpha(data_euDist,xy_euDist,range.alpha,k,graph=FALSE)
    
    alpha=0.1
    min=100
    for (i in 1:length(cr$Qnorm[,1])){
        q0=NULL;q1=NULL;
        q0=cr$Qnorm[i,1]
        q1=cr$Qnorm[i,2]
        if (q0==1 | q1==1 | q0 < 0.01 | q1<0.01){
            next
        }
        m_tmp=NULL
        m_tmp=2-q0-q1
        if(m_tmp<min){
            alpha=(i-1)/10
            min=m_tmp
        }
    }
    return(alpha)
}

resFrmClustGeo <- function(data_euDist,xy_euDist,alpha,k,labels_true){
    tree=NULL
    tree = hclustgeo(data_euDist,xy_euDist,alpha)
    labels_predicted=NULL
    labels_predicted=cutree(tree,k)
    randIdx=NULL;AdjustedRandIdx=NULL;
    randIdx=round(rand.index(labels_predicted,labels_true), digits = 2)
    AdjustedRandIdx=round(adj.rand.index(labels_predicted,labels_true),digits = 2)
    res=NULL
    res=list("randIdx"=randIdx,"AdjustedRandIdx"=AdjustedRandIdx)
    return(res)
}

main <- function(simulationData,savePath){
    finalRes4backup_clustGeo=list()
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
            
            
            print(paste("start test type:",parameter_value,'method:clustGeo',sep = ' '))
            tmpRes=c()
            
            data=NULL;xy=NULL;data_euDist=NULL;xy_euDist=NULL;
            data=simulationData[[i1]][[i2]][[1]]$mat
            data_euDist = dist(data)
            xy=simulationData[[i1]][[i2]][[1]]$xy
            xy_euDist=dist(xy)
            parameters_clustGeo=NULL
            parameters_clustGeo=getParamtersFrmClustGeo(data_euDist,xy_euDist,k)
            for(i3 in 1:length(simulationData[[i1]][[i2]])){
                print(paste("The data Id:",i3))
                data=NULL;data_euDist=NULL;xy=NULL;xy_euDist=NULL;labels_true=NULL;beta_true=NULL;
                data=simulationData[[i1]][[i2]][[i3]]$mat
                data_euDist = dist(data)
                xy=simulationData[[i1]][[i2]][[i3]]$xy
                xy_euDist=dist(xy)
                labels_true=simulationData[[i1]][[i2]][[i3]]$cl_new
                beta_true=simulationData[[i1]][[i2]][[i3]]$beta
                
                # get the results
                res_clustGeo=NULL
                res_clustGeo=resFrmClustGeo(data_euDist,xy_euDist,parameters_clustGeo,k,labels_true)
                
                # get the outputs rand index and ARandIndex
                randIdx=NULL;AdjustedRandIdx=NULL;
                randIdx=res_clustGeo$randIdx
                AdjustedRandIdx=res_clustGeo$AdjustedRandIdx
                
                # the euclidean distance of the beta_predicted and the beta_true
                outlierAccuracy="NaN";beta_euDist="NaN";
                tmpRes=append(tmpRes,c(randIdx,AdjustedRandIdx,outlierAccuracy,beta_euDist))
            }
            
            # convert list to matrix
            tmpRes=matrix(tmpRes,ncol=4,byrow = TRUE)
            colnames(tmpRes)=c("RI","ARI","outlierAcc","beta_euDist")
            oneTestTypeRes["res_clustGeo"]=list(tmpRes)
            print(paste("End test type:",parameter_value,'| method:clustGeo',sep = ' '))
            cat("\n")# new line
            
            # save the results of 4 methos to allTestTypeRes
            allTestTypeRes[parameter_value]=list(oneTestTypeRes)
        }
        
        # save the results of all types to finalRes4backup
        finalRes4backup_clustGeo[testType]=list(allTestTypeRes)
        save(finalRes4backup_clustGeo,file = paste(savePath,"finalRes4backup_clustGeo.rdata",sep = ''))
    }
    return(TRUE)
}
savePath <- "clustGeo/"
load("simulated_data_list.RData")
main(simu_listAndList,savePath)
