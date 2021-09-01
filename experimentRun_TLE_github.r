library("fossil")
library(RobMixReg)
library(robustbase)

main <- function(simulationData,savePath){
    
    finalRes4backup_TLE=list()
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
            
            print(paste("start test type:",parameter_value,'method:TLE',sep = ' '))
            tmpRes=c()
            for(i3 in 1:length(simulationData[[i1]][[i2]])){
                print(paste("The data Id:",i3))
                data=NULL;xy=NULL;labels_true=NULL;beta_true=NULL;
                data=simulationData[[i1]][[i2]][[i3]]$mat
                data=data.frame(data)
                xy=simulationData[[i1]][[i2]][[i3]]$xy
                labels_true=simulationData[[i1]][[i2]][[i3]]$cl_new
                beta_true=simulationData[[i1]][[i2]][[i3]]$beta
                
                errorInfo=NULL
                res_TLE=NULL
                errorInfo=try({
                    # get the results
                    res_TLE=TLE(formula= as.formula("y~x"),data=data, nc=k,tRatio=0.5,MaxIt=1000)
                },silent = TRUE)
                if ("try-error" %in% class(errorInfo)){
                    print("Error!")
                    tmpRes=append(tmpRes,c(-1,-1,-1,-1))
                    next
                }
                # get the outputs rand inex and ARandIndex
                labels_predicted=NULL;randIdx=NULL;AdjustedRandIdx=NULL;
                labels_predicted=res_TLE@ctleclusters
                
                randIdx=round(rand.index(labels_predicted,labels_true),digits = 2)
                AdjustedRandIdx=round(adj.rand.index(labels_predicted,labels_true),digits = 2)
                
                # the euclidean distance of the beta_predicted and the beta_true
                beta_predicted=c()
                for (i in 1:k){
                    coef_x=NULL
                    coef_x=round(res_TLE@compcoef["coef.x",i][[1]],digits = 2)
                    beta_predicted=append(beta_predicted,c(coef_x))
                }
                print("beta_true:")
                print(sort(beta_true))
                print("beta_predicted:")
                print(sort(beta_predicted))
                beta_euDist=round(dist(rbind(sort(beta_true),sort(beta_predicted)))[1],digits=2)
                
                # get the outlier Accuracy
                outlierAccuracy=NULL
                outlierAccuracy= round(length(intersect(res_TLE@indout , which(labels_true<0))) / length(which(labels_true<0)),digits = 2)
                tmpRes=append(tmpRes,c(randIdx,AdjustedRandIdx,outlierAccuracy,beta_euDist))
            }
            
            # convert list to matrix
            tmpRes=matrix(tmpRes,ncol=4,byrow = TRUE)
            colnames(tmpRes)=c("RI","ARI","outlierAcc","beta_euDist")
            oneTestTypeRes["res_TLE"]=list(tmpRes)
            print(paste("End test type:",parameter_value,'| method:TLE',sep = ' '))
            cat("\n")# new line
            # for res_robSpaReg end
            
            # save the results of 4 methos to allTestTypeRes
            allTestTypeRes[parameter_value]=list(oneTestTypeRes)
        }
        
        # save the results of all types to finalRes4backup
        finalRes4backup_TLE[testType]=list(allTestTypeRes)
        save(finalRes4backup_TLE,file = paste(savePath,"finalRes4backup_TLE.rdata",sep = ''))
    }
    return(TRUE)
}
savePath <- "TLE/"
load("simulated_data_list.RData")
main(simu_listAndList,savePath)
