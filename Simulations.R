########## SIMULATIONS FOR LARGE p THROUGH PARALLEL COMPUTING ##############

# parallel computation #
library(foreach)
library(doParallel)


#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

# Loading required functions on clusters
clusterCall(cl, function() source("Sim_MultFunData.R"))
clusterCall(cl, function() source("depthGram.R"))

#Simulation settings
n=85            #nb. of non-outlying curves
n.out=15        #nb. of outlying curves
type=7          #all magnitude, shape and joint outliers
N=100           #nb. of time points
replicates=200  #simulation runs


#Simulation
for (p in c(10000,50000)){   #Dimension of the multivariate functional data set
  
   for (model in 1:4){       #Simulation Model

     for(crate in c(0,0.25,0.5,0.75,1)){   #Contamination rate: proportion of dimension in which the outlier curve is indeed outlier

        resultsDG <- foreach(i=1:replicates) %dopar% { 
     
         #Data Generation
         Sim<-Sim_mfdata(n,N,p,n.out,c=crate,model=model,type.out=type)
  
         Data<-Sim$values
         w.outm<-Sim$w.outm  #Ids and outlying dimensions of outliers
         w.outs<-Sim$w.outs  
         w.outj<-Sim$w.outj
     
         rm(Sim)

         DG<-depthGram(Data, marg.out = T)
      
         DG$w.outm<-w.outm  #We register the "true" outliers together with the result of the DepthGram
         DG$w.outs<-w.outs
         DG$w.outj<-w.outj
  
         return(DG)
  
         rm(Data,DG)
        }
        
        save(resultsDG, file=paste0("Results_Sim_model_",model,"_crate_",crate,"_p_",p,".RData"))
        rm(resultsDG)
        gc()
     }
   }
}
        
stopCluster(cl)
