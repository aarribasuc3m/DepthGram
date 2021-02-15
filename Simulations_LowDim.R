########## SIMULATIONS FOR SMALL p THROUGH PARALLEL COMPUTING ##############
########## COMPARISON WITH FOM AND MSPLOT
# USE OF mrfDepth package to implement FOM
# USE OF functions and code in msplot_code folder to implement MSplot 
# (msplot_code folder as downloaded from  https://www.tandfonline.com/doi/suppl/10.1080/10618600.2018.1473781?scroll=top, supplement.rar)

# parallel computation #
library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

# Loading required functions and libraries on clusters
clusterCall(cl, function() source("Sim_MultFunData.R"))
clusterCall(cl, function() source("depthGram.R"))
clusterCall(cl, function() source("OutDetectionRule_DG.R"))
clusterCall(cl, function() source("other_functions.R"))
clusterCall(cl, function() source("msplot_code/loading.R"))  #auxiliary functions for msplot: PACKAGES INSTALLATION MIGHT BE REQUIRED!
clusterCall(cl, function() source("msplot_code/msplot.R"))   #msplot function
clusterCall(cl, function() library(mrfDepth))                #library for fom
 

#Simulation settings
n=85            #nb. of non-outlying curves
n.out=15        #nb. of outlying curves
n.tot=n+n.out   #total nb. of curves
type=7          #all magnitude, shape and joint outliers
N=100           #nb. of time points
replicates=200  #simulation runs


#Simulation
for (p in c(10,50)){   #Dimension of the multivariate functional data set
  
  for (model in 1:4){       #Simulation Model
    
    for(crate in c(0,0.25,0.5,0.75,1)){   #Contamination rate: proportion of dimension in which the outlier curve is indeed outlier
      
      results_sim <- foreach(i=1:replicates) %dopar% { 
        
        ## Data Generation
        Sim<-Sim_mfdata(n,N,p,n.out,c=crate,model=model,type.out=type)
        
        Data<-Sim$values
        w.outm<-Sim$w.outm  #Ids and outlying dimensions of outliers
        w.outs<-Sim$w.outs  
        w.outj<-Sim$w.outj
        
        rm(Sim)
        
        ## DepthGram
        DG<-depthGram(Data, marg.out = T)
        
        DG$w.outm<-w.outm  #We register the "true" outliers together with the result of the DepthGram
        DG$w.outs<-w.outs
        DG$w.outj<-w.outj
        
        out.DG<-out_det_DG(DG,n.tot)  #Outlier Identification for the DepthGram
        
        
        ## FOM
        #Functional p-variate version
        X3d=array(dim=c(N,n.tot,p))
        for (k in 1:p){
          X3d[,,k]=t(Data[[k]])
        }
        DO_X <-  mrfDepth::fOutl(x = X3d, type = "fDO", distOptions = list(rmZeroes = TRUE,maxRatio = 3),diagnostic = TRUE)
        set.seed(seed=NULL)  #because there's a set.seed call inside fOutl 
        FOM<-fom_ggplot(DO_X, type="fO",cutoff = TRUE,plot=FALSE)
        out.FOM=FOM$out
        
        #Functional univariate version
        X3d=array(dim=c(N*p,n.tot,1))
        for (k in 1:p){
          X3d[((k-1)*N+1):(k*N),,1]=t(Data[[k]])
        }
        DO_X2 <-  mrfDepth::fOutl(x = X3d, type = "fDO", distOptions = list(rmZeroes = TRUE,maxRatio = 3),diagnostic = TRUE)
        set.seed(seed=NULL)  #because there's a set.seed call inside fOutl 
        FOM2<-fom_ggplot(DO_X2, type="fO",cutoff = TRUE,plot=FALSE)
        out.FOM.1d=FOM2$out
        
        ## MSplot
        #Functional p-variate version
        X3d=array(dim=c(n.tot,N,p))  #Data reshaping into a 3d array
        for (k in 1:p){
          X3d[,,k]=Data[[k]]
        }
        MSP<-msplot(data = X3d, depth.dir="RP",plot=FALSE,dirout=T)
        out.MSplot<-MSP$out.dir
       
        #Functional univariate version
        X3d=array(dim=c(n.tot,N*p,1))
        for (k in 1:p){
          X3d[,((k-1)*N+1):(k*N),1]=Data[[k]]
        }
        MSP2<-msplot(data = X3d[,,1], depth.dir="RP",plot=FALSE,dirout=T)
        out.MSplot.1d<-MSP2$out.dir
  
        
        return(list(out.DG=out.DG,out.mag.DG=DG$mag.out.det,out.mag.true=DG$w.outm,out.sh.true=DG$w.outs,
                    out.jt.true=DG$w.outj,
                    out.FOM=out.FOM,out.FOM.1d=out.FOM.1d,out.MSplot=out.MSplot,out.MSplot.1d=out.MSplot.1d))
        
        }
      
      save(results_sim, file=paste0("Results_LowDimSim_model_",model,"_crate_",crate,"_p_",p,".RData"))
      rm(results_sim)
      gc()
    }
  }
}

stopCluster(cl)
