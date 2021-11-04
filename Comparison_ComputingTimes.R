########## COMPARISON OF COMPUTING TIMES OF DepthGram, MS-plot and FOM ##############
# USE OF mrfDepth package to implement FOM
# USE OF functions and code in msplot_code folder to implement MSplot 
# (msplot_code folder as downloaded from  https://www.tandfonline.com/doi/suppl/10.1080/10618600.2018.1473781?scroll=top, supplement.rar)

# parallel computation #
library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) # nb. of cores depending on the system specification
registerDoParallel(cl)

# Loading required functions and libraries on clusters
clusterCall(cl, function() {source("depthGram.R")})
clusterCall(cl, function() {source("msplot_code/loading.R")})  #functions for msplot: PACKAGES INSTALLATION MIGHT BE REQUIRED!
clusterCall(cl, function() {library(mrfDepth)})                #library for fom

###########################################
#LOW-DIM SETTING: DepthGram, FOM AND MSplot 
n = 100                    #nb. of curves
p = c(seq(10,90,10),99)    #vector of nb. of dimensions
N = seq(100,500,100)       #vector of nb. of time points

lp = length(p)
lN = length(N)

Times = array (dim=c(lp,5,lN))

#for (k in 1:lN){   
  for (k in 2:lN){   
      results_CompTimes <- foreach(j=1:lp) %dopar% { 
        
        # Data generation in form of a 3d array
        X = array(rnorm(n*p[j]*N[k]),dim=c(n,N[k],p[j]))
        
        # MSplot computing time
        MSplot = system.time(DirOut(data = X, depth.dir="RP",D_value=FALSE))[3]
        
        # Data transformation into 1-dim functional data set
        Xlong = matrix(nrow=n,ncol=N[k]*p[j])
        for (h in 1:p[j]){
          Xlong[,((h-1)*N[k]+1):(h*N[k])] = X[,,h]
        }
        
        # MSplot 1-dim version computing time
        MSplot1d = system.time(DirOut(data = Xlong, depth.dir="RP",D_value=FALSE))[3]  
        
        # Data transformation into FOM format for 1-dim functional data set
        tXlong = array(dim=c(N[k]*p[j],n,1))
        tXlong[,,1] = t(Xlong)
        
        # FOM 1-dim version computing time
        FOM1d <- system.time(fOutl(x = tXlong, type = "fDO", distOptions = list(rmZeroes = TRUE,maxRatio = 3),diagnostic = FALSE))[3]
        rm(Xlong,tXlong)
        
        # Data transformation into list format and 3d array for FOM p-variate version
        Data = list(length(p[j]))
        tX= array(dim = c(N[k],n,p[j]))
        for (h in 1:p[j]){
             Data[[h]] = X[,,h]
             tX[,,h] = t(Data[[h]])
        }
        rm(X)
        
        # DepthGram computing time
        DG = system.time(depthGram(Data, marg.out = F))[3]
        rm(Data)
        
        # FOM version computing time
        FOM = system.time(fOutl(x = tX, type = "fDO", distOptions = list(rmZeroes = TRUE,maxRatio = 3),diagnostic = FALSE))[3]
        rm(tX)
        
        results = c(DG,FOM,FOM1d,MSplot,MSplot1d)
        return(results)
        
      }
      
      Times[,,k] = matrix(unlist(results_CompTimes),ncol=5,byrow=T)
}      

Times = reshape2::melt(Times)
names(Times) = c("p", "Method", "N", "time")
Method = c("DepthGram", "FOM p-dim", "FOM 1-dim", "MSplot p-dim", "MS-plot 1-dim")
Times$Method =  Method[Times$Method]
Times$N = N[Times$N]
Times$p = p[Times$p]

Times$N = factor(Times$N)
levels(Times$N) = paste("N =", N)

### Plot
ggplot(Times, aes(x=p, y=time)) + geom_line(aes(color=Method))+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  annotation_logticks(sides="l")  +  # "l" stands for left
  facet_grid(~N)+
  theme_light()+ 
  theme(legend.position="bottom",strip.text = element_text(size=14, face = "bold"), axis.title= element_text(size=14),
        axis.text= element_text(size=12),title=element_text(size=16),legend.text=element_text(size=12),legend.title = element_text(size=14))+ 
  xlab("p") + ylab("Time (sec.) - Log scale") + ggtitle("Computation times") 


############################################
#HIGH-DIM SETTING: ONLY DepthGram AND MSplot 
n = 100                                                      #nb. of curves
p = c(100,seq(5000,100000,5000),seq(110000,200000,10000))    #vector of nb. of dimensions
N = 100                                                      #nb. of time points

lp = length(p)

results_CompTimes <- foreach(j=1:lp) %dopar% { 
    
    # Data generation in form of a 3d array
    X = array(rnorm(n*p[j]*N),dim=c(n,N,p[j]))
    
    # MSplot computing time
    MSplot = system.time(DirOut(data = X, depth.dir="RP",D_value=FALSE))[3]
    
    # Data transformation into 1-dim functional data set
    Xlong = matrix(nrow=n,ncol=N*p[j])
    for (h in 1:p[j]){
      Xlong[,((h-1)*N+1):(h*N)] = X[,,h]
    }
    
    # MSplot 1-dim version computing time
    MSplot1d = system.time(DirOut(data = Xlong, depth.dir="RP",D_value=FALSE))[3]  
    rm(Xlong)
    
    # Data transformation into list format 
    Data = list(length(p[j]))
    for (h in 1:p[j]){
      Data[[h]] = X[,,h]
    }
    rm(X)
    
    # DepthGram computing time
    DG = system.time(depthGram(Data, marg.out = F))[3]
    rm(Data)

    results = c(DG,MSplot,MSplot1d)
    return(results)
}
  
Times = matrix(unlist(results_CompTimes),ncol= 3,byrow=T)
Times = reshape2::melt(Times)
names(Times) = c("p", "Method", "time")
Times$p = p[Times$p]
Method = c("DepthGram", "MSplot p-dim", "MS-plot 1-dim")
Times$Method = Method[Times$Method]
### Plot  
ggplot(Times, aes(x=p, y=time)) + geom_line(aes(color=Method)) +
      theme_light()+ 
      theme(legend.position="bottom",strip.text = element_text(size=14, face = "bold"), axis.title= element_text(size=14),
            axis.text= element_text(size=12),title=element_text(size=16),legend.text=element_text(size=12),legend.title = element_text(size=14))+ 
      xlab("p") + ylab("Time (sec.)") + ggtitle("Computation times") 

    


stopCluster(cl)


