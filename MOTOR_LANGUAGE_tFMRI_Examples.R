library(roahd)
source("depthGramPlot.R") 

#The size of the datasets for the LANGUAGE and Motor tFMRI experiments doesn't allow to build a single list for each experiment
#containing the whole multivariate functional data set. Therefore, we can not use the function depthGram directly in this setting.
#In the following, we describe all the steps to conduct the Depthgram analysis on these data sets from individuals nifti files

###########################
### LANGUAGE EXPERIMENT ###
###########################

### 1. READING THE DATA FOR THE LANGUAGE EXPERIMENT (the files are NOT provided for size constraints)

## Individual data are contained in the files ID_tfMRI_LANGUAGE_RL.nii.gz (brain activity, ~250 MB) and ID_brainmask_fs.2.nii.gz (brain mask, 23 KB)
## were ID stands for the individual ID. There are 100 pairs of these files.

# library("oro.nifti") #Allows to read nii files
# names.lang<-list.files(pattern="*\\_tfMRI_LANGUAGE_RL.nii.gz")  # names of individual files for brain activity
# names.mask<-list.files(pattern="*\\_brainmask_fs.2.nii.gz")     # names of individual files for brain masks
# n<-length(names.lang) #nb. of individuals


## 1st individual exploration
# fmri<-readNIfTI(names.lang[1], verbose = FALSE, warn = -1, reorient = FALSE, call = NULL) # oro.nifti
# N=dim(fmri)[4]   #316    time points
# n1=dim(fmri)[1]  #91     dim1 of brain
# n2=dim(fmri)[2]  #109    dim2 of brain
# n3=dim(fmri)[3]  #91     dim3 of brain

## Extraction of intersection mask
# bin.min<-array(1,dim=c(n1,n2,n3))  
# 
# for (i in 1:n){  #Reading individual masks
#   aux<-readNIfTI(names.mask[i], verbose = FALSE, warn = -1, reorient = FALSE, call = NULL) # oro.nifti
#   aux<-aux[ dim1=1:n1, dim2=1:n2, dim3=1:n3]
#   bin.min<-bin.min*aux                            #1 if voxel takes value 1 for all individuals, 0 otherwise
# }
# 
# bin.vmin<-as.vector(bin.min)
# w.min<-which(bin.min==1)               
# vox.min<-arrayInd(w.min, c(n1,n2,n3))  # 3-dim indexes of voxels in the intersection mask 

# p.min=dim(vox.min)[1] #192631
# 
# vox.N<-vcbind(matrix(rep(vox.min,each=N),nrow=p.min*N),rep(1:N,p.min))

## We first read all nifti files and store selected voxels (intersection mask) in RData objets: it speeds lecture of files
# for (i in 1:n){  #Reading accross individulas
#   aux<-readNIfTI(names.lang[i], verbose = FALSE, warn = -1, reorient = FALSE, call = NULL) #oro.nifti
#   aux<-aux[dim1=1:n1,dim2=1:n2,dim3=1:n3,dim4=1:N]
#   aux=aux[vox.N]
#   save(aux,file=paste0("Ind",i,"_LANGUAGE.RData"))
#   rm(aux)
# }
# rm(vox.N)

### 2. COMPUTATIONS OF MBD AND MEI OVER VOXELS (DIMENSIONS) AND MARGINAL MAGNITUDE AND SHAPE OUTLIER DETECTION

## Reading individual files, extraction of signal in the intersection mask, storing in different files by blocks of 10000 voxels

# p.div<-10000    
# p=p.min
# p.times=ceiling(p/p.div)

# mbd.d <-array(0,dim=c(n,p))   #array for mbd over voxels
# mei.d <-array(0,dim=c(n,p))   #array for mei over voxels
# mag.out.det<-list(length=p)   #list for marginal magnitude outliers detection
# shp.out.det<-list(length=p)   #list for marginal shape outliers detection
# n2=ceiling(n*0.5)
# a2=a0=-2/(n*(n-1))            #values for marginal shape outlier detection
# a1=2*(n+1)/(n-1)
# FO = FB = 1.5                 #values for marginal shape and magnitude outlier detection 

## For each block of 10000 voxels, we read all individuals data and store into a list (each list element corresponds to a voxel)

# for (j in 1:p.times){ #j stands for each block of p.div voxels
#   
#   b=min(p,(p.div*j))
#   a=p.div*(j-1)+1
#   r=b-a+1
#   
#   mat.ind<-c() 
#   
#   for (i in 1:n){  #Reading across individulas
#     load(paste0("Ind",i,"_LANGUAGE.RData"))
#     mat.ind<-cbind(mat.ind,matrix(aux[aN:bN],nrow=p.div,ncol=N,byrow=TRUE))
#     rm(aux)
#   }
#   
#   fmri.vx.list<-vector("list",b-a+1) #this is a list of length equal to the nb. of voxels (usually p.div, but not always). 
#   for (k in 1:r){                    #Each element is an N*n matrix
#     fmri.vx.list[[k]]<-matrix(mat.ind[k,],nrow=N)
#   
#     ### mbd and mei calculations over voxels (we don't use roadh::MBD and roah::MEI to share common computations and speed up the process)
#     i=a+k-1
#     x<-t(fmri.vx.list[[k]])  #transpose because fmri.vx.list[[k]] is N*n and not n*N 
#     rk = apply(x, 2, function(v) (rank(v)))
#     up=n-rk
#     mbd.d[,i] = (rowSums(up * (rk - 1))/N + n - 1)/(n * (n - 1)/2)
#     mei.d[,i] = rowSums(up+1)/(n*N)
#
#    ### Marginal magnitude and shape outlier detection (we don't use roahd::fbplot and roahd::outliergram to share common computations and speed up the process)
# 
#    index<-order(mbd.d[,i],decreasing=T)  
#    center=x[index[1:n2],]
#    inf=apply(center,2,min)
#    sup=apply(center,2,max)
#    dist=FB*(sup-inf)
#    upper=sup+dist
#    lower=inf-dist
#    mag.out.det[[i]]<-which(colSums((t(x) <= lower) + (t(x) >= upper))>0)
#    dist=(a0+a1*mei.d[,i]+a2*n^2*mei.d[,i]^2)-mbd.d[,i]
#    q=quantile(dist,probs = c(0.25,0.75))
#    lim=FO*(q[2]-q[1])+q[2]
#    shp.out.det[[i]]<- which(dist>lim) 
#    rm(index,center,inf,sup,dist,upper,lower,q,lim)
#  }
#   
#   save(fmri.vx.list,file=paste0("datos_fmriLANGUAGE_voxels_",a,"_to_",b ,".RData"))
#   rm(fmri.vx.list,mat.ind)
# }
# 
# p1 = round(p/2)   # Due to file size constraints in GitHub repositories, mbd.d and mei.d matrices are split and saved in 2 different files each
# mbd.d1=mbd.d[,1:p1]
# mbd.d2=mbd.d[,(p1+1):p]
# mei.d1=mei.d[,1:p1]
# mei.d2=mei.d[,(p1+1):p]
# save(mbd.d1,file="LANGUAGE_mbd_by_voxel1.RData")
# save(mbd.d2,file="LANGUAGE_mbd_by_voxel2.RData")
# save(mei.d1,file="LANGUAGE_mei_by_voxel1.RData")
# save(mei.d2,file="LANGUAGE_mei_by_voxel2.RData")
# save(mag.out.det, shp.out.det, file="LANGUAGE_marginal_outliers_by_voxel.RData")


### 3. COMPUTATIONS OF MBD AND MEI OVER TIME (ORIGINAL DATA SET AND "CORRELATION CORRECTED" DATA SET)

# corr = rep(1,p) #vector to store correlation of mei.d across voxels
# p_coord=c() #vector to store the coordinates of the dimensions that need to be "unfolded" for the time/correlation DepthGram
# for (j in 2:p){
#   corr[j]=corr[j-1]*sign(cor(mei.d[,j-1],mei.d[,j]))
#   if(corr[j]==-1){
#     p_coord=c(p_coord,j)
#   }
# }

# mbd.t <-array(0,dim=c(n,N))   #array for mbd over time
# mei.t <-array(0,dim=c(n,N))   #array for mei over time
# mbd.t2 <-array(0,dim=c(n,N))  #array for mbd over time/corr
# mei.t2 <-array(0,dim=c(n,N))  #array for mei over time/corr

# for (j in 1:N){ #j stands for each possible time point
  
#   mat.ind<-c() 
  
  # for (i in 1:n){  #For each time point j we read across individuals and get all signals for t=j
    # load(paste0("Ind",i,"_LANGUAGE_maskMIN.RData")) #contains aux
    # aux2=matrix(aux,nrow=p)
    # rm(aux)
    # mat.ind<-rbind(mat.ind,aux2[,j])
    # rm(aux2)
  #}
  
  ### mbd and mei calculations over time (we don't use roadh::MBD and roah::MEI to share common computations and speed the process)
  # rk = apply(mat.ind, 2, function(v) (rank(v)))
  # up=n-rk
  # mbd.t[,j] = (rowSums(up * (rk - 1))/p + n - 1)/(n * (n - 1)/2)
  # mei.t[,j]<-rowSums(up+1)/(n*p)

  ### mbd and mei calculations on the "correlation corrected" sample

  # mat.ind2=mat.ind
  # mat.ind2[,p_coord]=-1*mat.ind[,p_coord]
  # rk = apply(mat.ind2, 2, function(v) (rank(v)))  
  # up=n-rk
  # mbd.t2[,j] = (rowSums(up * (rk - 1))/p + n - 1)/(n * (n - 1)/2)
  # epi.t2[,j]<-rowSums(up+1)/(n*p)
  
  # save(mat.ind,file=paste0("datos_fmriLANGUAGE_time_point_",j,".RData"))
  # save(mat.ind2,file=paste0("datos_fmriLANGUAGE_time_point2_",j,".RData"))
  # rm(mat.ind,rk,up)
# }

# save(mbd.t,file="LANGUAGE_mbd_by_time.RData")
# save(mei.t,file="LANGUAGE_mei_by_time.RData")
# save(mbd.t2,file="LANGUAGE_mbd_by_time_corr.RData")
# save(mei.t2,file="LANGUAGE_mei_by_time_corr.RData")


### 4. FINAL CALCULATIONS FOR THE LANGUAGE EXPERIMENT (the files required from this step are provided)

# Load the data (due to file size constraints in GitHub repositories, mbd.d and mei.d matrices were split and saved in 2 different files each)
load(file="tFMRI_experiment_files/LANGUAGE_mbd_by_voxel1.RData") #contains mbd.d1, 1st part of mbd.d, the mbd matrix on voxels (dimensions)
load(file="tFMRI_experiment_files/LANGUAGE_mbd_by_voxel2.RData") #contains mbd.d2, 2nd part of mbd.d, the mbd matrix on voxels (dimensions)
mbd.d = cbind(mbd.d1,mbd.d2)
load(file="tFMRI_experiment_files/LANGUAGE_mei_by_voxel1.RData") #contains mei.d1, 1st part of mei matrix on voxels (dimensions)
load(file="tFMRI_experiment_files/LANGUAGE_mei_by_voxel2.RData") #contains mei.d2, 2nd part of mei matrix on voxels (dimensions)
mei.d = cbind(mei.d1, mei.d2)
load(file="tFMRI_experiment_files/LANGUAGE_mbd_by_time.RData") #contains mbd.t, the mbd matrix for time
load(file="tFMRI_experiment_files/LANGUAGE_mei_by_time.RData") #contains mei.t, the mei matrix for time
load(file="tFMRI_experiment_files/LANGUAGE_mbd_by_time_corr.RData") #contains mbd.t2, the mbd matrix for time/corr
load(file="tFMRI_experiment_files/LANGUAGE_mei_by_time_corr.RData") #contains mei.t2, the mei matrix for time/corr

# Build the object output of the DepthGram function by calculating mbd of mei and
# mei of mbd for voxels, time and time/correlation (MBD and MEI functions in roahd package)
DG <- data.frame(mbd.mei.d=MBD(mei.d), mei.mbd.d=MEI(mbd.d), mbd.mei.t=MBD(mei.t),
                 mei.mbd.t=MEI(mbd.t), mbd.mei.t2=MBD(mei.t2), mei.mbd.t2=MEI(mbd.t2))


### 5. DEPTHGRAM PLOT FOR THE LANGUAGE EXPERIMENT

n=nrow(mbd.d)
ids=rep("",n)
ids[c(81,84,86)] <- as.character(c(81,84,86))

DGplot <- depthGramPlot(DG, print=T, text.labels=ids, ax.lims=F)


### 6. MARGINAL OUTLIERS VOXEL-WISE FOR THE LANGUAGE EXPERIMENT

load(file="tFMRI_experiment_files/LANGUAGE_marginal_outliers_by_voxel.RData") #contains lists mag.out.det and shp.out.det with the ids of potential magnitude and shape outliers on each voxel

sort(table(unlist(shp.out.det)),decreasing = T)  #To see the number of voxels in which each individual has been detected as shape outlier
sort(table(unlist(mag.out.det)),decreasing = T)  #To see the number of voxels in which each individual has been detected as magnitude outlier

###########################
#### MOTOR EXPERIMENT #####
###########################

### 1 to 3 AS ABOVE (we provide the resulting data sets required for points 4 and 5)


### 4. FINAL CALCULATIONS FOR THE MOTOR EXPERIMENT (the files required from this step are provided)

# Load the data (due to file size constraints in GitHub repositories, mbd.d and mei.d matrices were split and saved in 2 different files each)
load(file="tFMRI_experiment_files/MOTOR_mbd_by_voxel1.RData") #contains mbd.d1, 1st part of mbd.d, the mbd matrix on voxels (dimensions)
load(file="tFMRI_experiment_files/MOTOR_mbd_by_voxel2.RData") #contains mbd.d2, 2nd part of mbd.d, the mbd matrix on voxels (dimensions)
mbd.d = cbind(mbd.d1,mbd.d2)
load(file="tFMRI_experiment_files/MOTOR_mei_by_voxel1.RData") #contains mei.d1, 1st part of mei matrix on voxels (dimensions)
load(file="tFMRI_experiment_files/MOTOR_mei_by_voxel2.RData") #contains mei.d2, 2nd part of mei matrix on voxels (dimensions)
mei.d = cbind(mei.d1, mei.d2)
load(file="tFMRI_experiment_files/MOTOR_mbd_by_time.RData") #contains mbd.t, the mbd matrix for time
load(file="tFMRI_experiment_files/MOTOR_mei_by_time.RData") #contains mei.t, the mei matrix for time
load(file="tFMRI_experiment_files/MOTOR_mbd_by_time_corr.RData") #contains mbd.t2, the mbd matrix for time/corr
load(file="tFMRI_experiment_files/MOTOR_mei_by_time_corr.RData") #contains mei.t2, the mei matrix for time/corr


# Build the object output of the DepthGram function by calculating mbd of mei and
# mei of mbd for voxels, time and time/correlation (MBD and MEI functions in roahd package)
DG <- data.frame(mbd.mei.d=MBD(mei.d), mei.mbd.d=MEI(mbd.d), mbd.mei.t=MBD(mei.t),
                 mei.mbd.t=MEI(mbd.t), mbd.mei.t2=MBD(mei.t2), mei.mbd.t2=MEI(mbd.t2))


### 5. DEPTHGRAM PLOT FOR THE MOTOR EXPERIMENT

n=nrow(mbd.d)
ids=rep("",n)
ids[c(20,28,30,39,66,81,89)] <- as.character(c(20,28,30,39,66,81,89))

DGplot <- depthGramPlot(DG, print=T, text.labels=ids, ax.lims=F)


### 6. MARGINAL OUTLIERS VOXEL-WISE FOR THE MOTOR EXPERIMENT

load(file="tFMRI_experiment_files/MOTOR_marginal_outliers_by_voxel.RData") #contains lists mag.out.det and shp.out.det with the ids of potential magnitude and shape outliers on each voxel

sort(table(unlist(shp.out.det)),decreasing = T)  #To see the number of voxels in which each individual has been detected as shape outlier
sort(table(unlist(mag.out.det)),decreasing = T)  #To see the number of voxels in which each individual has been detected as magnitude outlier
