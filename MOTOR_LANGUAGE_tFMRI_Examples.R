library(roahd)
source("depthGramPlot.R") 

#The size of the datasets for the LANGUAGE and Motor tFMRI experiments doesn't allow to build a single list for each experiment
#containing the whole multivariate functional data set. Therefore, we can not use the function depthGram directly in this setting.

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

### 2. COMPUTATIONS OF MBD AND MEI OVER VOXELS (DIMENSIONS)

## Reading individual files, extraction of signal in the intersection mask, storing in different files by blocks of 10000 voxels

# p.div<-10000    
# p=p.min
# p.times=ceiling(p/p.div)

# mbd.d <-array(0,dim=c(n,p))
# mei.d <-array(0,dim=c(n,p))

## For each block of 10000 voxels, we read all individuals data and store into a list (each list element corresponds to a voxel)

# for (j in 1:p.times){ #j stands for each block of p.div voxels
#   
#   b=min(p,(p.div*j))
#   a=p.div*(j-1)+1
#   r=b-a+1
#   aN=(a-1)*N+1
#   bN=b*N
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
#     ### mbd and mei calculations over voxels (we don't use roadh::MBD and roah::MEI to share common computations and speed the process)
#     i=a+k-1
#     x<-t(fmri.vx.list[[k]])  #transpose because fmri.vx.list[[k]] is N*n and not n*N 
#     rk = apply(x, 2, function(v) (rank(v)))
#     up=n-rk
#     mbd.d[,i] = (rowSums(up * (rk - 1))/N + n - 1)/(n * (n - 1)/2)
#     mei.d[,i] = rowSums(up+1)/(n*N)
#   }
#   
#   save(fmri.vx.list,file=paste0("datos_fmriLANGUAGE_voxels_",a,"_to_",b ,".RData"))
#   rm(fmri.vx.list,mat.ind)
# }
# 
# save(mbd.d,file="LANGUAGE_maskMAX_mbd_vxl_by_vxl.RData")
# save(mei.d,file="LANGUAGE_maskMAX_epi_vxl_by_vxl.RData")


### 3. COMPUTATIONS OF MBD AND MEI OVER TIME (ORIGINAL DATA SET AND "CORRELATION CORRECTED" DATA SET)

# corr = rep(1,p) #vector to store correlation of mei.d across voxels
# p_coord=c() #vector to store the coordinates of the dimensions that need to be "unfolded" for the time/correlation DepthGram
# for (j in 2:p){
#   corr[j]=corr[j-1]*sign(cor(mei.d[,j-1],mei.d[,j]))
#   if(corr[j]==-1){
#     p_coord=c(p_coord,j)
#   }
# }

# mbd.t <-array(0,dim=c(n,N))
# mei.t <-array(0,dim=c(n,N))
# mbd.t2 <-array(0,dim=c(n,N))
# mei.t2 <-array(0,dim=c(n,N))

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

# Load the data
load(file="tFMRI_experiment_files/LANGUAGE_mbd_by_voxel.RData") #contains mbd.d, the mbd matrix on voxels (dimensions)
load(file="tFMRI_experiment_files/LANGUAGE_mei_by_voxel.RData") #contains mei.d, the mei matrix on voxels (dimensions)
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


###########################
#### MOTOR EXPERIMENT #####
###########################

### 1 to 3 AS ABOVE (we provide the resulting data sets requiered for points 5 and 6)


### 4. FINAL CALCULATIONS FOR THE MOTOR EXPERIMENT (the files required from this step are provided)

# Load the data
load(file="tFMRI_experiment_files/MOTOR_mbd_by_voxel.RData") #contains mbd.d, the mbd matrix on voxels (dimensions)
load(file="tFMRI_experiment_files/MOTOR_mei_by_voxel.RData") #contains mei.d, the mei matrix on voxels (dimensions)
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

