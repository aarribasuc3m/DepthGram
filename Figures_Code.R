#### Code to reproduce the Figures in the paper

library(roahd)
library(ggplot2)
source("depthGram.R")
source("depthGramPlot.R")
source("Sim_MultFunData.R")
####### Figure 1 #######
########################

### Data generation
n=85
N=100

time_grid = seq( 0, 1, length.out = N )

C1 <- exp_cov_function( time_grid, alpha = 0.3, beta = 0.2 )
cholC1=chol(C1)

Data<-list()
Data[[1]] <- matrix(rep(cos(2*pi*(time_grid)-0.25),n),nrow=n,byrow=T)+matrix(rep(seq(0,0.05*n,,n),each=N),nrow=n,byrow = T)+
           0.4*matrix(rnorm(n*N), nrow = n, ncol = N) %*% cholC1
 
Data[[2]] <- exp(-(matrix(rep(time_grid,n),nrow=n,byrow=T)-0.5)^2/0.1)*matrix(rep(seq(0,2,,n),each=N),nrow=n,byrow = T)+
          0.08*matrix(rnorm(n*N), nrow = n, ncol = N) %*% cholC1

### Outliers
joint.out.1 <- Data[[1]][10,]+0.1*matrix(rnorm(N), nrow = 1) %*% cholC1
joint.out.2 <- exp(-(time_grid-0.5)^2/0.1)*1.85+0.08*matrix(rnorm(N), nrow = 1) %*% cholC1
sh.out.1 <- cos(2*pi*(time_grid)-0.25)+n*0.03*0.5+0.4*matrix(rnorm(N), nrow = 1) %*% cholC1
sh.out.2 <- Data[[1]][1,]*0.3+0.375+0.08*matrix(rnorm(N), nrow = 1) %*% cholC1
mg.out.1 <- Data[[1]][n,]+1.75+0.4*matrix(rnorm(N), nrow = 1) %*% cholC1
mg.out.2 <- exp(-(time_grid-0.5)^2/0.1)*1.95+0.08*matrix(rnorm(N), nrow = 1) %*% cholC1

id.central <- sort(depthGram(Data)$mei.mbd.d,index.return =T)$ix[1] #a central curve

Data[[1]] <- rbind(Data[[1]] , joint.out.1, sh.out.1, mg.out.1)
Data[[2]] <- 3*rbind(Data[[2]], joint.out.2, sh.out.2, mg.out.2)

### Plot
color <- c(gray.colors(n, start = 0, end = 0.9, gamma = 1.5, alpha = NULL),"red","green","blue")
color[id.central] <- "magenta"
lw=c(rep(1,n),3,3,3)
lw[id.central] <- 3
plot(mfData(time_grid,Data),col=color,lwd=lw,
     xlab = 't', ylab = list( expression("x"[1]*"(t)"), expression("x"[2]*"(t)") ),
     main = list( 'First Dimension', 'Second Dimension' ) )


####### Figure 2 #######
########################

### DepthGram
DG<-depthGram(Data)
### Plot
lbs=rep("",85)
lbs[id.central]=paste0(id.central)
DGplot<-depthGramPlot(DG,col=color,text.labels=c(lbs,"86","87","88"))


####### Figure 3 #######
########################

#### Top panel - Curves

### Data generation (same data set as before, but negative association between components)
Data[[2]] <- Data[[2]][n:1,]

### Outliers
joint.out.2 <-exp(-(time_grid-0.5)^2/0.1)*0.22+0.08*matrix(rnorm(N), nrow = 1) %*% cholC1
mg.out.2 <- exp(-(time_grid-0.5)^2/0.1)*0.1-0.05+0.08*matrix(rnorm(N), nrow = 1) %*% cholC1

Data[[2]] <- rbind(Data[[2]], 3*joint.out.2, 3*sh.out.2, 3*mg.out.2)

### Plot
plot(mfData(time_grid,Data),col=color,lwd=lw,
     xlab = 't', ylab = list( expression("x"[1]*"(t)"), expression("x"[2]*"(t)") ),
     main = list( 'First Dimension', 'Second Dimension' ) )


#### Mid panel - Parallel coordinates
ntot=n+3
Timepoints=c(20,50)
p=2 #Nb of dimensions

Data.t<-list()      #Data for time points 20 and 50
for (j in 1:2){
  Data.t[[j]]<-matrix(0,nrow=ntot,ncol=p)
  for (i in 1:p){
    Data.t[[j]][,i]<-Data[[i]][,Timepoints[j]]
  }
}  

### Plot
plot(mfData(1:2,Data.t),col=color,lwd=lw,
     xlab = 't', ylab = list( expression("x"[20]*"(t)"), expression("x"[50]*"(t)") ),
     main = list( expression("x(t"[20]*") over dimensions"), expression("x(t"[50]*") over dimensions")) )

#### Bottom Panel - DepthGram
### DepthGram
DG<-depthGram(Data)
### Plot
DGplot<-depthGramPlot(DG,col=color,text.labels=c(lbs,"86","87","88"))


####### Figure 4 #######
########################

#### Top panel - Curves
### Data generation (same data set as before + a third dimension, independence among components)
Data[[2]][1:n,] <- Data[[2]][sample(1:n),]
Data[[2]][ntot,] <- 3*(exp(-(time_grid-0.5)^2/0.1)*1.95+0.08*matrix(rnorm(N), nrow = 1) %*% cholC1 )  #blue curve

Data[[3]] <- matrix(rep(4*time_grid,n),nrow=n,byrow=T)+matrix(rep(runif(n,0,0.05*n),each=N),nrow=n,byrow = T)+
             0.4*matrix(rnorm(n*N), nrow = n, ncol = N) %*% cholC1

Data[[3]][id.central,] <- 4*time_grid+2.1+0.4*matrix(rnorm(N), nrow = 1) %*% cholC1 #central curve

### Outliers
joint.out.3 <- 4*time_grid+0.2+0.2*matrix(rnorm(N), nrow = 1) %*% cholC1
sh.out.3 <- 4*time_grid+Data[[1]][1,]*2+2+0.08*matrix(rnorm(N), nrow = 1) %*% cholC1
mg.out.3 <- 4*time_grid+4+0.2*matrix(rnorm(N), nrow = 1) %*% cholC1

Data[[3]] <- rbind(Data[[3]], joint.out.3, sh.out.3, mg.out.3)

### Plot
plot(mfData(time_grid,Data),col=color,lwd=lw,
     xlab = 't', ylab = list( expression("x"[1]*"(t)"), expression("x"[2]*"(t)"), expression("x"[3]*"(t)") ),
     main = list( 'First Dimension', 'Second Dimension', 'Third Dimension' ) )

#### Bottom Panel - DepthGram
### DepthGram
DG<-depthGram(Data)
### Plot
DGplot<-depthGramPlot(DG,col=color,text.labels=c(lbs,"86","87","88"))


####### Figure 5 #######
########################

n=85       #nb. of non-outlying curves
n.out=15   #nb. of outlying curves
type=7     #all magnitude, shape and joint outliers
N=100      #nb. of time points
p=10000    #dimension
c=1        #contamination proportion

time_grid = seq( 0, 1, length.out = N )

#colors for plots

color.mg <- hcl(h=210, l=c(40,50,60,70,80), c=c(20,40,60,80,100)) #blue scale for magnitude outliers
color.sh <- c("coral","indianred1","firebrick1","orangered","red2") #pink/red scale for shape outliers
color.jt <-c("#4A047D", "#650AA8", "#8D20DD", "#AF56F1", "#D7A0FF") #purple scale for magnitude outliers

color<-c(rep("grey80",n-2),rep(1,2),color.mg,color.sh,color.jt)
lw=c(rep(0.5,n-2),rep(0.5,2),rep(0.8,n.out-5),rep(1.2,5))
lt=rep(1,n+n.out)
lt[n-1] <-2


for (i in 1:4) {
  
    #### i-th Row - Model i
    Data<-Sim_mfdata(n,N,p,n.out,c,model=i,type.out=type)$values
    ### Plot
    
    plot(mfData(time_grid,list(Data[[1]],Data[[3000]],Data[[7001]],Data[[10000]])),col=color,lwd=lw,lty=lt,
      xlab = 't', ylab = list( expression("x"[1]*"(t)"), expression("x"[3000]*"(t)"), expression("x"[7001]*"(t)") , expression("x"[10000]*"(t)") ),
      main = list( '1st dimension', '3000th dimension', '7001st dimension', '10000th dimension'  ) )

}


####### Figure 6 ####### Figures 1-7 in the Supplementary material
########################

#colors for plots
gray="grey50"
blue="#00BFC4" 
purple="purple4"
red="#F8766D"
colors=c(blue,red,purple,gray)
n=100

#Repeat for each model and value of p
model=1
p=50

depths=c()

for(crate in c(0,0.25,0.5,0.75,1)){
    
      #LOAD HERE THE OUTPUT FILES OF Simulations.R #Sample files are provided for model=1, p=50, crate in [0,1]
      filename=paste0("Sample_Results_Files/Results_Sim_model_",model,"_crate_",crate,"_p_",p,".RData") 
      load(file=filename) #loads resultsDG
      
      A<-lapply(resultsDG,function(x){data.frame(mbd=c(x$mbd.mei.d,x$mbd.mei.t,x$mbd.mei.t2), mei=c(x$mei.mbd.d,x$mei.mbd.t,x$mei.mbd.t2))})
      depths<-rbind(depths,do.call(rbind,A))
      
      if (crate==1){  l=length(resultsDG) }
      
      rm(resultsDG,A) 
}
    
ids=c(rep("Rest of the sample",85),rep("Mag. outliers", 5), rep("Shape outliers",5), rep("Joint outliers",5))
depths$label=rep(ids,3*l*5)
depths$type=rep(rep(c("DepthGram on Dimensions","DepthGram on Time","DepthGram on Time/Correlation"),each=n),l*5)
depths$crate=rep(c(0,0.25,0.5,0.75,1),each=3*n*l)

    
depths$label=factor(depths$label, levels = c("Mag. outliers", "Shape outliers", "Joint outliers", "Rest of the sample"))
depths$crate=factor(depths$crate)
levels(depths$crate)= c("c = 0", "c = 0.25", "c = 0.5", "c = 0.75", "c = 1")
    

ggplot(depths, aes(x=1-mei, y=mbd) ) + stat_density_2d(aes(alpha = stat(nlevel),fill=label), geom = "polygon") +
    geom_density_2d(aes(colour=label)) +facet_grid(crate~ type)+ scale_fill_manual(name="Observations",values=colors) +
    scale_colour_manual(guide=FALSE,values=colors)+scale_alpha_continuous(name="Frequency on each group",range=c(0,1))+
    theme_light()+ theme(legend.position="bottom",strip.text = element_text(size=14, face = "bold"), axis.title= element_text(size=14),
                         axis.text= element_text(size=12),title=element_text(size=16),legend.text=element_text(size=12),legend.title = element_text(size=14))+ 
    xlab(expression(atop(1-MEI~(paste(bold('MBD') ) ) ) ) ) +
    ylab(expression(atop(MBD~(paste(bold('MEI') ) ) ) ) ) +
    ggtitle(paste("Simulation Summary - Model",model,", p=",p))
    
rm(depths)


###### Figures 7,8  #######            model=1, 2 for Figures 7, 8
###########################            model=3, 4 for Figures 8 and 9 in the Supplementary Materials
n=85       #nb. of non-outlying curves
n.out=15   #nb. of outlying curves
type=7     #all magnitude, shape and joint outliers
N=100      #nb. of time points
p=50       #dimension
c=1        #contamination proportion
model=1    #simulation model

### Data generation
Sim<-Sim_mfdata(n,N,p,n.out,c=c,model=model,type.out=type)

Data<-Sim$values

### colors for plots
color <- c(rep(8,n),color.mg,color.sh,color.jt) #colors for magnitude, shape and joint outliers defined above
ids=c(rep("",n),as.character((n+1):(n+n.out)))

#### Top Row - DepthGram 

### DepthGram
DG<-depthGram(Data)
### Plot
DGplot<-depthGramPlot(DG,col=color, limits=T)

#### Bottom Row - FOM and MSplot 
# USE OF mrfDepth package to implement FOM

### FOM
library(mrfDepth)               #library for FOM
source("other_functions.R")     #auxiliary functions for FOM and MSplot graphic display
#Functional p-variate version
X3d=array(dim=c(N,n+n.out,p))
for (k in 1:p){
  X3d[,,k]=t(Data[[k]])
}

FOM_X <-  mrfDepth::fOutl(x = X3d, type = "fDO", distOptions = list(rmZeroes = TRUE,maxRatio = 3),diagnostic = TRUE)
### Plot
fom_ggplot(FOM_X,col=color,cutoff = TRUE,subtit="p-dim",ids=ids,sp=3,st=5)

#Functional univariate version
X3d=array(dim=c(N*p,n+n.out,1))
for (k in 1:p){
  X3d[((k-1)*N+1):(k*N),,1]=t(Data[[k]])
}

FOM_X2 <-  mrfDepth::fOutl(x = X3d, type = "fDO", distOptions = list(rmZeroes = TRUE,maxRatio = 3),diagnostic = TRUE)
### Plot
fom_ggplot(FOM_X2, col=color,cutoff = TRUE,ids=ids,subtit="1-dim",sp=3,st=5)

### MSplot 
# USE OF functions and code in msplot_code folder to implement MSplot 
# (msplot_code folder as downloaded from  https://www.tandfonline.com/doi/suppl/10.1080/10618600.2018.1473781?scroll=top, supplement.rar)

source("msplot_code/loading.R")   #auxiliary functions for msplot: PACKAGES INSTALLATION MIGHT BE REQUIRED!
source("msplot_code/msplot.R")    #msplot function
source("msplot_code/DirOut.R")    
#Functional p-variate version
X3d=array(dim=c(n+n.out,N,p))  #Data reshaping into a 3d array
for (k in 1:p){
  X3d[,,k]=Data[[k]]
}
MS_X<-DirOut(data = X3d, depth.dir="RP",D_value=TRUE)
### Plot
msplot_ggplot(MS_X,col=color,ids=ids,fmdata=X3d,subtit="p-dim",sp=3,st=5)
#msplot(X3d,depth.dir="RP",dirout=T,col.normal=8)

#Functional univariate version
X3d=array(dim=c(n+n.out,N*p,1))
for (k in 1:p){
  X3d[,((k-1)*N+1):(k*N),1]=Data[[k]]
}
MS_X2<-DirOut(data = X3d[,,1], depth.dir="RP",D_value=TRUE)  #Outlyingness calculation
### Plot
msplot_ggplot(MS_X2,col=color,ids=ids,subtit="1-dim",sp=3,st=5)
#msplot(data = X3d[,,1], depth.dir="RP",dirout=T,col.normal=8)


######################## 
# See MOTOR_LANGUAGE_tFMRI_Examples.R for
# Figures 10 and 11 in the Supplementary Materials
########################

######################## 
# See Comparison_ComputingTimes.R for
# Figures 10 and 11 in the Supplementary Materials
########################
