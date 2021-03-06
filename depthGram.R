depthGram<- function(Data,marg.out=F,FB=1.5,FO=1.5,ids=NULL){  #Function to compute (not plot) the 3 DepthGram representations 
                                                               #from a p-variate funcional data set
  
  ### Data is a list of length p (nb. of coordinates). Each element is an n*N matrix, n individuals, N time points
  ### marg.out: if TRUE, the function returns shape and magnitude outliers over each dimension
  ### FB: factor for boxplot rule in functional boxplot for marginal magnitude outlier detection (if marg.out==T)
  ### FO: factor for boxplot rule in outliergram for marginal magnitude outlier detection (if marg.out==T)
  ### ids: labels for individual observations
  
  p=length(Data)   
  n=nrow(Data[[1]])
  N=ncol(Data[[1]])
  
  if (marg.out==T){
    a2=a0=-2/(n*(n-1))
    a1=2*(n+1)/(n-1)
  }
  
  #########################
  # Dimension-wise
  #########################
  mbd.d<-array(0,dim=c(n,p))    #n*p matrix with mbd's on each dimension 
  mei.d<-array(0,dim=c(n,p))    #n*p matrix with mbd's on each dimension 
  mag.out.det<-list(length=p)   #list for marginal magnitude outliers detection
  shp.out.det<-list(length=p)   #list for marginal shape outliers detection
  n2=ceiling(n*0.5)
  
  corr.mei<-vector("numeric",length=p)  #vector containing the sign of correlation between mei in each dimension and the next one
  corr.mei[1]=1
  
  rmat.mat<-array(0,dim=c(n,N,p)) #Array with all observation ranks on each time-point/dimension
  wp<-c()  #storing components in which -1 transformation will be applied (since negative mei correlation exists)
  
  for (i in 1:p) { ### Over dimensions    

    x<-Data[[i]]
    #MBD and MEI computation on dimension i
    rmat=apply(t(x),1,rank) 
    rmat.mat[,,i] <- rmat
    down=rmat-1
    up=n-rmat
    mbd.d[,i]<-(rowSums(up*down)/N+n-1)/(n*(n-1)/2)  
    mei.d[,i]<-rowSums(up+1)/(n*N)
    #MEI correlation between i and i-1 dimensions
    if(i>1){corr.mei[i]=corr.mei[i-1]*sign(cor(mei.d[,i],mei.d[,i-1]))}
    if(corr.mei[i]==-1){
      wp=c(wp,i)
    }
    
    #Marginal outlier detection on dimension i: functional boxplot and outliergram with factors FB and FO respectively
    if (marg.out==T){
      
      index<-order(mbd.d[,i],decreasing=T)  
      center=x[index[1:n2],]
      inf=apply(center,2,min)
      sup=apply(center,2,max)
      dist=FB*(sup-inf)
      upper=sup+dist
      lower=inf-dist
      mag.out.det[[i]]<-which(colSums((t(x) <= lower) + (t(x) >= upper))>0)
      dist=(a0+a1*mei.d[,i]+a2*n^2*mei.d[,i]^2)-mbd.d[,i]
      q=quantile(dist,probs = c(0.25,0.75))
      lim=FO*(q[2]-q[1])+q[2]
      shp.out.det[[i]]<- which(dist>lim) 
      rm(index,center,inf,sup,dist,upper,lower,q,lim)
    }
    
    rm(x,rmat,down,up)
  }
  rm(Data)
  

  ##### MEI of MBD and MBD of MEI across dimensions
  
  mei.mbd.d<-roahd::MEI(mbd.d) #managing ties
  mbd.mei.d<-roahd::MBD(mei.d)
  
  ########################
  # Time-wise
  ########################
  mbd.t<-array(0,dim=c(n,N))    #n*N matrix with mbd's on each time point 
  mei.t<-array(0,dim=c(n,N))    #n*N matrix with mei's on each time point 
  mbd.t2<-array(0,dim=c(n,N))   #n*N matrix with mbd's on each time point for the "corrected" data set
  mei.t2<-array(0,dim=c(n,N))   #n*N matrix with mbd's on each time point for the "corrected" data set
  
  
  for (i in 1:N) { ### Over time points
    
    rmat=rmat.mat[,i,]          #Getting observation ranks at time-point i
    #MBD and MEI computation on time-point i
    down=rmat-1
    up=n-rmat
    mbd.t[,i]<-(rowSums(up*down)/p+n-1)/(n*(n-1)/2)  
    mei.t[,i]<-rowSums(up+1)/(n*p)
    
    #MBD and MEI computation on dimension i for the "corrected" data set
    if(length(wp)>0){
      down[,wp]=n-rmat[,wp]
      up[,wp]=rmat[,wp]-1
      mbd.t2[,i]<-(rowSums(up*down)/p+n-1)/(n*(n-1)/2)  
      mei.t2[,i]<-rowSums(up+1)/(n*p)
    }else{
      mbd.t2[,i]<-mbd.t[,i]
      mei.t2[,i]<-mei.t[,i]
    }
    
    rm(down,rmat,up)
  }
  
  
  ##### MEI of MBD and MBD of MEI across time-points (original and corrected data sets)
  mei.mbd.t<-roahd::MEI(mbd.t)
  mei.mbd.t2<-roahd::MEI(mbd.t2)
  mbd.mei.t<-roahd::MBD(mei.t)
  mbd.mei.t2<-roahd::MBD(mei.t2)
  
  
  ##### RETURN
  list.return<- list(mbd.mei.d=mbd.mei.d,mei.mbd.d=mei.mbd.d,
                     mbd.mei.t=mbd.mei.t,mei.mbd.t=mei.mbd.t,
                     mbd.mei.t2=mbd.mei.t2,mei.mbd.t2=mei.mbd.t2,
                     shp.out.det=shp.out.det,mag.out.det=mag.out.det,
                     mbd.d=mbd.d, mei.d=mei.d,mbd.t=mbd.t,mei.t=mei.t,mbd.t2=mbd.t2,mei.t2=mei.t2)
  
  if(!is.null(ids)){
    list.return<- lapply(list.return, function(x) {
      if (is.list(x)){
          lapply(x,function(y){names(y)<-ids[y];y})
      }else{
      if (!is.null(dim(x))){row.names(x)<-ids}else{if(length(x)<n){names(x)<-ids[x]}else{names(x)<-ids}}}
      return(x) 
      })
  } 
  
  list.return$corr.mei<-corr.mei
  
  return(list.return)
}

