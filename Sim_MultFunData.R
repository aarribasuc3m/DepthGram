Sim_mfdata<- function(n,N,p,n.out,c=1,model=1,type.out=7){ #Generation of multivariate functional data 
                                                           #with dependence across dimensions and outliers
  
  #n: number of observations
  #N: number of Time points
  #p: functional dimension
  #n.out: number of atypical observations
  #c: intensity of outlyingness, from 0 to 1 - proportion of dimensions in which the atypical observation is indeed an outlying curve
  #model: Simulation model: 1,2,3,4
  ##type.out: type of outliers 1 mag., 2 shape, 3 joint, 4 mag&shape, 5 mag&joint, 6 shape&joint, 7 mag&shape&joint
  
  ###########################
  #### 1. Number and ID of outliers
  ###########################
  if (type.out==7){
    n3=round(n.out/3)
    n.out=c(n3,n3,n.out-2*n3)
  }else{
    if (type.out>=4 ){
      n2=round(n.out/2)
      n.out=c(n2,n.out-n2)
    }
  }
  n.tot=n+sum(n.out)
  
  mag.out.id=c()
  sha.out.id=c()
  jt.out.id=c()
  n.out.m=0
  n.out.s=0
  n.out.j=0
  
  if (type.out %in% c(1,4,5,7)){
    n.out.m=n.out[1]    #nb. of magnitude outliers
    mag.out.id=(n+1):(n+n.out.m)  #indexes of magnitude outliers
  }
  if (type.out %in% c(2,4,6,7)){
    if(type.out==2){
      n.out.s=n.out[1]   #nb. of shape outliers
    }else{
      n.out.s=n.out[2]   #nb. of shape outliers
    }
    sha.out.id= (n+n.out.m+1):(n+n.out.m+n.out.s)  #indexes of shape outliers
  }
  if (type.out %in% c(3,5,6,7)){
    if(type.out==3 | type.out==6){
      n.out.j=n.out[1]                #nb. of joint outliers
    }else{
      if(type.out==5){
        n.out.j=n.out[2]
      }else{n.out.j=n.out[3]}
    }
    jt.out.id= (n+n.out.m+n.out.s+1):(n+n.out.m+n.out.s+n.out.j)  #indexes of joint outliers
  }
  
  
  ###########################
  #### 2. Simulation settings
  ###########################
  time_grid = seq( 0, 1, length.out = N )
  Data<-list(length=p)
  nc=floor(c*p)

  C1 =   outer(time_grid, time_grid, function(s, t) (0.3 * exp(-abs(s - t)/0.3)))
  CholCov = chol(C1)
  
  if (model==1 | model==3){
    f<-function(x){
      sin( 4 * pi * x )
   }
    fout<-function(x){
      cos( 4 * pi * x +pi/2)
    }
    h<-function(x,i,p){
      alpha=i/p
      1+2*x^(1+alpha)*(1-x)^(2-alpha)
    }
  }else{
    f<-function(x){
      4*x
   }
    fout<-function(x){
      4*x+2*sin(4*(x+0.5)*pi)
    }
    h<-function(x,i,p){
      alpha=i/p
      if (i%%2==0){
        return(-1-2*x^(1+alpha)*(1-x)^(2-alpha))
      }else{return(1+2*x^(1+alpha)*(1-x)^(2-alpha))}  
    }
  }
  
  centerline0 =t(matrix(rep(rnorm(n.tot),N),nrow=N,byrow=T)+f(time_grid))
  w.outs=c()
  w.outj=c()
  w.outm=c()
  if (c==1){
    if (n.out.s>0){
      centerline_out =t(matrix(rep(rnorm(n.out.s),N),nrow=N,byrow=T)+fout(time_grid))
      w.outs=1:p
    }
    if (n.out.m>0){w.outm=1:p}
    if (n.out.j>0){w.outj=1:p}
  }else{
    if(n.out.s>0 & c>0 ){
      centerline_out =t(matrix(rep(rnorm(n.out.s),N),nrow=N,byrow=T)+fout(time_grid))
      w.outs=matrix(0,nrow=n.out.s,ncol=nc)
      for (i in 1:n.out.s){
        w.outs[i,]=sample(1:p,nc,replace=FALSE,prob=rep(1/p,p))   #which components contain shape outlying curves 
      }                                                            
    }
    if(n.out.m>0 & c>0 ){
      w.outm=matrix(0,nrow=n.out.m,ncol=nc)
      for (i in 1:n.out.m){
        w.outm[i,]=sample(1:p,nc,replace=FALSE,prob=rep(1/p,p))   #which components contain magnitude outlying curves 
      }
    }
    if(n.out.j>0 & c>0 ){
      w.outj=matrix(0,nrow=n.out.j,ncol=nc)
      for (i in 1:n.out.j){
        w.outj[i,]=sample(1:p,nc,replace=FALSE,prob=rep(1/p,p))   #which components contain joint outlying curves 
      }
    }
  }
  
  if (model==3 | model==4){ #joint outliers move across dimensions "simmetrically around the center of the curve cloud"
    if(n.out.j>0 & c>0){
      m=apply(centerline0[1:n,],1,mean)
      ms=sort(m,decreasing = F,index.return=T)
      u=runif(1)
      nr2=ceiling(n.out.j/2)*(u<0.5)+floor(n.out.j/2)*(u>=0.5)
      q=c(runif(nr2,0.05,0.25),runif(n.out.j-nr2,0.75,0.95))
      ind<-round(n*q)
      index.outj<-ms$ix[ind] 
      aux<-centerline0[jt.out.id,]                                                    
      centerline0[jt.out.id,]<-centerline0[index.outj,]
      centerline0[index.outj,]<-aux
      m=apply(centerline0[1:n,],1,mean)
      ms=sort(m,decreasing = F,index.return=T)
      ind<-round(n*(1-q))
      inv.index.outj<-ms$ix[ind] 
    }
  }
  
  
  ###########################
  #### 3. Process Generation
  ###########################
  
  
  for (j in 1:p){ #Generation of the process dimension-wise
    
    Data[[j]]=t(t(centerline0)*h(time_grid,j,p))+ matrix(rnorm(n.tot*N), nrow = n.tot, ncol = N) %*% CholCov
    
    #Magnitude outliers
    if(n.out.m>0 & c==1){  
      Data[[j]][mag.out.id,]=Data[[j]][mag.out.id,]+10
    }else{
      if (n.out.m>0 & c>0){
        i=which(w.outm==j,arr.ind = T)[,1]
        if(length(i)>0){
            k=mag.out.id[i]
            Data[[j]][k,]=Data[[j]][k,]+10
        }
      }
    }
    
    #Shape outliers
    if (n.out.s>0 & c==1){  
      Data[[j]][sha.out.id,]=t(t(centerline_out)*h(time_grid,j,p))+ matrix(rnorm(n.out.s*N), nrow = n.out.s, ncol = N) %*% CholCov
    }else{
      if (n.out.s>0& c>0){
        i=which(w.outs==j,arr.ind = T)[,1]
        if(length(i)>0){
            k=sha.out.id[i]
            Data[[j]][k,]=centerline_out[i,]*h(time_grid,j,p)+ matrix(rnorm(N*length(i)), nrow = length(i), ncol = N) %*% CholCov
        }
      }
    }
    
    #Joint outliers
    if (model==1 | model==2){ #joint outliers component vary randomly across dimensions
       if (n.out.j>0& c==1){
          w=sample(1:n,n.out.j,replace=FALSE)
          Data[[j]][jt.out.id,]=t(t(centerline0[w,])*h(time_grid,j,p))+ matrix(rnorm(n.out.j*N), nrow = n.out.j, ncol = N) %*% CholCov
       }else{
          if (n.out.j>0& c>0){
             i=which(w.outj==j,arr.ind = T)[,1]
             if(length(i)>0){
               k=jt.out.id[i]
               w=sample(1:n,length(k),replace=FALSE)
               Data[[j]][k,]=centerline0[w,]*h(time_grid,j,p)+ matrix(rnorm(N*length(i)), nrow = length(i), ncol = N) %*% CholCov
             }
          }
        }
    }else{  #models 3 and 4: joint outliers components vary "simmetrically around the center of the curve cloud" across dimensions 
       if (n.out.j>0& c==1 & j%%2==1){
          Data[[j]][jt.out.id,]=t(t(centerline0[inv.index.outj,])*h(time_grid,j,p))+ matrix(rnorm(n.out.j*N), nrow = n.out.j, ncol = N) %*% CholCov
        }else{
           if (n.out.j>0 & c>0 & j%%2==1){
               i=which(w.outj==j,arr.ind = T)[,1]
               if(length(i)>0){
                 k=jt.out.id[i]
                 Data[[j]][k,]=centerline0[inv.index.outj[i],]*h(time_grid,j,p)+ matrix(rnorm(N*length(i)), nrow = length(i), ncol = N) %*% CholCov
               }
           }
        }
    }
    
  }
  
  return(list(values=Data,mag.out.id=mag.out.id,shape.out.id=sha.out.id,jt.out.id=jt.out.id,
              w.outm=w.outm,w.outs=w.outs,w.outj=w.outj))
  
}

  
