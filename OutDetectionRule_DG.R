P2<-function(x,n){  #Equation of the theoretical parabola in DepthGram representations
  a0=2/n
  a2=-n/(2*(n-1))
  return(a0+x+a2*x^2)
}

out_det_DG<-function(DG,n){  #Basic outlier detection rule for the DepthGram
  
  dist.d=DG$mbd.mei.d-P2(1-DG$mei.mbd.d,n)
  q3=quantile(dist.d,0.75)
  q1=quantile(dist.d,0.25)
  out.d=which(dist.d>q3+1.5*(q3-q1))
  
  dist.t2=DG$mbd.mei.t2-P2(1-DG$mei.mbd.t2,n)
  q3=quantile(dist.t2,0.75)
  q1=quantile(dist.t2,0.25)
  out.t2=which(dist.t2>q3+1.5*(q3-q1))
  
  out=c(out.d,out.t2)
  
  return(out)
}
