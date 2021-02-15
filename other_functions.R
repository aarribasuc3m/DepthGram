#### Other functions used in plots and through simulation settings

fom_ggplot<-function (fOutlResult, col=NULL,cutoff = FALSE,main="",log=F,plot=TRUE,ids=c(),subtit="",sp=3,st=5) 
{
  ###Own version of the mrfDepth::fom function (to return outliers ID's and match aesthetics with DepthGram)
  #fOutlResult: result of fOutl function 
  #cutoff: whether to include the "fence" on plots and restun outliers id's
  #log: whether to return the plot on log scale
  #subtit: subtitle text
  

    if (missing(fOutlResult)) {
      stop("Input argument fOutlResult is required.")
    }
    if (!is.list(fOutlResult)) {
      stop("fOutlResult must be a list returned from a call to fOutl.")
    }
    InputNames <- names(fOutlResult)
    if (!("fOutl" %in% class(fOutlResult))) {
      stop("fOutlResult must be a list returned from a call to fOutl.")
    }
    if (!("distType" %in% InputNames)) {
      stop("fOutlResult must be a list returned from a call to fOutl.")
    }
    if (!("weights" %in% InputNames)) {
      stop("fOutlResult must be a list returned from a call to fOutl.")
    }
    if (!("crossDistsX" %in% InputNames)) {
      stop(paste("fOutlResult must be a list returned from a call to fOutl", 
                 "with option diagnostic = TRUE"))
    }
    if (!("locOutlX" %in% InputNames)) {
      stop(paste("fOutlResult must be a list returned from a call to fOutl", 
                 "with option diagnostic = TRUE"))
    }
    NFunc <- nrow(fOutlResult$crossDistsX)
    NTObs <- ncol(fOutlResult$crossDistsX)
    AOValues <- fOutlResult$crossDistsX
    fAO <- AOValues %*% fOutlResult$weights
    sdAO = apply(AOValues, 1, FUN = function(y) sqrt(sum(fOutlResult$weights * 
                                                           (y - sum(fOutlResult$weights * y, na.rm = TRUE))^2, na.rm = TRUE)/(1 - 
                                                                                                                                1/length(fOutlResult$weights))))
    LocOutl <- rowSums(fOutlResult$locOutlX)
    PlotData <- data.frame(row.names = 1:NFunc)
    PlotData$fAO <- fAO
    PlotData$DispMeasure <- sdAO/(1 + fAO)
    xlabel<-fOutlResult$distType
    ylabel<-paste("sd(", fOutlResult$distType, ") / (1+", 
                  fOutlResult$distType, ")", sep = "")
    tit="FOM"
    
    
    if (is.null(col)){
      PlotData$colorvec <- rep("black", NFunc)
    }else{
      PlotData$colorvec <-col}
    if (cutoff) {
      CAO <- log(0.1 + sqrt((PlotData$fAO/median(PlotData$fAO))^2 + 
                              (PlotData$DispMeasure/median(PlotData$DispMeasure))^2))
      Fence <- qnorm(0.995) * mad(CAO) + median(CAO)
      theta <- seq(0, pi/2, length = (100))
      FenceData <- matrix(0, nrow = length(theta), ncol = 2)
      colnames(FenceData) <- c("x", "y")
      FenceData <- data.frame(FenceData)
      FenceData$x <- median(PlotData$fAO) * (exp(Fence) - 0.1) * 
        cos(theta)
      FenceData$y <- median(PlotData$DispMeasure) * (exp(Fence) - 
                                                       0.1) * sin(theta)
      if (is.null(col)){
        PlotData$colorvec <- rep("black", length(CAO))
        PlotData$colorvec[which(CAO > Fence)] <- "red"
      }
      out=which(CAO > Fence)
    }

  names(PlotData)<-c("x","y","col")
  
  
  if (log==T){
    xlabel=paste("log 1+",xlabel)
    ylabel=paste("log 1+",ylabel)
    tit=paste(tit,", log-scale")
    PlotData$x=log(1+PlotData$x)
    PlotData$y=log(1+PlotData$y)
  }
  
  pch=rep(1,NFunc) #empty circle; 16 solid circle
  pch[out]=19 #solid circle ;8 asterisk
  
  
  Plot<- ggplot(data=PlotData,aes(x,y))+geom_point(aes(colour=as.factor(col)),shape=pch,size=sp,show.legend=FALSE)+
    geom_text(data=PlotData,aes(x,y,colour=as.factor(col)),label=ids,hjust=-0.15, vjust=-0.15,size=st,show.legend=FALSE)+
    scale_colour_manual(values=sort(unique(col)))+xlab("fAO")+ylab("DispMeasure")+#labs(title=tit)+theme(plot.title = element_text(hjust = 0.5))
    theme_light()+ggtitle(paste("FOM",subtit)) +theme(axis.title= element_text(size=14),axis.text= element_text(size=12),title=element_text(size=14))
  
  
  if (cutoff) {
    if (log==T){
      FenceData$x=log(1+FenceData$x)
      FenceData$y=log(1+FenceData$y)
    }
    Plot <- Plot + geom_path(mapping = aes_string(x = "x", 
                                                  y = "y"), data = FenceData, color = "black", linetype = 2, size = 1)
  }
  Plot <- Plot + xlab(xlabel)
  Plot <- Plot + ylab(ylabel)
  # Plot <- Plot + mrfDepth_theme() + guides(shape = guide_legend(title = "z"))
  if (plot==TRUE){
    Plot
  }
  
  
  return(list(Plot=Plot,PlotData=PlotData,out=out))
  
}
source("msplot_code/facCal_num.R")
library(abind)

msplot_ggplot<-function(result,col=NULL,log=F,plot=TRUE,ids=c(),fmdata=NULL,subtit="",sp=3,st=5){
  ###Own version of the msplot function (to match aesthetics with DepthGram and provide log-scale plots)
  #result: result of DirOut function 
  #log: whether to return the plot on log scale
  #subtit: subtitle text
  
  
  mo=result$out_avr  ###MS-plot (norm of MO vs VO; extracted from msplot.R lines 218 to 227 and 239 to 244)
  vo=result$out_var
  n=length(vo)
  if (is.null(dim(mo))){
    
    
    D=result$D
    factor=facCal_num(n,2) 
    fac1=factor$fac1
    cutoff1=factor$fac2 #cut off value for testing/outlier detection#
    num=sum(fac1*D>cutoff1) #number of outliers#
    cutoff=cutoff1/fac1
    
    mo <-result$out_avr
    vo <-result$out_var
    out.dir=which(result$D>cutoff)
    medcurve=which.min(result$D)
    
    M=cbind(mo,vo)
    ans=cov.rob(M,method="mcd")
    L=solve(chol.default(solve(ans$cov)))
    theta=1:200/199
    circle=abind(sin(theta*2*pi),cos(theta*2*pi),along=2)
    x=ans$center[1]+(cutoff^(1/2)*L%*%t(circle))[1,]
    y=ans$center[2]+(cutoff^(1/2)*L%*%t(circle))[2,]
    elip.data=data.frame(x=x,y=y,a=rep(1,length(x)))
    

    if (is.null(col)){
      col <- rep("black", n)
      col[out.dir] <- "red"
    }
      
    col.point=col
    pch=rep(1,n) #empty circle; 16 solid circle
    pch[out.dir]=19 #solid circle ;8 asterisk
    
    ms.data=data.frame(x=mo,y=vo,out=col.point,pch=pch)
    
    Plot<-ggplot(data=ms.data,aes(x=x,y=y))+geom_point(col=col.point,shape=pch,size=sp)+
      geom_text(data=ms.data,aes(x=x,y=y),col=col.point,label=ids,hjust=-0.15, vjust=-0.15,size=st,show.legend=FALSE)+
      geom_path(data=elip.data,aes(x=x,y=y),show.legend=FALSE,colour="black")+
      xlab("MO")+ylab("VO")+#labs(title=paste("MS-Plot",subtit))+theme(plot.title = element_text(hjust = 0.5))
      theme_light()+ggtitle(paste("MS-Plot",subtit)) +theme(axis.title= element_text(size=14),axis.text= element_text(size=12),title=element_text(size=14))
    
    
    
    
    
    MSP=c()
  }else{
    MSP<-msplot(data = fmdata, depth.dir="RP",dirout=T,plot=FALSE)
    out.dir<-MSP$out.dir
    
    mo=MSP$mo  
    vo=MSP$vo
    MO<-(apply(mo^2,1,sum))^(1/2)
    #ids=as.character(1:length(MO))
    if (is.null(col)){
      col <- rep("black", n)
      col[out.dir] <- "red"
    }
    col.point=col
    pch=rep(1,n) #empty circle; 16 solid circle
    pch[out.dir]=19 #solid circle ;8 asterisk
    
    ms.data=data.frame(x=MO,y=vo) 
    if (log==T){
      ms.data=data.frame(x=log(MO+1),y=log(vo+1)) 
      
      Plot<-ggplot(data=ms.data,aes(x=x,y=y))+geom_point(col=col.point,shape=pch,size=sp)+
        geom_text(data=ms.data,aes(x=x,y=y),col=col.point,label=ids,hjust=-0.15, vjust=-0.15,size=st,show.legend=FALSE)+
        xlab("log 1+||MO||")+ylab("log 1+VO")+#labs(paste("MS-Plot",subtit))+theme(plot.title = element_text(hjust = 0.5))
        theme_light()+ggtitle(paste("MS-Plot",subtit,"(log-scale)")) +theme(axis.title= element_text(size=14),axis.text= element_text(size=12),title=element_text(size=14))
      
    }else{
      Plot<-ggplot(data=ms.data,aes(x=x,y=y))+geom_point(col=col.point,shape=pch)+
        geom_text(data=ms.data,aes(x=x,y=y),col=col.point,label=ids,hjust=-0.15, vjust=-0.15,size=st,show.legend=FALSE)+
        xlab("||MO||")+ylab("VO")+    theme_light()+ggtitle(paste("MS-Plot",subtit)) +theme(axis.title= element_text(size=14),axis.text= element_text(size=12),title=element_text(size=14))
      
      
      
    }
  }
  
  if (plot==TRUE){
    Plot
  }
  
  return(list(Plot=Plot,PlotData=ms.data,MSP=MSP))
}

