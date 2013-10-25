plotA0D0=function(p,dd,save="pdf",live=FALSE) {
  ## save="pdf"; live=TRUE; load("CHECKPOINT.Rdata"); source("funcs.R")
	suppressMessages(require(IRanges))
	suppressMessages(require(plotrix))
	if(length(which(save%in%c("pdf","png")))==0) {
		stop("Unsupported plot type, only 'pdf' or 'png' are allowed!")
	}
	load(paste(p$Dir,"/",p$Project,"/",p$Project,".Rdata",sep=""))
	pos=dd$anno$plottable[dd$anno$chrom<=getSpeciesAutosome(p)]
	rem=which(dd$anno$polymorphic[dd$anno$chrom<=getSpeciesAutosome(p)])
	ii=apply(cohort$breaks[,6:7],1,function(idx) idx[1]:idx[2])
	qG=unlist(lapply(1:length(ii),function(idx,cohort,ii) rep(cohort$scored[idx,"A0"],length(ii[[idx]])),cohort,ii))
	qL=unlist(lapply(1:length(ii),function(idx,cohort,ii) rep(cohort$scored[idx,"D0"],length(ii[[idx]])),cohort,ii))
	qG[rem]=NA
	qL[rem]=NA

	dG=approx(pos,qG,n=length(pos))
	dG$y[c(1,length(dG$y))]=0
	dL=approx(pos,qL,n=length(pos))
	dL$y[c(1,length(dL$y))]=0

	if(!live) {
		if(save=="pdf") {
			pdf(file=paste(p$Dir,"/",p$Project,"/",p$Project,"Fraction",".pdf",sep=""),width=11,height=8.5)
		} else {
			png(file=paste(p$Dir,"/",p$Project,"/",p$Project,"Fraction",".png",sep=""),width=900,height=600)
		}
	} 

	makeChromPolygons=function(pos,dd,ymin,ymax) { 
          for(i in seq(1,22,2)) { 
            x=range(pos[dd$anno$chrom==i])
            polygon(c(x[1],x[1],x[2],x[2]),c(ymin,ymax,ymax,ymin),col="gray89",border="gray89")
          } 
	}
        
        labelChromPos <- function(pos,dd) {
          for(i in seq(1,22,1)) { 
            x=range(pos[dd$anno$chrom==i])
            text(mean(x),1,i)
          } 
	}
        

	layout(matrix(c(1,1,2,2,3,3),3,2,byrow=T),widths=rep(1,4),heights=c(0.33,0.02,0.33))
	par(xaxs="i",yaxs="i")

        ymax=0.25
	ybase=(seq(11)-1)/10 
	yhigh=ybase
	ylab=ybase
	par(mar=c(0,5,0.5,0.1))
	plot(1,1,type="n",xlim=range(dG$x),ylim=c(0,ymax),axes=FALSE,ylab="A0 (Gain)",xlab="")
	makeChromPolygons(pos,dd,0,ymax)
	idx=T
	lines(dG$x[idx],dG$y[idx],type="h",col="red3")
	axis(2,at=yhigh,labels=ylab,las=2,cex.axis=0.8)

        plot(1,1,type='n',xlab="",ylab="",axes=FALSE,xlim=range(dG$x))
        labelChromPos(pos,dd)
        
	par(mar=c(0.5,5,0,0.1))

        ybase=rev(-ybase)
	ylab=ybase
	yhigh=ybase
        
	plot(1,1,type="n",xlim=range(dL$x),ylim=c(-ymax,0),axes=FALSE,ylab="D0 (Loss)",xlab="")
	makeChromPolygons(pos,dd,-ymax,0)
	idx=T
	lines(dL$x[idx],-(dL$y[idx]),type="h",col="darkblue")
	axis(2,at=yhigh,labels=ylab,las=2,cex.axis=0.8)

	if(!live) dev.off()
}

loglog <- function(x) {
  ifelse(x>.25,NA,log(-log(x)))
}

makeChromPolygons=function(pos,dd,ymin,ymax) { 
  for(i in seq(1,22,2)) { 
    x=range(pos[dd$anno$chrom==i])
    polygon(c(x[1],x[1],x[2],x[2]),c(ymin,ymax,ymax,ymin),col="gray89",border="gray95")
  } 
} 


plotQ=function(p,dd,save="pdf",live=FALSE) {
  ## save="pdf"; live=TRUE; load("CHECKPOINT.Rdata"); source("funcs.R")
  suppressMessages(require(IRanges))
  suppressMessages(require(plotrix))
  if(length(which(save%in%c("pdf","png")))==0) {
    stop("Unsupported plot type, only 'pdf' or 'png' are allowed!")
  }
  load(paste(p$Dir,"/",p$Project,"/",p$Project,".Rdata",sep=""))
  pos=dd$anno$plottable[dd$anno$chrom<=getSpeciesAutosome(p)]
  rem=which(dd$anno$polymorphic[dd$anno$chrom<=getSpeciesAutosome(p)])
  ii=apply(cohort$breaks[,6:7],1,function(idx) idx[1]:idx[2])
  
  qG=unlist(lapply(1:length(ii),function(idx,cohort,ii) rep(cohort$q[idx,"Gain"],length(ii[[idx]])),cohort,ii))
  qL=unlist(lapply(1:length(ii),function(idx,cohort,ii) rep(cohort$q[idx,"Loss"],length(ii[[idx]])),cohort,ii))
  qG[rem]=NA
  qL[rem]=NA
  
  dG=approx(pos,qG,n=length(pos))
  dG$y[c(1,length(dG$y))]=1
  dL=approx(pos,qL,n=length(pos))
  dL$y[c(1,length(dL$y))]=1
  
  xRange=c(0,max(max(dL$x),max(dG$x)))
  
  scale.loglog=FALSE
  if(scale.loglog) {
    plot(dG$x,loglog(dG$y),type='h',col="#880000",
         ylim=c(-6,6),xlim=xRange,axes=F,xlab="",ylab="q-value")
    par(new=T)
    plot(dL$x,-loglog(dL$y),type='h',col="#000088",
         ylim=c(-6,6),xlim=xRange,axes=F,xlab="",ylab="")
    
    yAxis=2*(-3:3)
    axis(2,at=yAxis,labels=F)
    for(y in yAxis) {
      yExp=round(exp(abs(y))/log(10),0)
      if(y!=0) {
        text(par("usr")[1],y,call("^",10,-yExp),xpd=TRUE,pos=2,offset=1)
      } else {
        text(par("usr")[1],y,call("^",10,0),xpd=TRUE,pos=2,offset=1)
      }
    }
    box()
    abline(h=0,lwd=3,col=8)
  } else {
    plot(dG$x,-log10(dG$y),type='h',col="#880000", xlim=xRange,axes=F,xlab="",ylab="",ylim=c(-16,16))
    makeChromPolygons(pos,dd,-20,20)
    
    lines(dG$x,-log10(dG$y),type='h',col="#880000")
    lines(dL$x,log10(dL$y),type='h',col="#000088")
    yAxis=4*(-4:4)
    axis(2,at=yAxis,labels=F,line=-1.5)
    for(y in yAxis) {
      yExp=round(y,0)
      if(y!=0) {
        text(par("usr")[1],y,call("^",10,-yExp),xpd=TRUE)
      } else {
        text(par("usr")[1],y,call("^",10,0),xpd=TRUE)
      }
    }

    xbox=range(pos)
    ybox=par()$usr[3:4]
    polygon(c(xbox[1],xbox[1],xbox[2],xbox[2]),c(ybox[1],ybox[2],ybox[2],ybox[1]),lwd=2)
    mtext("q-value",2,line=2)
  }
}

##load("CHECKPOINT.Rdata")
source("funcs.R")

halt("INTERACTIVE")
plotA0D0(p,dd,save="png",live=F)
plotQ(p,dd)
dev.copy2pdf(file="test.pdf",width=11,height=8.5)
