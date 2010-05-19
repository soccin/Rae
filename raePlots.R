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


load("CHECKPOINT.Rdata")
source("funcs.R")

halt("INTERACTIVE")
plotA0D0(p,dd,save="png",live=F)
