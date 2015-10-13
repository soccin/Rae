#####################################################################
###
### IO functions
###
#####################################################################

options(error=traceback)

# --------------------------------------------------------------------
canWrite=function(path) {
	res=(file.access(path,mode=2)==0)
	names(res)=NULL
	return(res);
}

# --------------------------------------------------------------------
isDirectory=function(path) {
	if(!file.exists(path)) {
		return(FALSE)
	}
	path=gsub("[/\\\\]$","",path)
	info=file.info(path)
	if(identical(info$isdir,TRUE)) {
		return(TRUE)
	} else {
		return(FALSE)
	}
}

# --------------------------------------------------------------------
isExisting=function(path) file.exists(path)

# --------------------------------------------------------------------
isFile=function(path) {
	if(!file.exists(path)) return(FALSE)
	info=file.info(path)
	!info$isdir;
}

# --------------------------------------------------------------------
mkdir=function(path) {
	res=dir.create(path)
	return(res)
}

# --------------------------------------------------------------------
remove=function(path) {
	if(!isExisting(path)) return(FALSE)
	if(!isDirectory(path)) return(file.remove(path))
	if(length(list.files(path))!=0) return(FALSE)
	(unlink(path,recursive=TRUE)==0)
}

#####################################################################
###
### Analysis modules
###
#####################################################################

# -------------------------------------------------------------------
# gatherProjectSegmentation(params)
#
# Given individual-sample segmentation by CBS and subsequent single-
# sample Rdata files with a single variable named 'd' storing the
# output of DNAcopy segment comment, this will gather and normalize
# segmentation into a final data structure used throughout the
# balance of RAE.
# -------------------------------------------------------------------

gatherProjectSegmentation=function(p) {
	loc=paste(p$Dir,p$Project,sep="/")
	if(!isDirectory(paste(p$Dir,p$Project,sep="/")) & canWrite(p$Dir)) mkdir(loc)
	cbs=dir(p$Data,pattern="Rdata",full.names=TRUE)

	cat("   >> Generating platform annotation from randomly chosen sample...\n")
	chr=seqlengths(Genome)
	chr=chr[-grep("(random|hap|M)",names(chr))]
	load(cbs[sample(1:length(cbs),1)])

	anno=list()
	anno$chrom=d$data[,"chrom"]
	anno$maploc=d$data[,"maploc"]
	anno$plottable=rep(NA,length(anno$chrom))
	anno$polymorphic=rep(FALSE,length(anno$chrom))
	for(i in unique(anno$chrom)) {
		pidx=which(anno$chrom==i)
		anno$plottable[pidx]=(anno$maploc[pidx]/1000000)+ifelse(i==1,0,sum(chr[1:(i-1)]/1000000))
		if(grepl("hg",p$Build)) {; # Don't currently have CNV repository for mouse, could build one, but later...
			probes=IRanges(start=anno$maploc[pidx],end=anno$maploc[pidx])
			poly=IRanges(start=cnv$Start[cnv$Chr==i],end=cnv$End[cnv$Chr==i])
			anno$polymorphic[pidx[unique(matchMatrix(findOverlaps(poly,probes))[,1])]]=TRUE
		}
	}

	dd=NULL
	dd$build=p$Build
	dd$anno=as.data.frame(anno)
	dd$data=dd$anno[!(dd$anno[,1]>getSpeciesAutosome(p)),1:2]
	dd$segs=NULL
	dd$output=NULL
	dd$dn=NULL
	dd$offset=NULL
	dd$mass_median=NULL
	dd$diploid_peak_status=NULL

	if(p$SNP6) bounds=getStartsAndEnds(dd$anno)
	pbound=getStartsAndEndsP(dd$anno)

	out=paste(p$Dir,"/",p$Project,"/cbs-gathered-",p$Project,"-",timeStamp(),sep="")
	pdf(file=paste(out,".pdf",sep=""),width=8,height=11)
	par(mfrow=c(2,2))

	cat("   >> Processing individual sample-segmentation for",length(cbs),"samples...\n")
	for(i in 1:length(cbs)) {
		cat("        ...",gsub("^.*[/]","",cbs[i]),":",sep="")
		load(cbs[i])
		d$output$ID=colnames(d$data)[3]; # Verify segments labeled by sample ID
		if(p$SNP6) d=extendStartAndEnd(d,bounds); # Handle NA probes when segment boundaries
		if(nrow(d$data)!=nrow(dd$anno)) {
			cat("Samples not all run on the same platform processed similarly!\n")
			stop("Probe count mismatch!")
		}
		ds=processCBS(d,p,pbound)
		dd$data=cbind(dd$data,ds$data[,3,drop=FALSE])
		dd$segs=cbind(dd$segs,ds$segs[,3,drop=FALSE])
		dd$output=rbind(dd$output,ds$output)
		dd$offset=c(dd$offset,ds$offset)
		dd$dn=c(dd$dn,ds$dn)
		dd$mass_median=c(dd$mass_median,ds$mass_median)
		dd$diploid_peak_status=c(dd$diploid_peak_status,ds$diploid_peak_status)
	}
	dev.off()
	cat("\n")
	cat("   >> Writing/saving output...\n")
	tmp=dd$output[,1:6]
	tmp$loc.start=as.integer(tmp$loc.start)
	tmp$loc.end=as.integer(tmp$loc.end)
	write.table(tmp,file=paste(gsub("/cbs","/IGV-cbs",out),".seg",sep=""),sep="\t",eol="\n",quote=F,row.names=F)
	write.table(
		cbind(p),
		file=paste(p$Dir,p$Project,"param.input.log",sep="/"),
		sep="=",eol="\n",quote=F,row.names=T,col.names=F
	)
	save(dd,file=paste(out,".Rdata",sep=""),compress=T)
	rm(anno,d,ds,cbs)
	gc()
	return(dd)
}

# -------------------------------------------------------------------
# parameterizeMultiComponentModel(params,dd)
#
# Uses segmentation to derive parameters (E and Beta) for each of the
# four discriminators in each sample based on the density distribution
# of total autosomal segmentation
# -------------------------------------------------------------------

parameterizeMultiComponentModel=function(p,dd) {
	ds=dd
	rm(dd)
	gc()

	ds$segs=cbind(ds$data[,1:2],ds$segs)
	ds$data=ds$data[!(ds$anno[,1]>getSpeciesAutosome(p)),]
	ds$segs=ds$segs[!(ds$anno[,1]>getSpeciesAutosome(p)),]
	ds$output=ds$output[!(ds$output$chrom>getSpeciesAutosome(p)),]
	ds=compute.rsegs(ds)

	out=paste(p$Dir,"/",p$Project,"/params-",p$Project,"-",timeStamp(),sep="")
	pdf(file=paste(out,".pdf",sep=""),height=8,width=8)
	par(mfrow=c(2,2))
	parm=list()

	for(i in 3:ncol(ds$segs)) {
		dn=dlrs(ds$data[,i])
		pd=density(ds$rsegs[,i],na.rm=T,from=-1,to=1,n=N.density,adj=2)
		plot(pd,main=colnames(ds$rsegs)[i],xlim=c(-1,1))
		abline(v=c(-dn,dn),col=2)
		abline(v=0,col=8)
		fwhm=get.fwhm(pd)
		abline(v=fwhm$px,col=8,lty=2)
		abline(v=fwhm$nx,col=8,lty=2)

		### New exclusion rules based on TCGA data
		if(!(fwhm$nx<0 & 0<fwhm$px)) {
			# Excluding if the full-width at half max does not span zero
			# This can be considered a fatal error
			parm$Ed0=c(parm$Ed0,NA)
			parm$Ea0=c(parm$Ea0,NA)
			parm$betaD0=c(parm$betaD0,NA)
			parm$betaA0=c(parm$betaA0,NA)
			parm$Ed1=c(parm$Ed1,NA)
			parm$Ea1=c(parm$Ea1,NA)
			parm$betaD1=c(parm$betaD1,NA)
			parm$betaA1=c(parm$betaA1,NA)
			parm$names=c(parm$names,colnames(ds$data)[i])
			parm$fwhm.nx=c(parm$fwhm.nx,NA)
			parm$fwhm.px=c(parm$fwhm.px,NA)
			parm$SN.fwhm=c(parm$SN.fwhm,NA)
			parm$SEGN=c(parm$SEGN,NA)
			next
		}

		# Calculate measure of segment noise, looking at adjacent jumps
		segn=unlist(lapply(
			1:getSpeciesAutosome(p),
			function(ch) diff(ds$output$seg.mean[(ds$output$ID==colnames(ds$segs)[i] & ds$output$chrom==ch)])
		))
		segn=median(abs(segn))/dn

		N=length(pd$y)
		pi=peaksign(pd$y,13)
		px=which(pi!=0)
		ps=pi[px]
		pm=px[which(ps>0)]
		np=pm[which(pm<fwhm$dpi)]
		np1=np[which(pd$x[np]> -2*dn & pd$x[np]< -dn/2 & pd$x[np] < fwhm$nx)]

		if(length(np1)>0) {
			bnp=np1[which.max(pd$y[np1])]
			points(pd$x[bnp],pd$y[bnp],col=3,pch=3)
			Ed0=(pd$x[bnp]+0)/2
			f=function(x,E0,tol) fer(fwhm$nx,E0,x)-tol
			while(f(-100,Ed0,BETA.tol)>0) Ed0=Ed0*1.1
			betaD0=uniroot(f,c(-0.001,-100),E0=Ed0,BETA.tol)$root
			if(abs(betaD0)<10) betaD0=-10
			col="#993333"
		} else {
			beta.fix=-15
			g=function(x,beta,tol) fer(fwhm$nx,x,beta)-tol
			Ed0=uniroot(g,c(-10,10),beta=beta.fix,BETA.tol)$root
			betaD0=beta.fix
			col="#333399"
		}
		betaD1=betaD0
		Ed1=quantile(ds$data[,i],c(0.025),na.rm=TRUE)
		if(Ed1>1.4*Ed0) Ed1=quantile(ds$data[,i],c(0.01),na.rm=TRUE)
		if(Ed1>1.4*Ed0) Ed1=1.4*Ed0

		xx=(-100:100)/100
		par(new=TRUE)
		plot(xx,fer(xx,Ed0,betaD0),xlim=c(-1,1),ylim=c(0,1),type='l',col=col,yaxt='n',xlab='',ylab='')
		text(-1,0.95,pos=4,paste("bD0=",round(betaD0,1)))
		text(-1,0.90,pos=4,paste("Ed0=",round(Ed0,3)))
		text(-1,0.70,pos=4,paste("fwhm.nx",round(fwhm$nx,4)))
		par(new=TRUE)
		plot(xx,fer(xx,Ed1,betaD1),xlim=c(-1,1),ylim=c(0,1),type='l',col=col,yaxt='n',xlab='',ylab='',lty=3)
		text(-1,0.55,pos=4,paste("bD1=",round(betaD1,1)))
		text(-1,0.50,pos=4,paste("Ed1=",round(Ed1,3)))

		pp=pm[which(pm>fwhm$dpi)]
		pp1=pp[which(pd$x[pp] < 2*dn & pd$x[pp] > dn/2 & pd$x[pp] > fwhm$px)]
		if(length(pp1)>0) {
			bpp=pp1[which.max(pd$y[pp1])]
			points(pd$x[bpp],pd$y[bpp]/max(pd$y),col=3,pch=3)
			Ea0=(pd$x[bpp]+fwhm$px)/2
			f=function(x,E0,tol)fer(fwhm$px,E0,x)-tol
			while(f(100,Ea0,BETA.tol)>0) Ea0=Ea0*1.1
			betaA0=uniroot(f,c(0.001,100),E0=Ea0,BETA.tol)$root
			if(abs(betaA0)<10) betaA0=10
			col="#993333"
		} else {
			beta.fix=20
			g=function(x,beta,tol)fer(fwhm$px,x,beta)-tol
			Ea0=uniroot(g,c(-10,10),beta=beta.fix,BETA.tol)$root
			betaA0=beta.fix
			col="#333399"
		}
		Ea1=(2*Ea0)-fwhm$px
		if(Ea1<1.4*Ea0) Ea1=1.4*Ea0
		betaA1=1/log(2)

		par(new=TRUE)
		plot(xx,fer(xx,Ea0,betaA0),xlim=c(-1,1),ylim=c(0,1),type='l',col=col,yaxt='n',xlab='',ylab='')
		text(0.2,0.95,pos=4,paste("bA0=",round(betaA0,2)))
		text(0.2,0.90,pos=4,paste("Ea0=",round(Ea0,2)))
		par(new=TRUE)
		plot(xx,ferA1(xx,Ea1,betaA1),xlim=c(-1,1),ylim=c(0,1),type='l',col=col,yaxt='n',xlab='',ylab='',lty=2)
		text(0.2,0.75,pos=4,paste("bA1=",round(betaA1,2)))
		text(0.2,0.70,pos=4,paste("Ea1=",round(Ea1,2)))

		pdN=density(ds$rsegs[,i],na.rm=T,from=-2,to=2,n=N.density,adj=2)
		xx=pdN$x
		yy=pdN$y
		yd=pdN$y*fer(xx,Ed0,betaD0)
		ydd=pdN$y*fer(xx,Ed1,betaD1)
		yaa=pdN$y*ferA1(xx,Ea1,betaA1)
		ya=pdN$y*fer(xx,Ea0,betaA0)
		ym=pmax(yd,ya)
		SN.fwhm=sum(ym)/(abs(fwhm$nx)+fwhm$px)
		plot(xx,ym,main=colnames(ds$rsegs)[i],xlim=c(-2,2),type='l')
		lines(xx,ydd,col=2,lty=2)
		lines(xx,yaa,col=2,lty=2)

		parm$Ed0=c(parm$Ed0,Ed0)
		parm$Ea0=c(parm$Ea0,Ea0)
		parm$betaD0=c(parm$betaD0,betaD0)
		parm$betaA0=c(parm$betaA0,betaA0)
		parm$Ed1=c(parm$Ed1,Ed1)
		parm$Ea1=c(parm$Ea1,Ea1)
		parm$betaD1=c(parm$betaD1,betaD1)
		parm$betaA1=c(parm$betaA1,betaA1)
		parm$names=c(parm$names,colnames(ds$data)[i])
		parm$fwhm.nx=c(parm$fwhm.nx,fwhm$nx)
		parm$fwhm.px=c(parm$fwhm.px,fwhm$px)
		parm$SN.fwhm=c(parm$SN.fwhm,SN.fwhm)
		parm$SEGN=c(parm$SEGN,segn)
	}
	dev.off()

	parm=data.frame(parm)
	pp=parm
	save(pp,file=paste(out,".Rdata",sep=""))
	rm(ds)
	gc()
	return(pp)
}

# -------------------------------------------------------------------
# scoreProjectAnalysis(params,dd,pp,anno)
#
# Use segmentation and parameterization to generate the unified
# breakpoint profile, score the genome in each sample, combine those
# scores across all tumors and assign the tumor/region lesion map.
# -------------------------------------------------------------------

scoreProjectAnalysis=function(p,dd,pp) {

	cohort=NULL
	cohort$build=dd$build
	cohort$desc=paste("RAE analysis for",p$Project,"on",timeStamp())
	cohort$breaks=NULL
	cohort$scored=NULL

	segments=dd$output
	segments=segments[!(segments[,2]>getSpeciesAutosome(p)),]
	segments=segments[!(segments[,5]<3),]
	segments=cbind(segments[,1:4],NA,segments[,5:9])
	colnames(segments)[5]="Empty"

	###
	### Exclude samples for QC reasons
	###

	cat("   >> Checking for excluded samples...")
	samples=colnames(dd$segs)
	excl=getExcludedSamples(pp,p,LOG=FALSE)

	if(excl$status) {
		midx=match(excl$samples,samples)
		segments=segments[-which(segments$ID%in%excl$samples),]
		samples=samples[-midx]
		cat("excluding ",length(midx)," samples",sep="")
	} else {
		cat("no exclusions")
	}
	cat(", done\n")

	cat("   >> Generating breakpoint profile...")
  bps=getOutlierCorrectedSampleList(samples,segments)
  cohort$breaks=getProcessedBreakpoints(p,bps,segments,dd$anno)
	indices=apply(cohort$breaks[,6:7],1,function(b) b[1]:b[2])
  cat("found ",nrow(cohort$breaks)," breaks, done\n",sep="")
	rm(segments)
	gc()

	cat("   >> Transforming segmentation and scoring breakpoint profile...\n")
	mS=match(samples,colnames(dd$segs))
	mP=match(colnames(dd$segs)[mS],pp$names)
	lesions=NULL
	lesions$breaks=cohort$breaks

	cat("      >> Scoring low-level gains...")
	lesions$a0=lapply(indices,function(ii) apply(t(fer(t(dd$segs[ii,mS]),pp$Ea0[mP],pp$betaA0[mP])),2,getWindowWeightedMean))
	lesions$a0=as.matrix(do.call(rbind,lesions$a0))
	cat("done!\n")

	cat("      >> Scoring highly-level gains...")
	lesions$a1=lapply(indices,function(ii) apply(t(ferA1(t(dd$segs[ii,mS]),pp$Ea1[mP],pp$betaA1[mP])),2,getWindowWeightedMean))
	lesions$a1=as.matrix(do.call(rbind,lesions$a1))
	cat("done!\n")

	cat("      >> Scoring likely heterozygous losses...")
	lesions$d0=lapply(indices,function(ii) apply(abs(t(fer(t(dd$segs[ii,mS]),pp$Ed0[mP],pp$betaD0[mP]))),2,getWindowWeightedMean))
	lesions$d0=as.matrix(do.call(rbind,lesions$d0))
	cat("done!\n")

	cat("      >> Scoring likely homozygous deletions...")
	lesions$d1=lapply(indices,function(ii) apply(abs(t(fer(t(dd$segs[ii,mS]),pp$Ed1[mP],pp$betaD1[mP]))),2,getWindowWeightedMean))
	lesions$d1=as.matrix(do.call(rbind,lesions$d1))
	cat("done!\n")

	cat("   >> Summarizing scores across breakpoints...")
	colnames(lesions$a0)=colnames(lesions$a1)=colnames(lesions$d0)=colnames(lesions$d1)=samples
	save(lesions,file=paste(p$Dir,"/",p$Project,"/",p$Project,"-lesions.Rdata",sep=""))
	cohort$scored=cbind(
		rowMeans(lesions$a0,na.rm=T),
		rowMeans(lesions$a1,na.rm=T),
		rowMeans(lesions$d0,na.rm=T),
		rowMeans(lesions$d1,na.rm=T)
	)
	colnames(cohort$scored)=c("A0","A1","D0","D1")
	cat("done\n")

	cat("   >> Calculating analytical p-values...")
	A=sqrt(lesions$a0^2 + lesions$a1^2)
	D=sqrt(lesions$d0^2 + lesions$d1^2)
	cohort$p=cbind(
		calc.pv(rowSums(A,na.rm=T),fit.pv.params(A),logT=TRUE),
		calc.pv(rowSums(D,na.rm=T),fit.pv.params(D),logT=TRUE),
		calc.pv(rowSums(lesions$a1,na.rm=T),fit.pv.params(lesions$a1),logT=TRUE),
		calc.pv(rowSums(lesions$d1,na.rm=T),fit.pv.params(lesions$d1),logT=TRUE)
	)
	cohort$p=apply(cohort$p,2,function(p) ifelse(exp(p)==0,MINDBL,exp(p)))
	cohort$q=apply(cohort$p,2,p.adjust,method="BH")
	colnames(cohort$p)=colnames(cohort$q)=c("Gain","Loss","Amplification","Deletion")
	cat("done\n")

	muA=mean(A,na.rm=T)
	muD=mean(D,na.rm=T)
	convergeparamA=(0.71*mean((A-muA)^3,na.rm=T))/((sqrt(mean((A-muA)^2,na.rm=T))^3)*sqrt(length(samples)))
	convergeparamD=(0.71*mean((D-muD)^3,na.rm=T))/((sqrt(mean((D-muD)^2,na.rm=T))^3)*sqrt(length(samples)))
	if(max(convergeparamA,convergeparamD)>1) {
		cat("   >> WARNING: Gaussian is poor approximation for analytical p-value calculation!\n")
	}
	if(length(samples)<20) {
		cat("   >> WARNING: Likely too few samples in analysis to trust analytical p- and q-values!\n")
	}

	cat("   >> Bootstrapping for error on analytical p-values...")
	bs=NULL
	bs$A=NULL
	bs$D=NULL
	N=ncol(A)

	rand=lapply(1:100,function(x) sample(1:N,floor(N*0.75)))
	for(i in 1:length(rand)) {
		rA=A[,rand[[i]]]
		rD=D[,rand[[i]]]
		bs$A=cbind(bs$A,-log10(calc.pv(rowSums(rA,na.rm=T),fit.pv.params(rA))))
		bs$D=cbind(bs$D,-log10(calc.pv(rowSums(rD,na.rm=T),fit.pv.params(rD))))
	}
	cohort$scored=cbind(cohort$scored,apply(bs$A,1,sd,na.rm=T),apply(bs$D,1,sd,na.rm=T))
	colnames(cohort$scored)[5:6]=c("DeltaGain","DeltaLoss")
	cohort$scored[is.na(cohort$scored[,"DeltaGain"]),"DeltaGain"]=MAXDBL
	cohort$scored[is.na(cohort$scored[,"DeltaLoss"]),"DeltaLoss"]=MAXDBL
	cat("done\n")

	save(cohort,file=paste(p$Dir,"/",p$Project,"/",p$Project,".Rdata",sep=""),compress=T)

	dd$cohort=cohort
	dd$lesions=lesions

	return(dd)
}

# -------------------------------------------------------------------
getRegionsOfInterest=function(p,dd,pp,errorFactor=1) {
	load(paste(p$Dir,"/",p$Project,"/",p$Project,"-lesions.Rdata",sep=""))
	load(paste(p$Dir,"/",p$Project,"/",p$Project,".Rdata",sep=""))

	# Only quiet human polymorphisms
	if(grepl("hg",p$Build)) cohort=redactCNV(cohort,p)

	# Get significant genomic events
	gain=getSignificantRegionsGain(cohort,p)
	loss=getSignificantRegionsLoss(cohort,p)
	if(is.null(gain) & is.null(loss)) {
		cat("There exists no regions of statistically significant aberrations!\n")
		return(NULL)
	}

	# Process the events into regions and peaks
	model=list()
	if(!is.null(gain)) {
		results=getStageOneRegions(gain,p)
		model=processStageResults(p,results,model,gain,errorFactor,"Amp")
	}
	if(!is.null(loss)) {
		results=getStageOneRegions(loss,p)
		model=processStageResults(p,results,model,loss,errorFactor,"Del")
	}
	model=as.data.frame(model)
	model=model[order(model$chromosome,model$region_start,model$peak_start),]

	# Annotate the events
	model=getEventFrequency(model,cohort,lesions,0.9,0.25)
	if(grepl("hg",p$Build)) {
		model=getNoncodingGeneContent(model)
		model=getProteinCodingGeneContent(model)
	}
	cols=colnames(model)
	cols=c(cols[-11],"SamplesWithSpanningAlteration")
	model=model[,cols]
	rownames(model)=NULL

	# Annotate sample membership
	model=getSampleMembership(model,dd,pp)

	# Write final, minimally annotated output
	write.table(
		model,file=paste(p$Dir,"/",p$Project,"/",p$Project,"-summary-lesions.txt",sep=""),
		sep="\t",eol="\n",quote=F,row.names=F
	)
}

# -------------------------------------------------------------------
getSampleMembership=function(model,dd,pp,perc_arm_len=0.75,min_overlap_len=0.05) {
  cyto=feat[feat$Type=="Cytoband",]
  cen=sapply(unique(cyto$Chr),function(ch) {
    idx=which(cyto$Chr==ch)
    m=idx[max(grep("p",cyto$Name[idx]))]
    as.integer(cyto$End[m])
  })
  names(cen)=unique(cyto$Chr)
  model$SamplesWithQualifiedSpanningAlteration=rep(NA,nrow(model))

									# Get segmentation
  seg=dd$output
  seg=seg[!(seg$chrom>getSpeciesAutosome(p)),]
  seg=seg[!(seg$num.mark<3),]
  samples=unique(as.character(as.vector(seg$ID)))

	  								# Transform segmentation
  pidx=match(seg$ID,pp$names)
  seg$A0=fer(seg$seg.mean,pp$Ea0[pidx],pp$betaA0[pidx])
  seg$A1=ferA1(seg$seg.mean,pp$Ea1[pidx],pp$betaA1[pidx])
  seg$D0=abs(fer(seg$seg.mean,pp$Ed0[pidx],pp$betaD0[pidx]))
  seg$D1=abs(fer(seg$seg.mean,pp$Ed1[pidx],pp$betaD1[pidx]))

  rem=NULL
  for(i in 1:nrow(model)) {

                                        # Use peak boundaries if they exist, otherwise regional
    st=ifelse(is.na(model$peak_start[i]),model$region_start[i],model$peak_start[i])
    ed=ifelse(is.na(model$peak_end[i]),model$region_end[i],model$peak_end[i])
    ch=model$chromosome[i]

                                        # Segments and samples contributing signal to region
    this=NULL
    if(model$model[i]=="Amp") {
      this=seg[which(seg$chrom==ch & seg$loc.end>=st & seg$loc.start<=ed & seg$A0>=0.9),]
    } else {
      this=seg[which(seg$chrom==ch & seg$loc.end>=st & seg$loc.start<=ed & seg$D0>=0.9),]
    }
	size=this$loc.end-this$loc.start+1
  	arm=ifelse(st>=cen[ch],length(Genome[[ch]])-cen[ch],cen[ch])
  	this=this[!(size/arm)>perc_arm_len,]
  	this=this[order(this$ID),]
  	if(nrow(this)==0) {
      model$SamplesWithQualifiedSpanningAlteration[i]=NA
      rem=c(rem,i)
    } else {
      ov=sapply(1:nrow(this),function(ii,this,st,ed) overlap2(st,ed,this$loc.start[ii],this$loc.end[ii]),this,st,ed)
      samp=unique(this$ID)
      fin=unlist(lapply(samp,function(s,this,ov) sum(ov[this$ID==s]),this,ov))
      model$SamplesWithQualifiedSpanningAlteration[i]=paste(samp[(fin/(ed-st+1))>min_overlap_len],collapse=",")
    }
  }
  if(length(rem)>0) model=model[-rem,]
  return(model)
}

# -------------------------------------------------------------------
getProteinCodingGeneContent=function(model) {
	isoforms=do.call(RangesList,lapply(unique(model$chromosome),function(ch) {
		idx=which(feat$Type=="Gene" & feat$Chr==ch)
		IRanges(start=feat$Start[idx],end=feat$End[idx],names=feat$Name[idx])
	}))
	model=getGenicContent(model,isoforms)
	N=ncol(model)
	colnames(model)[c(N-1,N)]=c("n_genes","genes")
	return(model)
}

# -------------------------------------------------------------------
getNoncodingGeneContent=function(model) {
	isoforms=do.call(RangesList,lapply(unique(model$chromosome),function(ch) {
		idx=which(feat$Type=="MicroRNA" & feat$Chr==ch)
		IRanges(start=feat$Start[idx],end=feat$End[idx],names=feat$Name[idx])
	}))
	model=getGenicContent(model,isoforms)
	model$content[grep("Nearest",model$content)]="None"
	N=ncol(model)
	colnames(model)[c(N-1,N)]=c("n_microRNAs","microRNAs")
	return(model)
}

# -------------------------------------------------------------------
getGenicContent=function(model,isoforms) {
	chr_idx=match(model$chromosome,unique(model$chromosome))
	for(i in 1:nrow(model)) {
		region=IRanges(
			start=ifelse(is.na(model$peak_start[i]),model$region_start[i],model$peak_start[i]),
			end=ifelse(is.na(model$peak_end[i]),model$region_end[i],model$peak_end[i])
		)
		found=matchMatrix(findOverlaps(isoforms[[chr_idx[i]]],region))[,2]
		if(length(found)>0) {
			content=unique(gsub("^.*:","",names(isoforms[[chr_idx[i]]])[found]))
			model$n_content[i]=length(content)
			# Arbitrary gene count to redact to prevent large gene lists in output
			model$content[i]=ifelse(length(content)>200,"Redacted for size",paste(content,collapse=","))
		} else {
			closest=gsub("^.*:","",names(isoforms[[chr_idx[i]]])[nearest(region,isoforms[[chr_idx[i]]])])
			model$n_content[i]=0
			model$content[i]=paste("[Nearest:",closest,"]",sep="")
		}
	}
	return(model)
}

# -------------------------------------------------------------------
getEventFrequency=function(model,cohort,lesions,thr0,thr1) {
	N=ncol(lesions$a0)
	for(i in 1:nrow(model)) {
		bidx=which(
			cohort$breaks[,1]==model$chromosome[i] &
			cohort$breaks[,2]>=ifelse(is.na(model$peak_start[i]),model$region_start[i],model$peak_start[i]) &
			cohort$breaks[,3]<=ifelse(is.na(model$peak_end[i]),model$region_end[i],model$peak_end[i])
		)
		status0=NULL
		status1=NULL
		if(model$model[i]=="Amp") {
			status0=apply(lesions$a0[bidx,,drop=F],2,max)>=thr0
			status1=apply(lesions$a1[bidx,,drop=F],2,max)>=thr1
		} else {
			status0=apply(lesions$d0[bidx,,drop=F],2,max)>=thr0
			status1=apply(lesions$d1[bidx,,drop=F],2,max)>=thr1
		}
		model$freqX0[i]=sprintf("%.1f",(sum(status0,na.rm=T)/N)*100)
		model$freqX1[i]=sprintf("%.1f",(sum(status1,na.rm=T)/N)*100)
		model$SamplesWithSpanningAlteration[i]=paste(unique(names(which(status0)),names(which(status1))),collapse=",")
	}
	return(model)
}

# -------------------------------------------------------------------
processStageResults=function(p,results,model,base,errorFactor,flag) {
	for(j in 1:nrow(results)) {
		model$model=c(model$model,flag)
		model$chromosome=c(model$chromosome,results[j,1])
		#rev=isCentromeric(cen[results[j,1],2],results[j,1],results[j,2],results[j,3],cohort$breaks)
		#model$region_start=c(model$region_start,rev$REVSTART)
		#model$region_end=c(model$region_end,rev$REVEND)
		model$region_start=c(model$region_start,results[j,2])
		model$region_end=c(model$region_end,results[j,3])
		model$peak_start=c(model$peak_start,NA)
		model$peak_end=c(model$peak_end,NA)
		model$size=c(model$size,results[j,3]-results[j,2]+1)
		model$q_value=c(model$q_value,results[j,10])
		if(results[j,11]==1) {; # This region may contain peaks
			pks=getStageTwoPeaks(p,rownames(results)[j],base,errorFactor=errorFactor)
			if (length(pks)==0) next;
			for(k in 1:nrow(pks)) {
				model$model=c(model$model,flag)
				model$chromosome=c(model$chromosome,results[j,1])
				model$region_start=c(model$region_start,results[j,2])
				model$region_end=c(model$region_end,results[j,3])
				model$peak_start=c(model$peak_start,pks[k,2])
				model$peak_end=c(model$peak_end,pks[k,3])
				model$size=c(model$size,pks[k,3]-pks[k,2]+1)
				model$q_value=c(model$q_value,pks[k,10])
			}
		}
	}
	return(model)
}

# -------------------------------------------------------------------
# Map 'likely' alteration types from breakpoints to genes in
# individual tumors. Done on the BSgenome (TBD) version of genes unless
# a absolute file path is provided to the geneSet parameter. This
# file must have have columns: name, chr, start, end.
# (add microRNAs??)
# -------------------------------------------------------------------

writeGeneMap=function(p,geneSet=NULL,genes=NULL,removePolymorhic=FALSE) {
	# Load gene data.
	if(is.null(geneSet)) {
		genes=feat[feat$Type=="Gene" | feat$Type=="MicroRNA",]
	} else {
		genes=try(read.table(geneSet,header=T,as.is=T,colClasses=c("character",rep("numeric",3))),silent=TRUE)
		if(class(genes)=="try-error") {
			cat("   >> Provided gene set is not available or formatted properly...skipping gene map!\n")
			return(NULL)
		}
	}
	genes=genes[order(genes$Chr,genes$Start),]
	genes=genes[!(genes$Chr>getSpeciesAutosome(p)),]

	# Load project data
	load(paste(p$Dir,"/",p$Project,"/",p$Project,"-lesions.Rdata",sep=""))
	load(paste(p$Dir,"/",p$Project,"/",p$Project,".Rdata",sep=""))
	rownames(cohort$breaks)=paste("chr",cohort$breaks[,1],":",cohort$breaks[,2],"-",cohort$breaks[,3],sep="")
	p$A0=p$D0=p$D1=0.9
	p$A1=0.25

	# Find breakpoints overlapping isoforms
	breaks=do.call(RangesList,lapply(1:getSpeciesAutosome(p),function(ch) {
		idx=which(cohort$breaks[,1]==ch)
		IRanges(start=cohort$breaks[idx,2],end=cohort$breaks[idx,3],names=rownames(cohort$breaks)[idx])
	}))
	isoforms=do.call(RangesList,lapply(1:getSpeciesAutosome(p),function(ch) {
		IRanges(start=genes$Start[genes$Chr==ch],end=genes$End[genes$Chr==ch],names=genes$Name[genes$Chr==ch])
	}))
	cov=lapply(1:getSpeciesAutosome(p),function(ch) matchMatrix(findOverlaps(isoforms[[ch]],breaks[[ch]])))

	# Fill in for genes with missing coverage by nearest neighbor
	cov=lapply(1:getSpeciesAutosome(p),function(ch) {
		missing=which(1:length(isoforms[[ch]])%in%cov[[ch]][,2]==FALSE)
		cov[[ch]]=rbind(cov[[ch]],cbind(nearest(isoforms[[ch]][missing,],breaks[[ch]]),missing))
		cov[[ch]][order(cov[[ch]][,2]),]
	})

	# All genes indexed by overlapping or nearest breakpoints
	map=lapply(1:getSpeciesAutosome(p),function(ch) {
		cbind(names(breaks[[ch]])[cov[[ch]][,1]],names(isoforms[[ch]])[cov[[ch]][,2]])
	})
	map=as.matrix(do.call(rbind,map))
	map=split(
		match(map[,1],rownames(cohort$breaks)),
		match(map[,2],genes$Name)
	)

	# Assign gene status, extremal signal for multiple overlapping (intragenic) breakpoints
	cna=lapply(map,assignGeneStatus,p,lesions)
	cna=t(as.data.frame(cna))
	cna=cbind(cna,genes[,c("Chr","Start","End")])
	rownames(cna)=genes$Name
	colnames(cna)=c(colnames(lesions$a0),"Chr","Start","End")
	save(cna,file=paste(p$Dir,p$Project,"GeneMap.Rdata",sep="/"))
}

# -------------------------------------------------------------------

# -------------------------------------------------------------------
plotGenome=function(p,dd,save="pdf",live=FALSE) {
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

	if(!live) {
		if(save=="pdf") {
			pdf(file=paste(p$Dir,"/",p$Project,"/",p$Project,".pdf",sep=""),width=6,height=5)
		} else {
			png(file=paste(p$Dir,"/",p$Project,"/",p$Project,".png",sep=""),width=450,height=375)
		}
	}

	makeChromPolygons=function(pos,dd,ymin,ymax) {
		for(i in seq(1,22,2)) {
			x=range(pos[dd$anno$chrom==i])
			polygon(c(x[1],x[1],x[2],x[2]),c(ymin,ymax,ymax,ymin),col="gray89",border="gray89")
		}
	}

	layout(matrix(c(1,1,2,2,3,3,4,4),4,2,byrow=T),widths=rep(1,4),heights=c(0.33,0.165,0.165,0.33))
	par(xaxs="i",yaxs="i")
	ybase=10^-seq(0,20,by=5)
	yhigh=10^-round(seq(20,ceiling(-log10(as.numeric(sprintf("%.1g",min(dG$y))))),length=5))
	ylab=yhigh
	ylab[1]=NA
	par(mar=c(0,5,0.5,0.1))
	plot(1,1,type="n",xlim=range(dG$x),ylim=rev(range(yhigh)),log="y",axes=FALSE,ylab="Q-value (gain)",xlab="")
	makeChromPolygons(pos,dd,1e-20,min(dG$y))
	idx=which(dG$y<1e-20)
	lines(dG$x[idx],dG$y[idx],type="h",col="red3")
	axis(2,at=yhigh,labels=ylab,las=2,cex.axis=0.8)
	#axis.break(2,1e-40)
	par(mar=c(0,5,0,0.1))
	plot(1,1,type="n",xlim=range(dG$x),ylim=rev(range(ybase)),log="y",axes=FALSE,ylab="",xlab="")
	makeChromPolygons(pos,dd,1,1e-20)
	idx=which(dG$y>=0.9)
	lines(dG$x[-idx],dG$y[-idx],type="h",col="red3")
	axis(2,at=ybase,las=2,cex.axis=0.8)

	par(mar=c(0,5,0,0.1))
	plot(1,1,type="n",xlim=range(dL$x),ylim=-log10(range(ybase)),axes=FALSE,ylab="",xlab="")
	makeChromPolygons(pos,dd,0,20)
	idx=which(dL$y>=0.9)
	lines(dL$x[-idx],-log10(dL$y[-idx]),type="h",col="darkblue")
	ylab=ybase
	ylab[1]=NA
	axis(2,at=-log10(ybase),labels=ylab,las=2,cex.axis=0.8)
	abline(h=0)
	par(mar=c(0.5,5,0,0.1))
	yhigh=10^-round(seq(20,ceiling(-log10(as.numeric(sprintf("%.1g",min(dL$y))))),length=5))
	ylab=yhigh
	ylab[1]=NA
	plot(1,1,type="n",xlim=range(dL$x),ylim=-log10(range(yhigh)),axes=FALSE,ylab="Q-value (loss)",xlab="")
	makeChromPolygons(pos,dd,20,max(-log10(yhigh)))
	idx=which(dL$y<1e-20)
	lines(dL$x[idx],-log10(dL$y[idx]),type="h",col="darkblue")
	axis(2,at=-log10(yhigh),labels=ylab,las=2,cex.axis=0.8)
	#axis.break(2,40)
	if(!live) dev.off()
}

# -------------------------------------------------------------------
plotGenome2=function(p,dd,save="pdf",live=FALSE) {
	suppressMessages(require(IRanges))
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

	if(!live) {
		if(save=="pdf") {
			pdf(file=paste(p$Dir,"/",p$Project,"/",p$Project,".pdf",sep=""),width=6,height=5)
		} else {
			png(file=paste(p$Dir,"/",p$Project,"/",p$Project,".png",sep=""),width=450,height=375)
		}
	}
	layout(matrix(c(1,1,2,2),2,2,byrow=T),widths=rep(1,2),heights=c(0.5,0.5))
	par(mar=c(0.1,5,0.1,0.1))
	ylim=c(1,min(dG$y))
	plot(1,1,type="n",xlim=range(dG$x),ylim=ylim,log="y",axes=FALSE,ylab="Q-value (gain)",xlab="")
	for(i in seq(1,22,2)) {
		x=range(pos[dd$anno$chrom==i])
		polygon(c(x[1],x[1],x[2],x[2]),c(1,ylim[2],ylim[2],1),col="gray89",border="gray89")
	}
	polygon(dG,col="red3",border="red3")
	#ax=axTicks(2)
	#idx=c(max(which(ax<ylim[2])),which(ax==1))
	#axis(2,at=ax[round(seq(idx[1],idx[2],length=5))],las=2,cex.axis=0.8)
	ax=10^-round(seq(0,ceiling(-log10(as.numeric(sprintf("%.1g",min(dG$y))))),length=5))
	axis(2,at=ax,las=2,cex.axis=0.8)
	lines(range(dG$x),c(1,1),lwd=2)

	ylim=c(min(dL$y),1)
	plot(1,1,type="n",xlim=range(dL$x),ylim=ylim,log="y",axes=FALSE,ylab="Q-value (loss)",xlab="")
	for(i in seq(1,22,2)) {
		x=range(pos[dd$anno$chrom==i])
		polygon(c(x[1],x[1],x[2],x[2]),c(ylim[1],1,1,ylim[1]),col="gray89",border="gray89")
	}
	polygon(dL,col="darkblue",border="darkblue")
	#ax=axTicks(2)
	#idx=c(max(which(ax<ylim[1])),which(ax==1))
	#axis(2,at=ax[round(seq(idx[1],idx[2],length=5))],las=2,cex.axis=0.8)
	ax=10^-round(seq(0,ceiling(-log10(as.numeric(sprintf("%.1g",min(dL$y))))),length=5))
	axis(2,at=ax,las=2,cex.axis=0.8)
	lines(range(dL$x),c(1,1),lwd=2)
	if(!live) dev.off()
}

#####################################################################
###
### Miscellaneous
###
#####################################################################

# -------------------------------------------------------------------
assignGeneStatus=function(bidx,p,lesions) {
	N=ncol(lesions$a0)
	if(length(bidx)>1)  {; # Gene is spanned by multiple breakpoints
		gidx=bidx[apply(lesions$a0[bidx,],2,which.max)]; # In these regions, the index of the maximum value of A0
		a0=sapply(1:N,function(i) lesions$a0[gidx[i],i]); # Use max index of A0 to get single A0
		a1=sapply(1:N,function(i) lesions$a1[gidx[i],i]); # Use max index of A0 to get single A1
		lidx=bidx[apply(lesions$d1[bidx,],2,which.max)]; # In these regions, the index of the maximum value of D0
		d0=sapply(1:N,function(i) lesions$d0[lidx[i],i]); # Use max index of D0 to get single D0
		d1=sapply(1:N,function(i) lesions$d1[lidx[i],i]); # Use max index of D0 to get single D1
		gain=getDiscretizedValue(a0,a1,p$A0,p$A1)
		loss=getDiscretizedValue(d0,d1,p$D0,p$D1)
		return(ifelse(gain>loss,gain,-loss))
	} else {; # Gene is spanned by a single breakpoint
		gain=getDiscretizedValue(lesions$a0[bidx,],lesions$a1[bidx,],p$A0,p$A1)
		loss=getDiscretizedValue(lesions$d0[bidx,],lesions$d1[bidx,],p$D0,p$D1)
		return(ifelse(gain>loss,gain,-loss))
	}
}

# -------------------------------------------------------------------
calc.pv=function(statistic,params,lower=FALSE,logT=FALSE) {
	return(pnorm(statistic,mean=params$mean,sd=params$sd,lower=lower,log=logT))
}

# -------------------------------------------------------------------
compute.rsegs=function(d) {
	d$rsegs=d$segs[,1:2]
	for(i in 3:ncol(d$segs)) {
		rsegs=rep(NA,nrow(d$segs))
		id=colnames(d$segs)[i]
		out=d$output[d$output$ID==id,]
		for(j in seq(nrow(out))) {
			st=out$start.idx[j]
			ed=out$end.idx[j]
			n=ed-st+1
			rsegs[st:ed]=rnorm(n,mean=out$seg.mean[j],sd=out$seg.sd[j]/sqrt(out$num.mark[j]))
		}
		d$rsegs=cbind(d$rsegs,rsegs)
	}
	colnames(d$rsegs)=colnames(d$segs)
	return(d)
}

# -------------------------------------------------------------------
density.Norm=function(x,bw="nrd0",adjust=1,n=N.density,from=min(x,na.rm=T),to=max(x,na.rm=T),...) {
	pi=density(x,bw=bw,adjust=adjust,from=from,to=to,n=n,na.rm=TRUE,...)
	pi$y=pi$y/sum(pi$y)
	return(pi)
}

# -------------------------------------------------------------------
detectCopyNumberPeaks=function(trace,delta) {
	lookformax=1
	maxtab=NULL
	mxpos=NA
	mnpos=NA
	mx=-Inf
	mn=Inf
	for(i in 1:length(trace)) {
		this=trace[i]
		if(this>mx) {
			mx=this; mxpos=i
		}
		if(this<mn) {
			mn=this; mnpos=i
		}
		if (lookformax) {
			if (this<(mx-delta)) {
				maxtab=c(maxtab,mxpos); mn=this; mnpos=i; lookformax=0
			}
		} else {
			if (this>(mn+delta)) {
				mx=this; mxpos=i; lookformax=1
			}
		}
	}
	return(maxtab)
}

# -------------------------------------------------------------------
dlrs=function(x) median(abs(diff(x)),na.rm=T)

# -------------------------------------------------------------------
fer=function(E,Eo=0.5,b=25) return(1/(1+exp(-b*(E-Eo))))

# -------------------------------------------------------------------
ferA1=function(x,Eo,B=1/log(2)) return(theta(x,Eo)*((2*fer(x,Eo,B))-1))

# -------------------------------------------------------------------
ferDeriv=function(E,Eo,b) return((abs(b)*exp(b*(E-Eo)))/((exp(b*(E-Eo))+1)^2))

# -------------------------------------------------------------------
#
# Mstats must be the matrix of statistics, it must be breaks in rows
# and samples in cols (ie: ncol(Mstats)==N.samps). We want to compute
# the p-value of the sum over samples Mstats for each breakpoint.
# This function is called once to get the parameters used by the
# calc.pv function.
#
# Nrs == NULL no re-sampling
#
# Set Nrs to a number greater than 0 to resample the Mstats matrix.
# -------------------------------------------------------------------

fit.pv.params=function(Mstats,Nrs=NULL) {
	N.samp=ncol(Mstats)
	if(is.null(Nrs) || Nrs<1) {
		mm=N.samp*mean(as.numeric(Mstats),na.rm=T)
		ss=sqrt(N.samp*var(as.numeric(Mstats),na.rm=T))
	} else {
		rii=double(Nrs*nrow(Mstats))
		for(i in seq(Nrs)) {
			print(i)
			rii[seq(nrow(Mstats))+(i-1)*nrow(Mstats)]=rowSums(
				matrix(sample(Mstats[!is.na(Mstats)],prod(dim(Mstats)),replace=T),ncol=ncol(Mstats))
			)
		}
		mm=mean(rii)
		ss=sd(rii)
	}
	return(list(N.samp=N.samp,mean=mm,sd=ss))
}

# -------------------------------------------------------------------
getBreakpointCoverage=function(pos,breaks) {
	chr=pos[1]
	st=pos[2]
	ed=pos[3]

	if(ed<cohort$breaks[min(which(breaks[,1]==chr)),2]) {
		return(min(which(breaks[,1]==chr))); # Is 5' of coverage (telomeric)
	}
	if(st>cohort$breaks[max(which(cohort$breaks[,1]==chr)),3]) {
		return(max(which(breaks[,1]==chr))); # Is 3' of coverage (telomeric)
	}
	ii=which(breaks[,1]==chr & ((breaks[,2]>=st & breaks[,2]<=ed) | (breaks[,3]>=st & breaks[,3]<=ed)))
	ii=c(ii,which(breaks[,1]==chr & breaks[,2]<=st & breaks[,3]>=ed))
  if(length(ii)==0) {; # Get nearest breakpoint
		five=max(which(breaks[,1]==chr & breaks[,3]<=st))
		three=min(which(breaks[,1]==chr & breaks[,2]>=ed))
		ii=as.numeric(ifelse((st-breaks[five,3])<(breaks[three,2]-ed),five,three))
	}
	return(ii)
}

# -------------------------------------------------------------------
get.diploid.peak=function(pd,Nf=4,LOG=TRUE) {
	ii=which(pd$x>=-Nf*pd$bw & pd$x<=Nf*pd$bw)
	pk=try(peaksign(pd$y[ii],11),silent=TRUE)
	if(class(pk)=="try-error") return(ii[which.max(pd$y[ii])])
	pc=sum(pk==1)
	if(pc==1) {
		if(LOG) cat(" Diploid peak detected!\n")
		return(list(peak=ii[which(pk==1)],status=1))
	} else if (pc>1) {
		if(LOG) cat(" Diploid peak ambiguity, review!\n")
		return(list(peak=ii[which(pk==1)[which.max(pd$y[ii[which(pk==1)]])]],status=0))
	} else {
		if(LOG) cat(" No peaks; investigate!\n")
		return(list(peak=ii[which.max(pd$y[ii])],status=-1))
	}
}

# -------------------------------------------------------------------
getDiscretizedValue=function(x0,x1,thresh1,thresh2) {
	discretized=rep(0,length(x0))
	discretized[which(x0>=thresh1 & x1<thresh2)]=1; # Samples considered to be only single-copy loss or gain
	discretized[which(x0>=thresh1 & x1>=thresh2)]=2; # Samples considered to be either amplified or homozygous
	return(discretized)
}

# -------------------------------------------------------------------
getEvaluatedBreakpoint=function(x,segs) {
	idx=x[6]:x[7]
	a0=mean(apply(segs$a0[idx,,drop=FALSE],2,getWindowWeightedMean),na.rm=T)
	a1=mean(apply(segs$a1[idx,,drop=FALSE],2,getWindowWeightedMean),na.rm=T)
	d0=mean(apply(segs$d0[idx,,drop=FALSE],2,getWindowWeightedMean),na.rm=T)
	d1=mean(apply(segs$d1[idx,,drop=FALSE],2,getWindowWeightedMean),na.rm=T)
	return(c(a0,a1,d0,d1))
}

# -------------------------------------------------------------------
getEvaluatedCleavedBreakpoint=function(x,cn) return(apply(cn[x[6]:x[7],,drop=FALSE],2,getWindowWeightedMean))

# -------------------------------------------------------------------
getExcludedSamples=function(pp,p,LOG=FALSE) {
	if(!p$Exclude) {
		# Not using QC parameters, but checking for un-parameterized samples...
		rem=which(is.na(pp$Ea0))
		if(length(rem)>0) {
			return(list(status=TRUE,samples=as.character(as.vector(pp$names[rem]))))
		} else {
			return(list(status=FALSE,samples=NA))
		}
	}

	# Calculate thresholds for both values of signal-to-noise, based on FWHM and also
	# based on segmentation jumps, values < 1.5*IQR for log10 of the values triggers
	# removal of the sample for putative quality reasons
	stats1=stats::fivenum(log10(pp$SN.fwhm),na.rm=T)
	sn1=stats1[2]-1.5*diff(stats1[c(2,4)])

	stats2=stats::fivenum(log10(pp$SEGN),na.rm=T)
	sn2=stats2[2]-1.5*diff(stats2[c(2,4)])

	### Find samples to exclude
	rem=which(is.na(pp$Ed0) | (log10(pp$SN.fwhm)<sn1) | (log10(pp$SEGN)<sn2))

	### Log the exclusions
	if(LOG & length(rem)>0) {
		cat("----------------------------------------------------------------\n")
		cat("--- Found ",length(rem)," exclusions per QC criteria...\n",sep="")
		cat("----------------------------------------------------------------\n")

		ex=pp[rem,]
		reasons=rep(NA,length(rem))
		reasons[(is.na(ex$Ea0))]="Catastrophic failure, FWHM doesn't span zero"
		reasons[(log10(ex$SN.fwhm)<sn1)]="Low signal-to-noise or poorly positioned discriminators"
		reasons[(log10(ex$SEGN)<sn2)]="High segmentation noise (low adjacent signal jumps amid high DN)"
		for(i in 1:nrow(ex)) cat(as.character(ex$names[i]),"\t",reasons[i],"\n",sep="")
	}

	### Exclude
	if(length(rem)) {
		return(list(status=TRUE,samples=as.character(as.vector(pp$names[rem]))))
	} else {
		return(list(status=FALSE,samples=NA))
	}
}

# -------------------------------------------------------------------
get.fwhm=function(pd) {
	pi=peaksign(pd$y,13)
	px=which(pi!=0)
	ps=pi[px]
	pm=px[which(ps>0)]
	dpii=which.min(abs(which(pi>0)-length(pi)/2))
	dpi=which(pi>0)[dpii]
	N=length(pd$x)
	fwp=(dpi:N)[which(pd$y[dpi:N]-pd$y[dpi]/2<0)[1]]
	fwn=(dpi:1)[which(pd$y[dpi:1]-pd$y[dpi]/2<0)[1]]
	fwp.x=pd$x[fwp]
	fwn.x=pd$x[fwn]
	fwhm=abs(diff(pd$x[c(fwn,fwp)]))
	if(abs(fwn.x)<0.025) fwn.x=-0.025
	if(fwp.x<0.025) fwp.x=0.025
	return(list(fwhm=fwhm,pi=fwp,ni=fwn,dpi=dpi,px=fwp.x,nx=fwn.x))
}

# -------------------------------------------------------------------
getProcessedBreakpoints=function(p,samples,segments,anno,cleave=FALSE) {
	positions=segments[segments[,1]%in%samples,c(2:4,9:10)]
	breaks=NULL
	for(chr in 1:getSpeciesAutosome(p)) {
		coord=as.vector(anno[anno[,1]==chr,2])
		incr=ifelse(chr>1,min(which(anno[,1]==chr))-1,0)
		windows=getUniqueBreakpoints(positions[positions[,1]==chr,2:5],coord,incr)
		windows=cbind(windows,apply(windows,1,function(x) x[4]-x[3]+1))
		idx=which(windows[,5]==1)
		if(length(idx)>0) {
			consecutive=which(abs(diff(c(FALSE,diff(idx)==1)))==1)
			if(length(consecutive)%%2) consecutive=c(consecutive,length(idx))
			if(length(consecutive)>0) {
				cons=matrix(idx[consecutive],ncol=2,nrow=length(consecutive)/2,byrow=T)
				removals=NULL
				for(i in 1:nrow(cons)) {
					windows[cons[i,1],2]=windows[cons[i,2],2]
					windows[cons[i,1],4]=windows[cons[i,2],4]
					windows[cons[i,1],5]=length(cons[i,1]:cons[i,2])
					rx=cons[i,1]:cons[i,2]
					removals=c(removals,rx[2:length(rx)])
				}
				windows=windows[-removals,]
			}
			isolated=which(windows[,5]==1)
			if(length(isolated)>0) {
				nidx=isolated+1
				windows[nidx,1]=windows[isolated,1]
				windows[nidx,5]=windows[nidx,5]+1
				windows[nidx,3]=windows[nidx,3]-1
				windows=windows[-isolated,]
			}
		}
		breaks=rbind(
			breaks,
			cbind(
				chr,windows[,1:2,drop=F],
				windows[,2]-windows[,1]+1,
				windows[,5],
				windows[,3:4,drop=F]
			)
		)
		rownames(breaks)=NULL
	}
	colnames(breaks)=c("Chr","Start","End","Size","Length","StartIdx","EndIdx")
	return(breaks)
}

# -------------------------------------------------------------------
getOutlierCorrectedSampleList=function(samples,segments) {
	ssc=sapply(samples,function(x,segments) sum(segments[,1]==x),segments)
	box=boxplot(ssc,plot=FALSE)
	ridx=NULL
	if(length(box$out)>0) {
		for(k in 1:length(box$out)) ridx=c(ridx,which(ssc==box$out[k]))
	}
	ridx=unique(ridx)
	if(length(ridx)>0) {
		return(samples[-ridx])
	} else {
		return(samples)
	}
}

# -------------------------------------------------------------------
getSignificantRegions=function(breaks,q,score) {
	sig=list()
	sig$Chr=breaks[,1]
	sig$Start=breaks[,2]
	sig$End=breaks[,3]
	sig$Size=breaks[,4]
	sig$Length=breaks[,5]
	sig$StartIdx=breaks[,6]
	sig$EndIdx=breaks[,7]
	sig$X0=score[,1]
	sig$X1=score[,2]
	sig$Delta=score[,3]
	sig$Q=q
	return(as.data.frame(sig))
}
getSignificantRegionsGain=function(cohort,p) {
	sidx=which(cohort$q[,1]<=p$Qval)
	if(length(sidx)==0) return(NULL)
	return(getSignificantRegions(cohort$breaks[sidx,],cohort$q[sidx,1],cohort$scored[sidx,c(3:4,5)]))
}
getSignificantRegionsLoss=function(cohort,p) {
	sidx=which(cohort$q[,2]<=p$Qval)
	if(length(sidx)==0) return(NULL)
	return(getSignificantRegions(cohort$breaks[sidx,],cohort$q[sidx,2],cohort$scored[sidx,c(1:2,6)]))
}

# -------------------------------------------------------------------
redactCNV=function(cohort,p) {
	suppressMessages(require(IRanges))

	# Find rgions of the breakpoint profile that overlap regions
	# of known polymorphism. This requires that a region's coverage
	# by one or more CNVs is >50% before its considered polymorhic
	polymorphic=NULL
	for(ch in unique(cohort$breaks[,1])) {
		idx=which(cohort$breaks[,1]==ch)
		brks=IRanges(start=cohort$breaks[idx,2],end=cohort$breaks[idx,3])
		poly=IRanges(start=cnv$Start[cnv$Chr==ch],end=cnv$End[cnv$Chr==ch])
		cov1=lapply(1:length(brks),function(b) coverage(poly,shift=-start(brks[b,])+1,width=width(brks[b,])))
		cov2=unlist(lapply(cov1,function(cc) sum(width(slice(cc,lower=1,includeLower=TRUE)))/length(cc)))
		polymorphic=c(polymorphic,idx[cov2>0.5])
	}

	# For each polymorphic region, find the nearest somatic region
	# from the remaining breakpoints
	som=(1:nrow(cohort$breaks))[-polymorphic]
	near=som[nearest(
		IRanges(start=polymorphic,end=polymorphic),
		IRanges(start=som,end=som)
	)]

	# Reset the q-values for the polymorphic regions to equal the
	# significance of the nearest somatic region.
	cohort$q[polymorphic,1]=cohort$q[near,1]
	cohort$q[polymorphic,2]=cohort$q[near,2]
	cohort$q[polymorphic,3]=cohort$q[near,3]
	cohort$q[polymorphic,4]=cohort$q[near,4]
	return(cohort)
}

# -------------------------------------------------------------------
getSpeciesAutosome=function(p) {
	if(length(grep("hg",p$Build))>0) return(22)
	if(length(grep("mm",p$Build))>0) return(19)
}

# -------------------------------------------------------------------
getStageOneRegions=function(sig,p) {
	chr=unique(sig$Chr)
	results=NULL
	for(i in chr) {
		this=sig[sig$Chr==i,,drop=F]
		if(nrow(this)==1) {
			results=rbind(results,c(
				as.numeric(this$Chr),
				as.numeric(this$Start),
				as.numeric(this$End),
				as.numeric(this$Size),
				as.numeric(this$Length),
				as.numeric(this$StartIdx),
				as.numeric(this$EndIdx),
				as.numeric(this$X0),
				as.numeric(this$X1),
				as.numeric(this$Q),
				0
			))
			next
		}
		wide.midx=(this$EndIdx[1:nrow(this)-1]+1)==this$StartIdx[2:nrow(this)]
		wide.ends=c(which(wide.midx==FALSE),nrow(this))
		wide.starts=c(1,wide.ends[-length(wide.ends)]+1)
		wide.size=this$End[wide.ends]-this$Start[wide.starts]+1
		wide.snps=this$EndIdx[wide.ends]-this$StartIdx[wide.starts]+1
		for(j in 1:length(wide.starts)) {
			jj=wide.starts[j]:wide.ends[j]
			doPeaks=ifelse(sum(this$Q[jj]<=(as.numeric(p$Qval)/10))>0,1,0)
			results=rbind(results,c(
				this[wide.starts[j],1],
				this[wide.starts[j],2],
				this[wide.ends[j],3],
				wide.size[j],
				wide.snps[j],
				this[wide.starts[j],6],
				this[wide.ends[j],7],
				mean(this$X0[jj]),
				mean(this$X1[jj]),
				min(this$Q[jj]),
				doPeaks
			))
		}
	}
	colnames(results)=c("Chr","Start","End","Size","Length","StartIdx","EndIdx","<X0>","<X1>","Q","doPeaks")
	rownames(results)=paste("chr",results[,1],":",results[,2],"-",results[,3],sep="")
	return(results)
}

# -------------------------------------------------------------------
getStageTwoPeaks=function(p,locus,sig,errorFactor=1,shoulder=2) {
	options(warn=-1000)
	hasPeaks=FALSE
	loc=splitLocus(locus)
	id=locus

	ii=which(sig$Chr==loc[1] & sig$Start>=loc[2] & sig$End<=loc[3])
	pad=-log10(max(sig$Q[which(sig$Q<=p$Qval/10)]))
	sig$Q=-log10(sig$Q)
	#sig$Delta=-log10(sig$Delta)

	trace=sig$Q[ii]
	trace=c(pad,trace,pad)
	trace[trace<pad]=pad
	delta=c(0,sig$Delta[ii],0)/errorFactor
	pks=detectCopyNumberPeaks(trace,mean(delta[2:(length(delta)-1)])*shoulder)

	if(FALSE) {
		require(plotrix)
		plot(trace,ylim=c(median(trace-delta)/2,max(trace+delta)))
		plotCI(trace,uiw=delta,add=T)
		points(pks,trace[pks],pch=20,col=2)
	}

	if(length(pks)==0) {
		return(NULL)
	}
	pks=pks-1
	peaks=NULL
	highlight=NULL

	for(pp in pks) {
		if((sig$Q[ii[pp]]-((sig$Delta[ii[pp]])/errorFactor))<pad) {
			next
		}
		skip=pks%in%highlight
		if(sum(skip)>0) {
			if(pp==pks[skip]) {
				next
			}
		}
		peak=sig[ii[pp],]
		peak$Delta=peak$Delta/errorFactor
		merged=NULL
		if(pp<length(ii)) {
			prev=peak$EndIdx
			for(n in (pp+1):length(ii)) {
				val=sig$Q[ii[n]]
				if(((sig$StartIdx[ii[n]]-1)==prev) && (val<=(peak$Q+peak$Delta) && val>=(peak$Q-peak$Delta))) {
					peak$EndIdx=sig$EndIdx[ii[n]]
					peak$End=sig$End[ii[n]]
					prev=sig$EndIdx[ii[n]]
					merged=c(merged,ii[n])
					highlight=c(highlight,n)
				}
			}
		}
		if(pp>1) {
			prev=peak$StartIdx
			for(n in (pp-1):1) {
				val=-log10(sig$Q[ii[n]])
				if(((sig$EndIdx[ii[n]]+1)==prev) && (val<=(peak$Q+peak$Delta) && val>=(peak$Q-peak$Delta))) {
					peak$StartIdx=sig$StartIdx[ii[n]]
					peak$Start=sig$Start[ii[n]]
					prev=sig$StartIdx[ii[n]]
					merged=c(merged,ii[n])
					highlight=c(highlight,n)
				}
			}
		}
		merged=c(merged,ii[pp])
		peak$Size=peak$End-peak$Start+1
		peak$Length=peak$EndIdx-peak$StartIdx+1
		peak$X0=mean(sig$X0[merged])
		peak$X1=mean(sig$X1[merged])
		peak$Q=min(10^-sig$Q[merged])
		peak$L2=NULL
		peak$Delta=NULL
		rownames(peak)=paste(id,paste("chr",peak$Chr,":",peak$Start,"-",peak$End,sep=""),sep="|")
		peaks=rbind(peaks,peak)
		hasPeaks=TRUE
	}
	options(warn=0)
	return(peaks)
}

# -------------------------------------------------------------------
getStartsAndEnds=function(anno) {
	chr=unique(anno[,1])
	st=sapply(chr,function(ch) min(anno[anno[,1]==ch,2]))
	ed=sapply(chr,function(ch) max(anno[anno[,1]==ch,2]))
	return(list(Chr=chr,Starts=st,Ends=ed))
}
getStartsAndEndsP=function(anno) {
  chr=unique(anno[,1])
  idx=sapply(chr,function(ch) range(which(anno[,1]==ch)))
  return(list(Chr=chr,Starts=idx[1,],Ends=idx[2,]))
}

# -------------------------------------------------------------------
getUniqueBreakpoints=function(segments,coord,incr) {
	starts=unique(segments[,1])
	starts=starts[order(starts)]
	nas=which(is.na(match(starts,coord)))
	if(length(nas)>0) starts=starts[-nas]
	if(length(starts)>1) {
		ends=coord[c(match(starts[2:length(starts)],coord)-1,which(coord==max(coord)))]
	} else {
		ends=max(coord)
	}
	return(cbind(starts,ends,match(starts,coord)+incr,match(ends,coord)+incr))
}

# -------------------------------------------------------------------
getWindowWeightedMean=function(samp) {
	nas=which(is.na(samp))
	if(length(nas)>0) samp=samp[-nas]
	n=length(samp)
	g=unique(samp)
	if(length(g)==1) {
		w=rep(1,n)
	} else {
		w=NULL
		for(i in 1:length(g)) {
			gi=length(which(samp==g[i]))
			w=c(w,rep(gi/n,gi))
		}
	}
	return(weighted.mean(samp,w))
}

# -------------------------------------------------------------------
isCentromeric=function(centro,chr,st,ed,breaks) {
	if(st<centro & centro<ed) {; # Region spans the centromere
		cidx=which(breaks[,1]==chr & breaks[,2]<=centro & breaks[,3]>=centro)
		if(length(cidx)==0) {
			if((ed-centro+1)>(centro-st+1)) {
				cidx=min(which(breaks[,1]==chr & breaks[,2]>=centro))
			} else {
				cidx=max(which(breaks[,1]==chr & breaks[,3]<=centro))
			}
		}
		psize=centro-st+1
		qsize=ed-centro+1
		if(qsize>psize) {
			return(list(REVSTART=breaks[cidx+1,2],REVEND=ed))
		} else {
			return(list(REVSTART=st,REVEND=breaks[cidx,3]))
		}
	}
	return(list(REVSTART=st,REVEND=ed))
}

# -------------------------------------------------------------------
# SNP6 data comes with many short segments that are artifact and need
# to be removed. This will merge segments with probes <=limit to the
# closest adajacent segment.
# -------------------------------------------------------------------
mergeSegments=function(d,limit) {
	ds=d
	rem=which(d$output$num.mark<limit); # Candidate segments to merge
	chr=unique(d$output$chrom[rem]); # Chromsomes on which mergeable segments appear

	for(i in chr) {; # Merge is done on a per-chromosome basis
		ridx=rem[which(d$output$chrom[rem]==i)]
		st=d$output$loc.start[ridx]
		ed=d$output$loc.end[ridx]
		prime5seg=d$output$loc.end[(ridx-1)]; # The end breakpoint of the 5' segment
		if(ridx[1]==1) prime5seg=c(-99999999,prime5seg); # Ensure first segment maps to 3'
		prime3seg=d$output$loc.start[(ridx+1)]; # The start breakpoint of the 3' segment
		near=(st-prime5seg)<(prime3seg-ed)

		# If last sement in sample is <limit, merge to 5' adjacent segment
		if(max(ridx)==nrow(d$output)) near[length(near)]=TRUE

		# For those closer to the 5' adjacent segment, reset 3' breakpoint and probe counts
		nidx=ridx[near==T]-1
		ds$output$loc.end[nidx]=d$output$loc.end[ridx[near==T]]
		ds$output$num.mark[nidx]=ds$output$num.mark[nidx]+d$output$num.mark[ridx[near==T]]
		probes=lapply(ridx[near==T],function(rr) d$output$start.idx[rr]:d$output$end.idx[rr])
		ds$segs[unlist(probes),3]=rep(ds$output$seg.mean[nidx],lapply(probes,length))

		# For those closer to the 3' adjacent segment, reset 5' breakpoint and probe counts
		nidx=ridx[near==F]+1
		ds$output$loc.start[nidx]=d$output$loc.start[ridx[near==F]]
		ds$output$num.mark[nidx]=ds$output$num.mark[nidx]+d$output$num.mark[ridx[near==F]]
		probes=lapply(ridx[near==F],function(rr) d$output$start.idx[rr]:d$output$end.idx[rr])
		ds$segs[unlist(probes),3]=rep(ds$output$seg.mean[nidx],lapply(probes,length))
	}
	ds$output=ds$output[-rem,]
	return(ds)
}

# -------------------------------------------------------------------
P=function(N,x,bw=DN) return(density.orm(x,bw=bw,adj=(1/N)))

# -------------------------------------------------------------------
peaks=function(series,span=3,do.pad=TRUE) {
	if((span<-as.integer(span))%%2!= 1) stop("'span' must be odd")
	s1=1:1+(s=span %/% 2)
	if(span==1) return(rep.int(TRUE,length(series)))
	z=embed(series,span)
	v=apply(z[,s1]>z[,-s1,drop=FALSE],1,all)
	if(do.pad) {
		pad=rep.int(FALSE,s)
		c(pad,v,pad)
	} else v
}

# -------------------------------------------------------------------
peaksign=function(series,span=3,do.pad=TRUE) {
	# Purpose: return (-1 / 0 / 1) if series[i] is ( trough / "normal" / peak )
	# Author: Martin Maechler, Date: 25 Nov 2005
	if((span<-as.integer(span))%%2!=1 || span==1) stop("'span' must be odd and >= 3")
	s1=1+(s<-span%/%2)
	z=embed(series,span)
	d=z[,s1]-z[,-s1]
	ans=rep.int(0,nrow(d))
	ans[apply(d>0,1,all)]=as.integer(1)
	ans[apply(d<0,1,all)]=as.integer(-1)
	if(do.pad) {
		pad=rep.int(0,s)
		c(pad,ans,pad)
	} else ans
}

# -------------------------------------------------------------------
extendStartAndEnd=function(d,bounds) {
	sidx=sapply(bounds$Chr,function(x) which(d$output$chrom==x)[1]); # Start segment on each chromosome
	eidx=c(sidx[2:length(sidx)]-1,nrow(d$output)); # Terminal segment on each chromosome
	d$output$loc.start[sidx]=bounds$Starts; # Ensure 5' breakpoint starts on the first probe of chromosome
	d$output$loc.end[eidx]=bounds$Ends; # Ensure 3' breakpoint ends on the last probe of chromosome
	return(d)
}

# -------------------------------------------------------------------
processCBS=function(d,p,pbound) {
	d$segs=as.matrix(d$data)
	if(is.null(d$map) | ncol(d$output)!=9) {; # Testing of our pipline output
		seg=rep(NA,nrow(d$data))
		for(i in seq(nrow(d$output))) {
			rr=(
				d$data$chrom==d$output$chrom[i] &
				d$data$maploc>=d$output$loc.start[i] &
				d$data$maploc<=d$output$loc.end[i]
			)
			seg[rr]=d$output$seg.mean[i]
			d$output$start.idx[i]=range(which(rr))[1]
			d$output$end.idx[i]=range(which(rr))[2]
			d$output$seg.sd[i]=sd(d$data[rr,3],na.rm=T)
		}
		#d$output$start.idx=match(
		#	paste(d$output$chrom,d$output$loc.start,sep="-"),
		#	paste(d$data[,1],d$data[,2],sep="-")
		#)
		#d$output$end.idx=c(d$output$start.idx[2:nrow(d$output)]-1,nrow(d$data))
		d$segs[,3]=seg
	} else {
		d$segs[,3]=d$map$seg
	}
	d=verifyProbeIndices(d,pbound)
	if(p$SNP6) d=mergeSegments(d,8)
	d$dn=median(abs(diff(d$data[,3])),na.rm=TRUE); # Derivative noise
	ds=compute.rsegs(d); # Gaussian-noised segmentation
	q=quantile(ds$rsegs[,3],probs=c(0.2,0.8),na.rm=TRUE); # Mass of the data
	m=median(ds$rsegs[(ds$rsegs[,3]>=q[1] & ds$rsegs[,3]<=q[2]),3],na.rm=TRUE); # Trimmed median
	d$mass_median=m
	d$segs[,3]=d$segs[,3]-m; # Center by the mass of the data
	d$data[,3]=d$data[,3]-m; # Center by the mass of the data
	P.norm=density(d$segs[,3],bw=d$dn/8,from=-d$dn/2,to=d$dn/2,na.rm=T); # Segmentation density
	kk=get.diploid.peak(P.norm); # Detect position of the diploid peak
	d$diploid_peak_status=kk$status; # 1=detected, 0=ambiguous, -1=none found
	d$offset=-P.norm$x[kk$peak]; # Offset from diploid of the peak of segmentation
	plot(P.norm,xlim=c(-2,2),col="royalblue",ylim=c(0,max(P.norm$y)),main=colnames(d$segs)[3])
	points(P.norm$x[kk$peak],P.norm$y[kk$peak],cex=2,col="red3")
	d$segs[,3]=d$segs[,3]+d$offset; # Shift by the offset from diploid
	d$data[,3]=d$data[,3]+d$offset; # Shift by the offset from diploid
	d$output$seg.mean=d$output$seg.mean-m+d$offset; # Shift by the offset from diploid
	lines(density(d$segs[,3],bw=d$dn/8,from=-2,to=2,na.rm=TRUE))
	abline(v=c(-d$dn,d$dn),col="forestgreen")
	text(-2,max(P.norm$y),paste("DN=",sprintf("%.3f",d$dn),sep=""),pos=4)
	abline(v=0,col="darkgray",lty=2)
	colnames(d$segs)=colnames(d$data)
	d$data=d$data[!(d$data[,1]>getSpeciesAutosome(p)),]
	d$segs=d$segs[!(d$segs[,1]>getSpeciesAutosome(p)),]
	d$output=d$output[!(d$output$chrom>getSpeciesAutosome(p)),]
	class(d$data)=c("CNA","data.frame")
	class(d)="DNAcopy"
	return(d)
}

# -------------------------------------------------------------------
splitLocus=function(locus) {
	tmp=unlist(strsplit(locus,":"))
	chr=as.numeric(gsub("chr","",tmp[1]))
	st=as.numeric(unlist(strsplit(tmp[2],"-"))[1])
	ed=as.numeric(unlist(strsplit(tmp[2],"-"))[2])
	return(c(chr,st,ed))
}

# -------------------------------------------------------------------
theta=function(x,Eo) return(ifelse(x<Eo,0,1))

# -------------------------------------------------------------------
timeStamp=function() format(Sys.time(),"%Y%m%d")

# -------------------------------------------------------------------
verifyProbeIndices=function(d,bounds) {
	sidx=sapply(bounds$Chr,function(x) which(d$output$chrom==x)[1])
	eidx=c(sidx[2:length(sidx)]-1,nrow(d$output))
	d$output$start.idx[sidx]=bounds$Starts
	d$output$end.idx[eidx]=bounds$Ends
	return(d)
}

# -------------------------------------------------------------------
xcdf=function(x,v) return(sum(v<x)/length(v))

# -------------------------------------------------------------------
overlap2=function(w1,w2,v1,v2) {
	#print(c(w1,w2,v1,v2))
	if ((w1>w2) || (v1>v2)) {
		stop("Your intervals are improperly defined.")
	}
	if ((w2<v1) || (v2<w1)) {
		return(0)
	} else if ((w2<=v2) && (w1>=v1)) {
		return(w2-w1)
	} else if ((v2<=w2) && (v1>=w1)) {
		return(v2-v1)
	} else if ((v1<=w1) && (w1<=v2)) {
		return(v2-w1)
	} else if ((w2<=v2) && (w1>=v1)) {
		return(w2-w1)
	} else if ((w1<=v1) && (v1<=w2)) {
		return(w2-v1)
	}
}

