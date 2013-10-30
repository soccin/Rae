write.xls <- function(dd,filename,row.names=T,col.names=NA) {
  if (!is.data.frame(dd)) {
    dd <- data.frame(dd,check.names=F)
  }
  if(!row.names) {
    col.names=T
  }
  write.table(dd,file=filename,sep="\t",quote=FALSE,
              col.names=col.names,row.names=row.names)
}

fer=function(E,Eo=0.5,b=25) return(1/(1+exp(-b*(E-Eo))))
ferA1=function(x,Eo,B=1/log(2)) return(theta(x,Eo)*((2*fer(x,Eo,B))-1))
p=list()
p$A0=p$D0=p$D1=0.9
p$A1=0.25

library(IRanges)
load("raeImage.Rdata")
feat=read.table("FEAT.file",header=T,skip=1,as.is=T,colClasses=c(rep("character",2),rep("numeric",3)))

gAnnote=feat[feat$Type %in% c("Gene","MicroRNA"),]
gAnnote=gAnnote[gAnnote$Chr<23,]
genes=RangedData(IRanges(start=gAnnote$Start,end=gAnnote$End),
                    name=gAnnote$Name,space=gAnnote$Chr)

probeChr=factor(dd$data[,"chrom"])
probes=RangedData(IRanges(start=dd$data[,"maploc"],width=1),space=probeChr,value=seq(nrow(dd$data)))

breaks=dd$lesions$breaks
blocks=RangedData(IRanges(start=breaks[,"Start"],end=breaks[,"End"]),space=breaks[,"Chr"],value=seq(nrow(breaks)))

# fo=findOverlaps(blocks,probes)
# blockSegs=list()
# for(chrom in seq(fo)) {
#     cat(chrom,"\t")
#     bSegs=list()
#     for(i in seq(ncol(dd$segs))) {
#         cat(i,",")
#         bSegs[[i]]=sapply(tapply(fo[[chrom]]@subjectHits,fo[[chrom]]@queryHits,c),function(x){median(dd$segs[x,i],na.rm=T)})
#     }
#     cat("\n")
#     blockSegs[[chrom]]=do.call(cbind,bSegs)
# }
# segs=do.call(rbind,blockSegs)
# colnames(segs)=colnames(dd$segs)
# rownames(segs)=seq(nrow(segs))
# write.xls(cbind(breaks[,1:4],segs),"segmentsByBlocks.txt")

dl=dd$lesions
s0=ifelse(dl$a0>dl$d0,ifelse(dl$a1>.25,2,ifelse(dl$a0>.9,1,0)),ifelse(dl$d1>.9,-2,ifelse(dl$d0>.9,-1,0)))

rownames(s0)=seq(nrow(s0))

#fo=findOverlaps(genes,blocks)
fp=findOverlaps(genes,probes)



raeCalls=list()
for(chrom in seq(fp)) {

    cat("Chrom=",chrom,"\t")
    raeChrom=matrix(NA,nrow=fp[[chrom]]@queryLength,ncol=ncol(dd$segs))
    rownames(raeChrom)=values(genes[chrom])[[1]][,1]
    pMap=values(probes[chrom])[[1]][,1]
    bMap=values(blocks[chrom])[[1]][,1]

    # First assign calls for genes that span multiple probes

    multiProbes=names(which(table(fp[[chrom]]@queryHits)>1))
    if(len(multiProbes)>0){
        idx=lapply(multiProbes,function(x){fp[[chrom]]@subjectHits[fp[[chrom]]@queryHits==as.numeric(x)]})
        si=list()
        for(i in seq(ncol(dd$segs))) {
            cat(i,",")
            si[[i]]=sapply(idx,function(x){
                    median(dd$segs[pMap[x],i],na.rm=T)
                    })
        }
        sx=do.call(cbind,si)

        a0=t(apply(sx,1,function(x){fer(x,pp$Ea0,pp$betaA0)}))
        a1=t(apply(sx,1,function(x){ferA1(x,pp$Ea1,pp$betaA1)}))
        d0=t(apply(sx,1,function(x){abs(fer(x,pp$Ed0,pp$betaD0))}))
        d1=t(apply(sx,1,function(x){abs(fer(x,pp$Ed1,pp$betaD1))}))
        raeTmp=ifelse(sx>0,ifelse(a1>.25,2,ifelse(a0>.9,1,0)),ifelse(d1>.9,-2,ifelse(d0>.9,-1,0)))

        raeChrom[as.numeric(multiProbes),]=raeTmp
    }

    # Then genes that span one or less probes take block value (s0) for call

    blkVal=which(is.na(raeChrom[,1]))

    bRange=blocks[chrom]@ranges[[1]]
    gRange=genes[chrom]@ranges[[1]]
    idx=nearest(gRange[blkVal,],bRange)

    bi=list()
    for(i in seq(ncol(s0))){
        bi[[i]]=unlist(sapply(idx,function(x){if(length(x)>0){s0[bMap[x],i]}else{NA}}))
    }
    bx=do.call(cbind,bi)

    raeChrom[blkVal,]=bx

    raeCalls[[chrom]]=raeChrom
    cat("\n")
}

temp=do.call(rbind,raeCalls)
raeCalls=temp
colnames(raeCalls)=colnames(dd$segs)

write.xls(raeCalls,"GeneMatrix.txt")
