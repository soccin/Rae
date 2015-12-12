##
##

library(DNAcopy)

file=commandArgs(trailingOnly=TRUE)
print(file)

pos=regexpr("_A1_",file)
nameTmpl=paste(substr(file,1,pos-1),"A1_pSubset",sep="_")


##################
##################

alpha=0.01
nperm=10000
undo.SD=1.0
#undo.SD=0.5

####
####

path=strsplit(file,"/")[[1]]
load(file)
dd=as.data.frame(d$data)
rm(d)
load("probeSubset.Rdata")
ii=sort(match(probeSubset,rownames(dd)))
dd=dd[ii,]

cna=CNA(dd[,3],dd$chrom,dd$maploc,"logratio",colnames(dd)[3])
rownames(cna)=rownames(dd)


smCNA=smooth.CNA(cna)

smOut=as.matrix(smCNA)
rownames(smOut)=rownames(dd)

d=segment(smCNA,verbose=T,alpha=alpha,nperm=nperm,
  undo.splits='sdundo',
  undo.SD=undo.SD)

d$params=list()
d$params$revision = '$Id: doCBSV3.R 72fc2c193ef1 2010/01/07 23:31:53 socci $'
d$params$segment.opts=list(alpha=alpha,nperm=nperm,undo.splits='sdundo',undo.SD=undo.SD)


##
## Map segments onto genome
##
seg = rep(NA,nrow(d$data))
ssd = rep(NA,nrow(d$data))
nn  = rep(NA,nrow(d$data))
d$output$seg.sd=rep(NA,nrow(d$output))
for(j in seq(nrow(d$output))) {
  rr = (
        d$data$chrom==d$output$chrom[j]
        & d$data$maploc>=d$output$loc.start[j]
        & d$data$maploc<=d$output$loc.end[j]
        )
  seg[rr] = d$output$seg.mean[j]
  ssd[rr] = sd(d$data[rr,3],na.rm=T)
  Srr=sum(rr)
  nn[rr] = Srr
  d$output$seg.sd[j]=sd(d$data[rr,3],na.rm=T)
}

d$output$start.idx=match(
  paste(d$output$chrom,d$output$loc.start),
  paste(d$data[,1],d$data[,2]))
d$output$end.idx=c(d$output$start.idx[2:nrow(d$output)]-1,nrow(d$data))

d$map=data.frame(seg=seg,seg.sd=ssd,num=nn)

d$QC=list()
d$QC$dn=median(abs(diff(d$data[,3])),na.rm=TRUE)

rfile=paste(nameTmpl,"__CBS_out.Rdata",sep="__")
save(d,file=cc(rfile),compress=T)

