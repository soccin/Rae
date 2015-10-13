outDir="CragoProgressionC"
lesionFile=file.path(outDir,paste(outDir,"summary-lesions.txt",sep="-"))
ubmFile=file.path(outDir,paste(outDir,"lesions.Rdata",sep="-"))

dd=read.delim(lesionFile)
load(ubmFile)

## Using SpanningAlteration as this correspoinds to thresholding A0,D0
## Not sure what QualifiedSpanning means

samples=colnames(lesions$a0)

peaks=!is.na(dd$peak_start)
dPeaks=dd[peaks,]
dPeaks$Event=apply(dPeaks[,c(1,2,5,6)],1,function(x){gsub(" ","",paste(unlist(x),collapse=":"))})
rownames(dPeaks)=dPeaks$Event

dRegions=dd[!peaks,]
dRegions$Event=apply(dRegions[,1:4],1,function(x){gsub(" ","",paste(unlist(x),collapse=":"))})
rownames(dRegions)=dRegions$Event

ans=dRegions[,c(1,2,3,4,7,8,9,10)]

write.xls(ans,cc(outDir,"SignificantRegions.txt"))

eventMatrix=matrix("",nrow=len(samples),ncol=nrow(dRegions))
rownames(eventMatrix)=samples
colnames(eventMatrix)=rownames(dRegions)

for(i in seq(nrow(dRegions))){
	event=rownames(dRegions)[i]
	eventMatrix[strsplit(dRegions$SamplesWithSpanningAlteration[i],",")[[1]],event]="X"
}

write.xls(eventMatrix,cc(outDir,"EventMatrix.txt"))

