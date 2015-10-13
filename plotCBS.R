options(error=traceback)
library(DNAcopy)
args=commandArgs(trailing=T)

data=args[1:length(args)]

for(fname in data){
    bname=gsub("_S01_.*","",basename(fname))
    print(bname)
    load(fname)
    png(file=paste("cbsPlot",bname,".png",sep="___"),width=900,height=600)
    plot(d)
    dev.off()
}

