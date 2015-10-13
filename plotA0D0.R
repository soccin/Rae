load("raeImage.Rdata")

options(error=traceback)


breaks=data.frame(dd$cohort$breaks)

offset=cumsum(c(0,tapply(breaks$End,breaks$Chr,max)))
width=tapply(breaks$End,breaks$Chr,max)

breaks$Begin=breaks$Start+offset[breaks$Chr]
breaks$Stop=breaks$End+offset[breaks$Chr]

a0=dd$lesions$a0>.8
d0=dd$lesions$d0>.8


xpos=(breaks$Begin+breaks$Stop)/2
meanA0=apply(a0,1,mean,na.rm=T)
meanD0=apply(-d0,1,mean,na.rm=T)

plot(xpos,meanA0,type='n',ylim=c(-.6,.6),axes=F,xlab="",ylab="Gain/Loss %")
polygon(c(xpos,rev(xpos)),c(meanA0,rep(0,length(meanA0))),col="#880000",border=NA)
polygon(c(xpos,rev(xpos)),c(meanD0,rep(0,length(meanD0))),col="#000088",border=NA)

axis(2)
abline(v=offset,col=8,lty=2)
text(width+offset[-length(offset)]-width/2,-.55+.05*(seq(width)%%2),seq(width))
abline(h=0)
box()
