analyse=function(){

resultfiles=list.files(path="Output_Data/",pattern="*.out")
length(resultfiles)

PDs=c()
stds=c()
means=c()

pdf("Results.pdf")
resno=1
op=par(mfrow=c(2,2),mai=c(0.7,0.65,0.3,0.15))
for(result in resultfiles){
	# Read in Starting population mean and standard deviation from filename
	start=strsplit(result,c(".out"))[[1]]
	start=strsplit(start,"_")
	smean=start[[1]][3]
	sstd=start[[1]][4]
	smean=as.numeric(substring(smean,2))
	sstd=as.numeric(substring(sstd,2))
	# Read file contents
	dat=read.delim(paste("Output_Data",result,sep="/"),sep=" ",header=FALSE,stringsAsFactors=FALSE)
	colnames(dat)=c("time", "PD", "passno", "num", "dead", "Afrac", "Bfrac", "deadfrac","commfrac")
	# Get cell passage times
	L=length(dat$time)
	pass=dat$passno[2:L]-dat$passno[1:L-1]
	pass=c(0,pass)
	dat$pass=pass
	pt=dat[dat$pass==1,]
	# Get the time at which the last uncommitted cell becomes committed
	lst=dat[dat$commfrac==0,]
	lst=lst[1,]

	PDs=c(PDs,max(dat$PD))
	means=c(means,smean)
	stds=c(stds,sstd)

	# Make some plots

		plot(dat$PD,dat$deadfrac,type="l",xlab="Population doublings",ylab="Fraction of senescent cells",lwd=2)
		abline(v=lst$PD,col="blue")
		abline(v=pt$PD,col="red",lty=2)
		
		#plot(dat$PD,dat$Afrac,type="l",xlab="Population doublings",ylab="Fraction of cells labelled A",ylim=c(0,1),lwd=2)
		#abline(v=pt$PD,col="red",lty=2)
		plot(dat$PD,dat$commfrac,type="l",xlab="Population doublings",ylab="Fraction of uncommitted cells",ylim=c(0,1),lwd=2)
		abline(v=lst$PD,col="blue")
		abline(v=pt$PD,col="red",lty=2)
		
		plot(dat$time,dat$deadfrac,type="l",xlab="Time (d)",ylab="Fraction of senescent cells",lwd=2)
		abline(v=lst$time,col="blue")
		abline(v=pt$time,col="red",lty=2)
		
		plot(dat$time,dat$PD,type="l",xlab="Time (d)",ylab="Population doublings",lwd=2,ylim=c(0,100))
		abline(v=lst$time,col="blue")
		abline(v=pt$time,col="red",lty=2)
	
	title(sprintf("%04d",resno),outer=TRUE,line=-1)
	resno=resno+1
}
par(op)

op=par(mai=c(1,1,1,1))
h<-hist(PDs,100,main="Observed population lifespans",xlab="PD")
xfit<-seq(min(PDs),max(PDs),length=400) 
yfit<-dnorm(xfit,mean=mean(means,na.rm=TRUE),sd=sd(means,na.rm=TRUE)) 
#yfit<-dnorm(xfit,mean=mean(means,na.rm=TRUE),sd=mean(stds,na.rm=TRUE)) 
yfit <- yfit*diff(h$mids[1:2])*length(PDs) 
#lines(xfit, yfit, col="blue", lwd=2)
par(op)
mean(PDs)
mean(means)
mean(stds)
sd(PDs)

#plot(means,PDs)
#cor(means,PDs)

dev.off()

}
