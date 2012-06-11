readfiles=function(folder="."){
	resultfiles=list.files(path=folder,pattern="*.txt")
	bigdat=data.frame()
	for(result in resultfiles){
		# Read file contents
		dat=read.delim(paste(folder,result,sep="/"),sep="\t",header=TRUE,stringsAsFactors=FALSE)
		#colnames(dat)=c("time", "PD", "passno", "num", "dead", "Afrac", "Bfrac", "deadfrac","commfrac")
		dat$filename=substring(result,1,nchar(result)-4)
		dat$X=NULL
		bigdat=rbind(bigdat,dat)
	}
	return(bigdat)
	}

makeplot=function(dat){
	# Get cell passage times
	L=length(dat$Time)
	pass=dat$Passage[2:L]-dat$Passage[1:L-1]
	pass=c(0,pass)
	dat$pass=pass
	pt=dat[dat$pass==1,]
	# Get the time at which the last uncommitted cell becomes committed
	lst=dat[dat$FracCommit==0,]
	lst=lst[1,]

	op=par(mfrow=c(2,2),mai=c(0.7,0.65,0.3,0.15))

	plot(dat$PD,dat$FracDead,type="l",xlab="Population doublings",ylab="Fraction of senescent cells",lwd=2)
	abline(v=lst$PD,col="blue")
	abline(v=pt$PD,col="red",lty=2)
	
	#plot(dat$PD,dat$Afrac,type="l",xlab="Population doublings",ylab="Fraction of cells labelled A",ylim=c(0,1),lwd=2)
	#abline(v=pt$PD,col="red",lty=2)
	plot(dat$PD,dat$FracCommit,type="l",xlab="Population doublings",ylab="Fraction of uncommitted cells",ylim=c(0,1),lwd=2)
	abline(v=lst$PD,col="blue")
	abline(v=pt$PD,col="red",lty=2)
	
	plot(dat$Time,dat$FracDead,type="l",xlab="Time (d)",ylab="Fraction of senescent cells",lwd=2)
	abline(v=lst$Time,col="blue")
	abline(v=pt$Time,col="red",lty=2)
	
	plot(dat$Time,dat$PD,type="l",xlab="Time (d)",ylab="Population doublings",lwd=2,ylim=c(0,100))
	abline(v=lst$Time,col="blue")
	abline(v=pt$Time,col="red",lty=2)
	
	title(dat$filename[1],outer=TRUE,line=-1)
	par(op)
}

histPD=function(bigdat){
	PDs=tapply(bigdat$PD,bigdat$filename,max)
	op=par(mai=c(1,1,1,1))
	h<-hist(PDs,100,main="Observed population lifespans",xlab="PD")
	xfit<-seq(min(PDs),max(PDs),length=400) 
	par(op)
	return(h)
}



