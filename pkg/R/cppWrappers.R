library(RcppGSL)

simsen <- function(Capacity=1000,Prat=0.25,Divrate=0.5,Deltarep=1.0,Pc=0.275,Meandiv=43,Stddiv=0,Inoc=5,Froot="simsen"){
	.Call( "simsen", Capacity, Prat, Divrate, Deltarep, Pc, Meandiv, Stddiv, Inoc, Froot, PACKAGE = "senesceR" )
}

iPop=data.frame(labels=c(1,2,1,2),divPots=c(10,20,10,20),commStates=c(1,1,1,1))

simsendist <- function(Capacity=1000,Prat=0.25,Divrate=0.5,Deltarep=1.0,Pc=0.275,Meandiv=43,Stddiv=0,InocPop=iPop,Froot="simsendist"){
	.Call( "simsendist", Capacity, Prat, Divrate, Deltarep, Pc, Meandiv, Stddiv, InocPop, Froot, PACKAGE = "senesceR" )
}
