library(RcppGSL)

simsen <- function(Capacity=1000,Prat=0.25,Divrate=0.5,Deltarep=1.0,Pc=0.275,Meandiv=43,Stddiv=0){
	.Call( "simsen", Capacity, Prat, Divrate, Deltarep, Pc, Meandiv, Stddiv, PACKAGE = "senesceR" )
}

