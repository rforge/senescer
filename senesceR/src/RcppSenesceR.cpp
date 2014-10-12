#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "senesce.h"
#include <Rcpp.h>

using namespace Rcpp;

RcppExport SEXP simsen(SEXP Capacity,SEXP Prat,SEXP Divrate,SEXP Deltarep,SEXP Pc,SEXP Meandiv,SEXP Stddiv, SEXP Inoc, SEXP Froot)
{
	int capacity=as<int>(Capacity);
	float prat=as<float>(Prat);
	float divrate=as<float>(Divrate);
	float deltarep=as<float>(Deltarep);
	float pc=as<float>(Pc);
	int meandiv=as<int>(Meandiv);
	int stddiv=as<int>(Stddiv);
	int inoc=as<int>(Inoc);
	string fRoot=as<string>(Froot);

	Senesce cells(capacity,prat,divrate,deltarep,pc,meandiv,stddiv);
	cells.inoculate(inoc);
	cells.simulate();
	cells.writeTable(fRoot);
	
	return R_NilValue; 
}

RcppExport SEXP simsendist(SEXP Capacity,SEXP Prat,SEXP Divrate,SEXP Deltarep,SEXP Pc,SEXP Meandiv,SEXP Stddiv, SEXP InocPop, SEXP Froot)
{
	int capacity=as<int>(Capacity);
	float prat=as<float>(Prat);
	float divrate=as<float>(Divrate);
	float deltarep=as<float>(Deltarep);
	float pc=as<float>(Pc);
	int meandiv=as<int>(Meandiv);
	int stddiv=as<int>(Stddiv);
	string fRoot=as<string>(Froot);
	
	DataFrame inocPop = DataFrame(InocPop);

	IntegerVector labels=inocPop["labels"];
	IntegerVector divPots=inocPop["divPots"];
	IntegerVector commStates=inocPop["commStates"];
	int popsize=labels.size();

	Senesce cells(capacity,prat,divrate,deltarep,pc,meandiv,stddiv);
	cells.inoculate(0);
	for (int i=0;i<min(popsize,capacity);i++){
		cells.flask[i][0]=labels[i];
		cells.flask[i][1]=divPots[i];
		cells.flask[i][2]=commStates[i];
	}
	cells.num=popsize;
	cells.num0=popsize;
	
	cells.simulate();
	cells.writeTable(fRoot);
	
	return R_NilValue; 
}
