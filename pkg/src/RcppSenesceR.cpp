#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "senesce.h"
#include <Rcpp.h>

using namespace Rcpp;

RcppExport SEXP simsen(SEXP Capacity,SEXP Prat,SEXP Divrate,SEXP Deltarep,SEXP Pc,SEXP Meandiv,SEXP Stddiv)
{
	int capacity=as<int>(Capacity);
	float prat=as<float>(Prat);
	float divrate=as<float>(Divrate);
	float deltarep=as<float>(Deltarep);
	float pc=as<float>(Pc);
	int meandiv=as<int>(Meandiv);
	int stddiv=as<int>(Stddiv);

	Senesce cells(capacity,prat,divrate,deltarep,pc,meandiv,stddiv);
	cells.inoculate(5);
	cells.simulate();
	cells.writeTable();
	
	return R_NilValue; 
}


