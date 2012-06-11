#include "senesce.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>

#define round(x) (x<0?ceil((x)-0.5):floor((x)+0.5))

Senesce::Senesce(int Capacity,float Prat,float Divrate,float Deltarep,float Pc,float Meandiv,float Stddiv){

	capacity=Capacity;
	prat=Prat;
	divrate=Divrate;
	deltarep=Deltarep;
	P=Pc;
	meandiv=Meandiv;
	stddiv=Stddiv;

	// Prepare GSL random number generator (Mersenne twister)
	r = gsl_rng_alloc(gsl_rng_mt19937);
	// Seed RNG with current system time (in seconds)
	gsl_rng_set(r,time(0));

	passage=int(round(float(capacity)*prat));

	// 0: Label, 1: Division potential, 2: Commitment state
	nfeat=3; 	      // Number of features a cell has
	nreps = 2000;	      // Maximum number of reports (ie max simulation time = deltarep*nreps)
	repcols = 9;	      // Number of descriptions in rep (report) object
	histbins = 201;       // Number of bins for remaining division potential +1

	// Allocate Global recording variables
	flask = new int*[capacity];
	for(int i=0; i<capacity; i++) flask[i]=new int[nfeat];
	transfer = new int*[passage];
	for(int i=0; i<passage; i++) transfer[i]=new int[nfeat];
	ind = new int[capacity];
	rep = new float*[nreps];
	for(int i=0; i<nreps; i++) rep[i]=new float[repcols];
	hist = new int*[nreps];
	for(int i=0; i<nreps; i++) hist[i]=new int[histbins];

	// Initialise recording variables
	for (int i=0; i<capacity; i++){
		for (int j=0; j<nfeat; j++) flask[i][j]=0;	
	}
	for (int i=0; i<passage; i++){
		for (int j=0; j<nfeat; j++) transfer[i][j]=0;	
	}
	for (int i=0; i<capacity; i++) ind[i]=0;
	for (int i=0; i<nreps; i++){
		for (int j=0; j<repcols; j++) rep[i][j]=0;	
	}
	for (int i=0; i<nreps; i++){
		for (int j=0; j<histbins; j++) hist[i][j]=0;	
	}
	std::cout << starting << std::endl;
}

Senesce::~Senesce(){
	for (int i = 0; i < capacity; ++i) delete[] flask[i];
	for (int i = 0; i < passage; ++i) delete[] transfer[i];
	delete[] ind;
	for (int i = 0; i < nreps; ++i) delete[] rep[i];
	for (int i = 0; i < nreps; ++i) delete[] hist[i];
}


void Senesce::report(){
	// Calculate proportion of typeA typeB (and dead) cells
	int typeA=0;
	int typeB=0;
	int dead=0;
	int uncomm=0;

	for(int i=0; i<num; i++){
		if(flask[i][0]==1) typeA+=1;
		if((flask[i][1]==0)&&(flask[i][2]==1)) dead+=1;
		if((flask[i][0]>0)&&(flask[i][2]==0)){
			uncomm+=1;
			hist[reprow][histbins-1]+=1;
		}else{
			hist[reprow][flask[i][1]]+=1;
		}
	}

	typeB=num-typeA;
	float Afrac=float(typeA)/float(num);
	float Bfrac=float(typeB)/float(num);
	float deadfrac=float(dead)/float(num);
	float uncommfrac=float(uncomm)/float(num);
	PD=PD0+log(float(num)/float(num0))/log(2.0);	
	rep[reprow][0]=simtime;
	rep[reprow][1]=PD;
	rep[reprow][2]=passno;
	rep[reprow][3]=num;
	rep[reprow][4]=dead;
	rep[reprow][5]=Afrac;
	rep[reprow][6]=Bfrac;
	rep[reprow][7]=deadfrac;
	rep[reprow][8]=uncommfrac;
	reprow=reprow+1;
}

void Senesce::inoculate(int Starting){
	starting=Starting;
	// Initialise counters and cell index
	simtime = 0.0;
	PD = 0.0;
	PD0 = 0.0;
	num = starting;
	num0 = starting;
	dividing = starting;
	trep = 0.0;
	passno = 1;
	reprow = 0;

	// Empty flask
	for(int i=0;i<capacity;i++){
		for(int j=0;j<3;j++){
			flask[i][j]=0;
		}
	}
	// Inoculate cells
	for(int i=0;i<starting;i++){
		flask[i][0]=1+gsl_rng_uniform_int (r, 2);
		flask[i][1]=meandiv;
		flask[i][2]=0;
	}
	std::cout << starting << std::endl;
}

void Senesce::growingCulture(){
	while (num<capacity){
		for(int i=0;i<capacity;i++) ind[i]=0;
		int numdiv=0;
		for(int i=0;i<num;i++){
			bool condition=((flask[i][1]>0)&&(flask[i][2]==1))|((flask[i][0]>0)&&(flask[i][2]==0));
			//bool condition=(flask[i][1]>0);
			if(condition){
				ind[numdiv]=i;
				numdiv+=1;
			}
		}
		if(numdiv==0) {
			PD=PD0+log(float(num)/float(num0))/log(2.0);			
			break;
		}
		if(simtime>=trep){
            		report();
			trep=trep+deltarep;
		}
		// Assume all cells have same probability of division
		float h0 = float(numdiv)*divrate;
		float deltat=gsl_ran_exponential(r,1.0/h0);
		simtime=simtime+deltat;
		// Pick a cell to divide (all equally probable)
		int candidate=gsl_rng_uniform_int(r,numdiv);
		// If cell uncommitted, roll 2 dice to define commitment state of daughters
		if (flask[ind[candidate]][2]==0){
			float rm = gsl_rng_uniform(r);
			float rd = gsl_rng_uniform(r);
            		// Copy label to daughter cell
            		flask[num][0]=flask[ind[candidate]][0];
            		if(rm<=P){
				// Mother commits
                    		flask[ind[candidate]][1]=meandiv;
                		flask[ind[candidate]][2]=1;
			}
           		if(rd<=P){
				// Daughter commits
                    		flask[num][1]=meandiv;
                		flask[num][2]=1;
			}	
		}else{
			// Create new daughter cell with same label, division potential and commitment state
			for(int j=0;j<nfeat;j++) flask[num][j]=flask[ind[candidate]][j];
			// Reduce cell proliferative potential of committed cells
			flask[ind[candidate]][1]-=1;
			flask[num][1]-=1;
		}
		// Update cell number
		num+=1;		
	}
	
}

void Senesce::passagingCells(){
	dividing=0;
	report();
	// Empty transfer flask
	for(int i=0; i<passage; i++) {
		for(int j=0; j<nfeat; j++) transfer[i][j]=0;
	}
	// Create and shuffle list of cells
	int * ind;
	ind = new int [num];
	for(int i=0; i<num;i++) ind[i]=i;
	gsl_ran_shuffle (r, ind, num, sizeof (int));
	// Transfer selected cells and count those that can divide
	for(int i=0; i<std::min(num,passage);i++){
		for(int j=0; j<nfeat; j++) transfer[i][j]=flask[ind[i]][j];
		bool condition=((flask[i][1]>0)&&(flask[i][2]==1))|((flask[i][0]>0)&&(flask[i][2]==0));
		if(condition) dividing+=1;
	}
	delete [] ind;
	// Create empty flask
	for(int i=0; i<capacity; i++){
		for(int j=0; j<nfeat; j++) flask[i][j]=0;
	}
	// Load empty flask with selected cells	
	for(int i=0; i<std::min(num,passage);i++){
		for(int j=0; j<nfeat; j++) flask[i][j]=transfer[i][j];
	}
	num=std::min(num,passage);
	num0=std::min(num,passage);
	passno=passno+1;
	PD0=PD;
	if(dividing>0) report();
}

void Senesce::simulate(){

	std::cout << "Inoculated population." << std::endl;
	for(int i=0;i<10;i++){
		std::cout << flask[i][0] << "\t" << flask[i][1] << "\t" << flask[i][2] << std::endl;
	}
	// Keep growing cells and passaging until population is exhausted
	while(dividing>0){
		int pstart=time(0);
		growingCulture();
		int ptime=time(0)-pstart;
		std::cout << "Confluent population." << std::endl;
		for(int i=0;i<10;i++){
			std::cout << flask[i][0] << "\t" << flask[i][1] << "\t" << flask[i][2] << std::endl;
		}
		passagingCells();
		if (dividing==0) report();
		std::cout << "Simulation time: " << ptime << std::endl;
	}

}

void Senesce::writeTable(){
	// Write results to file
	std::ofstream outfile;
	int resint = gsl_rng_uniform_int(r,999999999);
	char filename[50];
	sprintf (filename, "ResBig_%09d.txt", resint);
	outfile.open(filename);
	// Make a header
	outfile << "Time" << "\t" << "PD" << "\t" << "Passage" << "\t" << "Cells" << "\t" << "Dead" << "\t" << "FracA" << "\t" << "FracB" << "\t" << "FracDead" << "\t" << "FracCommit" << "\t";
	for (int i=0;i<histbins-1;i++) outfile << "D" << i << "\t";
	outfile << "Uncommitted" << "\t" << std::endl;

	for (int i=0;i<reprow;i++){
		for (int j=0;j<repcols;j++){
			outfile << rep[i][j] << "\t";
		}
		for (int j=0;j<histbins;j++){
			outfile << hist[i][j] << "\t";
		}
		outfile << std::endl;
	}
	outfile.close();
}
