#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <string>
using namespace std;

class Senesce{

public:
	int capacity;   // Size of flask (or size of 50% of flask if not allowing confluence)
	int starting;   // Starting population size
	float prat; 	// Fraction of cells to transfer at each passage
	int passage; 	// Number of cells to transfer at each passage
	float divrate;  // Divisions per cell per unit time
	float deltarep; // Report about cell population every deltarep time units
	float P;        // Probability of uncomitted cell becoming committed upon division
	int meandiv;    // Mean number of divisions a (committed) cell undergoes before senescence
	int stddiv;     // Standard deviation for number divisions a committed cell undergoes before senescence

	// 0: Label, 1: Division potential, 2: Commitment state
	int nfeat; 	// Number of features a cell has
	int nreps;	// Maximum number of reports (ie max simulation time = deltarep*nreps)
	int repcols;	// Number of descriptions in rep (report) object
	int histbins;   // Number of bins for remaining division potential +1

	// Global recording variables
	int ** flask;
	int ** transfer;
	int * ind;
	float simtime, PD, PD0, trep;
	int num, num0, passno, reprow, dividing;
	float ** rep;
	int ** hist;

	gsl_rng *r; 
	
	Senesce(int Capacity,float Prat,float Divrate,float Deltarep,float Pc,float Meandiv,float Stddiv);
    	~Senesce();  

	void report();
	void inoculate(int Starting);
	void growingCulture();
	void passagingCells();
	void simulate();
	void writeTable(string fileRoot);
};
