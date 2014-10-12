#include "senesce.h"

int main()
{
	for(int i=0; i<200; i++){
		Senesce cells(1000000,0.25,0.5,1.0,0.275,43,0);
		cells.inoculate(1);
		cells.simulate();
		cells.writeTable("FullLifecycle");
	}
	return 0;
}


