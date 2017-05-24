
#include <time.h>				// Get clock time, to measure total run time
#include <iomanip>
#include "bond.h"
#include "site.h"
#include "siteRP.h"
#include "SiteRP_cont.h"


int main(int argc, char* argv[])
{
/*
	if (argc < 4) { // We expect 3 arguments: the program name, the source path and the destination path
		std::cerr << "Usage: " << argv[0] << " <CORRELATION> <TRIAL-COUNTING> <RANDOMSEED>" << std::endl;
		return 1;
	}

	clock_t t1, t2;													// Creates clock variables to track the runtime
	t1 = clock();
	srand(atoi(argv[3]));

	float cfor = atof(argv[1]);
	int trial = atoi(argv[2]);

	SiteRP TriLattice;
	TriLattice.OneTrialTest(cfor,trial);

	t2 = clock();
	float clocktime((float)t2 - (float)t1);
	std::cout << "\n The total run time was " << clocktime / CLOCKS_PER_SEC << std::endl;
	return 0;
*/

    SiteRP_cont a;
    a.ReadVerticeInfo();
    a.BuildNetwork();
    a.CoordNumber();
    //a.test();
    //a.RigidAtomWriteBack();
}
