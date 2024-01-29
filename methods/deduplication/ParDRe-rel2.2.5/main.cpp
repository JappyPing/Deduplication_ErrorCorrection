#include "SingleEnd.h"
#include "PairedEnd.h"
#include "Options.h"

int main(int argc, char* argv[]) {

	// Initialize MPI
	MPI::Init(argc,argv);

	Options* options = new Options();

	double stime, etime;

	int numP=MPI::COMM_WORLD.Get_size();
	int myRank=MPI::COMM_WORLD.Get_rank();
	Utils::log("Process %d/%d: Initialized\n", myRank, numP);

	/*parse the arguments*/
	if (!options->parse(argc, argv)) {
		options->printUsage();
		return 0;
	}

	MPI::COMM_WORLD.Barrier();

	/*get the startup time*/
	stime = MPI::Wtime();

	/*run the engines*/
	if (options->isPaired() == false) {
		/*for single-end alignment*/
		SingleEnd* single = new SingleEnd(options);

		/*run the alignment*/
		single->execute();
		/*release aligner*/
		delete single;
	} else {
		/*for single-end alignment*/
		PairedEnd* paired = new PairedEnd(options);

		/*run the alignment*/
		paired->execute();
		/*release aligner*/
		delete paired;
	}

	/*report the wall clock time*/
	etime = MPI::Wtime();
	fprintf(stderr, "Process %d/%d: Overall time: %.2f seconds\n", myRank, numP, etime - stime);

	// Terminate MPI
	MPI::Finalize();
	return 0;
}
