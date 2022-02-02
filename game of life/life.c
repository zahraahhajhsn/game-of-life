

#include "Life.h"		// For function's definitions and instructions.
#include "Defaults.h"
#include <mpi.h> 	// For Life's constants

int main(int argc, char ** argv) {
int rank;
int count;
struct life_t life;

init(&life, &argc, &argv);


//MPI_Init(&argc,&argv);

//MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	for (count = 0; count < life.generations; count++) {
	//	MPI_Init(&argc,&argv);

       // MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	   if (life.do_display)
		//do_draw(&life);
			
		callEvalRules(life,&argc,&argv);
		//eval_rules(&life,rank);	//	Finds out how many neighbor's each cell has then updates 
							//the cell's state to DEAD or ALIVE.
		//reverse_corners(&life);
	//	update_grid(&life);	//	Updates the old grid with the state value of each cell.

		//throttle(&life);	//	Slows down the simulation to make X display easier to watch.
	}

   
	cleanup(&life);			// Writes the coordinates of all the ALIVE cells into the output files
							//and frees all the memory used by the program
	 
	exit(EXIT_SUCCESS);



//MPI_Finalize();

return 0;
}

