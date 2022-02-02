#ifndef BCCD_LIFE_H
#define BCCD_LIFE_H
//#define Nthreads 5

#include "XLife.h"    // For display routines
#include "Defaults.h" // For Life's constants

#include <time.h>     // For seeding random
#include <stdlib.h>   // For malloc et al.
#include <stdbool.h>  // For true/false
#include <getopt.h>   // For argument processing
#include <string.h>
#include <unistd.h>
//#include <mpi.h>
#include "mpi.h"
#include <omp.h>
#include <stdio.h>    // For file i/o

int         init (struct life_t * life, int * c, char *** v);
void        callEvalRules(struct life_t  life,int *argc, char *** argv);
void        reverse_corners(struct life_t * life);
void        eval_rules (struct life_t * life,int rank);
void       copy_bounds (struct life_t * life);
void       update_grid (struct life_t * life);
void          throttle (struct life_t * life);
void        allocate_grids (struct life_t * life);
void        init_grids (struct life_t * life);
void        write_grid (struct life_t * life);
void        free_grids (struct life_t * life);
double     rand_double ();
void    randomize_grid (struct life_t * life, double prob);
void           cleanup (struct life_t * life);
void        parse_args (struct life_t * life, int argc, char ** argv);
void             usage ();


/*
	init_env()
		Initialize runtime environment.
*/
int init (struct life_t * life, int * c, char *** v) {
	int argc          = *c;
	char ** argv      = *v;
	life->size        = 1;
	life->throttle    = -1;
	life->ncols       = DEFAULT_SIZE;
	life->nrows       = DEFAULT_SIZE;
	life->generations = DEFAULT_GENS;
	life->do_display  = DEFAULT_DISP;
	life->infile      = NULL;
	life->outfile     = NULL;

	srandom(time(NULL));

	parse_args(life, argc, argv);

	init_grids(life);

	if (life->do_display) {
		setupWindow(life);
		moveWindow(life);
	}
}

/*
	eval_rules()
		Evaluate the rules of Life for each cell; count
		neighbors and update current state accordingly.
*/

void callEvalRules(struct life_t  life,int *argc, char *** argv){
	int rank;
       MPI_Init(argc,argv);
       MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
	   do_draw(&life);
	   eval_rules(&life,rank);
	   MPI_Finalize();
	   reverse_corners(&life);
	   update_grid(&life);	//	Updates the old grid with the state value of each cell.

		throttle(&life);	//	Slows down the simulation to make X display easier to watch.
	   

}

void reverse_corners(struct life_t * life){
   int ** grid      = life->grid;
   int ** ngrid  =life->next_grid;
   int ncols = life->ncols;
	int nrows = life->nrows;
   
   ngrid[0][0]=grid[ncols][nrows];
   ngrid[ncols][nrows]=grid[0][0];
  
   ngrid[0][nrows]=grid[ncols][0];
   ngrid[ncols][0]=grid[0][nrows];
  
}
void sendrows(struct life_t * life,int rank){
	int *ghost_upper,*ghost_lower;
	int ncols=life->ncols;
	ghost_lower=(int *) malloc(sizeof(int) * (ncols));
	ghost_upper=(int *) malloc(sizeof(int) * (ncols));
	int jmax,i;

	int ** grid      = life->grid;
    MPI_Request request,request1;
    if(rank==9){
        jmax=rank*10+15;
    }
    else{
        jmax=rank*10+10;
    }


    for(i=0 ; i<ncols ; i++){
        ghost_upper[i]=grid[i][0];
    }
    for(i=0 ; i<ncols ; i++){
		if (rank==9){
         ghost_lower[i]=grid[i][rank*10+14];
		}
		else
		{
			 ghost_lower[i]=grid[i][rank*10+9];
		}
		
    }
    if(rank!=0){
    MPI_Isend(&ghost_upper,DEFAULT_SIZE,MPI_INT,rank-1,0,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
       if(rank!=9){
         MPI_Isend(&ghost_lower,DEFAULT_SIZE,MPI_INT,rank+1,1,MPI_COMM_WORLD,&request1);
         MPI_Wait(&request1,MPI_STATUS_IGNORE);
        }
         else{
          MPI_Isend(&ghost_lower,DEFAULT_SIZE,MPI_INT,0,1,MPI_COMM_WORLD,&request1);
          MPI_Wait(&request1,MPI_STATUS_IGNORE);
        }

    } 
    else{
    MPI_Isend(&ghost_upper,DEFAULT_SIZE,MPI_INT,9,0,MPI_COMM_WORLD,&request);
    MPI_Wait(&request,MPI_STATUS_IGNORE);
    MPI_Isend(&ghost_lower,DEFAULT_SIZE,MPI_INT,rank+1,1,MPI_COMM_WORLD,&request1);
         MPI_Wait(&request1,MPI_STATUS_IGNORE);
    }

}

void eval_rules (struct life_t * life,int rank) {
    int i,j,k,l,neighbors,jmax;

	int ncols = life->ncols;
	int nrows = life->nrows;

	int ** grid      = life->grid;
	int ** next_grid = life->next_grid;
    int  *ghost_upper_row , *ghost_lower_row;

	ghost_upper_row=(int *) malloc(sizeof(int) * (ncols));
	ghost_lower_row=(int *) malloc(sizeof(int) * (ncols));
      if(rank==9){
        jmax=rank*10+15;
    }
    else{
        jmax=rank*10+10;
    }

     sendrows(life,rank);
   
    if(rank!=0){
        MPI_Recv(&ghost_upper_row,DEFAULT_SIZE,MPI_INT,rank-1,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        if(rank==9)
        MPI_Recv(&ghost_lower_row,DEFAULT_SIZE,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        else{
            MPI_Recv(&ghost_lower_row,DEFAULT_SIZE,MPI_INT,rank+1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
    }
    else{
        MPI_Recv(&ghost_lower_row,DEFAULT_SIZE,MPI_INT,rank+1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Recv(&ghost_upper_row,DEFAULT_SIZE,MPI_INT,9,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
 //omp_set_dynamic(1);

#pragma omp parallel for num_threads(5) private(j)
 for(i=0; i< ncols; i++){ 
     for(j=rank*10; j<jmax ; j++){
        
		 neighbors=0;

		  if(j==rank*10 || j==jmax-1){

			   if(j==jmax-1){
                 //the lower left corner
                if(i==0 && j!= nrows-1){
					 

                 if(grid[i][j-1]!=DEAD){
					 #pragma omp critical
                   neighbors++;
				 }
                 if(grid[i+1][j-1]!=DEAD){
					 #pragma omp critical
                   neighbors++;
				 }
                 if(grid[i+1][j]!=DEAD){
					 #pragma omp critical
                   neighbors++;   
				 }
                 if(grid[ncols][j]!=DEAD){
					 #pragma omp critical
                   neighbors++;
				 }
                 if(grid[ncols][j-1]!=DEAD){
					 #pragma omp critical
                   neighbors++;
				 }
                 if(ghost_lower_row[0]!=DEAD){
					 #pragma omp critical
                    neighbors++;
				 }
                 if(ghost_lower_row[1]!=DEAD){
					 #pragma omp critical
                    neighbors++;
				 }
                 if(ghost_lower_row[ncols]!=DEAD){
					 #pragma omp critical
                    neighbors++;
				 }
				
             }//i==0
             //the lower right corner
             else if(i==ncols-1 && j!=nrows-1){
				
                 if(grid[i-1][j]!=DEAD){
					 #pragma omp critical
                   neighbors++;
				 }
                 if(grid[i][j-1]!=DEAD){
					 #pragma omp critical
                   neighbors++;
				 }
                 if(grid[i-1][j-1]!=DEAD){
					 #pragma omp critical
                   neighbors++;  
				 }
                 if(grid[0][j]!=DEAD){
					 #pragma omp critical
                   neighbors++;
				 }
                 if(grid[0][j-1]!=DEAD){
					 #pragma omp critical
                   neighbors++;
				 }
                 if(ghost_lower_row[0]!=DEAD){
					 #pragma omp critical
                    neighbors++;
				 }
                 if(ghost_lower_row[ncols-2]!=DEAD){
					 #pragma omp critical
                    neighbors++;
				 }
                 if(ghost_lower_row[ncols-1]!=DEAD){
					 #pragma omp critical
                    neighbors++;
				 }
				
 
             }//ncols-1

             // all cells in the last row except the corners
             else{
                   for(k=i-1;k<i+1;k++){
                     for(l=j;l<j-1;l++){
                         if(!(k==i&&l==j)&& grid[k][l]!=DEAD){
							 #pragma omp critical
                             neighbors++;
							 
                         }
                     }
                  if(ghost_lower_row[k]!=DEAD){
				     #pragma omp critical
                     neighbors++;                        
				     
                   }
				 
               
				   }
              }//end of else   
    
             }//jmax-1
			 if(j==rank*10){
                 //the upper right corner
                if(i==ncols-1&& j!=0){
			     

                 if(grid[i-1][j]!=DEAD){
					 #pragma omp critical
					 neighbors++;
				 }
                 if(grid[i][j-1]!=DEAD){
					 #pragma omp critical
                   neighbors++;
				 }
                 if(grid[i-1][j-1]!=DEAD){
					 #pragma omp critical
                   neighbors++;   
				 }
                 if(grid[0][j]!=DEAD){
					 #pragma omp critical
                   neighbors++;
				 }
                 if(grid[0][j-1]!=DEAD){
					 #pragma omp critical
                   neighbors++;
				 }
                 if(ghost_lower_row[0]!=DEAD){
					 #pragma omp critical
                    neighbors++;
				 }
                 if(ghost_lower_row[ncols-2]!=DEAD){
					 #pragma omp critical
                    neighbors++;
				 }
                 if(ghost_lower_row[ncols-1]!=DEAD){
					 #pragma omp critical
                    neighbors++;
				 }
				  
                }//end of if
               //the upper left corner
                if(i==0 && j!=0){
				 
                if(grid[i+1][j]!=DEAD){
					#pragma omp critical
                   neighbors++;
				}
                 if(grid[i][j+1]!=DEAD){
					 #pragma omp critical
                   neighbors++;
				 }
                 if(grid[i+1][j+1]!=DEAD){
					 #pragma omp critical
                   neighbors++;   
				 }
                 if(grid[ncols][j+1]!=DEAD){
					 #pragma omp critical
                   neighbors++;
				 }
                 if(grid[ncols][j]!=DEAD){
					 #pragma omp critical
                   neighbors++;
				 }
                 if(ghost_upper_row[0]!=DEAD){
					 #pragma omp critical
                    neighbors++;
				 }
                 if(ghost_upper_row[1]!=DEAD){
					 #pragma omp critical
                    neighbors++;
				 }
                 if(ghost_upper_row[ncols]!=DEAD){
					 #pragma omp critical
                    neighbors++;
				 }
                   
			    }//end of if
             //all cells in the upper row except the corners
             else{
                   for(k=i-1;k<i+1;k++){
                     for(l=j;l<j+1;l++){
                         if(!(k==i&&l==j)&& grid[k][l]!=DEAD){
							 #pragma omp critical
                             neighbors++;
							 
                         }
                     }
                  if(ghost_upper_row[k]!=DEAD){
				  #pragma omp critical
                     neighbors++;                        
				  
                 }
				   }
                }//end of else

             }//rank*10 

		  }
		   else if(i==0){
            // all cells in the first colomn except the corners
             if(j!=rank*10 || j!=jmax-1){
                 for(l=j-1;l<j+1;l++){
                     for(k=i;k<i+1;k++){
                         if(!(k==i&&l==j)&& grid[k][l]!=DEAD){
						 #pragma omp critical
                             neighbors++; 
						 
					 }                       
                     }
                     if(grid[ncols][j]!=DEAD){
					  #pragma omp critical
                       neighbors++;
						 
					 }
                 }
             }//end of else if

			   //all cells in last colomn except the corners
        else if(i==ncols-1){
             if(j!=rank*10 || j!=jmax-1){
                 for(l=j-1;l<j+1;l++){
                     for(k=i;k<i-11;k++){
                         if(!(k==i&&l==j)&& grid[k][l]!=DEAD){
							 #pragma omp critical
                             neighbors++;  
							 
						 }                      
                     }
                     if(grid[0][j]!=DEAD){
						#pragma omp critical
                       neighbors++;
						
					 }
                 }
             }
         }//end of else if 

		 // all cells inside the grid except edges and corners
         else{
             for (k = i-1; k <= i+1; k++) {
				for (l = j-1; l <= j+1; l++) {
					if (!(k == i && l == j) && grid[k][l] != DEAD){
						#pragma omp critical
						neighbors++;
						
					}
				}
			}
         }//end of else 

		 // update state
         if (neighbors < LOWER_THRESH || neighbors > UPPER_THRESH)
				next_grid[i][j] = DEAD;

		 else if (grid[i][j] != DEAD || neighbors == SPAWN_THRESH){
			#pragma omp critical
				next_grid[i][j] = grid[i][j]+1;
			
          }
         
	 }
 }
}

}//end of function 

/*
	update_grid()
		Copies temporary values from next_grid into grid.
*/
void update_grid (struct life_t * life) {
	int i,j;
	int ncols = life->ncols;
	int nrows = life->nrows;
	int ** grid      = life->grid;
	int ** next_grid = life->next_grid;

	for (i = 0; i < ncols; i++)
		for (j = 0; j < nrows; j++)
			grid[i][j] = next_grid[i][j];
}

/*
	throttle()
		Slows down the simulation to make X display easier to watch.
		Has no effect when run with --no-display.
*/
void throttle (struct life_t * life) {
	unsigned int delay;
	int t = life->throttle;

	if (life->do_display && t != -1) {
		delay = 1000000 * 1/t;
		usleep(delay);
	}
}

/*
	allocate_grids()
		Allocates memory for a 2D array of integers.
*/
void allocate_grids (struct life_t * life) {
	int i,j;
	int ncols = life->ncols;
	int nrows = life->nrows;

	life->grid      = (int **) malloc(sizeof(int *) * (ncols));
	life->next_grid = (int **) malloc(sizeof(int *) * (ncols));

	for (i = 0; i < ncols; i++) {
		life->grid[i]      = (int *) malloc(sizeof(int) * (nrows));
		life->next_grid[i] = (int *) malloc(sizeof(int) * (nrows));
	}
}
void init_grids (struct life_t * life) {
	FILE * fd;
	int i,j;

	if (life->infile != NULL) {
		if ((fd = fopen(life->infile, "r")) == NULL) {
			perror("Failed to open file for input");
			exit(EXIT_FAILURE);
		}

		if (fscanf(fd, "%d %d\n", &life->ncols, &life->nrows) == EOF) {
			printf("File must at least define grid dimensions!\nExiting.\n");
			exit(EXIT_FAILURE);
		}
	}

	allocate_grids(life);

	for (i = 0; i < life->ncols; i++) {
		for (j = 0; j < life->nrows; j++) {
			life->grid[i][j]      = DEAD;
			life->next_grid[i][j] = DEAD;
		}
	}

	if (life->infile != NULL) {
		while (fscanf(fd, "%d %d\n", &i, &j) != EOF) {
			life->grid[i][j]      = ALIVE;
			life->next_grid[i][j] = ALIVE;
		}
		
		fclose(fd);
	} else {
		randomize_grid(life, INIT_PROB);
	}
}
/*
	write_grid()
		Dumps the current state of life.grid to life.outfile.
		Only outputs the coordinates of !DEAD cells.
*/
void write_grid (struct life_t * life) {
	FILE * fd;
	int i,j;
	int ncols   = life->ncols;
	int nrows   = life->nrows;
	int ** grid = life->grid;

	if (life->outfile != NULL) {
		if ((fd = fopen(life->outfile, "w")) == NULL) {
			perror("Failed to open file for output");
			exit(EXIT_FAILURE);
		}

		fprintf(fd, "%d %d\n", ncols, nrows);

		for (i = 0; i < ncols; i++) {
			for (j = 0; j < nrows; j++) {
				if (grid[i][j] != DEAD)
					fprintf(fd, "%d %d\n", i, j);
			}
		}

		fclose(fd);
	}
}
/*
	free_grids()
		Frees memory used by an array that was allocated 
		with allocate_grids().
*/
void free_grids (struct life_t * life) {
	int i;
	int ncols = life->ncols;

	for (i = 0; i < ncols; i++) {
		free(life->grid[i]);
		free(life->next_grid[i]);
	}

	free(life->grid);
	free(life->next_grid);
}

/*
	rand_double()
		Generate a random double between 0 and 1.
*/
double rand_double() {
	return (double)random()/(double)RAND_MAX;
}

/*
	randomize_grid()
		Initialize a Life grid. Each cell has a [prob] chance
		of starting alive.
*/
void randomize_grid (struct life_t * life, double prob) {
	int i,j;
	int ncols = life->ncols;
	int nrows = life->nrows;

	for (i = 0; i < ncols; i++) {
		for (j = 0; j < nrows; j++) {
			if (rand_double() < prob)
				life->grid[i][j] = ALIVE;
		}
	}
}

/*
	cleanup()
		Prepare process for a clean termination.
*/
void cleanup (struct life_t * life) {
	write_grid(life);
	free_grids(life);

	if (life->do_display)
		free_video(life);

}

/*
	usage()
		Describes Life's command line option
*/
void usage () {
	printf("\nUsage: Life [options]\n");
	printf("  -c|--columns number   Number of columns in grid. Default: %d\n", DEFAULT_SIZE);
	printf("  -r|--rows number      Number of rows in grid. Default: %d\n", DEFAULT_SIZE);
	printf("  -g|--gens number      Number of generations to run. Default: %d\n", DEFAULT_GENS);
	printf("  -i|--input filename   Input file. See README for format. Default: none.\n");
	printf("  -o|--output filename  Output file. Default: none.\n");
	printf("  -h|--help             This help page.\n");
	printf("  -t[N]|--throttle[=N]  Throttle display to N generations/second. Default: %d\n",
		DEFAULT_THROTTLE);
	printf("  -x|--display          Use a graphical display.\n");
	printf("  --no-display          Do not use a graphical display.\n"); 
	printf("     Default: %s\n",
		(DEFAULT_DISP ? "do display" : "no display"));
	printf("\nSee README for more information.\n\n");

	exit(EXIT_FAILURE);
}
/*
	parse_args()
		Make command line arguments useful
*/
void parse_args (struct life_t * life, int argc, char ** argv) {
	int opt       = 0;
	int opt_index = 0;
	int i;

	for (;;) {
		opt = getopt_long(argc, argv, opts, long_opts, &opt_index);

		if (opt == -1) break;

		switch (opt) {
			case 0:
				if (strcmp("no-display", long_opts[opt_index].name) == 0)
					life->do_display = false;
				break;
			case 'c':
				life->ncols = strtol(optarg, (char**) NULL, 10);
				break;
			case 'r':
				life->nrows = strtol(optarg, (char**) NULL, 10);
				break;
			case 'g':
				life->generations = strtol(optarg, (char**) NULL, 10);
				break;
			case 'x':
				life->do_display = true;
				break;
			case 'i':
				life->infile = optarg;
				break;
			case 'o':
				life->outfile = optarg;
				break;
			case 't':
				if (optarg != NULL)
					life->throttle = strtol(optarg, (char**) NULL, 10);
				else
					life->throttle = DEFAULT_THROTTLE;
				break;
			case 'h':
			case '?':
				usage();
				break;

			default:
				break;
		}
	}

	// Backwards compatible argument parsing
	if (optind == 1) {
		if (argc > 1)
			life->nrows       = strtol(argv[1], (char**) NULL, 10);
		if (argc > 2)
			life->ncols       = strtol(argv[2], (char**) NULL, 10);
		if (argc > 3)
			life->generations = strtol(argv[3], (char**) NULL, 10);
		if (argc > 4)
			// 0 interpreted as false, all other values true
			life->do_display  = strtol(argv[4], (char**) NULL, 10);
	}
}



#endif
