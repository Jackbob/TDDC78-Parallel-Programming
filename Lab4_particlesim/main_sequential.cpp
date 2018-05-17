#include <cstdlib>
#include <ctime>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <vector>
#include <mpi.h>
#include "coordinate.h"
#include "definitions.h"
#include "physics.h"


//Feel free to change this program to facilitate parallelization.

float rand1(){
	return (rand()/(float) RAND_MAX);
}

void init_collisions(std::vector<bool> collisions, unsigned int max){
	for(unsigned int i=0;i<max;++i)
		collisions[i]= false;
}

void calculateRowCol(int& row, int& col, int nProc){
	row = (int)std::ceil(std::sqrt(nProc));
	col = (int)std::floor(std::sqrt(nProc));
}


int main(int argc, char** argv){
	/* Define variables */
	unsigned int time_stamp = 0, time_max;
	float pressure = 0;
	int rank{}, world{}, root{0}, gridcol, gridrow;
	int dims[2], period[2], coord[2];
	std::vector<Particle> Particles;
	std::vector<bool> collisions;
	cord_t wall;


	/* Initialize MPI environment */
	MPI_Init(nullptr,nullptr);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world);
	MPI_Comm grid_comm;

    /* Create MPI data type for particles */
	MPI_Datatype MPI_Particle, oldtypes[1];
	int blockcounts[1];
	MPI_Aint offsets[1];
	offsets[0] = 0;
	oldtypes[0] = MPI_FLOAT;
	blockcounts[0] = 4;
	MPI_Type_create_struct(1, blockcounts, offsets, oldtypes, &MPI_Particle);
	MPI_Type_commit(&MPI_Particle);

	calculateRowCol(gridrow, gridcol, world);
	std::cout << gridrow << std::endl;
	std::cout << gridcol << std::endl;
	dims[0] = gridrow;
	dims[1] = gridcol;
	period[0] = period[1] = 0;
	MPI_Dims_create(world, 2, dims); //Creates a division of processors in 2-D cartesian grid.
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, 2  ,&grid_comm);


	// parse arguments
	if(argc != 2) {
		fprintf(stderr, "Usage: %s simulation_time\n", argv[0]);
		fprintf(stderr, "For example: %s 10\n", argv[0]);
		exit(1);
	}

	time_max = atoi(argv[1]);

    std::cout << "Rank " << rank << " out of " << world << "\n";
	MPI_Cart_coords(grid_comm, rank, 2, coord);
	printf("Rank %d coordinates are %d %d\n", rank, coord[0], coord[1]);

    //Initializes the wall and Particle on the root/master processor.

	// 1. set the walls
	wall.y0 = wall.x0 = 0;
	wall.x1 = BOX_HORIZ_SIZE;
	wall.y1 = BOX_VERT_SIZE;

	// 2. allocate particle buffer and initialize the Particle
	Particles = std::vector<Particle>(INIT_NO_PARTICLES);
	collisions = std::vector<bool>(INIT_NO_PARTICLES);



	srand( time(nullptr) + 1234 );

	float r, a;
	for(int i=0; i<INIT_NO_PARTICLES; i++){
		// initialize random position
		Particles[i].x = static_cast<float>(wall.x0 + rand1()*BOX_HORIZ_SIZE);
		Particles[i].y = static_cast<float>(wall.y0 + rand1()*BOX_VERT_SIZE);

		// initialize random velocity
		r = rand1()*MAX_INITIAL_VELOCITY;
		if(r < 50) r = 49;
		a = static_cast<float>(rand1()*2*PI);
		Particles[i].vx = r*cos(a);
		Particles[i].vy = r*sin(a);
	}


	unsigned int p, pp;

	/* Main loop */
	for (time_stamp=0; time_stamp<time_max; time_stamp++) { // for each time stamp

		init_collisions(collisions, INIT_NO_PARTICLES);

		for(p=0; p<INIT_NO_PARTICLES; p++) { // for all Particle
			if(collisions[p]) continue;

			/* check for collisions */
			for(pp=p+1; pp<INIT_NO_PARTICLES; pp++){
				if(collisions[pp]) continue;
				float t=collide(&Particles[p], &Particles[pp]);
				if(t!=-1){ // collision
					collisions[p]=collisions[pp]=1;
					interact(&Particles[p], &Particles[pp], t);
					break; // only check collision of two Particle
				}
			}

		}

		// move Particle that has not collided with another
		for(p=0; p<INIT_NO_PARTICLES; p++)
			if(!collisions[p]){
				feuler(&Particles[p], 1);

				/* check for wall interaction and add the momentum */
				pressure += wall_collide(&Particles[p], wall);
			}


	}

	printf("Average pressure = %f\n", pressure / (WALL_LENGTH*time_max));
    MPI_Finalize();
	return 0;

}

