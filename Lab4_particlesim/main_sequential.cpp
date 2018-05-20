#include <cstdlib>
#include <ctime>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <vector>
#include <mpi.h>
#include <algorithm>
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

int calcRowRank(float y){
	int n_proc;
	MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
	float rowSplit = (float)BOX_VERT_SIZE / n_proc;
	return (int)(y/rowSplit);
}

bool boundaryCheck(Particle &p, cord_t wall){

	if((p.y < wall.y0 || p.y > wall.y1) && (p.y > 0 && p.y < BOX_VERT_SIZE) )
		return true;
	else
		return false;
}


int main(int argc, char** argv){
	/* Define variables */
	unsigned int time_stamp = 0, time_max;
	float pressure = 0;
	int rank{}, world{}, root{0};
	std::vector<Particle> Particles;
	std::vector<Particle> sendParticlesUp;
	std::vector<Particle> sendParticlesDown;
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


	// parse arguments
	if(argc != 2) {
		fprintf(stderr, "Usage: %s simulation_time\n", argv[0]);
		fprintf(stderr, "For example: %s 10\n", argv[0]);
		exit(1);
	}

	time_max = atoi(argv[1]);

    std::cout << "Rank " << rank << " out of " << world << "\n";

    //Initializes the wall and Particle on the root/master processor.

	// 1. set the walls
	float rowSplit = (float)(BOX_VERT_SIZE/world);
	wall.y0 = rank * rowSplit;
	wall.y1 = wall.y0 + rowSplit;
	wall.x0 = 0;
	wall.x1 = BOX_HORIZ_SIZE;


	// 2. allocate particle buffer and initialize the Particle
	Particles = std::vector<Particle>(INIT_NO_PARTICLES);
	collisions = std::vector<bool>(INIT_NO_PARTICLES);


	srand( time(nullptr) + 1234 );

	float r, a;
	for(int i=0; i<INIT_NO_PARTICLES; i++){
		// initialize random position
		Particles[i].x = static_cast<float>(wall.x0 + rand1()*BOX_HORIZ_SIZE);
		Particles[i].y = (wall.y0 + rand1()*rowSplit);

		// initialize random velocity
		r = rand1()*MAX_INITIAL_VELOCITY;
		if(r < 50) r = 49;
		a = static_cast<float>(rand1()*2*PI);
		Particles[i].vx = r*cos(a);
		Particles[i].vy = r*sin(a);
	}


	unsigned int p, pp;

	MPI_Request reqUp, reqDown;
	MPI_Status statUp, statDown, statRec;
	Particle *sendBufDown, *sendBufUp;
	int flag = 0;

	/* Main loop */
	for (time_stamp=0; time_stamp<time_max; time_stamp++) { // for each time stamp

		if(flag){
			int recCount{};
			MPI_Get_count(&statRec, MPI_Particle, &recCount);
			Particle *recbuf = new Particle[recCount];
			MPI_Recv(recbuf, recCount, MPI_Particle, statRec.MPI_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for(int i=0; i<recCount; i++)
				Particles.emplace_back(recbuf[i]);

			delete[] recbuf;
			std::cout << "fisk \n";
		}

		init_collisions(collisions, collisions.size());
		sendParticlesUp.clear();
		sendParticlesDown.clear();

		for(p=0; p<Particles.size(); p++) { // for all Particle
			if(collisions[p]) continue;

			/* check for collisions */
			for(pp=p+1; pp<Particles.size(); pp++){
				if(collisions[pp]) continue;

				float t=collide(&Particles[p], &Particles[pp]);
				if(t!=-1){ // collision
					collisions[p] = collisions[pp] = true;
					interact(&Particles[p], &Particles[pp], t);


					if(boundaryCheck(Particles[p], wall)) {
						if(calcRowRank(Particles[p].y) < rank)
							sendParticlesUp.emplace_back(Particles[p]);
						else
							sendParticlesDown.emplace_back(Particles[p]);

						std::cout << "Erasing... \n" ;
						//Particles.erase( Particles.begin() + p );
						//collisions.erase( collisions.begin() + p );
					}
					if(boundaryCheck(Particles[pp], wall)) {
						if(calcRowRank(Particles[pp].y) < rank)
							sendParticlesUp.emplace_back(Particles[pp]);
						else
							sendParticlesDown.emplace_back(Particles[pp]);

						std::cout << "Erasing... \n" ;
						//Particles.erase( Particles.begin() + pp );
						//collisions.erase( collisions.begin() + pp );
					}



					break; // only check collision of two Particle
				}
			}

		}

		// move Particle that has not collided with another
		for(p=0; p<Particles.size(); p++)
			if(!collisions[p]){
				feuler(&Particles[p], 1);


				if(boundaryCheck(Particles[p], wall)) {

					if(calcRowRank(Particles[p].y) < rank)
						sendParticlesUp.emplace_back(Particles[p]);
					else
						sendParticlesDown.emplace_back(Particles[p]);

					//Particles.erase( Particles.begin() + p );
					//collisions.erase( collisions.begin() + p );
				}
				/* check for wall interaction and add the momentum */
				pressure += wall_collide(&Particles[p], wall);
			}

		if(!sendParticlesUp.empty()) {
			std::cout << "Sending... \n" ;
			//MPI_Wait(&reqUp, &statUp);
			delete[] sendBufUp;
			int sendcount = static_cast<int>(sendParticlesUp.size());
			sendBufUp = new Particle[sendcount];
			for(int i=0; i<sendParticlesUp.size(); i++)
				sendBufUp[i] = sendParticlesUp[i];
			//MPI_Isend(sendBufUp, sendcount, MPI_Particle, rank-1, 0, MPI_COMM_WORLD, &reqUp);
		}

		if(!sendParticlesDown.empty()) {
			std::cout << "Sending... \n" ;
			//MPI_Wait(&reqDown, &statDown);
			delete[] sendBufDown;
			int sendcount = static_cast<int>(sendParticlesDown.size());
			sendBufDown = new Particle[sendcount];
			for(int i=0; i<sendParticlesDown.size(); i++)
				sendBufDown[i] = sendParticlesDown[i];
			//MPI_Isend(sendBufDown, sendcount, MPI_Particle, rank+1, 0, MPI_COMM_WORLD, &reqDown);
		}


		//MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &statRec);

	}

	printf("Average pressure = %f\n", pressure / (WALL_LENGTH*time_max));
    MPI_Finalize();
	return 0;

}

