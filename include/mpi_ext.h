#ifndef _MPI_EXT_H
#define _MPI_EXT_H
#ifdef WITH_MPI
#include "mpi.h"
#endif

#define ROOT 0

#ifdef WITH_MPI
	#define BEGIN_ROOT if(rank == ROOT){
	#define END_ROOT }
#else
	#define BEGIN_ROOT
	#define END_ROOT
#endif 

extern int rank;

#endif