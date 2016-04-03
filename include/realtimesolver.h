#ifndef _REALTIMESOLVER_H
#define _REALTIMESOLVER_H

#include <sys_types.h>
#include <quad_types.h>
#include <riccati.h>

#define GETLAMBDA(t, lambda) (lambda + (NLAMBDA * t))

typedef struct realtimesolver_res_struct
{
	realtype s[NSTATE];
	realtype q[NCONTR];
	realtype lambda[NLAMBDA];
	rqci mu[NADDCONSTR];
}rtsolres;

typedef struct realtimesolver_struct
{
	realtype *vec;
	realtype *lambda; //NLAMBDA
	rqci* activeSet;
	rqci* muSet;
	riccati* Ric;
	QuadConfig * conf;
	rtsolres *res;
	rqci ntime;
	rqci nhorizon;
	realtype mesh_h;
	realtype* data;
}realtimesolver;


realtimesolver* initialize_rtsolver(rqci nhorizon, realtype (*init_state_contr) (rqci indx, void *data), realtype (*init_lambda) (rqci indx, void *data), void* data);
realtimesolver* free_rtsolver(realtimesolver* solverRT);
rqci calculateSolution(rqci t, realtimesolver *rtsol);
void performNewtonAndShift(realtimesolver *rtsol);
void estimateNewHorizonPoint(realtimesolver *rtsol);
void storeFirstIteration(rqci t, realtimesolver *rtsol);

//void calculateSolution(realtype* vec);

#endif