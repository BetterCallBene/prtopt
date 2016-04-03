
#ifndef _RICCATI_H
#define _RICCATI_H

#include <sys_base.h>
#include <cs_types.h>
#include <quad_types.h>

typedef struct cslu_struct
{
	csn* N;
	css* S;
}cslu;

typedef struct riccati_step_tmp_struct
{
	cslu* z3; //z3
	cs* z4;
	realtype* z5;
}riccati_step_tmp;

typedef struct riccati_step_struct
{
	csi nvar;
	csi naddConstr;
	cs *P;
	cs *Q;
	cs *M;
	cs *R;
	cs *A;
	cs *B;
	cs *D;
	realtype nabla_lambda[NLAMBDA]; //in R^NLAMBDA
	realtype nabla_s_star[NSTATE];
	realtype delta_s[NSTATE]; //R^NLAMBDA
	realtype delta_lambda[NLAMBDA];
	realtype delta_q[NCONTR];
	realtype* delta_mu;
	riccati_step_tmp *tmp;
}riccati_step;


typedef struct riccati_struct
{
	csi nhorizon;
	csi nbytesDeltaMu;
	riccati_step** steps;
}riccati;

riccati* initialize_ric(csi nhorizon);
riccati *free_riccati(riccati* ric);
void doStep_prepare(csi i, csi nactive, cs* Q, cs*M, cs*R, cs* A, cs *B, cs* D, riccati* Ric);
void doStep(csi i, realtype* LD_i, riccati* Ric);
void solveStep(csi i, riccati* Ric);
csi assembleDelta(riccati* Ric, realtype** delta);
realtype * assembleMu(csi k, riccati* Ric, rqci* active);

#endif