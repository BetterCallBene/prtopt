#ifndef _QUAD_TYPES_H
#define _QUAD_TYPES_H

#include <sys_types.h>

#define INDR 	0
#define INDQ	3
#define INDV	7
#define INDW	10
#define INDU    13

#define NCONTR  4
#define NSTATE 13
#define NVAR (NCONTR + NSTATE)
#define NQUAD 4
#define NLAMBDA NSTATE
#define NADDCONSTR 8

#define NASIZE (NSTATE * NSTATE)
#define NBSIZE (NSTATE * NCONTR)
#define NJSIZE (NASIZE + NBSIZE)
#define NEQ (NSTATE + NJSIZE)
//#define NEQ NSTATE

#define NJACOBI (NSTATE * NVAR)
#define NSUBHESSE (NVAR * NVAR)
#define NHESSE (NSTATE * NSUBHESSE)

#define IMth(y, i, j) (y[NSTATE + 		   (i-1) + (j-1)*NSTATE])
#define INth(y, i, j) (y[NSTATE + NASIZE + (i-1) + (j-1)*NSTATE])
#define IJsol(y, i, j) (y[NSTATE + (i-1) + (j-1)*NSTATE])

#define GETSTATE(t, vec) (vec + (NVAR * t))
#define GETCONTR(t, vec) (vec + (NVAR * t + NSTATE))
#define GETVAR(t, vec) (vec + (NVAR * t))

typedef struct quadConfig_struct
{
	realtype g;
	realtype m;
	realtype kT;
	realtype kQ;
	realtype d;
	realtype IM;
	realtype Iges[3];
	realtype u[NCONTR];
	realtype umax;
	realtype umin;
} QuadConfig;

#endif