#ifndef _LAGRANGE_H
#define _LAGRANGE_H
#include <sys_types.h>
#include <cs_types.h>
#include <realtimesolver.h>
#include <riccati.h>

realtype* getLD(rqci k, realtimesolver* solverRT, realtype *eq_constr, cs* eqConD, realtype *ineqh, cs* ineqConD);
rqci getLDD(rqci k, realtimesolver* solverRT, cs* eqConD, cs* ineqConD);

#endif