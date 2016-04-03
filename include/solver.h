#ifndef _SOLVER_H
#define _SOLVER_H
#include <sys_types.h>
#include <quad_types.h>

/* Problem Constants */

realtype* helperCreateInitialConditions(realtype* y0);
rqci numDiff_nDvec(realtype* yold, QuadConfig* data);
rqci integrate(realtype mesh_h, const realtype* y0, QuadConfig* data, realtype* yOut);
rqci integrate1(realtype mesh_h, realtype* y0, QuadConfig* data, realtype* F, realtype* J);
rqci integrate2(realtype mesh_h, realtype* y0, QuadConfig* data, realtype* F, realtype* A, realtype*B);
#endif