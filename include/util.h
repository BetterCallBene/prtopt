#ifndef _UTIL_H
#define _UTIL_H
#include <stdlib.h>
#include <time.h>
#include <sys_types.h>
#include <quad_types.h>
#include <cs_types.h>


void init_gen();
realtype* gen(rqci n);
void gen2(realtype* data, rqci n);

QuadConfig* getQuadConfig();
QuadConfig* getQuadConfigP();

realtype* getSteadyPoint();
realtype* getSteadyPointDyn(realtype *pos);

rqci numDiff_nD(rqci n, realtype* yold, QuadConfig* quad, void* pdOut, rqci (*func)(realtype*, QuadConfig*, realtype*),realtype (*Jac)(void*,rqci, rqci) );
/*Section Normen*/
void normQ(realtype* y);
rqci checkNorm(realtype* ydata);

realtype tic (void) ;
realtype toc (realtype t);

void print_vec(realtype* vec, rqci n, char* Name);
void iprint_vec(rqci* vec, rqci n, char* Name);
realtype *load_vec(FILE* f, csi n);

void *wrp_malloc (rqci n, size_t size);
void *wrp_calloc (rqci n, size_t size);
void *wrp_free (void *p);
void *wrp_realloc (void *p, rqci n, size_t size, rqci *ok);

cs *cs_copy(cs* A);
cs* cs_sub (cs *A, csi starti, csi endi, csi startj, csi endj, csi triplet);
csi cs_print_ext (const cs *A, csi brief, char* Name);
cs *emptyMatrix(rqci m, rqci n, rqci triplet);

realtype init_matlab_vec(rqci indx, void* data);
realtype init_matlab_lambda(rqci indx, void* data);

#endif