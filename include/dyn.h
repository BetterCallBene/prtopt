#ifndef _DYN_H
#define _DYN_H

#include <sys_types.h>
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <quad_types.h>

int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

int Jac(long int N, realtype t,
            N_Vector y, N_Vector fy, DlsMat J, void *user_data,
            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

rqci testJac(realtype* y, QuadConfig* quad);



#endif