#ifndef __RTOPT_GEN_H
#define __RTOPT_GEN_H

#include <sys_types.h>
#include <cs_types.h>
#include <quad_types.h>
#include <cvode/cvode_dense.h> 

/*
		| a d g |
	M:= | b e h | -> y:= [a b c d e f g h i]
		| c f i |
*/

#ifdef DENSE 
	typedef DlsMat MatrixType;
	#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */
	#define NewMat(n, m) NewDenseMat(n, m)
	#define DestrMat(M) DestroyMat(M)
#else
	typedef realtype** MatrixType;
	#define IJth(A,i,j) A[i-1][j-1]	
	#define NewMat(n, m) newDenseMat(n, m)
	#define DestrMat(M) destroyMat(M)
#endif

rqci gen_func(realtype* y,
		  QuadConfig* quad,
		  realtype *ydotOut
);

rqci gen_jac(realtype* y,
		 QuadConfig* quad,
		 MatrixType pd
);

realtype gen_cost(realtype* y, realtype* campos);
rqci gen_costD(realtype* y, realtype* campos, realtype* costD);
rqci gen_costDD(realtype* y, cs* costDD);
rqci gen_costDDQ(realtype* y, cs* costDDQ);
rqci gen_costDDR(realtype* y, cs* costDDR);
#endif