#include <dyn.h>
#include <util.h>
#include <quad_types.h>
#include <math.h>
#include <rtopt_gen.h>

realtype IJthfunc(void* m, rqci i, rqci j)
{
	MatrixType A = (MatrixType) m;
	return IJth(A, i, j);
}


rqci testJac(realtype* y, QuadConfig* quad)
{
	rqci errorCount = 0;
	realtype anaD;
	realtype ydot[NEQ];
	MatrixType pd = NewMat(NEQ, NEQ);

	gen_jac(y, quad, pd);
	errorCount = numDiff_nD(NEQ, y, quad, pd, gen_func, IJthfunc);	
	DestroyMat(pd);
	
	if(errorCount > 0)
	{
		printf("Failed: (testJac) -> Count of error %ld\n", errorCount);
		return (0);
	}
	printf("Success: (testJac)\n");
	return (1);
}

/*
void printOutput(realtype* y, realtype* ydot, MatrixType pd)
{
	rqci i, j;
	FILE * pFile;
	pFile = fopen ("outputC.m","w");
	fprintf(pFile, "y0C = zeros(%ld, 1);\n", NEQ);
	fprintf(pFile, "ydotC = zeros(%ld, 1);\n", NEQ);
	fprintf(pFile, "pdC = zeros(%ld, %ld);\n", NEQ, NEQ);
	fprintf(pFile, "M0C = zeros(%ld, %ld);\n", NSTATE, NSTATE);
	fprintf(pFile, "N0C = zeros(%ld, %ld);\n", NSTATE, NCONTR);

	for(i = 0; i < NEQ; i++)
		fprintf(pFile, "y0C(%ld)=%f;\n", i + 1, y[i]);
	for(i = 0; i < NEQ; i++)
		fprintf(pFile, "ydotC(%ld)=%f;\n", i + 1, ydot[i]);


	for(i = 0; i < NEQ; i++)
	{
		for(j = 0; j < NEQ; j++)
			fprintf(pFile, "pdC(%ld, %ld)=%f;\n", i + 1, j + 1, IJth(pd, i+1, j+1));
	}
	fclose(pFile);
}

*/

/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/*
 * f routine. Compute function f(t,y). 
 */

int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
	realtype* ydata = NV_DATA_S(y);
	realtype* ydotOut = NV_DATA_S(ydot);

	QuadConfig* data = (QuadConfig*) user_data;

/*
	void sys_func(realtype* y,
		  QuadConfig quad,
		  realtype *ydotOut
);
*/	
	gen_func(ydata, data, ydotOut);
	/*
	printf("Print F");
	for(int i = 0; i < NSTATE; i++)
		printf("%f\n", ydotOut[i]);
	printf("\n");
*/
	return(0);
}

int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  	realtype* ydata = NV_DATA_S(y);

  	QuadConfig* data = (QuadConfig*) user_data;
/*
  	void sys_jac(realtype* y,
		 QuadConfig quad,
		 MatrixType pdOut
	);
*/

	gen_jac(ydata,
		data,
		J
	);

	return(0);
}