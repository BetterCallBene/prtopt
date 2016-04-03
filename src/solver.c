#include <solver.h>
#include <dyn.h>
#include <util.h>
#include <quad_types.h>
//Workaround
#include <cvode_impl.h>
//End Workaround
#include <cvode/cvode.h> 
//#define DEB_SOLV

#define lfree     (cv_mem->cv_lfree)
#define IJsol2(y, i, j) (y[(i-1) + (j-1)*NSTATE])

/*
void PrintOutput(realtype* y)
{
	if(checkNorm(y))
		printf("Ok. Output norm is correct,\n");
	else
		printf("Failed: Output norm is not correct.\n");
	
	return;
}
*/

realtype* helperCreateInitialConditions(realtype* y0)
{
	rqci k;
	realtype* ydata;
	ydata = wrp_calloc(NEQ, sizeof(realtype));
	
	for(k = 0; k < NSTATE; k++)
		ydata[k] = y0[k];

	for(k = 1; k <= NSTATE; k++)
		IMth(ydata, k, k) = 1.0;
	
	return ydata;
}


rqci check_flag(void *flagvalue, char *funcname, rqci opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(0); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(0); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(0); }

  return(1);
}

void PrintFinalStats(void *cvode_mem, rqci t)
{
  long int nst = 0, nfe= 0, nsetups=0, nje = 0, nfeLS = 0, nni = 0, ncfn=0, netf=0, nge=0;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
  check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
  check_flag(&flag, "CVDlsGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals(cvode_mem, &nge);
  check_flag(&flag, "CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("Number Steps: %-6ld\n", nst);
  //printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
//	 nst, nfe, nsetups, nfeLS, nje);
 // printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
//	 nni, ncfn, netf, nge);
}


rqci integrate(realtype mesh_h, const realtype* y0, QuadConfig* data, realtype* yOut)
{
	rqci flag, flagr, iout, bNrm, mxsteps =50000;
	realtype reltol, abstol, t, t0, tout, tstop = 0;
	realtype *ptmp;
	N_Vector y, yfnormdot;
	CVodeMem cv_mem;
	void *cvode_mem;
	
  	y = NULL;
	cvode_mem = NULL;

	abstol = ATOL;
	reltol = RTOL;
	t0 = T0;
	iout = 0;  tout = mesh_h;
	#ifdef DEB_SOLV
	tstop = tic();
	#endif
	/* Create serial vector of length NEQ for I.C. and abstol */

	memcpy(yOut, y0, sizeof(realtype) * NEQ); 

	#ifdef DEB_SOLV
	print_vec(yOut, NSTATE, "yOut");
	#endif

	y=  N_VMake_Serial(NEQ, yOut);
	
	if (!check_flag((void *)y, "N_VMake_Serial", 0)) return(0);

	/* Call CVodeCreate to create the solver memory and specify the 
	 * Backward Differentiation Formula and the use of a Newton iteration */
	//cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	cvode_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
	if (!check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(0);

	/* Set the pointer to user-defined data */
	flag = CVodeSetUserData(cvode_mem, data);
	if(!check_flag(&flag, "CVodeSetUserData", 1)) return(0);
  
	/* Call CVodeInit to initialize the integrator memory and specify the
   	 * user's right hand side function in y'=f(t,y), the inital time T0, and
     * the initial dependent variable vector y. */
	flag = CVodeInit(cvode_mem, f, t0, y);
	if (!check_flag(&flag, "CVodeInit", 1)) return(0);

	//Call CVodeSVtolerances to specify the scalar relative tolerance
	// * and vector absolute tolerances 
	flag = CVodeSStolerances(cvode_mem, reltol, abstol);
	if (!check_flag(&flag, "CVodeSStolerances", 1)) return(0);

	flag = CVodeSetMaxNumSteps(cvode_mem, mxsteps);
	if (!check_flag(&flag, "CVodeSetMaxNumSteps", 1)) return(0);

	/* Call CVodeRootInit to specify the root function g with 2 components */
//	flag = CVodeRootInit(cvode_mem, 2, g);
//	if (check_flag(&flag, "CVodeRootInit", 1)) return(1);

	/* Call CVDense to specify the CVDENSE dense linear solver */
	flag = CVDense(cvode_mem, NEQ);
	if (!check_flag(&flag, "CVDense", 1)) return(0);

	/* Set the Jacobian routine to Jac (user-supplied) */
	flag = CVDlsSetDenseJacFn(cvode_mem, Jac);
	if (!check_flag(&flag, "CVDlsSetDenseJacFn", 1)) return(0);

	/* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */
  	//printf(" \n3-species kinetics problem\n\n");
  	
    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    if (!check_flag(&flag, "CVode", 1)) return(0);
    #ifdef DEB_SOLV
    tstop = toc(tstop);
    //printf("tstop %lf\n", tstop); 
    //PrintFinalStats(cvode_mem, t);
    print_vec(yOut, NSTATE, "solution solver1");
    #endif 
/* Free integrator memory */
    //PrintFinalStats(cvode_mem);

//Workaround
    cv_mem = (CVodeMem) (cvode_mem);
    if(lfree != NULL)
    	lfree(cv_mem);

	CVodeFree(&cvode_mem);
	N_VDestroy_Serial(y);
	
/* Print some final statistics */
	return(1);
}

rqci integrate1(realtype mesh_h, realtype* y0, QuadConfig* data, realtype* F, realtype* J)
{
	rqci flag = 0;
	realtype ysol[NEQ]={0};
	flag = integrate(mesh_h, y0, data, ysol);
	if(F)
	{
		memcpy(F, ysol			, NSTATE * sizeof(realtype));
		#ifdef DEB_SOLV
		print_vec(F, NSTATE, "solution solver2");
		#endif
	}
	if(J)
		memcpy(J, ysol + NSTATE	, NJSIZE * sizeof(realtype));
	return flag;
}

rqci integrate2(realtype mesh_h, realtype* y0, QuadConfig* data, realtype* F, realtype* A, realtype*B)
{
	rqci flag = 0;
	realtype ysol[NEQ]={0};
	flag = integrate(mesh_h, y0, data, ysol);
	if(F)
	{
		memcpy(F, ysol			, NSTATE * sizeof(realtype));
		#ifdef DEB_SOLV
		print_vec(F, NSTATE, "solution solver2");
		#endif
	}
	if(A)
		memcpy(A, ysol + NSTATE	, NASIZE * sizeof(realtype));
	if(B)
		memcpy(B, ysol + NSTATE	+ NASIZE, NBSIZE * sizeof(realtype));
	return flag; 
}
/*
rqci numDiff_nD(rqci n, realtype* yold, QuadConfig* quad, void* pdOut, rqci (*func)(realtype*, QuadConfig*, realtype*),realtype (*Jac)(void*,rqci, rqci) )
{
	rqci ind = 0;
	rqci i, j;
	realtype anaD, numD, error;
	realtype ywork[n], pfunc[n], nfunc[n];
	//Copy
	for(i = 0; i < n; i++)
		ywork[i] = yold[i];

	for(i = 0; i < n;i++)
	{
		
		ywork[i] = yold[i] + EPS;
		func(ywork, quad, pfunc);
		ywork[i] = yold[i] - EPS;
		func(ywork, quad, nfunc);
		ywork[i] = yold[i];

		
		for(j = 0; j < n; j++)
		{
			anaD = Jac(pdOut, j+1, i+1);
			numD = (pfunc[j] - nfunc[j])/2.0/EPS;
			error = fabs(anaD - numD);
			
			if(error > 1e-3)
			{
				printf("Ana: %f Num: %f, Error: (%ld, %ld) -> %f\n", anaD, numD, i+1, j+1, error);
				ind++;
			}
		}
	}
	return ind;
}
*/

rqci numDiff_nDvec(realtype* yold, QuadConfig* data)
{
	rqci i, j, ind= 0;
	realtype mesh_h = 1.0;
	realtype* u;
	realtype anaD, numD, error;
	realtype Jac[NEQ]={0};
	realtype* ywork, uold[NCONTR], pfunc[NEQ], nfunc[NEQ];
/*
	for(i = 0; i < n; i++)
		ywork[i] = yold[i];
*/
	u = data->u;
	memcpy(uold, u, NCONTR * sizeof(realtype));
	ywork = helperCreateInitialConditions(yold);

	integrate(mesh_h, ywork, data, Jac);

	for(i = 0; i < NSTATE; i++)
	{
		ywork[i] = yold[i] + EPS;
		integrate1(mesh_h, ywork, data, pfunc, NULL);
		ywork[i] = yold[i] - EPS;
		integrate1(mesh_h, ywork, data, nfunc, NULL);
		ywork[i] = yold[i];

		for(j = 0; j < NSTATE; j++)
		{
			anaD = IMth(Jac, j+1, i+1);

			numD = (pfunc[j] - nfunc[j])/2.0/EPS;
			error = fabs(anaD - numD);
			
			if(error > 1e-2)
			{
				printf("Ana: %f Num: %f, Error: (%ld, %ld) -> %f\n", anaD, numD, i+1, j+1, error);
				ind++;
			}
		}
	}

	
	
	for(i = 0; i < NCONTR; i++)
	{
		
		data->u[i] = uold[i] + EPS;
		integrate1(mesh_h, ywork, data, pfunc, NULL);
		data->u[i] = uold[i] - EPS;
		integrate1(mesh_h, ywork, data, nfunc, NULL);
		data->u[i] = uold[i];

		for(j = 0; j < NSTATE; j++)
		{
			anaD = INth(Jac, j+1, i+1);

			numD = (pfunc[j] - nfunc[j])/2.0/EPS;
			error = fabs(anaD - numD);
			error = error/fabs(anaD);
			
			if(error > 9e-1)
			{
				printf("NCONTR: Ana: %f Num: %f, Error: (%ld, %ld) -> %f\n", anaD, numD, i+1, j+1, error);
				ind++;
			}
		}

	}
	wrp_free(ywork);
}
/*
int main(int argc, char** argv)
{
	
	rqci bNrm;
	//realtype* test;
	//realtype* ysol = wrp_malloc(NEQ, sizeof(realtype));
	QuadConfig* config = getQuadConfig();
	
	init_gen();
	realtype* testdata = gen(NEQ);
	normQ(testdata);
	gen2(config->u, NCONTR);

	numDiff_nDvec(testdata, config);

	/*
	config->u[0] = 1000;
	config->u[1] = 1000;
	config->u[2] = 1000;
	config->u[3] = 1000;
	*/
	//print_vec(config->u, NCONTR, "u");
	
	//testJac(testdata, config);
	//realtype* sol = wrp_malloc(NEQ, sizeof(realtype));
	

	//dyn_func(testdata, config, sol);
	//testJac(testdata, config);
	//print_vec(testdata, NSTATE, "testdata");
	//normQ(testdata);
	//bNrm = checkNorm(testdata);
	//if(bNrm)
	//	printf("Norm: Check\n");
	//integrate(testdata, config, ysol);
	//testJsol(testdata, config);
	

	//print_vec(ysol, NSTATE, "ysol");

	//wrp_free(sol);
	//wrp_free(ysol);
//	wrp_free(config);
//	wrp_free(testdata);
//	return (0);
//}
