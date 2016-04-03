#include <util.h>
#include <test_types.h>

realtype init_matlab_vec(rqci indx, void* data)
{
	rqci nhorizon = 0, nintervals = 0, nvec = 0;
	static rqci bInit = 0;
	static realtype* vec;


	if(!bInit)
	{
		test* Test = (test*) data;
		nintervals = Test->nhorizon + 1;
		nvec = nintervals * NVAR;
		vec = load_vec(Test->fvec, nvec);
		bInit = 1;
	}

	return vec[indx];
}

realtype init_matlab_lambda(rqci indx, void* data)
{
	rqci nhorizon = 0, nintervals = 0, nlambda = 0;
	static rqci bInit = 0;
	static realtype* lambda;

	test* Test = (test*) data;
	if(!bInit)
	{
		test* Test = (test*) data;
		nintervals = Test->nhorizon + 1;
		nlambda = nintervals * NLAMBDA;
		lambda = load_vec(Test->flambda, nlambda);
		bInit = 1;
	}
	return lambda[indx];
}