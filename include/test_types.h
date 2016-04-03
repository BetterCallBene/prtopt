#ifndef _TEST_TYPES_H
#define _TEST_TYPES_H
	typedef struct test_struct
	{
		rqci nhorizon;
		FILE* fvec;
		FILE* flambda;
		FILE* fLDi;
		FILE* fconfig;
		FILE* fdelta_ric;
		realtype *delta;
		rqci* n_mu_i;
		rqci nLDi;
		cs** Q;
		cs** M;
		cs** R;
		cs** A;
		cs** B;
		cs** D;
	}test;
#endif

	 