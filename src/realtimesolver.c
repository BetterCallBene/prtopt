#include <realtimesolver.h>
#include <lagrange.h>
#include <constraints.h>
#include <cs.h>
#include <dyn.h>
#include <util.h>
#include <test_types.h>
#include <math.h>
#include <cost.h>
#include <multipleshooting.h>
#include <time.h>

#define NPOS 3
#define GETPOS(t, vec) (vec + (NPOS * t))
static realtype start_pos[NPOS] = {2.0, 0.0, 5.0};


realtype* skierCamPos_Short5(rqci t, realtype* cam_pos_cur)
{
	rqci fileid = 5, k = 0, npos = 0, N = 0;
	FILE * fcam_pos_file = 0, *fconfig = 0;
	realtype* pos = 0;
	char strfcam_pos[256] ="";
	char strfconfig[256] ="";
	static rqci bInit = 0;
	static realtype* skier_cam_pos_vec = 0;

	npos = NPOS;

	if(!bInit)
	{
		printf("prtopt: Init skierCamPos_Short5\n");
		sprintf(strfconfig, "Data/config%ld.dat", fileid);
		sprintf(strfcam_pos, "Data/campos%ld.dat", fileid);
    	fcam_pos_file = fopen(strfcam_pos, "r");
    	fconfig = fopen(strfconfig, "r");
    	fscanf(fconfig, "%ld", &N);
    	N *= NPOS;
    	printf("N: %ld\n", N);

    	skier_cam_pos_vec = load_vec(fcam_pos_file, N);
    	fclose(fcam_pos_file);
    	fclose(fconfig);
    	bInit  = 1;
	}
	pos = GETPOS(t, skier_cam_pos_vec);

	for(k = 0; k < npos; k++)
		cam_pos_cur[k] = pos[k];


	return skier_cam_pos_vec;
}

realtype init_with_steadypoint(rqci indx, void* data)
{
	rqci l = 0, t = 0, *nhorizon = 0;
	
	static rqci bInit = 0;
	static realtype* steadyPoint = 0;

	nhorizon = (rqci*) data;

	if(!bInit)
	{
		steadyPoint = getSteadyPointDyn(start_pos);
		bInit = 1;
	}
	l = indx % NVAR; 
	t = indx / NVAR;

	return ((t == *nhorizon && l >=NSTATE) ? 0.0 : steadyPoint[l]);
}
realtype init_lambda_runable(rqci indx, void *data)
{
	rqci l = 0, t = 0, *nhorizon = 0;

	nhorizon = (rqci*) data;

	l = indx % NVAR; 
	t = indx / NVAR;

	return (t == *nhorizon ? 1.0 : (t+1));
}

realtimesolver* free_rtsolver(realtimesolver* solverRT)
{
	rqci nTime = 0, k = 0;
	
	nTime = solverRT->ntime;

	for(k = 0; k < nTime; k++)
		wrp_free(solverRT->res[k]);

	wrp_free(solverRT->res);
	wrp_free(solverRT->conf);
	wrp_free(solverRT->muSet);
	wrp_free(solverRT->activeSet);
	wrp_free(solverRT->lambda);
	wrp_free(solverRT->vec);
	wrp_free(solverRT->data);
	printf("Finished: free_rtsolver\n");
	return wrp_free(solverRT);
}

realtimesolver* initialize_rtsolver(rqci nhorizon, realtype (*init_state_contr) (rqci indx, void *data), realtype (*init_lambda) (rqci indx, void *data), void* data)
{
	
	rqci nTime, l = 0, k = 0, nintervals = 0, 
			naddContr =0, nvec = 0, nlambda = 0, nmu =0,
			nactiveSet = 0, nhorizonvars = 0, *mu = 0;
	realtype *vec = 0, *lambda = 0;

	realtimesolver* solverRT = 0;
	
	nTime = NTIME;
	printf("nhorizon: Init %ld\n", nhorizon);
	nintervals = nhorizon + 1;
	naddContr = NADDCONSTR;

	nhorizonvars = nhorizon * NVAR;
	nvec = nintervals * NVAR;
	nlambda = nintervals * NLAMBDA;
	nmu = nintervals * naddContr;
	nactiveSet = nintervals * naddContr;

	solverRT = wrp_malloc(1, sizeof(realtimesolver));
	solverRT->data = 0;
	solverRT->mesh_h = 1.0;
	solverRT->nhorizon = nhorizon;
	solverRT->ntime = nTime;
	solverRT->res = wrp_malloc(nTime, sizeof(rtsolres*));
	solverRT->vec  = wrp_calloc(nvec, sizeof(realtype));
	solverRT->lambda = wrp_calloc(nlambda, sizeof(realtype));
	solverRT->activeSet = wrp_calloc(nactiveSet, sizeof(rqci));
	solverRT->muSet = wrp_calloc(nmu, sizeof(nmu));
	solverRT->conf = getQuadConfig();

	vec = solverRT->vec;
	lambda = solverRT->lambda;
	mu = solverRT->muSet;

	for(k = 0; k <nTime; k++)
		solverRT->res[k] = wrp_malloc(1, sizeof(rtsolres));

	for(k = 0; k < nvec; k++)
		vec[k] = init_state_contr(k, data);

	for(k = 0; k < nlambda; k++)
		lambda[k] = init_lambda(k, data);

	for(k = 0; k < nmu; k++)
		mu[k] = 1;

	printf("Finished: initialize_rtsolver\n");
	return solverRT;
}

// void calculateSolution(realtype* vec)
// {
	
// 	realtype umin, umax;
// 	rqci k = 0, nactive = 0, nhorizon =0, naddContr = 0, nalloc = 0, bRes = 0;
// 	rqci *active = 0, *mu = 0;
// 	realtype* vec = 0, *ineqh = 0, *h = 0, *LD = 0;
// 	cs *hD = 0, *ineqhD = 0, *hDRet = 0;


// 	cs *hDcsc =0, *ineqhDcsc = 0;
// 	realtimesolver* rtsol =0;
// 	riccati* Ric = 0;
// 	QuadConfig* conf;


// 	rtsol = initialize_rtsolver();
// 	vec = rtsol->vec;
// 	active = rtsol->activeSet;
// 	mu = rtsol->muSet;
// 	conf = rtsol->conf;

// 	nhorizon = NHORIZON;
// 	naddContr = NADDCONSTR;
// 	k = nhorizon;
// 	print_vec(vec, (nhorizon + 1) * NVAR , "vec");

// 	rtsol->Ric = initialize_ric(nhorizon);
// 	Ric = rtsol->Ric;
	
// 	//iprint_vec(active, (nhorizon + 1) * naddContr , "active");
// 	//iprint_vec(mu, (nhorizon + 1) * naddContr , "mu");

// 	for(k = nhorizon; k >= 0; k--)
//  	{
//  		ineqh = 0; ineqhD = 0;
// 		nactive = checkIfActive(k, vec, conf, active, mu);
// 		hDRet = get_eq_con_at_t(k, vec, conf, hDRet, &h, &hD);
		
// 		get_ineq_con_at_t_act(k, vec, conf, active, &ineqh, &ineqhD);

// 		hDcsc = hD ? cs_compress(hD) : NULL;
// 		ineqhDcsc = ineqhD ? cs_compress(ineqhD) : NULL;

// 		nalloc = NSTATE + NVAR + nactive;

// 		//cs_print_ext(hDcsc, 0, "hD");
// 		//cs_print_ext(ineqhDcsc, 0, "ineqhD");
// 		printf("t=%ld\n", k);
// 		LD = getLD(k, rtsol, h, hDcsc, ineqh, ineqhDcsc);
// 		bRes = getLDD(k, rtsol, hDcsc, ineqhD);

// 		if(!LD|| !bRes)
// 		{
// 			printf("Error:\n");
// 			return 0;
// 		}

// 		doStep(k, LD, Ric);
		
// 		cs_spfree(hDcsc);
// 		cs_spfree(ineqhDcsc);
// 		wrp_free(LD);
// 		wrp_free(h);
// 		wrp_free(ineqh);
// 		cs_spfree(hD);
// 		cs_spfree(ineqhD);
// 	}

// 	solveStep(0, Ric);

// 	for(k = 1; k <= nhorizon;k++)
// 	{
// 		solveStep(k, Ric);
// 	}
// 	cs_spfree(hDRet);
// 	free_riccati(Ric);	
// 	free_rtsolver(rtsol);
// }

rqci calculateSolution(rqci t, realtimesolver *rtsol)
{
	
	realtype umin, umax, mesh_h;
	rqci k = 0, nactive = 0, nhorizon =0, naddContr = 0, nalloc = 0, bRes = 0;
	rqci *active = 0, *mu = 0;
	realtype* vec = 0, *ineqh = 0, *h = 0, *LD = 0;
	cs *hD = 0, *ineqhD = 0, *hDRet = 0;

	cs *hDcsc =0, *ineqhDcsc = 0;
	riccati* Ric = 0;
	
	vec = rtsol->vec;
	active = rtsol->activeSet;
	mu = rtsol->muSet;

	mesh_h = rtsol->mesh_h;
	nhorizon = rtsol->nhorizon;
	naddContr = NADDCONSTR;
	k = nhorizon;
	//print_vec(vec, (nhorizon + 1) * NVAR , "vec");

//	rtsol->Ric = initialize_ric(nhorizon);
	Ric = rtsol->Ric;
	realtype* data = rtsol->data;
	
	//iprint_vec(active, (nhorizon + 1) * naddContr , "active");
	//iprint_vec(mu, (nhorizon + 1) * naddContr , "mu");

	for(k = nhorizon; k >= 0; k--)
 	{
 		ineqh = 0; ineqhD = 0;
		nactive = checkIfActive(k, vec, active, mu);
		hDRet = get_eq_con_at_t(k, vec, nhorizon, mesh_h, hDRet, &h, &hD, data);
		
		get_ineq_con_at_t_act(k, vec, nhorizon, active, &ineqh, &ineqhD);

		hDcsc = hD ? cs_compress(hD) : NULL;
		ineqhDcsc = ineqhD ? cs_compress(ineqhD) : NULL;

		nalloc = NSTATE + NVAR + nactive;

		//cs_print_ext(hDcsc, 0, "hD");
		//cs_print_ext(ineqhDcsc, 0, "ineqhD");
	//	printf("t=%ld\n", k);
		LD = getLD(k, t, rtsol, h, hDcsc, ineqh, ineqhDcsc);
		bRes = getLDD(k, t, rtsol, hDcsc, ineqhD);

		if(!LD|| !bRes)
		{
			printf("Error:\n");
			return 0;
		}

		doStep(k, LD, Ric);
		
		cs_spfree(hDcsc);
		cs_spfree(ineqhDcsc);
		wrp_free(LD);
		wrp_free(h);
		wrp_free(ineqh);
		cs_spfree(hD);
		cs_spfree(ineqhD);
	}
	cs_spfree(hDRet);

	solveStep(0, Ric);

	// Hier Signal an Steuerung ->

	storeFirstIteration(t, rtsol);

	for(k = 1; k <= nhorizon;k++)
		solveStep(k, Ric);
	
//	free_riccati(Ric);
	return (1);	
}

void performNewtonAndShift(realtimesolver *rtsol)
{
	rqci k = 0, nhorizon, kminus = 0, j = 0, nintervals = 0, nvec = 0;
	rqci *active, *mu;
	realtype* vec = 0, *lambda = 0, *delta_mu = 0;
	
	riccati* Ric = 0;
	riccati_step* RicStep =0;
	
	realtype *pMinus = 0, *pCur = 0;
	rqci *piMinus= 0, *piCur = 0;

	vec = rtsol->vec;
	lambda = rtsol->lambda;
	Ric = rtsol->Ric;
	active = rtsol->activeSet;
	mu = rtsol->muSet;
	nhorizon = rtsol->nhorizon;
	nintervals = nhorizon + 1;
	nvec = nintervals * NVAR;
	for(k = 1; k < nhorizon; k++)
	{
		RicStep = Ric->steps[k];
		kminus = k - 1; //Bug
		pCur   = GETSTATE(k, vec);
		pMinus = GETSTATE(kminus, vec);
		for(j = 0; j < NSTATE; j++)
			pMinus[j] = pCur[j] + RicStep->delta_s[j];

		pCur   = GETCONTR(k, vec);
		pMinus = GETCONTR(kminus, vec);

		for(j = 0; j < NCONTR; j++)
			pMinus[j] = pCur[j] + RicStep->delta_q[j];

		pCur   = GETLAMBDA(k, lambda);
		pMinus = GETLAMBDA(kminus, lambda);

		for(j = 0; j < NLAMBDA; j++)
			pMinus[j] = pCur[j] + RicStep->delta_lambda[j];

		delta_mu = assembleMu(k, Ric, active);

		piCur   = GETMU(k, mu);
		piMinus = GETMU(kminus, mu);

		for(j = 0; j < NADDCONSTR; j++)
			piMinus[j] = piCur[j] + (rqci)delta_mu[j];
		delta_mu = wrp_free(delta_mu);
	}
	RicStep = Ric->steps[k];
	pCur   = GETSTATE(k, vec);
	pMinus = GETSTATE(kminus, vec);
	for(j = 0; j < NSTATE; j++)
		pMinus[j] = pCur[j] + RicStep->delta_s[j];

	pCur   = GETLAMBDA(k, lambda);
	pMinus = GETLAMBDA(kminus, lambda);

	for(j = 0; j < NLAMBDA; j++)
		pMinus[j] = pCur[j] + RicStep->delta_lambda[j];

}

void estimateNewHorizonPoint(realtimesolver *rtsol)
{
	rqci j= 0, k = rtsol->nhorizon, kminus = 0, kminusminus = 0;
	rqci *mu;
	realtype* vec, *lambda;
	realtype *pMinus = 0, *pCur = 0;
	rqci *piMinus= 0, *piCur = 0;

	kminus = k - 1;
	kminusminus = kminus - 1;
	vec = rtsol->vec;
	lambda = rtsol->lambda;
	mu = rtsol->muSet;

	pCur   = GETSTATE(k, vec);
	pMinus = GETSTATE(kminus, vec);

	for(j = 0; j < NSTATE; j++)
		 pCur[j] = pMinus[j];

	pCur   = GETCONTR(kminus, vec);
	pMinus = GETCONTR(kminusminus, vec);

	for(j = 0; j < NCONTR; j++)
		 pCur[j] = pMinus[j];

	pCur   = GETLAMBDA(k, lambda);
	pMinus = GETLAMBDA(kminus, lambda);

	for(j = 0; j < NLAMBDA; j++)
		 pCur[j] = pMinus[j];

	piCur   = GETMU(k, mu);
	piMinus = GETMU(kminus, mu);

	for(j = 0; j < NADDCONSTR; j++)
		piCur[j] = piMinus[j];

}

void storeFirstIteration(rqci t, realtimesolver *rtsol)
{
	rqci k = 0, nhorizon, kminus = 0, j = 0;
	rqci *active, *mu;
	realtype* vec = 0, *lambda = 0, *delta_mu = 0;
	
	riccati* Ric = 0;
	riccati_step* RicStep =0;
	rtsolres* res = 0;
	
	realtype *pMinus = 0, *pCur = 0;
	rqci *piMinus= 0, *piCur = 0;

	vec = rtsol->vec;
	lambda = rtsol->lambda;
	Ric = rtsol->Ric;
	active = rtsol->activeSet;
	mu = rtsol->muSet;
	nhorizon = rtsol->nhorizon;

	res = rtsol->res[t];

	RicStep = Ric->steps[k];
	kminus = k - 1; //Bug
	pCur   = GETSTATE(k, vec);
	
	for(j = 0; j < NSTATE; j++)
		res->s[j] = pCur[j] + RicStep->delta_s[j];
	

	pCur   = GETCONTR(k, vec);
	
	for(j = 0; j < NCONTR; j++)
		res->q[j] = pCur[j] + RicStep->delta_q[j];

	pCur   = GETLAMBDA(k, lambda);

	for(j = 0; j < NLAMBDA; j++)
		res->lambda[j] = pCur[j] + RicStep->delta_lambda[j];

	delta_mu = assembleMu(k, Ric, active);

	piCur   = GETMU(k, mu);
	
	for(j = 0; j < NADDCONSTR; j++)
		res->mu[j] = piCur[j] + (rqci)delta_mu[j];
	delta_mu = wrp_free(delta_mu);
}

void test_calculation_solution()
{

	realtimesolver* rtsol = 0;
	rqci t = 0, k = 0, fileid = 4, nhorizon = 0, nintervals = 0, ndelta = 0, nvec = 0, ndelta_matlab = 0;
	test* Test = wrp_malloc(1, sizeof(test));

	realtype *vec = 0;
	realtype diff_delta, delta_matlab, delta_rel, relT;


	realtype * delta = 0;
	char filename1[256] ={0};
    char filename2[256] ={0};
    char filename3[256] ={0};
    char filename4[256] ={0};

    sprintf(filename1, "Data/config%ld.dat", fileid);
    sprintf(filename2, "Data/vec%ld.dat", fileid);
    sprintf(filename3, "Data/lambda%ld.dat", fileid);
    sprintf(filename4, "Data/delta_ric%ld.dat", fileid);
	
	Test->fconfig = fopen(filename1, "r"); 
	Test->fvec = fopen(filename2, "r");
	Test->flambda = fopen(filename3, "r"); 
	Test->fdelta_ric = fopen(filename4, "r"); 
	
	fscanf (Test->fconfig, "%ld %ld %ld\n", &nhorizon, &t, &ndelta_matlab);
	Test->nhorizon = nhorizon;
	Test->delta = load_vec(Test->fdelta_ric, ndelta_matlab);

	nintervals = nhorizon + 1;
	nvec = nintervals * NVAR;

	

	rtsol = initialize_rtsolver(nhorizon,
		init_matlab_vec, init_matlab_lambda, 
		Test
	);
/*
	rtsol = initialize_rtsolver(nhorizon,
		init_with_steadypoint, init_lambda_runable, 
		&nhorizon
	);
*/
	print_vec(rtsol->lambda, NLAMBDA, "lambda");
	vec = rtsol->vec;
	printf("t: %ld\n", t);
	rtsol->Ric = initialize_ric(nhorizon);
	relT = tic();
	calculateSolution(t, rtsol);
	relT = toc(relT);
	printf("calc t: %f\n", relT);
 	ndelta = assembleDelta(rtsol->Ric, &delta);
 	printf("ndelta: %ld\n", ndelta);
 	printf("ndelta_matlab: %ld\n", ndelta_matlab);

 	for(k = 0; k < ndelta; k++)
 	{
 		delta_matlab = Test->delta[k];
 		diff_delta = fabs(delta_matlab - delta[k]);
 		delta_rel = fabs(diff_delta / delta_matlab);
 		if(diff_delta > 0.1 && delta_rel > 0.001)
 			printf("[%ld]: Matlab: %f, C: %f, Diff: %f\n", k, delta_matlab, delta[k], diff_delta);
 	}

 	performNewtonAndShift(rtsol);
 	estimateNewHorizonPoint(rtsol);

 	rtsol->Ric = free_riccati(rtsol->Ric);


	rtsol = free_rtsolver(rtsol);

	fclose(Test->fconfig);
 	fclose(Test->fvec);
 	fclose(Test->flambda);
 	fclose(Test->fdelta_ric);
 	wrp_free(Test->delta);
 	wrp_free(delta);
 	wrp_free(Test);
}

void saveData(realtimesolver* rtsol)
{
	time_t t;
    struct tm *ts;
    char str[13]; 
    char filename[30] ="solfolder/prtopt_%s_%s.dat";
    char filename_s[FILENAME_MAX] = "";
    char filename_q[FILENAME_MAX] = "";
    char filename_lambda[FILENAME_MAX] ="";
    char filename_mu[FILENAME_MAX] ="";

    rqci nTime, k = 0, l =0;
    FILE *fstate = 0, *fq=0, *flambda=0, *fmu= 0;

    realtype *s, *q, *lambda, *mu;
    rtsolres **res = 0, *kres = 0;
    t = time(NULL);
    ts = localtime(&t);
    
    strftime(str, 13, "%d%m%y%H%M%S", ts);

    sprintf(filename_s, filename, str, "s");
    sprintf(filename_q, filename, str, "q");
    sprintf(filename_lambda, filename, str, "lambda");
    sprintf(filename_mu, filename, str, "mu");
   
    fstate = fopen(filename_s, "w");
    fq = fopen(filename_q, "w");
    flambda = fopen(filename_lambda, "w");
    fmu = fopen(filename_lambda, "w");

    nTime = rtsol->ntime;

    res = rtsol->res;
    printf("prtopt: save data\n");

    for(k = 0; k < nTime; k++)
    {
    	kres = res[k];
    	for(l = 0; l < NSTATE; l++)
    		fprintf(fstate, "%f\n", kres->s[l]);
    	for(l = 0; l < NCONTR; l++)
    		fprintf(fq, "%f\n", kres->q[l]);
    	for(l = 0; l < NLAMBDA; l++)
    		fprintf(flambda, "%f\n", kres->lambda[l]);
    	for(l = 0; l < NADDCONSTR; l++)
    		fprintf(fmu, "%ld\n", kres->mu[l]);
    }

    fclose(fmu);
    fclose(flambda);
    fclose(fq);
    fclose(fstate);
}

#ifdef WITH_MPI
int rank = 0;
#endif

int main(int argc, char** argv)
{
	rqci nhorizon = 0;
	realtype costValue = 0;
	realtype *vec =0;
	realtimesolver* rtsol = 0;
	realtype *data = 0, *pos_vec = 0;
	rqci t, tend =NTIME;
#ifdef WITH_MPI
	rqci mesh_h = 0;
#endif 
#ifndef WITH_MPI
	if(argc == 2)
	{
		nhorizon = atol(argv[1]);
		printf("nhorizon: %ld\n", nhorizon);
	}
	else
	{
		printf("prtopt: Falsche Anzahl von Eingabeparamtern: %ld\n", argc);
		return 1;
	}
#else
	int numtasks = 0;
	int rc = 0;
	rc = MPI_Init(&argc,&argv);
    if (rc != MPI_SUCCESS) {
    	printf ("Error starting MPI program. Terminating.\n");
     	MPI_Abort(MPI_COMM_WORLD, rc);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if(numtasks < 3)
	{
		printf("numtasks must creater than > 2");
		printf("Error starting MPI program. Terminating.\n");
	}
	nhorizon = (rqci) (numtasks - 1);
#endif
BEGIN_ROOT
	f_cam_pos = skierCamPos_Short5;
	pos_vec = f_cam_pos(0, start_pos);// <-start_pos: global variable
	print_vec(start_pos, NPOS, "prtopt: starting position");
	printf("MPI: Init successed, nhorizon: %ld\n", nhorizon);
 	rtsol = initialize_rtsolver(nhorizon,
		init_with_steadypoint, init_lambda_runable, 
		&nhorizon
	);
#ifdef WITH_MPI
 	mesh_h = rtsol->mesh_h;
#endif
	vec = rtsol->vec;
END_ROOT 	
 	
 	for(t = 0; t < tend; t++)
 	{
#ifdef WITH_MPI
 		data = shoot(nhorizon, mesh_h, vec);
#endif 		
BEGIN_ROOT
		rtsol->data = data;
		
 		rtsol->Ric = initialize_ric(nhorizon);
 		calculateSolution(t, rtsol);

 		if(t > 0)
 		{
 			costValue = costAll(nhorizon, vec);
 			printf("Cost function: %f\n", costValue);
 		}
 		
 		performNewtonAndShift(rtsol);
 		estimateNewHorizonPoint(rtsol);
 		
 		rtsol->Ric = free_riccati(rtsol->Ric);
 		rtsol->data = wrp_free(rtsol->data);
END_ROOT
 	}
BEGIN_ROOT
	saveData(rtsol);
 	rtsol = free_rtsolver(rtsol);
 	wrp_free(pos_vec);
END_ROOT


#ifdef WITH_MPI
 	MPI_Finalize();
BEGIN_ROOT 	
 	printf("MPI: Finalized\n");
END_ROOT
#endif
	return 0;
}