#include <lagrange.h>
#include <cs.h>
#include <cost.h>
#include <quad_types.h>
#include <constraints.h>
#include <util.h>
#include <math.h>
#include <test_types.h>


cs* getDiag(rqci n, realtype alpha)
{
	rqci k = 0;
	cs*A; 

	A = cs_spalloc(n, n, n, 1, 0);
	for(k = 0; k < n; k++)
	{
		A->i[k] = k;
		A->p[k] = k;
		A->x[k] = alpha;
	}
	A->p[k] = n;
	A->nz = -1;
	return A;
}

realtype* getLD(rqci k, realtimesolver* solverRT, realtype *eq_constr, cs* eqConD, realtype *ineqh, cs* ineqConD)
{
	rqci kplus = 0, nactivei = 0, nalloc = 0, nhorizon = 0;
	rqci *activeSet, *muSet;

	realtype lambda[NVAR] = {0};
	realtype eqConDvec[NVAR]={0}, ineqConDvec[NVAR]={0};
	realtype *LD, *cD = 0, *RTlambdaNext = 0, *RTlambdaCurrent = 0;
	realtype *y, *vec = 0, *mu_t =0;

	realtype *u0;

	cs* eqConDT = 0, *ineqConDT = 0;
	nhorizon = solverRT->nhorizon;
	kplus = k + 1;

	if((eqConD && !CS_CSC(eqConD)) || (ineqConD && !CS_CSC(ineqConD)))
		return NULL;

	activeSet = solverRT->activeSet;
	muSet = solverRT->muSet;
	vec = solverRT->vec;

	y = GETVAR(k, vec);
	RTlambdaCurrent = GETLAMBDA(k, solverRT->lambda);
	RTlambdaNext = GETLAMBDA(kplus, solverRT->lambda);

	nactivei = checkIfActive(k, vec, activeSet, muSet);
	nalloc = NSTATE + NVAR + nactivei;
	
	LD = wrp_calloc(nalloc, sizeof(realtype));
	
	cD = costD(k, y);
	u0 = GETCONTR(k, vec);
	
//	print_vec(cD, NVAR, "costD");	
	
	if(k < nhorizon)
	{
		//cs_print_ext(eqConD, 0, "eqConD");
		eqConDT = cs_transpose(eqConD, 1);
		//cs_print_ext(eqConDT, 0, "eqConDT");
		cs_gaxpy(eqConDT, RTlambdaNext, eqConDvec);
		//cs_print_ext(eqConDT, 0, "eqConDT");
		
		cs_spfree(eqConDT);

		if(nactivei > 0)
		{
			mu_t = GetMu_t_act(k, activeSet, muSet);
			ineqConDT = cs_transpose(ineqConD, 1);
			
			cs_gaxpy(ineqConDT, mu_t, ineqConDvec);
			cs_spfree(ineqConDT);
			//print_vec(ineqConDvec, NVAR, "ineqConD");
			wrp_free(mu_t);
		}
		
	}

	if (k > 0)
	{
		for(k=0; k < NSTATE; k++)
			lambda[k] = (-1.0) * RTlambdaCurrent[k];
	}

	for(k = 0; k < NSTATE; k++)
		LD[k] = (-1.0)* eq_constr[k];

	for(k = 0; k < NVAR; k++)
		LD[k+NSTATE] = (-1.0)*(cD[k] + lambda[k] + eqConDvec[k] + ineqConDvec[k]);

	//print_vec(eq_constr, NSTATE, "eq_constr");
	//print_vec(lambda, NVAR, "lambda");
	//print_vec(eqConDvec, NVAR, "eqConDvec");
	//print_vec(ineqConDvec, NVAR, "ineqConDvec");

	for(k = 0; k < nactivei && ineqh; k++)
		LD[k+NSTATE+NVAR] = (-1.0) * ineqh[k];
	
	wrp_free(cD);
	return LD;
}

rqci getLDD(rqci k, realtimesolver* solverRT, cs* eqConD, cs* ineqConD)
{
	rqci i = 0, n_mu_i =0, nactivei = 0, nz = 0, nstate = 0, ncontr = 0, col =0, row = 0;
	rqci* activeSet, *muSet;

	realtype alpha  = 0, halpha = 1.0, hbeta = 1.0;
	realtype *y, *vec = 0;
	realtype x;

	riccati* Ric = 0;

	cs *M, *R, *A, *B, *Q, *D, *ineqConDresize = 0;
	cs *diagIAlpha;
	cs *DQalpha = 0, *Qall = 0;
	cs *DRalpha = 0, *Rall = 0;

	cs *MT, *tmp;


	alpha = ALPHA_LDD; 
	nstate = NSTATE;
	ncontr = NCONTR;

	activeSet = solverRT->activeSet;
	muSet = solverRT->muSet;
	vec = solverRT->vec;
	Ric = solverRT->Ric;

	if((eqConD && !CS_CSC(eqConD)) || (ineqConD && !CS_TRIPLET(ineqConD)))
		return (0);

	y = GETVAR(k, vec);
	nactivei = checkIfActive(k, vec, activeSet, muSet);
	
	if(ineqConD)
	{
		nz = ineqConD->nz;
		ineqConDresize = cs_spalloc(nactivei, ncontr, nz, 1, 1);
		for(i = 0; i < nz; i++)
		{
			//ineqConD->p[i] -= nstate;
			row = ineqConD->i[i];
			col = ineqConD->p[i] - nstate;
			x = ineqConD->x[i];
			cs_entry(ineqConDresize, row, col, x);
		}
	}


	A = cs_sub(eqConD, zero_ind(1), zero_ind(NSTATE), zero_ind(1), zero_ind(NSTATE), 0);
	B = cs_sub(eqConD, zero_ind(1), zero_ind(NSTATE), zero_ind(NSTATE+1), zero_ind(NVAR), 0);
	Q = costDDQ(k, y);
	R = costDDR(k, y);
	M = emptyMatrix(NSTATE, NCONTR, 0);
	D = ineqConDresize;
	DQalpha = getDiag(NSTATE, alpha);
	DRalpha = getDiag(NCONTR, alpha);

	Qall = cs_add(Q, DQalpha, halpha, hbeta);
	Rall = cs_add(R, DRalpha, halpha, hbeta);

//	cs_print_ext(Qall, 0, "Q");

	cs_spfree(Q);
	cs_spfree(R);

	cs_spfree(DQalpha);
	cs_spfree(DRalpha);

//	printf("time: %ld\n", k);
/*
	cs_print_ext(A, 0, "A");
	cs_print_ext(B, 0, "B");
	cs_print_ext(Qall, 0, "Q");
	cs_print_ext(Rall, 0, "R");
	cs_print_ext(M, 0, "M");
	cs_print_ext(D, 0, "D");
*/
  	doStep_prepare(k, nactivei, Qall, M, Rall, A, B, D, Ric);  

  	return (1);
}
/*
realtype *lambda; //NLAMBDA
	realtype *vec;
	realtype *mu;
	Quadconfig * conf;
	rqci* activeSet;
	rqci* muSet;
	riccati* Ric;
*/


void init_matrizen(rqci nhorizon, test* Test)
{
	rqci k, nintervals = nhorizon + 1;

	Test->Q = wrp_malloc(nintervals, sizeof(cs*));
	Test->M = wrp_malloc(nhorizon, sizeof(cs*));
	Test->R = wrp_malloc(nhorizon, sizeof(cs*));
	Test->A = wrp_malloc(nhorizon, sizeof(cs*));
	Test->B = wrp_malloc(nhorizon, sizeof(cs*));
	Test->D = wrp_malloc(nhorizon, sizeof(cs*));
	Test->n_mu_i = wrp_malloc(nhorizon, sizeof(rqci));

	for(k = 0; k < nintervals; k++)
	{
		Test->Q[k] = 0;
		if(k < nhorizon)
		{
			Test->M[k] = 0;
			Test->R[k] = 0;
			Test->A[k] = 0;
			Test->B[k] = 0;
			Test->D[k] = 0;
		}
	}
}

void free_matrizen(rqci nhorizon, test* Test)
{
	rqci k, nintervals = nhorizon + 1;
	for(k = 0; k < nintervals; k++)
	{
		cs_spfree(Test->Q[k]);
		if(k < nhorizon)
		{
			cs_spfree(Test->M[k]);
			cs_spfree(Test->R[k]);
			cs_spfree(Test->A[k]);
			cs_spfree(Test->B[k]);
			cs_spfree(Test->D[k]);
		}
	}
	wrp_free(Test->n_mu_i);
	wrp_free(Test->Q);
	wrp_free(Test->M);
	wrp_free(Test->R);
	wrp_free(Test->A);
	wrp_free(Test->B);
	wrp_free(Test->D);
}
void load_matrizen(rqci k, rqci i, rqci nhorizon, test* Test)
{
	char tmp[256];
	const char filename_LDDi_Qi[256] ={"Data/LDDi%ld_Q%ld.dat"};
    const char filename_LDDi_Mi[256] ={"Data/LDDi%ld_M%ld.dat"};
    const char filename_LDDi_Ri[256] ={"Data/LDDi%ld_R%ld.dat"};
    const char filename_LDDi_Ai[256] ={"Data/LDDi%ld_A%ld.dat"};
    const char filename_LDDi_Bi[256] ={"Data/LDDi%ld_B%ld.dat"};
    const char filename_LDDi_Di[256] ={"Data/LDDi%ld_D%ld.dat"};
    const char filename_LDDi_n_mu_i[256] ={"Data/LDDi%ld_n_mu_i%ld.dat"};
    
    rqci n_mu_i =0;
    cs* A = 0, *csc =0;

    memset(tmp, 0, sizeof(tmp));
	sprintf(tmp, filename_LDDi_Qi, k, i+1);
	FILE* ftmp = fopen(tmp, "r");
	printf("i->%ld\n", i);
	printf("Filename: %s\n", tmp);
	Test->Q[i] = cs_load(ftmp);
	Test->Q[i]->m = NSTATE;
	Test->Q[i]->n = NSTATE;
	fclose(ftmp);


	if (i < nhorizon)
	{
		memset(tmp, 0, sizeof(tmp));
		sprintf(tmp, filename_LDDi_n_mu_i, k, i+1);
		FILE* ftmp = fopen(tmp, "r");
		fscanf(ftmp, "%ld\n", &n_mu_i);
		Test->n_mu_i[i] = n_mu_i;
		fclose(ftmp);


		memset(tmp, 0, sizeof(tmp));
	    sprintf(tmp, filename_LDDi_Mi, k, i+1);
	    ftmp = fopen(tmp, "r");
		Test->M[i] = cs_load(ftmp);
		
		if(Test->M[i]->m != NSTATE || Test->M[i]->n != NCONTR)
			cs_entry(Test->M[i], NSTATE-1, NCONTR-1, 0);

		fclose(ftmp);

		memset(tmp, 0, sizeof(tmp));
	    sprintf(tmp, filename_LDDi_Ri, k, i+1);
	    ftmp = fopen(tmp, "r");
		Test->R[i] = cs_load(ftmp);

		//Test->R[i]->m = NCONTR;
		//Test->R[i]->n =	NCONTR;	
		if(Test->R[i]->m != NCONTR || Test->R[i]->n != NCONTR)
			cs_entry(Test->R[i], NCONTR-1, NCONTR-1, 0);
		fclose(ftmp);

	    memset(tmp, 0, sizeof(tmp));
	    sprintf(tmp, filename_LDDi_Ai, k, i+1);
	    ftmp = fopen(tmp, "r");
		Test->A[i] = cs_load(ftmp);
//		Test->A[i]->m =  NSTATE;
//		Test->A[i]->n =  NSTATE;

		if(Test->A[i]->m != NSTATE || Test->A[i]->n != NSTATE)
			cs_entry(Test->A[i], NSTATE-1, NSTATE-1, 0);

		fclose(ftmp);

		memset(tmp, 0, sizeof(tmp));
		sprintf(tmp, filename_LDDi_Bi, k, i+1);
		ftmp = fopen(tmp, "r");
		Test->B[i] = cs_load(ftmp);
//		Test->B[i]->m =  NSTATE;
//		Test->B[i]->n =  NCONTR;

		if(Test->B[i]->m != NSTATE || Test->B[i]->n != NCONTR)
			cs_entry(Test->B[i], NSTATE-1, NCONTR-1, 0);

		fclose(ftmp);

		if(n_mu_i > 0)
		{
			memset(tmp, 0, sizeof(tmp));
			sprintf(tmp, filename_LDDi_Di, k, i+1);
			ftmp = fopen(tmp, "r");
			Test->D[i] = cs_load(ftmp);

			if(Test->D[i]->m != n_mu_i || Test->D[i]->n != NCONTR)
				cs_entry(Test->D[i], n_mu_i-1, NCONTR-1, 0);

		//Test->D[i]->m =  
		//Test->D[i]->n =
			fclose(ftmp);
		}
	}
	
}

void compare(rqci k, rqci nhorizon, realtimesolver* rtsol, cs* Qmat, cs* Mmat, cs* Rmat, cs* Amat, cs* Bmat, cs* Dmat, rqci kactiveMatlab)
{
	riccati *Ric = 0;
	riccati_step* RicStep = 0;

	realtype nrm1 = 0, alpha = 1.0, beta = -1.0;
	rqci kactive = 0;
	cs* R;
	cs *QcscMat, *McscMat, *RcscMat, *AcscMat, *BcscMat, *DcscMat, *Dcsc = 0;
	cs *diffQ=0, *diffM=0, *diffR=0, *diffA=0, *diffB=0, *diffD=0;

	Ric = rtsol->Ric;
	RicStep = Ric->steps[k];

	//cs_print_ext(Qmat, 0, "Qmat");
	QcscMat = cs_compress(Qmat);
	//cs_print_ext(QcscMat, 0, "QcscMat");
	diffQ = cs_add(RicStep->Q, QcscMat, alpha, beta);
	nrm1 = cs_norm(diffQ);
	printf("error norm Q: %f\n", nrm1);
	
	cs_spfree(diffQ);
	cs_spfree(QcscMat);
	kactive = RicStep->naddConstr;
	
	if(k < nhorizon)
	{
		R = RicStep->R;
		McscMat = cs_compress(Mmat);
		RcscMat = cs_compress(Rmat);
		AcscMat = cs_compress(Amat);
		BcscMat = cs_compress(Bmat);
		
		
		if(kactive >0)
		{
			DcscMat = cs_compress(Dmat);
			Dcsc= cs_compress(RicStep->D);
			diffD = cs_add(Dcsc, DcscMat, alpha, beta);

			nrm1 = cs_norm(diffD);
			printf("error norm D: %f\n", nrm1);
			if(nrm1 > 0.01)
				cs_print_ext(diffD, 0, "Diff D");
			

			cs_spfree(DcscMat);
			cs_spfree(Dcsc);
		}
	
		diffM = cs_add(RicStep->M, McscMat, alpha, beta);
		diffR = cs_add(R, RcscMat, alpha, beta);
		diffA = cs_add(RicStep->A, AcscMat, alpha, beta);
		diffB = cs_add(RicStep->B, BcscMat, alpha, beta);
		//diffD = cs_add(Dcsc, DcscMat, alpha, beta);

		

		nrm1 = cs_norm(diffM);
		printf("error norm M: %f\n", nrm1);
		if(nrm1 > 0.0)
			cs_print_ext(diffM, 0, "Diff M");

		nrm1 = cs_norm(diffR);
		printf("error norm R: %f\n", nrm1);
		if(nrm1 > 0.0)
			cs_print_ext(diffR, 0, "Diff R");

		nrm1 = cs_norm(diffA);
		printf("error norm A: %f\n", nrm1);
		if(nrm1 > 0.01)
			cs_print_ext(diffA, 0, "Diff A");

		nrm1 = cs_norm(diffB);
		printf("error norm B: %f\n", nrm1);
		if(nrm1 > 0.9)
			cs_print_ext(diffB, 0, "Diff B");

		printf("nactive Matlab %ld\n", kactiveMatlab);
		printf("nactive %ld\n", kactive);

		
		cs_spfree(diffM);
		cs_spfree(diffR);
		cs_spfree(diffA);
		cs_spfree(diffB);
		cs_spfree(diffD);
		cs_spfree(McscMat);
		cs_spfree(RcscMat);
		cs_spfree(AcscMat);
		cs_spfree(BcscMat);
		
	}
}

void test_lagrange()
{
	rqci m =0, k = 0, fileid = 3, bRes = 0, nalloc = 0, nhorizon = 0, nLDi =0, nintervals = 0, nvec = 0, nlambda =0, nactive = 0;
	rqci *active = 0, *mu = 0, n_mu_i =0;
	cs *Mmat, *Qmat, *Rmat, *Amat, *Bmat, *Dmat; 
	realtype error = 0, mesh_h = 0;
	realtype *vec, *lambda, *ineqh = 0, *h = 0, *LD = 0, *LD_matlab0 = 0, *LD_matlab =0;

	cs *hD = 0, *ineqhD = 0, *hDRet = 0;
	cs *hDcsc =0, *ineqhDcsc = 0;

	realtimesolver* rtsol = 0;
	riccati* Ric =0;
	QuadConfig* conf 
		= getQuadConfig();

	test* Test = wrp_malloc(1, sizeof(test));

	char filename1[256] ={0};
    char filename2[256] ={0};
    char filename3[256] ={0};
    char filename4[256] ={0};

    sprintf(filename1, "Data/config%ld.dat", fileid);
    sprintf(filename2, "Data/LDi%ld.dat", fileid);
    sprintf(filename3, "Data/vec%ld.dat", fileid);
    sprintf(filename4, "Data/lambda%ld.dat", fileid);
	
	Test->fconfig = fopen(filename1, "r"); 
	Test->fLDi = fopen(filename2, "r");
	Test->fvec = fopen(filename3, "r");
	Test->flambda = fopen(filename4, "r"); 

	fscanf (Test->fconfig, "%ld %ld\n", &nhorizon, &nLDi);
	Test->nhorizon = nhorizon;
	Test->nLDi = nLDi;

	printf("horizon: %ld\n", nhorizon);
	printf("nLDi: %ld\n", nLDi);


	printf("Version 0.0.1\n");
	rtsol = initialize_rtsolver(nhorizon,
		init_matlab_vec, init_matlab_lambda, 
		Test
	);
	
	LD_matlab0 = load_vec(Test->fLDi, nLDi);
	LD_matlab = LD_matlab0 + nLDi;
 	vec = rtsol->vec;
 	lambda = rtsol->lambda;
 	active = rtsol->activeSet;
 	mu = rtsol->muSet;
 	nintervals = nhorizon + 1;

 	nvec = nintervals * NVAR;
 	nlambda = nintervals * NLAMBDA;

 	//print_vec(vec, nvec, "vec");
 	//print_vec(lambda, nlambda, "lambda");

 	rtsol->Ric = initialize_ric(nhorizon);
	Ric = rtsol->Ric;
	rtsol->mesh_h = 1.0/nhorizon;
	mesh_h = rtsol->mesh_h;
		
	init_matrizen(nhorizon, Test);

	for(k = nhorizon; k >= 0; k--)
 	{
 		ineqh = 0; ineqhD = 0;
		nactive = checkIfActive(k, vec, active, mu);
		hDRet = get_eq_con_at_t(k, vec, nhorizon, mesh_h, hDRet, &h, &hD, NULL);
		
		get_ineq_con_at_t_act(k, vec, nhorizon, active, &ineqh, &ineqhD);

		hDcsc = hD ? cs_compress(hD) : NULL;
		ineqhDcsc = ineqhD ? cs_compress(ineqhD) : NULL;

		nalloc = NSTATE + NVAR + nactive;
		

		//cs_print_ext(hDcsc, 0, "hD");
		//cs_print_ext(ineqhDcsc, 0, "ineqhD");
		printf("t=%ld\n", k);
		LD = getLD(k, rtsol, h, hDcsc, ineqh, ineqhDcsc);
		bRes = getLDD(k, rtsol, hDcsc, ineqhD);
//		print_vec(LD, nalloc, "LD" );

		Mmat =0; Qmat=0; Rmat=0; Amat=0; Bmat=0; Dmat=0;

		//if(k < nhorizon)
		{
			LD_matlab -= nalloc;
			for(m = 0; m < nalloc; m++)
			{
				error = fabs(LD[m] - LD_matlab[m]);
				if(error > 0.01)
					printf("k: %ld, m: %ld, error: %f, C: %f, Matlab: %f\n", k, m, error, LD[m], LD_matlab[m]);
			}
		}

		load_matrizen(fileid, k, nhorizon, Test);
		Qmat = Test->Q[k];
		n_mu_i = 0;
		if(k < nhorizon)
		{
			Mmat = Test->M[k];
			Rmat = Test->R[k];
			Amat = Test->A[k];
			Bmat = Test->B[k];
			Dmat = Test->D[k];
			n_mu_i = Test->n_mu_i[k];
		}

		compare(k, nhorizon, rtsol, Qmat, Mmat, Rmat, Amat, Bmat, Dmat, n_mu_i);
		/*
		cs_print_ext(Test->Q[k], 0, "Q");
		if(k < nhorizon)
		{
			cs_print_ext(Test->M[k], 0, "M");
			cs_print_ext(Test->R[k], 0, "R");
			cs_print_ext(Test->A[k], 0, "A");
			cs_print_ext(Test->B[k], 0, "B");
			cs_print_ext(Test->D[k], 0, "D");
		}
		*/
		//doStep(k, LD, Ric);
		cs_spfree(hDcsc);
		cs_spfree(ineqhDcsc);
		wrp_free(LD);
		wrp_free(h);
		wrp_free(ineqh);
		cs_spfree(hD);
		cs_spfree(ineqhD);
	}
	free_matrizen(nhorizon, Test);
	wrp_free(LD_matlab0);
 	cs_spfree(hDRet);
 	free_riccati(Ric);
 	rtsol = free_rtsolver(rtsol);

 	//
 	fclose(Test->fconfig);
 	fclose(Test->fLDi);
 	fclose(Test->fvec);
 	fclose(Test->flambda);
 	wrp_free(Test);
}
/*
int main(int argc, char** argv)
{
	test_lagrange();
	return 0;
}
*/

// int main(int argc, char**argv)
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
// 	return 0;
// }
