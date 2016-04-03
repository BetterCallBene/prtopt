#include <riccati.h>
#include <constraints.h>
#include <cs.h>
#include <util.h>
#include <quad_types.h>

//static realtype tic (void) { return (clock () / (realtype) CLOCKS_PER_SEC) ; }
//static realtype toc (realtype t) { realtype s = tic () ; return (CS_MAX (0, s-t)) ; }


//#define zero_ind(i) (i - 1)

typedef struct test_struct
{
	csi nhorizon;
	csi naddConstr;
	csi dimLDDi;
	csi dimHesse;
	cs *LDD_i;
	cs *hesse_L;
	realtype *LD_i;
	realtype *sol;
}test;

/*
void print_vec2(csi* vec, int n)
{
	csi i;
	printf("Vector length: %ld\n", n);
	for(i = 0; i < n; i++)
		printf("%ld : %ld\n", i, vec[i]);
	printf("\n");	
}
void print_vec(realtype* vec, int n)
{
	csi i;
	printf("Vector length: %ld\n", n);
	for(i = 0; i < n; i++)
		printf("%ld : %f\n", i, vec[i]);
	printf("\n");	
}

void print_vec1(realtype* vec, int n, char* Name)
{
	csi i;
	printf("Name: %s, Vector length: %ld\n", Name, n);
	for(i = 0; i < n; i++)
		printf("%ld : %f\n", i, vec[i]);
	printf("\n");	
}
*/


csi cs_chgPosMatrix(cs* A, csi newm, csi newn, csi starti, csi startj)
{
	csi i, j, k, nz, *Ai, *Aj;
	if(!CS_TRIPLET(A)) return (0);

	Ai = A->i;
	Aj = A->p;
	A->m = newm;
	A->n = newn;
	nz = A->nz;

	for(k = 0; k < nz; k++)
	{
		Ai[k]+=starti; 
		Aj[k]+=startj;
	}
	return (1);
}

cs* cs_addsubmatrix(cs* A, cs* Asub, csi newm, csi newn)
{
	csi* Ai, *Aj, *Asubi, *Asubj;
	csi nzmax, nzA, nzAsub, nz;
	realtype *Ax, *Asubx;
	if(!CS_TRIPLET(A) || !CS_TRIPLET(Asub)) return NULL;

	nzA = A->nz;
	nzAsub = Asub->nz;
	A->m = newm;
	A->n = newn;
	nz = nzA + nzAsub;
	
	A->nz = nz;
	cs_sprealloc (A, 0);

	Ai = A->i; Aj = A->p; Ax = A->x;
	Asubi = Asub->i; Asubj = Asub->p; Asubx = Asub->x;
	
	memcpy(Ai + nzA, Asubi, nzAsub * sizeof(csi));
	memcpy(Aj + nzA, Asubj, nzAsub * sizeof(csi)); 
	memcpy(Ax + nzA, Asubx, nzAsub * sizeof(realtype));

	return A;
}

cs* cs_cootranspose(cs* A)
{
	csi n, m, nzmax;
	cs* coo;
	m = A->m; n = A->n; nzmax = A->nzmax;
	
	if(!CS_TRIPLET(A)) return NULL;
	coo = cs_spalloc(n, m, nzmax, 1, 1);
	
	memcpy(coo->i, A->p, nzmax * sizeof(csi));
	memcpy(coo->p, A->i, nzmax * sizeof(csi));
	memcpy(coo->x, A->x, nzmax * sizeof(realtype));

	coo->nz = A->nz;
	
	cs_sprealloc(coo, 0);
	return coo;
}

csi cs_coosort(cs*A)
{
	csi nz, i, col, colnext, row, rownext;
	csi* Ap, *Ai;
	realtype xtmp;
	realtype *Ax;
	if(!CS_TRIPLET(A)) return (0);

	Ap = A->p; Ai = A->i; Ax = A->x;
	nz = A->nz;
	
	for (nz; nz>1; nz--){
    	for (i=0; i<nz-1; i++){
    		
    		col = Ap[i];
    		colnext = Ap[i+1];
    		row = Ai[i];
    		rownext = Ai[i+1];
    		if(col > colnext || (col == colnext && row > rownext))
    		{
    			xtmp = Ax[i];

    			Ai[i] = rownext;
    			Ai[i+1] = row;

    			Ax[i] = Ax[i+1];
    			Ax[i+1] = xtmp;

    			Ap[i] = colnext;
    			Ap[i+1] = col;
    		}
    	} // ende innere for-Schleife
  	}
  	return (1);
}

cs* cs_csccoo(cs *A)
{
	cs* coo;
	csi i, j, p, n, m, nzmax, *Ap, *Ai;
	realtype *Ax;

	if (!CS_CSC (A)) return (NULL) ;    /* check inputs */
	Ap = A->p ; Ai = A->i ; Ax = A->x ; n = A->n; m = A->m;
	nzmax = A->nzmax;
	coo = cs_spalloc(m, n, nzmax, 1, 1);
	for (j = 0 ; j < n ; j++)
    {
        p = Ap [j] ;                        /* get current location of col j */
        for ( ; p < Ap [j+1] ; p++)
	    {
	        i = Ai [p];
	        if(Ax) cs_entry(coo, i, j, Ax[p]);
	    }
	}
	cs_sprealloc (coo, 0);
	return coo;
}

cs* cs_ipmatrix(const csi *pivo, const cs* B)
{
	cs *PB, *CPB;
	csi p, i, j, n, m, *Bp, *Bi;
	realtype *Bx;
	if (!B || !CS_CSC (B)) return (NULL);

	m = B->m; n = B->n; Bp = B->p ; Bi = B->i ; Bx = B->x;
	PB = cs_spalloc (m, n, m * n, 1, 1) ;
	for (j = 0 ; j < n ; j++)
    {
        p = Bp [j] ;
        for ( ; p < Bp [j+1] ; p++)
        {
        	i = Bi [p];
        	cs_entry (PB, pivo ? pivo[i] : i, j, Bx[p]);
        }
    }
    cs_sprealloc (PB, 0); 
    CPB = cs_compress(PB);
    cs_spfree(PB);
    return CPB;
}

cs* cs_lsolve2 (const cs *L, const cs *B)
{
    csi p, i, j, n, m, *Bp, *Bi;
    cs *T, *CT;
    realtype *Bx ;
    realtype *x;
    if (!CS_CSC (B) || !CS_CSC(L)) return (0) ;                     
    n = B->n ; Bp = B->p ; Bi = B->i ; Bx = B->x ;
    m = B->m;
    x = cs_calloc(m, sizeof(realtype));
    T = cs_spalloc (0, 0, 1, 1, 1); 
    for (j = 0 ; j < n ; j++)
    {
    	p = Bp[j];

        for (;p < Bp[j+1]; p++)
        	x[Bi[p]] = Bx[p];

        cs_lsolve(L, x);
        for(i = 0; i < m; i++)
        {
        	if(x[i] != 0)
        	{
        		cs_entry(T, i, j, x[i]);
        		x[i] = 0;
        	}
        }
    }
    cs_sprealloc (T, 0) ;
    CT = cs_compress(T);
    cs_free(x);
    cs_spfree(T);
    return CT;
}

cs* cs_usolve2 (const cs *U, const cs *Y)
{
    csi p, i, j, n, m, *Yp, *Yi;
    cs *T, *CT;
    realtype *Yx ;
    realtype *x;
    if (!CS_CSC(U) || !CS_CSC (Y)) return (0) ;                     
    n = Y->n ; Yp = Y->p ; Yi = Y->i ; Yx = Y->x ;
    m = Y->m;
    x = cs_calloc(m, sizeof(realtype));
    T = cs_spalloc (0, 0, 1, 1, 1); 

    for (j = 0 ; j < n ; j++)
    {
    	p = Yp[j];

        for (;p < Yp[j+1]; p++)
        	x[Yi[p]] = Yx[p];
        cs_usolve(U, x);

        for(i = 0; i < m; i++)
        {
        	if(x[i] != 0)
        	{
        		cs_entry(T, i, j, x[i]);
        		x[i] = 0;
        	}
        }
    }
    
    cs_sprealloc (T, 0) ;
	CT = cs_compress(T) ;
	cs_spfree(T);
	cs_free(x);
    return CT; 
}

cslu* lu(cs* A)
{
	css* S;
	cslu* lu_struct;
	lu_struct = cs_malloc(1, sizeof(cslu));

	lu_struct->S = S = cs_sqr(0, A, 0);
 	lu_struct->N = cs_lu(A, S, TOL_RICCATI);

 	return lu_struct;
}

cs* lusolM(cslu* SN, cs* B)
{
	csn *N;
	css *S;

	cs *PB, *LS2, *US, *X;
	
	N = SN->N;
	S = SN->S;

	PB = cs_ipmatrix(N->pinv, B);
	LS2 = cs_lsolve2 (N->L, PB);
	cs_spfree(PB);
	
	US = cs_usolve2(N->U, LS2);
	cs_spfree(LS2);
	X = cs_ipmatrix(S->q, US);
	cs_spfree(US);

	return X;
}

cslu* sn_free(cslu* SN)
{
	if(!SN) return NULL;
	cs_nfree(SN->N);
	cs_sfree(SN->S);
	return (cs_free(SN));
}

csi lusol(cslu* SN, realtype* b, csi n)
{
	csi ok;
	realtype *x;
	csn *N;
	css *S; 
	

	if (!SN || !b) return (0) ; 

	N = SN->N;
	S = SN->S;

	x = cs_malloc (n, sizeof (realtype)) ;
	ok = (S && N && x) ;
	if(ok)
	{
		cs_ipvec (N->pinv, b, x, n) ;       /* x = b(p) */
    	cs_lsolve (N->L, x) ;               /* x = L\x */
    	cs_usolve (N->U, x) ;               /* x = U\x */
    	cs_ipvec (S->q, x, b, n);          /* b(q) = x */
	}
    /*rhs*/

	cs_free(x);

   	return(ok);
}

test* free_test(test* t)
{
	if(!t) return NULL;
	cs_spfree(t->hesse_L);
	cs_spfree(t->LDD_i);
	return (cs_free(t));
}

void getTestInfo(test* Test)
{
	csi m, n;
	m = Test->LDD_i->m; n = Test->LDD_i->n;
	printf("LDD_i (%ld, %ld)\n", m, n);
	m = Test->hesse_L->m; n = Test->hesse_L->n;
	printf("hesse_L (%ld, %ld)\n", m, n);
	printf("horizon(%ld), addConstr(%ld)\n", Test->nhorizon, Test->naddConstr);
}

/*
function doStep(o,i, LDD_i, LD_i, n_mu_i )
            o.n_mu{i} = n_mu_i;
            o.Q{i} = LDD_i(1:o.n_state,1:o.n_state);
            o.nabla_lambda{i} = LD_i(1:o.n_lambda);
            o.n_var{i} = o.n_lambda + o.n_state + o.n_contr + n_mu_i;
            
            if(i == o.horizon+1)
                %ldann ist n_mu_i = 0 da dann keine controls vorhanden sind
                o.P{i} = o.Q{i};
                o.nabla_s_star{i} = LD_i(o.n_lambda+1:o.n_lambda + o.n_state);
                
            else
                o.M{i} = LDD_i(1:o.n_state,o.n_state +1 :o.n_state + o.n_contr);
                o.R{i} = LDD_i(o.n_state +1 :o.n_state + o.n_contr, o.n_state+1 :o.n_state + o.n_contr);
                o.A{i} = LDD_i(o.n_state + o.n_contr + n_mu_i + 1: o.n_var{i}, 1:o.n_state);
                o.B{i} = LDD_i(o.n_state + o.n_contr + n_mu_i + 1: o.n_var{i}, o.n_state +1 : o.n_state + o.n_contr);
                o.nabla_q{i} = LD_i(o.n_lambda + o.n_state +1  : o.n_lambda + o.n_state + o.n_contr);
                if(n_mu_i == 0)
                    % Benutze Regel aus Riccati_woConstr
                    
                    
%                     cond = rcond(full( o.R{i} + o.B{i}' * o.P{i+1}  * o.B{i} ));
%                     
%                     if (cond > 1e4 || cond < 1e-4  )
%                         ges = full(o.R{i} + o.B{i}' * o.P{i+1}  * o.B{i}) ;
%                         Ri = full(o.R{i});
%                         Bi = full(o.B{i});
%                         Pip1 = full(o.P{i+1});
%                     end
                                          
                    
                    o.P{i} = o.Q{i} + (o.A{i}' * o.P{i+1} * o.A{i}) - ...
                        (o.M{i} + o.A{i}' * o.P{i+1} * o.B{i}) * ...
                        ((o.R{i} + o.B{i}' * o.P{i+1}  * o.B{i}) \ ...
                        (o.M{i}'  + o.B{i}' * o.P{i+1} * o.A{i}));
                    
                    o.nabla_s_star{i} = LD_i(o.n_lambda +1 : o.n_lambda + o.n_state) + ...
                        o.A{i}' * o.P{i+1} * o.nabla_lambda{i+1} + ...
                        o.A{i}' * o.nabla_s_star{i+1} - ...
                        (o.M{i} + o.A{i}' * o.P{i+1} * o.B{i}) * ...
                        (( o.R{i} + o.B{i}' * o.P{i+1} * o.B{i}) ...
                        \ ( o.nabla_q{i} + o.B{i}' * o.P{i+1} *  o.nabla_lambda{i+1} + o.B{i}' * o.nabla_s_star{i+1}));
                    %TODO: hier kann man o.B{i}' und auch o.A{i}' ausklammern, machen wir
                    %aber erst, wenn der Test geht.
                else
                    % Benutze neue Regel mit Constraints
                    o.D{i} = LDD_i(o.n_state + o.n_contr + 1 : o.n_state + o.n_contr + n_mu_i ,  o.n_state +1 : o.n_state + o.n_contr);
                    o.nabla_mu{i} = LD_i( o.n_lambda + o.n_state + o.n_contr +1 : o.n_lambda + o.n_state + o.n_contr + n_mu_i);
                    
                    o.z4{i} = [o.M{i}'  + o.B{i}' * o.P{i+1} * o.A{i}; zeros(n_mu_i, o.n_state)];
                    o.z3{i} = [ o.R{i} +  o.B{i}' * o.P{i+1} *   o.B{i}, o.D{i}' ;...
                        o.D{i} , zeros(n_mu_i)];
                    z2 = [o.M{i} + o.A{i}' * o.P{i+1} * o.B{i}, zeros(o.n_state, n_mu_i)];
                    z1 =  o.Q{i} + (o.A{i}' * o.P{i+1} * o.A{i});
                    
                    o.P{i} = z1 - z2 * (o.z3{i} \ o.z4{i});
                    
                    o.z5{i} = [ o.nabla_q{i} + o.B{i}' * ( o.P{i+1} * o.nabla_lambda{i+1} + o.nabla_s_star{i+1}) ; ...
                        o.nabla_mu{i}];
                    
                    z6 = LD_i(o.n_lambda +1 : o.n_lambda + o.n_state) + ...
                        o.A{i}' * (o.P{i+1} * o.nabla_lambda{i+1} + o.nabla_s_star{i+1});
                    
                    o.nabla_s_star{i} = z6 - z2 * (o.z3{i} \ o.z5{i});
                    
                end
            end
        end
*/

void doStep(csi i, realtype* LD_i, riccati* Ric)
{
	const realtype alpha = 1.0;
	const realtype beta = 1.0;
	const realtype negbeta = -1.0;


	csi k, nhorizon, n_mu_i, m_z4, nTmp, ok, n1, nzmax, *pwork;
	csn* N = NULL;
	css* S = NULL;
	cs *Plast, *glP, *P, *Q, *M, *MT, *R, *A, *B, *AT, *BT, *D, *DT;
	realtype x[NSTATE] ={0};
	realtype *nabla_s_star, *nabla_q, *nabla_lambda, *last_nabla_s_star, *last_nabla_lambda, *nabla_mu, *rhs;
	realtype tmpVec1[NSTATE]={0}, tmpVec2[NSTATE]={0}, tmpVec3[NCONTR]={0}, tmpVec4[NCONTR]={0},
		tmpVec5[NSTATE] = {0};
	
	cs *tmp, *tmp1, *T1, *T2 = NULL, *T3 = NULL, *T13 = NULL, *ATPlast, *BTPlast, *cooT1;
	cslu* SN = NULL;


	nhorizon = Ric->nhorizon;
	riccati_step* RicStep = Ric->steps[i];
	riccati_step_tmp* RicStepTmp = RicStep->tmp;
	
	//P = RicStep->P;
	Q = RicStep->Q;
	M = RicStep->M;
	R = RicStep->R;
	A = RicStep->A;
	B = RicStep->B;
	D = RicStep->D;
	n_mu_i = RicStep->naddConstr;

	nabla_lambda = RicStep->nabla_lambda;
	nabla_s_star = RicStep->nabla_s_star;

	memcpy(nabla_lambda, 
	 	LD_i,
	 	NLAMBDA * sizeof(realtype)
	);


	memcpy(nabla_s_star, 
	  		LD_i + zero_ind(NLAMBDA + 1), 
	  		(NSTATE) * sizeof(realtype) 
	);

	
	if(i == zero_ind(nhorizon +1))
	{
		//printf("i->%ld\n", i);
  		RicStep->P = cs_copy(RicStep->Q);
  	}
	else
	{
		//printf("Last: i->%ld\n", i+1);
		riccati_step* RicStepLast = Ric->steps[i+1];

 		last_nabla_lambda = RicStepLast->nabla_lambda;
 		last_nabla_s_star = RicStepLast->nabla_s_star;
 		Plast = RicStepLast->P;

 		AT = cs_transpose(A, 1);
 		BT = cs_transpose(B, 1);
 		MT = cs_transpose(M, 1);

 		rhs = cs_calloc(NCONTR + n_mu_i, sizeof(realtype));
 		RicStepTmp->z5 = cs_calloc(NCONTR + n_mu_i, sizeof(realtype));

 		nabla_q = cs_calloc(NCONTR, sizeof(realtype));
 		memcpy(nabla_q, LD_i + zero_ind(NLAMBDA + NSTATE + 1), NCONTR * sizeof(realtype));

 	 	/* T1 = R + BT * PLast * B  */
		BTPlast = cs_multiply(BT, Plast);
		tmp = cs_multiply(BTPlast, B);
		T1 = cs_add(R, tmp, alpha, beta);
		cs_spfree(tmp);


		/*T2 = M + AT * PLast * B  */
 		ATPlast = cs_multiply(AT, Plast);
 		tmp = cs_multiply(ATPlast, B);
 		T2 = cs_add(M, tmp, alpha, beta);
 		cs_spfree(tmp);
 		
 		/*T3 = MT + BT * Plast * A*/
 		tmp = cs_multiply(BTPlast, A);
 		T3 = cs_add(MT, tmp, alpha, beta);
 		cs_spfree(tmp);
 		
 		/* P = Q + (A' + Plast + A) */
		tmp = cs_multiply(AT, Plast);
		tmp1 = cs_multiply(tmp, A);
		glP = cs_add(Q, tmp1, alpha, beta);
		cs_spfree(tmp);
		cs_spfree(tmp1);

		cs_gaxpy (ATPlast, last_nabla_lambda, tmpVec1); //A{i}' * P{i+1} * nabla_lambda{i+1}
 		cs_gaxpy (AT, last_nabla_s_star, tmpVec2);		//A{i}' * nabla_s_star{i+1}
 		cs_gaxpy (BTPlast, last_nabla_lambda, tmpVec3); //B{i}' * P{i+1} *  nabla_lambda{i+1}
 		cs_gaxpy (BT, last_nabla_s_star, tmpVec4);		//B{i}' * nabla_s_star{i+1}

 		for(k = 0; k < NCONTR; k++)
 		{
			rhs[k] = nabla_q[k] + tmpVec3[k] + tmpVec4[k];
			RicStepTmp->z5[k] = rhs[k];
 		}

 		if(n_mu_i == 0)
 		{
 			SN = lu(T1);
 			T13 = lusolM(SN, T3);
 			//o.Q{i} + (o.A{i}' * o.P{i+1} * o.A{i})
 			
			tmp = cs_multiply(T2, T13);
			/* P = P - ... */
			/*Q + (A' + Plast + A)*/
			P = cs_add(glP, tmp, alpha, negbeta);
			//o.nabla_q{i} + o.B{i}' * o.P{i+1} *  
															   //o.nabla_lambda{i+1} + 
															   //o.B{i}' * o.nabla_s_star{i+1}
			lusol(SN, rhs, NCONTR);
			cs_gaxpy(T2, rhs, tmpVec5);

			for(k = 0; k < NSTATE; k++)
 				nabla_s_star[k] = nabla_s_star[k] + tmpVec1[k] + tmpVec2[k] 
 							- tmpVec5[k];
 			
 		}
 		else
 		{
 			m_z4 = NCONTR + n_mu_i;
 			DT = cs_cootranspose(D);

 			nabla_mu = cs_calloc(n_mu_i, sizeof(realtype));
 			memcpy(nabla_mu, LD_i + zero_ind(NLAMBDA + NVAR +1 ), n_mu_i * sizeof(realtype) ); 
 			
 			T3->m += n_mu_i; //z4
 			cooT1 = cs_csccoo(T1);
 			cs_spfree(T1);

 			cs_chgPosMatrix(DT, m_z4, m_z4, 0, zero_ind(NCONTR+1));
 			cs_chgPosMatrix(D, m_z4, m_z4, zero_ind(NCONTR+1), 0);

 			cs_addsubmatrix(cooT1, DT, m_z4, m_z4); //z3
 			cs_addsubmatrix(cooT1, D, m_z4, m_z4); //z3

 			T1 = cs_compress(cooT1);
 			 			
 			nTmp = T2->n + 1;
 			nzmax = T2->nzmax;
			T2->n +=n_mu_i; //z2
 			n1 = T2->n + 1;

 			/* cs_realloc (T2->p, n1, sizeof (csi), &ok) ; -> Speicherfehler*/
			pwork = cs_calloc(n1, sizeof(csi));
			memcpy(pwork, T2->p, sizeof(csi) * nTmp);

			for(k = nTmp; k < n1; k++)
 				pwork[k] = nzmax;

 			cs_free(T2->p);
 			T2->p = NULL;
 			T2->p = cs_malloc(n1, sizeof(csi));
 			memcpy(T2->p, pwork, sizeof(csi) * n1);

 			cs_free(pwork);
			
			//P -> z1
			SN = lu(T1);
			T13 = lusolM(SN, T3); //z3 / z4

			tmp = cs_multiply(T2, T13);
			P = cs_add(glP, tmp, alpha, negbeta); //z1 - z2 * (o.z3{i} \ o.z4{i});

			for(k = 0; k < n_mu_i; k++)
			{
				rhs[k+NCONTR] = nabla_mu[k];
				RicStepTmp->z5[k+NCONTR] = rhs[k+NCONTR];
			}

			lusol(SN, rhs, m_z4);
			cs_gaxpy(T2, rhs, tmpVec5);

			for(k = 0; k < NSTATE; k++)
 				nabla_s_star[k] = nabla_s_star[k] + tmpVec1[k] + tmpVec2[k] 
 							- tmpVec5[k];

 			cs_spfree(cooT1);
 			cs_free(nabla_mu);
 			cs_spfree(DT);
 			
 		}

 		RicStep->P = P;
 		RicStepTmp->z3 = SN;
 		RicStepTmp->z4 = T3;
 		cs_free(rhs);
 		cs_spfree(glP);
 		cs_spfree(T2);
 		cs_spfree(T13);
		cs_spfree(ATPlast);
		cs_spfree(T1);
		cs_spfree(tmp);
		cs_spfree(BTPlast);
		cs_free(nabla_q);
 		cs_spfree(AT);
 		cs_spfree(BT);
 		cs_spfree(MT);

	}
}

void doStep_prepare(csi i, csi nactivei, cs* Q, cs*M, cs*R, cs* A, cs *B, cs* D, riccati* Ric)
{
	csi n_var_all;

	riccati_step* RicStep;
	riccati_step_tmp* RicStepTmp;

	Ric->steps[i] = cs_malloc(1, sizeof(riccati_step));
	RicStep = Ric->steps[i];
	RicStep->tmp = cs_malloc(1, sizeof(riccati_step_tmp));
	RicStepTmp = RicStep->tmp;
	RicStep->naddConstr = nactivei;

	//printf("nactivei: %ld\n", nactivei);
	n_var_all = NLAMBDA + NSTATE + NCONTR + nactivei;
	
	RicStep->nvar = n_var_all;

	RicStep->P = NULL;
	RicStep->delta_mu = NULL;
	RicStep->Q = Q;
	RicStep->M = M;
    RicStep->R = R;
    RicStep->A = A;
    RicStep->B = B;
    RicStep->D = D;//cs_sub(LDD_i, zero_ind(NVAR + 1), zero_ind(NVAR + n_mu_i), zero_ind(NSTATE + 1),  zero_ind(NVAR), 1);

    RicStepTmp->z3 = NULL;
    RicStepTmp->z4 = NULL;
    RicStepTmp->z5 = NULL;
}


void doStep_prepare_test(csi i, cs* LDD_i, realtype* LD_i, csi n_mu_i, riccati* Ric)
{
	csi n_var_all;
	cs *Q, *M, *R, *A, *B, *D;

	n_var_all = NLAMBDA + NSTATE + NCONTR + n_mu_i;

	Q = cs_sub(LDD_i, zero_ind(1), zero_ind(NSTATE), zero_ind(1), zero_ind(NSTATE), 0);
	M = cs_sub(LDD_i, zero_ind(1), zero_ind(NSTATE), zero_ind(NSTATE +1),  zero_ind(NVAR), 0);
    R = cs_sub(LDD_i, zero_ind(NSTATE +1),  zero_ind(NVAR), zero_ind(NSTATE+1), zero_ind(NVAR), 0);
    A = cs_sub(LDD_i, zero_ind(NVAR + n_mu_i + 1), zero_ind(n_var_all), zero_ind(1), zero_ind(NSTATE), 0);
    B = cs_sub(LDD_i, zero_ind(NVAR + n_mu_i + 1), zero_ind(n_var_all), zero_ind(NSTATE + 1),  zero_ind(NVAR), 0);
    D = cs_sub(LDD_i, zero_ind(NVAR + 1), zero_ind(NVAR + n_mu_i), zero_ind(NSTATE + 1),  zero_ind(NVAR), 1);

    doStep_prepare(i, n_mu_i, Q, M, R, A, B, D, Ric);

	/*
	riccati_step* RicStep;
	riccati_step_tmp* RicStepTmp;
	
	
	Ric->steps[i] = cs_malloc(1, sizeof(riccati_step));
	RicStep = Ric->steps[i];
	RicStep->tmp = cs_malloc(1, sizeof(riccati_step_tmp));
	RicStepTmp = RicStep->tmp;
	RicStep->naddConstr = n_mu_i;
	n_var_all = NLAMBDA + NSTATE + NCONTR + n_mu_i;
	
	RicStep->nvar = n_var_all;
	
	RicStep->P = NULL;
	RicStep->delta_mu = NULL;
	RicStep->Q = cs_sub(LDD_i, zero_ind(1), zero_ind(NSTATE), zero_ind(1), zero_ind(NSTATE), 0);
	RicStep->M = cs_sub(LDD_i, zero_ind(1), zero_ind(NSTATE), zero_ind(NSTATE +1),  zero_ind(NVAR), 0);
    RicStep->R = cs_sub(LDD_i, zero_ind(NSTATE +1),  zero_ind(NVAR), zero_ind(NSTATE+1), zero_ind(NVAR), 0);
    RicStep->A = cs_sub(LDD_i, zero_ind(NVAR + n_mu_i + 1), zero_ind(n_var_all), zero_ind(1), zero_ind(NSTATE), 0);
    RicStep->B = cs_sub(LDD_i, zero_ind(NVAR + n_mu_i + 1), zero_ind(n_var_all), zero_ind(NSTATE + 1),  zero_ind(NVAR), 0);
    RicStep->D = cs_sub(LDD_i, zero_ind(NVAR + 1), zero_ind(NVAR + n_mu_i), zero_ind(NSTATE + 1),  zero_ind(NVAR), 1);

    RicStepTmp->z3 = NULL;
    RicStepTmp->z4 = NULL;
    RicStepTmp->z5 = NULL;
    */
    

    //RicStep->nabla_lambda = cs_calloc(NLAMBDA, sizeof(realtype));
    //RicStep->nabla_s_star = cs_calloc(NSTATE + NLAMBDA, sizeof(realtype));
    
    doStep(i, LD_i, Ric);
}

void solveStep(csi i, riccati* Ric)
{
	csi k = 0;
	csi n_mu_i, nhorizon, ncontrmu;
	realtype z1[NSTATE] ={0}, z2[NSTATE]={0};
	realtype Pnabla_lambda[NLAMBDA]={0};
	realtype *nabla_lambda, *nabla_s_star, *delta_s, *delta_lambda, *delta_q, *delta_mu;
	realtype *last_delta_s, *last_delta_q;
	realtype *rhs, *z5;

	cs *P, *Q, *M, *R, *A, *B;
	cs *lastP, *lastA, *lastB, *z4;
	cslu* z3; //z3

	char message[256]={0};
	
	nhorizon = Ric->nhorizon;

	riccati_step* RicStep = Ric->steps[i];
	riccati_step_tmp* RicStepTmp = RicStep->tmp;


	n_mu_i = RicStep->naddConstr;
	nabla_lambda = RicStep->nabla_lambda;
	nabla_s_star = RicStep->nabla_s_star;

	delta_s = RicStep->delta_s;
	delta_lambda = RicStep->delta_lambda;
	delta_q = RicStep->delta_q;
	

	P = RicStep->P;
	Q = RicStep->Q;
	M = RicStep->M;
	R = RicStep->R;
	A = RicStep->A;
	B = RicStep->B;

	z3 = RicStepTmp->z3;
	z4 = RicStepTmp->z4;
	z5 = RicStepTmp->z5;

	ncontrmu = NCONTR + n_mu_i;
	cs_gaxpy(P, nabla_lambda, Pnabla_lambda);

	if(i == zero_ind(1))
	{
		
		for(k = 0; k < NSTATE; k++)
		{
			//o.delta_lambda{i} = -o.P{i}* o.nabla_lambda{i} - o.nabla_s_star{i};
			delta_lambda[k] = -Pnabla_lambda[k] - nabla_s_star[k];
			//o.delta_s{i} = -o.nabla_lambda{i};
			delta_s[k] =  -nabla_lambda[k];
		}
	}
	else
	{
		riccati_step* RicStepLast = Ric->steps[i-1];
		last_delta_s = RicStepLast->delta_s;
		last_delta_q = RicStepLast->delta_q;

		lastP = RicStepLast->P;
		lastA = RicStepLast->A;
		lastB = RicStepLast->B;

		//z1 = (o.A{i-1} * o.delta_s{i-1} + o.B{i-1} * o.delta_q{i-1});
		cs_gaxpy(lastA, last_delta_s, z1);
		/* y = A*x+y */
		cs_gaxpy(lastB, last_delta_q, z1);
		cs_gaxpy(P, z1, z2);

		for(k = 0; k < NSTATE; k++)
		{
			//z2 = P{i} * z1
			//o.delta_lambda{i} = -o.P{i}* o.nabla_lambda{i} - o.nabla_s_star{i} + o.P{i} * z1 ;
			delta_lambda[k] = -Pnabla_lambda[k] - nabla_s_star[k] + z2[k];
			//o.delta_s{i} = - o.nabla_lambda{i}  + z1;
			delta_s[k] =  -nabla_lambda[k] + z1[k];
		}
	}

	sprintf(message, "delta_lambda[%ld]", i);
	//print_vec1(delta_lambda, NSTATE, message);
	memset(message, 0, 256);
	sprintf(message, "delta_s[%ld]", i);
	//print_vec1(delta_s, NSTATE, message);

	if( i != zero_ind(nhorizon+1))
	{
		Ric->nbytesDeltaMu +=n_mu_i; 
		RicStep->delta_mu = cs_calloc(n_mu_i, sizeof(realtype));

		delta_mu = RicStep->delta_mu;
		rhs = cs_calloc(ncontrmu, sizeof(realtype));

		cs_gaxpy(z4, delta_s, rhs);
		
		for(k = 0; k < ncontrmu; k++)
    		rhs[k] = z5[k] - rhs[k];

    	lusol(z3, rhs, ncontrmu);
    	
    	for(k = 0; k < NCONTR; k++)
    		delta_q[k] = rhs[k];
    	for(k = 0; k < n_mu_i; k++)
    		delta_mu[k] = rhs[NCONTR + k];

    	//print_vec1(delta_q, NCONTR, "delta_q");
    	//print_vec1(delta_mu, n_mu_i, "delta_mu");

    	cs_free(rhs);
	}
}
/*
function solveStep(o,i)
            if ( i ==1 )
                o.delta_lambda{i} = -o.P{i}* o.nabla_lambda{i} - o.nabla_s_star{i};
                o.delta_s{i} = -o.nabla_lambda{i};
            else
                z1 = (o.A{i-1} * o.delta_s{i-1} + o.B{i-1} * o.delta_q{i-1});
                
                % solve delta lambda
                o.delta_lambda{i} = -o.P{i}* o.nabla_lambda{i} - o.nabla_s_star{i} + o.P{i} * z1 ;
                % solve delta s
                o.delta_s{i} = - o.nabla_lambda{i}  + z1;
            end
            
            %solve for delta_q and delta_mu
            if( i ~= o.horizon +1)
                if ( o.n_mu{i} == 0)
%                     cond = rcond(full(o.R{i} + o.B{i}' * o.P{i+1} * o.B{i}));
%                     
%                     rhs = full( ...
%                         o.nabla_q{i} + ...
%                         o.B{i}' * o.P{i+1} * o.nabla_lambda{i+1} + ...
%                         o.B{i}' * o.nabla_s_star{i+1} - ...
%                         ( o.M{i}' + o.B{i}' * o.P{i+1} * o.A{i} ) * o.delta_s{i} );
%                     
%                     if( cond < 1e-4 || cond > 1e4)
%                         Bi = o.B{i};
%                         Pip1 = o.P{i+1};
%                         nlip1 = o.nabla_lambda{i+1};
%                         nsip1 = o.nabla_s_star{i+1};
%                         Mi = o.M{i};
%                         Ai = o.A{i};
%                         dsi = o.delta_s{i};
%                     end
                    
                    o.delta_q{i} = (o.R{i} + o.B{i}' * o.P{i+1} * o.B{i}) \ ...
                        ( ...
                        o.nabla_q{i} + ...
                        o.B{i}' * o.P{i+1} * o.nabla_lambda{i+1} + ...
                        o.B{i}' * o.nabla_s_star{i+1} - ...
                        ( o.M{i}' + o.B{i}' * o.P{i+1} * o.A{i} ) * o.delta_s{i} );
                    
                    % We to reset delta_mu, such that the previous solution
                    % is overwriten.
                    o.delta_mu{i} = [];
                    
                else
                    z6 = o.z3{i} \ (o.z5{i} - (o.z4{i} * o.delta_s{i}));
                    
                    o.delta_q{i} = z6(1:o.n_contr);
                    o.delta_mu{i} = z6(o.n_contr +1 : o.n_contr + o.n_mu{i});
                    
                end
            end
*/
riccati_step_tmp* free_riccati_step_tmp(riccati_step_tmp* RicStepTmp)
{
	if(!RicStepTmp) return NULL;
	sn_free(RicStepTmp->z3);
	cs_spfree(RicStepTmp->z4);
	cs_free(RicStepTmp->z5);
	return (cs_free(RicStepTmp));
}
riccati_step* free_riccati_step(riccati_step* RicStep)
{
	if(!RicStep) return NULL;
	free_riccati_step_tmp(RicStep->tmp);
//	cs_free(RicStep->nabla_q);
//	cs_free(RicStep->delta_lambda);
	cs_free(RicStep->delta_mu);
//	cs_free(RicStep->delta_s);
//	cs_free(RicStep->delta_q);
//	cs_free(RicStep->nabla_lambda);
//	cs_free(RicStep->nabla_s_star);
	cs_spfree(RicStep->P);
	cs_spfree(RicStep->D);
	cs_spfree(RicStep->B);
	cs_spfree(RicStep->A);
	cs_spfree(RicStep->R);
	cs_spfree(RicStep->M);
	cs_spfree(RicStep->Q);
	return (cs_free(RicStep));
}

riccati *free_riccati(riccati* ric)
{
	csi k, nhorizon;
	nhorizon = ric->nhorizon;

	//printf("nhorizon: %ld\n", nhorizon);

	for(k = 0;k < nhorizon + 1;k++)
		free_riccati_step(ric->steps[k]);

	cs_free(ric->steps);
	
	return cs_free(ric);
}

riccati* initialize_ric(csi nhorizon)
{
	csi k;
	riccati* Ric;
	Ric = cs_malloc(1, sizeof(riccati));
	Ric->steps = NULL;
	Ric->nhorizon = nhorizon;
	Ric->nbytesDeltaMu = 0;
	
	(Ric)->steps = cs_malloc(nhorizon+1, sizeof(riccati_step*));
	for(k = 0;k < nhorizon + 1;k++)
		Ric->steps[k] = NULL;

	return Ric;
}

realtype * assembleMu(csi k, riccati* Ric, rqci* active)
{
	rqci i = 0, j = 0, naddConstr = 0, active_t_k = 0;
	rqci *active_t = 0;
	realtype *delta_mu = 0;
	riccati_step *RicStep = 0;
	naddConstr = NADDCONSTR;

	RicStep = Ric->steps[k];
	if(!RicStep)
		return 0;

	active_t = GETACTIVE(k, active);
	delta_mu = wrp_calloc(naddConstr, sizeof(realtype));

	//for()

	for(i = 0; i < naddConstr; i++)
	{
		active_t_k = active_t[i];
		if(active_t_k)
			delta_mu[i] = RicStep->delta_mu[j++];

	}
	return delta_mu;
}
csi assembleDelta(riccati* Ric, realtype** delta)
{
	csi i;
	csi j;
	csi k = 0;
	csi n_mu_i;
	csi nbytes;
	csi nhorizon;
	csi nbytesDeltaMu = 0;
	
	realtype* delta_s; //R^NLAMBDA
	realtype* delta_lambda;
	realtype* delta_q;
	realtype* delta_mu;
	
	riccati_step* RicStep;

	nhorizon = Ric->nhorizon;
	nbytesDeltaMu = Ric->nbytesDeltaMu;

	nbytes = (NLAMBDA + NSTATE) * (nhorizon + 1) + (NCONTR) * nhorizon + nbytesDeltaMu;

	*delta = cs_calloc(nbytes, sizeof(realtype));

	for(i = 0; i < nhorizon; i++)
	{
		RicStep = Ric->steps[i];

		n_mu_i = RicStep->naddConstr;
		delta_lambda = RicStep->delta_lambda;
		delta_s = RicStep->delta_s;
		delta_q = RicStep->delta_q;
		delta_mu = RicStep->delta_mu;

		for(j = 0; j < NLAMBDA; j++)
			(*delta)[k++] = delta_lambda[j];

		for(j = 0; j < NSTATE; j++)
			(*delta)[k++] = delta_s[j];


		for(j = 0; j < NCONTR; j++)
			(*delta)[k++] = delta_q[j];

		for(j = 0; j < n_mu_i; j++)
			(*delta)[k++] = delta_mu[j];

	}
	//Don't forget horizon +1

	RicStep = Ric->steps[i];

	n_mu_i = RicStep->naddConstr;
	delta_lambda = RicStep->delta_lambda;
	delta_s = RicStep->delta_s;
	delta_q = RicStep->delta_q;
	delta_mu = RicStep->delta_mu;

	for(j = 0; j < NLAMBDA; j++)
		(*delta)[k++] = delta_lambda[j];

	for(j = 0; j < NSTATE; j++)
		(*delta)[k++] = delta_s[j];

	return nbytes;
}

void doRiccati(test* Test)
{
	csi ndelta;
	csi indLD_i;
	csi i = 0;
	csi nhorizon, naddConstr;
	realtype nrm = 0;
	realtype *LD_i, *delta, *sol; 

	cs* LDD_i; 

	riccati* Ric =NULL;
	riccati_step* RicStep;


	getTestInfo(Test);

	nhorizon = Test->nhorizon;
	naddConstr = Test->naddConstr;
	sol = Test->sol;
	LD_i = Test->LD_i;
	LDD_i = Test->LDD_i;

	Ric = initialize_ric(nhorizon);
	
	for(i = zero_ind(nhorizon+1); i >= zero_ind(1); i--)
	{
		indLD_i = zero_ind(i * (30+naddConstr) +1);
		printf("%ld\n",indLD_i);
		doStep_prepare_test(i, LDD_i, LD_i+indLD_i, naddConstr, Ric);
	}

	for(i = 0; i <= zero_ind(nhorizon+1); i++)
		solveStep(i, Ric);

	ndelta = assembleDelta(Ric, &delta);

	for(i = 0; i < ndelta; i++)
		nrm+=(delta[i] - sol[i])*(delta[i] - sol[i]);
	printf("Error: %f\n", nrm);

	cs_free(delta);
	free_riccati(Ric);
}




test* load_test_data(csi fileid)
{
	FILE *fLDD_i =NULL, *fhesse_L= NULL, *fB = NULL, *fConfig =NULL, *fsol = NULL;
	cs *T;
	csi m, n;

	test *Test ;
    Test = cs_calloc (1, sizeof (test)) ;

    char filename1[256] ={0};
    char filename2[256] ={0};
    char filename3[256] ={0};
    char filename4[256] ={0};
    char filename5[256] ={0};

    sprintf(filename1, "Data/config%ld.dat", fileid);
    sprintf(filename2, "Data/LDDi%ld.dat", fileid);
    sprintf(filename3, "Data/hesse_L%ld.dat", fileid);
    sprintf(filename4, "Data/grad_L%ld.dat", fileid);
    sprintf(filename5, "Data/sol%ld.dat", fileid);


    fConfig = fopen(filename1, "r");
	fLDD_i = fopen(filename2, "r");
	fhesse_L = fopen(filename3, "r");
	fB = fopen(filename4, "r");
	fsol = fopen(filename5, "r");

	if(!fConfig || !fLDD_i || !fhesse_L || !fB || !fsol )
		return NULL;

	fscanf (fConfig, "%ld %ld %ld %ld\n", &Test->nhorizon, &Test->naddConstr, &Test->dimLDDi, &Test->dimHesse);

	T = cs_load(fLDD_i);
	T->m = Test->dimLDDi;
	T->n = Test->dimLDDi;
	Test->LDD_i = cs_compress (T);
	cs_spfree (T) ;
	
	T = cs_load(fhesse_L);
	T->m = Test->dimHesse;
	T->n = Test->dimHesse;
	Test->hesse_L = cs_compress (T);
	cs_spfree (T) ; 

	m = Test->hesse_L->m;
	Test->LD_i = load_vec(fB, m);
	Test->sol = load_vec(fsol, m);

	fclose(fsol);
	fclose(fConfig);
	fclose(fLDD_i);
	fclose(fhesse_L);
	fclose(fB);

	return Test;
}

void free_test_data(test* data)
{
	cs_spfree(data->hesse_L);
	cs_spfree(data->LDD_i);
	cs_free(data->LD_i);
	cs_free(data->sol);
	cs_free(data);
}



cs* test_matrix()
{
	cs *T;
	T = cs_spalloc (0, 0, 1, 1, 1);
	cs_entry (T, 0, 0, 1.0);
	cs_entry (T, 0, 1, -1.0);
	cs_entry (T, 0, 2, -3.15);
	cs_entry (T, 1, 0, -2.0);
	cs_entry (T, 1, 1, 5.0);
	cs_entry (T, 1, 2, 0.0);
	cs_entry (T, 1, 3, 0.0);
	cs_entry (T, 2, 2, 4.0);
	cs_entry (T, 2, 3, 6.0);
	cs_entry (T, 2, 4, 4.0);
	cs_entry (T, 3, 0, -4.0);
	cs_entry (T, 3, 2, 2.0);
	cs_entry (T, 3, 3, 7.0);
	cs_entry (T, 4, 1, 8.0);
	cs_entry (T, 4, 4, -5.0);
	return (T);
}

cs* test_matrix1()
{
	cs *T;
	T = cs_spalloc (0, 0, 1, 1, 1);
	cs_entry (T, 0, 0, 1);
	cs_entry (T, 0, 1, 2);
	cs_entry (T, 0, 2, 3);
	cs_entry (T, 1, 0, 1);
	cs_entry (T, 1, 1, 1);
	cs_entry (T, 1, 2, 1);
	cs_entry (T, 2, 0, 3);
	cs_entry (T, 2, 1, 3);
	cs_entry (T, 2, 2, 1);
	
	return (T);
}

cs* test_matrixX()
{
	cs* X;
	X = cs_spalloc (0, 0, 1, 1, 1);

	cs_entry (X, 0, 0, 1);
	cs_entry (X, 0, 1, 6);
	cs_entry (X, 0, 2, -5);
	cs_entry (X, 0, 3, 1);
	cs_entry (X, 1, 0, 7);
	cs_entry (X, 1, 1, 5);
	cs_entry (X, 1, 2, 4);
	cs_entry (X, 1, 3, 0);
	cs_entry (X, 2, 0, -1);
	cs_entry (X, 2, 1, 9);
	cs_entry (X, 2, 2, 3);
	cs_entry (X, 2, 3, 0);

	return (X);
}


void test_lu()
{
	realtype t = 0, tout = 0;
	cslu *SN;
	
	cs *T, *A, *B, *X, *Test;

	T = test_matrix1();
	A = cs_compress(T);
	cs_spfree(T);
	
	
	T = test_matrixX();
	X = cs_compress(T);
	cs_spfree(T);

	B = cs_multiply(A, X);
	cs_spfree(X);
	SN = lu(A);

	X = lusolM(SN, B);
	
	sn_free(SN);
	cs_spfree(B);
	cs_spfree(X);
	cs_spfree(A);
}

void test_csccoo()
{
	cs *T, *A, *cooT1;
	T = test_matrix1();
	A = cs_compress(T);
	cs_spfree(T);
	cooT1 = cs_csccoo(A);
 	
 	cs_spfree(cooT1);
 	cs_spfree(A);
}
	

void doTestHorizon1()
{
	test* Test = load_test_data(1);
	if(Test)
		doRiccati(Test);
	else
		printf("Error in doTestHorizon1\n");
	free_test_data(Test);
}

void doTestHorizon2()
{
	test* Test = load_test_data(2);
	if(Test)
		doRiccati(Test);
	else
		printf("Error in doTestHorizon2\n");
	free_test_data(Test);
}
void doTestHorizon200()
{
	test* Test = load_test_data(200);
	if(Test)
		doRiccati(Test);
	else
		printf("Error in doTestHorizon200\n");
	free_test_data(Test);
}
/*
int main(int argv, char** argc)
{

	printf("Version 0.0.4\n");
	//test_csccoo();
	//test_lu();
	//doTestHorizon1();
	doTestHorizon200();

	
	return 0;
}
*/