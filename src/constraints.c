#include <constraints.h>
#include <cs.h>
#include <quad_types.h>
#include <solver.h>
#include <util.h>
//#define DEB_CONSTR

realtype* wind(realtype t, realtype* state)
{
	realtype* res;
	res = wrp_calloc(NSTATE, sizeof(realtype));
	memcpy(res, state, NSTATE * sizeof(realtype));
	return res;
}

realtype* GetMu_t_act(rqci t, rqci* active, rqci* mu)
{
	rqci i = 0, k=0;
	rqci nactive = 0, naddConstr, active_t_i = 0;
	rqci* mu_t, *active_t;
	realtype *mu_t_act;
	
	mu_t = GETMU(t, mu);

	naddConstr = NADDCONSTR;

	active_t = GETACTIVE(t, active);

	for(i = 0; i < naddConstr; i++)
		nactive += active_t[i];

	mu_t_act = wrp_calloc(nactive, sizeof(realtype));

	for(i = 0; i < naddConstr; i++)
	{
		active_t_i = active_t[i];
		if(active_t_i)
			mu_t_act[k++] = (realtype)mu_t[i];
	}

	return mu_t_act;
}


realtype* lockup_table_ineq_con(realtype* y)
{
	realtype umax, umin;
	realtype *table;

	QuadConfig* conf = getQuadConfig();
	umin = conf->umin;
	umax = conf->umax;

	table = wrp_calloc(NADDCONSTR, sizeof(realtype));
	table[0] = y[NSTATE] - umax;
	table[1] = y[NSTATE+1] - umax;
	table[2] = y[NSTATE+2] - umax;
	table[3] = y[NSTATE+3] - umax;

	table[4] = umin - y[NSTATE];
	table[5] = umin - y[NSTATE+1];
	table[6] = umin - y[NSTATE+2];
	table[7] = umin - y[NSTATE+3];
	return table;
}

cs* lockup_table_ineqD_con()
{
	cs* A = 0;
	cs* T = 0;

	T = cs_spalloc(NADDCONSTR, NVAR, NADDCONSTR, 1, 1);

	cs_entry(T, 0, NSTATE, 1);
	cs_entry(T, NCONTR, NSTATE, -1);

	cs_entry(T, 1, NSTATE + 1,  1);
	cs_entry(T, NCONTR + 1, NSTATE + 1,  -1);

	cs_entry(T, 2, NSTATE + 2,  1);
	cs_entry(T, NCONTR + 2, NSTATE + 2,  -1);

	cs_entry(T, 3, NSTATE + 3,  1);
	cs_entry(T, NCONTR + 3, NSTATE + 3,  -1);


	A = cs_compress(T);
	cs_spfree(T);

	return A;
}

rqci checkIfActive(rqci t, realtype* vec, rqci* active, rqci* mu)
{
	rqci k = 0, naddConstr, active_t_k, nactive = 0, *mu_t =0, *active_t = 0;
	realtype *y, *table;


	naddConstr = NADDCONSTR;

	y = GETSTATE(t, vec);
	mu_t = GETMU(t, mu);
	active_t = GETACTIVE(t, active);


	table = lockup_table_ineq_con(y);

	for(k=0; k < naddConstr;k++)
	{	
		active_t_k = (table[k] >= 0) & (mu[k] > 0);

		if(!active_t_k) mu[k] = 1;
		else nactive++;

		active_t[k] = active_t_k;
	}

	wrp_free(table);
	return nactive;
}
/*
void get_ineq_con_at_t(rqci t, realtype* vec, QuadConfig* conf, realtype** ineqh, cs** ineqhD)
{
	rqci nactive;
	rqci active_t[NADDCONSTR];
	nactive = NADDCONSTR;

	for(k = 0; k < nactive; k++)
		active_t[k] = 1;
	get_ineq_con_at_t_act(0, nactive, vec, conf, active_t, ineqh, ineqhD);
}
*/
void get_ineq_con_at_t_act(rqci t, realtype* vec, rqci nhorizon, rqci* active, realtype** ineqh, cs** ineqhD)
{
	rqci i = 0, j = 0, k = 0, naddConstr = 0, active_t_i = 0, bIneq = 0, nz, Ap, nactive =0, *active_t = 0;

	realtype Ax;
	realtype umax;
	realtype umin;
	realtype *y;

	realtype *table = 0;
	cs* tableD = 0, *Asub = 0;
	
	naddConstr = NADDCONSTR;

	y = GETSTATE(t, vec);
	active_t = GETACTIVE(t, active);

	for(i = 0; i < NADDCONSTR; i++)
		nactive += active_t[i];

	if(!nactive)
	{
		ineqh = 0;
		ineqhD =0;
		return;
	}

	(*ineqh) = wrp_calloc(nactive, sizeof(realtype));
	table = lockup_table_ineq_con(y);

	if(t == nhorizon)
	{
		//(*ineqhD) = NULL; //cs_spalloc(nactive, NVAR, 0, 1, 0);
		*ineqhD = NULL;
		
		for(i = 0; i < naddConstr; i++)
		{
			active_t_i = active_t[i];
			if(active_t_i)
				(*ineqh)[k++] = table[i];
		}
	}
	else
	{
		*ineqhD=cs_spalloc(nactive, NVAR, nactive * NVAR, 1, 1);
		tableD = lockup_table_ineqD_con();
		for(i = 0; i < naddConstr; i++)
		{
			active_t_i = active_t[i];
			if(active_t_i)
			{
				(*ineqh)[k] = table[i];
				Asub = cs_sub(tableD, i, i, zero_ind(1), zero_ind(NVAR), 1);
				nz = Asub->nz;
				for(j = 0; j < nz; j++)
				{
					Ap = Asub->p[j];
					Ax = Asub->x[j];
					cs_entry(*ineqhD, k++, Ap, Ax);
				}
				cs_spfree(Asub);
			}
		}
		cs_sprealloc(*ineqhD, 0);
		cs_spfree(tableD);
		
	}
	wrp_free(table);
}

void get_data_at_t(rqci k, realtype* data, realtype *F, realtype *J)
{
	realtype *vec_t = 0;
	int ind = 0;
	ind = (k * NEQ);
	vec_t = data + ind;
	
	if(F)
		memcpy(F, vec_t, NSTATE * sizeof(realtype));
		
	if(J)
		memcpy(J, vec_t + NSTATE, NJSIZE * sizeof(realtype));
	
}


cs* get_eq_con_at_t(rqci t, realtype* vec, rqci nhorizon, realtype mesh_h, cs* hDLast, realtype** h, cs** hD, realtype* data)
{
	rqci tminus;
	rqci k, i, j, n, m, nz;
	
	realtype *y0 = 0, *yI =0, *unext= 0, *ynext=0, *w= 0;
	cs* hDret = 0;
	QuadConfig* conf; 

	conf = getQuadConfig();

	y0 = GETSTATE(t, vec);
	
	n = NSTATE;
	m = NVAR;
	nz = n * m;
	*h = wrp_calloc(NSTATE, sizeof(realtype));
#ifdef DEB_CONSTR
	printf("mesh_h: %f\n", mesh_h);
#endif	
	if(t > 0)
	{
		tminus = t - 1;
		hDret = cs_spalloc(n, m, nz, 1, 1);
		//ynext = GETSTATE(t-1, vec); ->Fuck you
		//unext = GETCONTR(t-1, vec);
		ynext = GETSTATE(tminus, vec);
		unext = GETCONTR(tminus, vec);

		conf->u[0] = unext[0];
		conf->u[1] = unext[1];
		conf->u[2] = unext[2];
		conf->u[3] = unext[3];
		#ifdef DEB_CONSTR
		print_vec(unext, NCONTR, "unext");
		print_vec(ynext, NSTATE, "ynext");
		#endif
#ifdef WITH_MPI
		get_data_at_t(tminus, data, *h, hDret->x);
#else
		yI = helperCreateInitialConditions(ynext);
		integrate1(mesh_h, yI, conf, *h, hDret->x);
#endif
		for(k = 0; k < NJSIZE; k++)
		{
			i = k % NSTATE;
			j = k / NSTATE;
			hDret->i[k] = i;
			hDret->p[k] = j;
		}

		hDret->nz = nz;
		
		wrp_free(yI);
	}

	if(t == nhorizon)
		*hD = 0;//cs_spalloc(n, m, 0, 1, 1);
	else
	{
		*hD = cs_copy(hDLast);
		cs_spfree(hDLast);
	}

	if(t == 0)
	{
		w = wind(t, y0);
		for(k = 0; k < NSTATE;k++)
			(*h)[k] = w[k] - y0[k];
		wrp_free(w);
	}
	else
	{
		for(k = 0; k < NSTATE;k++)
			(*h)[k] = (*h)[k] - y0[k];
	}

	return hDret;
}


// int main(int argc, char** argv)
// {
// 	rqci kminus, nactive = 0;
// 	rqci n, k, l = 0, nelm = 0;
// 	realtype naddConstr, nhorizon;
// 	realtype *vec = NULL, *ineqh = 0, *h, *steadyPoint = NULL;
// 	QuadConfig* conf;
// 	cs *ineqhDT, *ineqhD = 0, *hD = NULL, *hDRet = NULL, *hDRetTmp = NULL, *tmp;
// 	realtype ineqhDvec[NVAR] ={0};
// 	rqci mu[(NHORIZON+1) * NADDCONSTR] ={1};
// 	rqci active[(NHORIZON+1) * NADDCONSTR] = {0};
// 	realtype* mu_t;
// 	rqci* active_t;

// 	for(k = 0; k <(NHORIZON+1) * NADDCONSTR; k++ )
// 		mu[k] = 1;
	
// 	nhorizon = NHORIZON;
//  	n = nhorizon +1;

//  	nelm = nhorizon * NVAR;

	

// 	//init_gen();
// 	//vec = gen(NVAR * n);
//  	conf = getQuadConfig();
// 	naddConstr = NADDCONSTR;

// 	steadyPoint = (realtype*) getSteadyPoint(conf);
// 	normQ(steadyPoint);
// 	steadyPoint[NVAR-4] = 2;
// 	steadyPoint[NVAR-1] = 0.0;
	
// 	vec =  wrp_calloc(n * NVAR, sizeof(realtype));

// 	steadyPoint[0] = 2;
// 	steadyPoint[2] = 5;

// 	for(k = 0; k < nelm ; k++)
// 	{
// 		l = k % NVAR;
// 		vec[k] = steadyPoint[l]; 
// 	}
// 	steadyPoint[2] = 10;
// 	for(k;k < nelm + NVAR; k++)
// 	{
// 		l = k % NVAR;
// 		vec[k] = steadyPoint[l];
// 	}
	
	

// 	for(k = nhorizon; k >= 0; k--)
// 	{
		
// 		printf("t:%ld\n", k);
// 		printf("\nInEquality->Constraints\n");
// 		nactive = checkIfActive(k, vec, conf, active, mu);
// 		active_t = GETACTIVE(k, active);
// 		iprint_vec(active_t, naddConstr, "active_t");
// 		printf("nactive: %ld\n", nactive);
// 		if(nactive > 0)
// 		{
// 			if(k < nhorizon)
// 			{
// 				memset(ineqhDvec, 0, sizeof(realtype) *NVAR );
// 				get_ineq_con_at_t_act(k, vec, conf, active, &ineqh, &ineqhD);
// 				mu_t = GetMu_t_act(k, active, mu);

// 				tmp = cs_compress(ineqhD);
// 				ineqhDT = cs_transpose(tmp, 1);
// 				//
// 				cs_print_ext(ineqhDT, 0, "ineqhDT");
// 				print_vec(mu_t, nactive, "mu_t");
// 				cs_gaxpy(ineqhDT, mu_t, ineqhDvec);
// 				print_vec(ineqhDvec, NVAR, "ineqhDvec");
// 				cs_spfree(tmp);
// 				cs_spfree(ineqhDT);
// 				cs_spfree(ineqhD);
// 				wrp_free(ineqh);
// 				wrp_free(mu_t);
// 			}	
// 		}
// /*
// 		printf("\nEquality->Constraints\n");
// 		hDRet = get_eq_con_at_t(k, vec, conf, hDRet, &h, &hD);
// 		print_vec(h, NSTATE, "h");
// 		cs_print(hD, 0);
// 		cs_spfree(hD);
// 		wrp_free(h);
// */
// 	}
		
// 	free(hDRet);



// 	wrp_free(steadyPoint);
// 	wrp_free(conf);
// 	wrp_free(vec);
// 	return 0;
// }
