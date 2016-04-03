#include <sys_base.h>
#include <rtopt_gen.h>
#include <quad_types.h>
#include <cs.h>
#include <util.h>

//realtype cost(realtype* y, QuadConfig* quad, realtype* campos);

realtype cost(rqci k, realtype* y)
{
	/*ToDo: */
	realtype campos[3];
	campos[0] = 2;
	campos[1] = 0;
	campos[2] = 5;
	return gen_cost(y, campos);
}

realtype costAll(rqci nhorizon, realtype* vec)
{
	rqci k = 0, ninterval = 0;
	ninterval = nhorizon + 1;
	realtype costValue = 0;
	realtype *y =0;

	for(k = 0; k < ninterval; k++)
	{
		y = GETVAR(k, vec);
		costValue += cost(k, y);
	}
	return costValue;
}

realtype* costD(rqci k, realtype* y)
{
	realtype *ret;
	realtype campos[3];
	campos[0] = 2;
	campos[1] = 0;
	campos[2] = 5;
	ret = wrp_calloc(NVAR, sizeof(realtype));
	gen_costD(y, campos, ret);
	return ret;
}

cs* costDD(rqci k, realtype* y)
{
	cs* T, *ret;
	T = cs_spalloc (NVAR, NVAR, NVAR, 1, 1) ;
	gen_costDD(y, T);
	ret = cs_compress(T);
	cs_spfree(T);
	return ret;
}

cs* costDDQ(rqci t, realtype* y)
{
	cs* T, *ret;
	T = cs_spalloc (NSTATE, NSTATE, NSTATE, 1, 1) ;
	gen_costDDQ(y, T);
	ret = cs_compress(T);
	cs_spfree(T);
	return ret;
}

cs* costDDR(rqci t, realtype* y)
{
	cs* T, *ret;
	T = cs_spalloc (NCONTR, NCONTR, NCONTR, 1, 1) ;
	gen_costDDR(y, T);
	ret = cs_compress(T);
	cs_spfree(T);
	return ret;
}
/*
int main(int argc, char* argv)
{
	QuadConfig* conf;
	realtype c;
	realtype* cD, *y;
	cs* cDD;

	init_gen();
	y = gen(NVAR);

	c = cost(y);
	cD = costD(y);
	cDD = costDD(y);

	printf("Cost %f\n", c);
	print_vec(cD, NVAR, "CostD");
	cs_print(cDD, 0);

	cs_spfree(cDD);
	wrp_free(cD);
	wrp_free(y);
	
	return 0;
}
*/