#include <util.h>

QuadConfig* getQuadConfig()
{
//	Iges=(0.0093886,0.0093886,0.018406)
//	g=9.81
//	m=1.022
//	kT=1.5e1
//	kQ=3e-01
//	d=0.22
//	IM=4.4466e-06
	static rqci bInit = 0;
	static QuadConfig* config = 0;
	if(!bInit)
	{
		config = 
			(QuadConfig*)wrp_malloc(1, sizeof(QuadConfig));
		config->g = 9.81;
		config->m = 1.022;
		config->kT = 1.5e1;
		config->kQ = 3e-01;
		config->d = 0.22;
		config->IM = 4.4466e-06;
		config->Iges[0] = 0.0093886;
		config->Iges[1] = 0.0093886;
		config->Iges[2] = 0.018406;
		config->umin = 0.1;
		config->umax = 1.0;
		config->u[0] = 0;
		config->u[1] = 0;
		config->u[2] = 0;
		config->u[3] = 0;
		bInit = 1;
	}
	return config;
}

QuadConfig* getQuadConfigP()
{
//	Iges=(0.0093886,0.0093886,0.018406)
//	g=9.81
//	m=1.022
//	kT=1.5e1
//	kQ=3e-01
//	d=0.22
//	IM=4.4466e-06
	QuadConfig* config = 0;
		config = 
			(QuadConfig*)wrp_malloc(1, sizeof(QuadConfig));
		config->g = 9.81;
		config->m = 1.022;
		config->kT = 1.5e1;
		config->kQ = 3e-01;
		config->d = 0.22;
		config->IM = 4.4466e-06;
		config->Iges[0] = 0.0093886;
		config->Iges[1] = 0.0093886;
		config->Iges[2] = 0.018406;
		config->umin = 0.1;
		config->umax = 1.0;
		config->u[0] = 0;
		config->u[1] = 0;
		config->u[2] = 0;
		config->u[3] = 0;
	return config;
}
realtype* getSteadyPoint()
{
	rqci k =0;
	realtype m, kT, g, u0;
	realtype* steadyPoint = NULL;

	QuadConfig* conf = getQuadConfig();

	
	m = conf->m;
	kT = conf->kT;
	g = conf->g;

	u0 = sqrt(1.0/4.0 * m * g * 1.0/kT);

	steadyPoint = wrp_calloc(NVAR, sizeof(realtype));
	steadyPoint[INDQ] = 1.0;

	steadyPoint[0] = 2;// Noch dynamisch machen
	steadyPoint[2] = 5;

	for(k = 0; k < NCONTR; k++)
		steadyPoint[NSTATE+k] = u0;
	return steadyPoint;
}