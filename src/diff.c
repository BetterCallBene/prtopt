#include <util.h>
#include <math.h>


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