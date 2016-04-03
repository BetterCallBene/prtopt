#include <util.h>

void init_gen()
{
	srand(time(0));
}
realtype* gen(rqci n)
{
	realtype* gdata;
	gdata = wrp_malloc(n, sizeof(realtype));
	rqci k;
	for(k = 0; k < n; k++)
	{
		gdata[k] = rand() / (realtype)RAND_MAX ;
	}
	return gdata;
}

void gen2(realtype* data, rqci n)
{
	rqci k;
	for(k = 0; k < n; k++)
		data[k] = rand() / (realtype)RAND_MAX ;
	
}