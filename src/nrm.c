#include <util.h>
#include <math.h>
void normQ(realtype* ydata)
{
	realtype norm;
	realtype q1, q2, q3, q4;
	
	q1 =ydata[INDQ];
	q2 =ydata[INDQ + 1];
	q3 =ydata[INDQ + 2];
	q4 =ydata[INDQ + 3];
	norm = sqrt(q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4);

	ydata[INDQ] = q1 / norm;
	ydata[INDQ + 1] = q2 / norm;
	ydata[INDQ + 2] = q3 / norm;
	ydata[INDQ + 3] = q4 / norm;
}

rqci checkNorm(realtype* ydata)
{
	realtype norm;
	realtype q1, q2, q3, q4;

	q1 =ydata[INDQ];
	q2 =ydata[INDQ + 1];
	q3 =ydata[INDQ + 2];
	q4 =ydata[INDQ + 3];
	norm = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

	if(fabs(norm - 1) < 1e-4)
		return (1);//->true
	return (0);//->false
}