#ifndef _COST_H
#define _COST_H
#include <sys_types.h>
#include <cs_types.h>
extern realtype* (*f_cam_pos) (rqci t, realtype*);

realtype costAll(rqci t, realtype* y);
realtype cost(rqci t, realtype* y);
realtype* costD(rqci t, realtype* y);
cs* costDD(rqci t, realtype* y);
cs* costDDQ(rqci t, realtype* y);
cs* costDDR(rqci t, realtype* y);

#endif
