#ifndef _CONSTRAINTS_H
#define _CONSTRAINTS_H

#include <sys_types.h>
#include <cs_types.h>
#include <quad_types.h>


#define GETACTIVE(t, active) (active + (NADDCONSTR * t))
#define GETMU(t, mu) (mu + (NADDCONSTR * t))

rqci checkIfActive(rqci t, realtype* vec, rqci* active, rqci* mu);
cs* get_eq_con_at_t(rqci t, realtype* vec, rqci nhorizon, realtype mesh_h, cs* hDLast, realtype** h, cs** hD, realtype *data);
void get_ineq_con_at_t_act(rqci t, realtype* vec, rqci nhorizon, rqci* active, realtype** ineqh, cs** ineqhD);
realtype* GetMu_t_act(rqci t, rqci* active, rqci* mu);
void get_data_at_t(rqci k, realtype* data, realtype *F, realtype *J);


#endif