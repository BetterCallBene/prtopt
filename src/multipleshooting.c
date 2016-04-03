#include <sys_types.h>
#include <quad_types.h>
#include <solver.h>
#include <util.h>


realtype* shoot(rqci nhorizon, realtype mesh_h, realtype *vec)
{
	int sendcnt = 0;
	int recvcnt = 0;
	
	rqci nintervall = 0, Nsol = 0;
	realtype yrecv[NVAR], ydot[NEQ];
	realtype *data = 0, *state = 0, *contr = 0, *stateadv=0;
	
	QuadConfig *config;

	nintervall = nhorizon + 1;
	sendcnt = NVAR;
	recvcnt = sendcnt;

	if(rank == ROOT)
		config = getQuadConfig();
	else
		config = getQuadConfigP();
	
	memset(yrecv, 0, sizeof(yrecv));
	memset(ydot, 0, sizeof(ydot));


BEGIN_ROOT
	Nsol = nintervall * NEQ;
	data = wrp_malloc(Nsol, sizeof(realtype));
END_ROOT
	MPI_Bcast(&mesh_h, 1, MPI_REALTYPE, ROOT, MPI_COMM_WORLD);
	MPI_Scatter (vec, sendcnt, MPI_REALTYPE, yrecv, 
			recvcnt,MPI_REALTYPE, ROOT, MPI_COMM_WORLD);

	state = GETSTATE(0, yrecv);
	contr = GETCONTR(0, yrecv);

	config->u[0] = contr[0];
	config->u[1] = contr[1];
	config->u[2] = contr[2];
	config->u[3] = contr[3];


	stateadv = helperCreateInitialConditions(state);

	integrate(mesh_h, stateadv, config, ydot);

	MPI_Gather (ydot, NEQ, MPI_REALTYPE, data, 
		NEQ,MPI_REALTYPE,ROOT,MPI_COMM_WORLD);
	
	wrp_free(stateadv);

	if(rank != ROOT)
		wrp_free(config);
	return data;
}
/*
int main(int argc, char *argv[])  {
	int numtasks, N, rank, root =0, sendcnt = 0, recvcnt = 0, rc = 0, Nsol =0;
	int i = 0;
	char message[256] = {0};
	rqci k = 0, m = 0;
	realtype mesh_h = 1.0;
	realtype *vec = 0, *y0adv = 0, *ydotall;

	realtype y0[NSTATE];
	realtype u0[NCONTR];
	realtype ydot[NEQ]; 

	QuadConfig* quadConfig;
	MPI_Status Stat;

	rc = MPI_Init(&argc,&argv);
    if (rc != MPI_SUCCESS) {
     printf ("Error starting MPI program. Terminating.\n");
     MPI_Abort(MPI_COMM_WORLD, rc);
    }
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	memset(u0, 0, sizeof(realtype) * NCONTR );

	quadConfig = getQuadConfig();

//Initblock auslagern
	if(rank == root)
	{
		N = numtasks * NSTATE;
		Nsol = numtasks * NEQ;
		ydotall = wrp_malloc(Nsol, sizeof(realtype));
		init_gen();
		vec = gen(N);

		for(i = 0; k < numtasks; k++)
			normQ(vec + k * NSTATE);
	
		gen2(u0, NCONTR);
		for(k = 0; k < numtasks; k++)
		{
			if(checkNorm(vec + k * NSTATE))
				printf("1->YES\n");
			else
				printf("1->NO\n");
		}
		
	}

	sendcnt = NSTATE;
	recvcnt = sendcnt;
	MPI_Scatter (vec, sendcnt, MPI_REALTYPE, &y0, 
			recvcnt,MPI_REALTYPE,root, MPI_COMM_WORLD);

	MPI_Bcast(u0, NCONTR, MPI_REALTYPE, root, MPI_COMM_WORLD);

	quadConfig->u[0] = u0[0];
	quadConfig->u[1] = u0[1];
	quadConfig->u[2] = u0[2];
	quadConfig->u[3] = u0[3];


	y0adv = helperCreateInitialConditions(y0);

	integrate(mesh_h, y0adv, quadConfig, ydot);
	
	wrp_free(y0adv);
	wrp_free(quadConfig);


	MPI_Gather (ydot, NEQ, MPI_REALTYPE, ydotall, 
		NEQ,MPI_REALTYPE,root,MPI_COMM_WORLD);

	if(rank == root)
	{
		print_vec(ydotall, Nsol, "ydotall");
		wrp_free(ydotall);
		wrp_free(vec);
	}

	MPI_Finalize();
	return 0;
}
*/