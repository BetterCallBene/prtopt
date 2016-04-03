#include <util.h>
#include <cs.h>


realtype *load_vec(FILE* f, csi n)
{
    csi i;
    realtype x;
    realtype* vec;
    
    vec = cs_calloc(n, sizeof(realtype));
    for ( i=0;fscanf (f, "%lg\n", &x) == 1;i++)
        vec[i] = x;

    printf("load_vec: i=%ld\n", i);
    return vec;
}

void print_vec(realtype* vec, rqci n, char* Name)
{
	rqci i;
	printf("Name: %s, Vector length: %ld\n", Name, n);
	for(i = 0; i < n; i++)
		printf("%ld : %f\n", i, vec[i]);
	printf("\n");	
}

void iprint_vec(rqci* vec, rqci n, char* Name)
{
	rqci i;
	printf("Name: %s, Vector length: %ld\n", Name, n);
	for(i = 0; i < n; i++)
		printf("%ld : %ld\n", i, vec[i]);
	printf("\n");	
}


/* print a sparse matrix; use %g for integers to avoid differences with csi */
csi cs_print_ext (const cs *A, csi brief, char* Name)
{
    csi p, j, m, n, nzmax, nz, *Ap, *Ai ;
    realtype *Ax ;

    printf("Name: %s\n", Name);
    if (!A) { printf ("(null)\n") ; return (0) ; }
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    nzmax = A->nzmax ; nz = A->nz ;
//    printf ("CSparse Version %ld.%ld.%ld, %s.  %s\n", CS_VER, CS_SUBVER,
//        CS_SUBSUB, CS_DATE, CS_COPYRIGHT) ;
    if (nz < 0)
    {
        printf ("%g-by-%g, nzmax: %g nnz: %g, 1-norm: %g\n", (realtype) m,
            (realtype) n, (realtype) nzmax, (realtype) (Ap [n]), cs_norm (A)) ;
        for (j = 0 ; j < n ; j++)
        {
            printf ("    col %g : locations %g to %g\n", (realtype) j, 
                (realtype) (Ap [j]), (realtype) (Ap [j+1]-1)) ;
            for (p = Ap [j] ; p < Ap [j+1] ; p++)
            {
                printf ("      %g : %g\n", (realtype) (Ai [p]), Ax ? Ax [p] : 1) ;
                if (brief && p > 20) { printf ("  ...\n") ; return (1) ; }
            }
        }
    }
    else
    {
        printf ("triplet: %g-by-%g, nzmax: %g nnz: %g\n", (realtype) m,
            (realtype) n, (realtype) nzmax, (realtype) nz) ;
        for (p = 0 ; p < nz ; p++)
        {
            printf ("    %g %g : %g\n", (realtype) (Ai [p]), (realtype) (Ap [p]),
                Ax ? Ax [p] : 1) ;
            if (brief && p > 20) { printf ("  ...\n") ; return (1) ; }
        }
    }
    return (1) ;
}
