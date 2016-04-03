#include <cs.h>
#include <cs_types.h>


cs *emptyMatrix(rqci m, rqci n, rqci triplet)
{
	cs *A, *tmp;
	tmp = cs_spalloc(m, n, 0, 1, 1);
	if(!triplet)
	{
		cs_entry(tmp, 0, 0, 0);//Bug
		A = cs_compress(tmp);
		cs_spfree(tmp);
	}
	else A = tmp;
	return A;
}

cs *cs_copy(cs* A)
{
	cs *B;
	csi m, n, nzmax, triplet;
	m = A->m;
	n = A->n;
	nzmax = A->nzmax;
	triplet = !(A && (A->nz == -1));

	B = cs_spalloc(m, n, nzmax, 1, triplet);
	
	memcpy(B->p, A->p, (triplet ? nzmax : n+1) * sizeof(csi));
	memcpy(B->i, A->i, nzmax * sizeof(csi));
	memcpy(B->x, A->x, nzmax * sizeof(realtype));
	B->nz = A->nz;

	cs_sprealloc(B, 0);
	
	return B;
}

cs* cs_sub (cs *A, csi starti, csi endi, csi startj, csi endj, csi triplet)
{
	cs* Asub;
    csi k = 0, i, j, p, nz = 0, n, msub, nsub, *Ap, *Ai, *Asubp, *Asubi;
    double *Ax, *Asubx;
    
	msub = endi - starti + 1;
	nsub = endj - startj + 1;

	
    if (!CS_CSC (A)) return (NULL) ;    /* check inputs */

	Asub = cs_spalloc (msub, nsub, msub*nsub, 1, triplet) ;

    Ap = A->p ; Ai = A->i ; Ax = A->x ; n = A->n;
    Asubp = Asub->p; Asubi = Asub->i; Asubx = Asub->x;
    for (j = 0 ; j < n ; j++)
    {
        p = Ap [j] ;                        /* get current location of col j */
        if(j >= startj && j <= endj)
        {	
        	if(!triplet)
        		Asubp [k++] = nz ;                    /* record new location of col j */
    		
	        for ( ; p < Ap [j+1] ; p++)
	        {
	        	i = Ai [p];
	        	if (i >= starti && i <=endi)
	            {
	            	if(!triplet)
	            	{
	                	if (Ax) Asubx [nz] = Ax [p] ;  /* keep A(i,j) */
	                	Asubi [nz++] = Ai [p] - starti;
	                }
	                else
	                {
	                	Asubi[nz] = i - starti;
	                	Asubp[nz] = j - startj;
	                	Asubx[nz++] = Ax[p];
	                	Asub->nz = nz;
	                }
	                	
	            }
	        }
    	}
    }
    if(!triplet)
    	Asubp [nsub] = nz ;	
   	/* finalize A */
    cs_sprealloc (Asub, 0) ; 
    return (Asub) ;
}
/*
csi is_sym (cs *A)
{
    csi is_upper, is_lower, j, p, n = A->n, m = A->m, *Ap = A->p, *Ai = A->i ;
    if (m != n) return (0) ;
    is_upper = 1 ;
    is_lower = 1 ;
    for (j = 0 ; j < n ; j++)
    {
        for (p = Ap [j] ; p < Ap [j+1] ; p++)
        {
            if (Ai [p] > j) is_upper = 0 ;
            if (Ai [p] < j) is_lower = 0 ;
        }
    }
    return (is_upper ? 1 : (is_lower ? -1 : 0)) ;
}*/