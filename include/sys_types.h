#ifndef __SYS_TYPES
#define __SYS_TYPES

#ifndef _SYS_CONFIG_H
#define _SYS_CONFIG_H
#include <sys_config.h>
#endif

#define WITH_MPI
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stddef.h>//->ptrdiff_t
#include <float.h>
#include <string.h>
#include <mpi_ext.h>

/*
	Realtime optimal control of a quadcopter
*/

#ifndef rqci //realtime optimal control of a quadcopter index
#define rqci signed long int//ptrdiff_t
#endif

#ifndef MPI_RQCI
#define MPI_RQCI MPI_LONG 
#endif

#ifndef csi
#define csi rqci
#endif

#ifndef booleantype
#define booleantype int
#endif

#define EPS 0.02
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define zero_ind(i) (i - 1)

/*
 *------------------------------------------------------------------
 * Type realtype
 * Macro RCONST
 * Constants BIG_REAL, SMALL_REAL, and UNIT_ROUNDOFF
 *------------------------------------------------------------------
 */

#if defined(SUNDIALS_SINGLE_PRECISION)

typedef float realtype;
#define MPI_REALTYPE MPI_FLOAT
# define RCONST(x) x##F
# define BIG_REAL FLT_MAX
# define SMALL_REAL FLT_MIN
# define UNIT_ROUNDOFF FLT_EPSILON

#elif defined(SUNDIALS_DOUBLE_PRECISION)

typedef double realtype;
#define MPI_REALTYPE MPI_DOUBLE
# define RCONST(x) x
# define BIG_REAL DBL_MAX
# define SMALL_REAL DBL_MIN
# define UNIT_ROUNDOFF DBL_EPSILON

#elif defined(SUNDIALS_EXTENDED_PRECISION)

typedef long double realtype;
#define MPI_REALTYPE MPI_LONG_DOUBLE
# define RCONST(x) x##L
# define BIG_REAL LDBL_MAX
# define SMALL_REAL LDBL_MIN
# define UNIT_ROUNDOFF LDBL_EPSILON

#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#endif