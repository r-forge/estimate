/*
 * $Id$
 *
 * .Call interface to inverse redundant discrete wavelet transform method
 *
 * Copyright (c) 2004 MD Anderson Cancer Center. All rights reserved.
 * Created by Paul Roebuck, Department of Bioinformatics, MDACC.
 */

#include <R.h>
#include <Rdefines.h>
#include "do_mirdwt.h"
#include "do_util.h"


/*
 * Macros
 */
#define isint(x) ((x - floor(x)) > 0.0 ? 0 : 1)
#define min(a,b) (((a) < (b)) ? (a) : (b))


/*
 * Public
 */
SEXP do_mirdwt(SEXP vntYl, SEXP vntYh, SEXP vntH, SEXP vntL)
{
    SEXP vntOut;
    SEXP vntX;
    SEXP vntLr;
    double *x, *h, *yl, *yh;
    int m, n, mh, nh, lh, L;

#ifdef DEBUG_RWT
    REprintf("In do_mirdwt(yl, yh, h, L)...\n");
#endif

    /*
     * Handle first parameter (numeric matrix)
     */
#ifdef DEBUG_RWT
    REprintf("\tfirst param 'yl'\n");
#endif
    if (GetMatrixDimen(vntYl, &m, &n) != 2)
    {
        error("'yl' is not a two dimensional matrix");
        /*NOTREACHED*/
    }
    PROTECT(vntYl = AS_NUMERIC(vntYl));
    yl = NUMERIC_POINTER(vntYl);
#ifdef DEBUG_RWT
    REprintf("yl[%d][%d] = 0x%p\n", m, n, yl);
#endif

    /*
     * Handle second parameter (numeric matrix)
     */
#ifdef DEBUG_RWT
    REprintf("\tsecond param 'yh'\n");
#endif
    if (GetMatrixDimen(vntYh, &mh, &nh) != 2)
    {
        error("'yh' is not a two dimensional matrix");
        /*NOTREACHED*/
    }
    PROTECT(vntYh = AS_NUMERIC(vntYh));
    yh = NUMERIC_POINTER(vntYh);
#ifdef DEBUG_RWT
    REprintf("yh[%d][%d] = 0x%p\n", mh, nh, yh);
#endif

    /*
     * Handle third parameter (numeric vector)
     */
#ifdef DEBUG_RWT
    REprintf("\tthird param 'h'\n");
#endif
    PROTECT(vntH = AS_NUMERIC(vntH));
    h = NUMERIC_POINTER(vntH);
    lh = GET_LENGTH(vntH);
#ifdef DEBUG_RWT
    REprintf("h[%d] = 0x%p\n", GET_LENGTH(vntH), h);
#endif

    /*
     * Handle fourth parameter (integer scalar)
     */
#ifdef DEBUG_RWT
    REprintf("\tfourth param 'L'\n");
#endif
    {
        PROTECT(vntL = AS_INTEGER(vntL));
        {
            int *piL = INTEGER_POINTER(vntL);
            L = piL[0];
        }
        UNPROTECT(1);
    }
#ifdef DEBUG_RWT
    REprintf("L = %d\n", L);
#endif

#ifdef DEBUG_RWT
    REprintf("\tcheck number of levels\n");
#endif
    if (L < 0)
    {
        error("The number of levels, L, must be a non-negative integer");
        /*NOTREACHED*/
    }

    /* check for consistency of rows and columns of yl, yh */
#ifdef DEBUG_RWT
    REprintf("\tcheck row/column consistency of signal components\n");
#endif
    if (min(m,n) > 1)
    {
        if (!((m == mh) && (3*n*L == nh)))
        {
            error("Dimensions of first two input matrices not consistent!");
            /*NOTREACHED*/
        }
    }
    else
    {
        if (!((m == mh) && (n*L == nh)))
        {
            error("Dimensions of first two input vectors not consistent!");
            /*NOTREACHED*/
        }
    }

#ifdef DEBUG_RWT
    REprintf("\tcheck dimen prereqs\n");
#endif
    /* Check the ROW dimension of input */
    if (m > 1)
    {
        double mtest = (double) m / pow(2.0, (double) L);
        if (!isint(mtest))
        {
            error("The matrix row dimension must be of size m*2^(L)");
            /*NOTREACHED*/
        }
    }

    /* Check the COLUMN dimension of input */
    if (n > 1)
    {
        double ntest = (double) n / pow(2.0, (double) L);
        if (!isint(ntest))
        {
            error("The matrix column dimension must be of size n*2^(L)");
            /*NOTREACHED*/
        }
    }

#ifdef DEBUG_RWT
    REprintf("\tcreate value objects\n");
#endif

    /* Create x value object */
    {
#ifdef DEBUG_RWT
        REprintf("\tcreate 'x' value object\n");
#endif
        PROTECT(vntX = NEW_NUMERIC(n*m));
        x = NUMERIC_POINTER(vntX);

        /* Add dimension attribute to value object */
#ifdef DEBUG_RWT
        REprintf("\tconvert 'x' value object to matrix\n");
#endif
        {
            SEXP vntDim;

            PROTECT(vntDim = NEW_INTEGER(2));
            INTEGER(vntDim)[0] = m;
            INTEGER(vntDim)[1] = n;
            SET_DIM(vntX, vntDim);
            UNPROTECT(1);
        }
    }

    /* Create Lr value object */
    {
#ifdef DEBUG_RWT
        REprintf("\tcreating 'Lr' value object\n");
#endif
        PROTECT(vntLr = NEW_INTEGER(1));
        INTEGER_POINTER(vntLr)[0] = L;
    }

#ifdef DEBUG_RWT
    REprintf("\tcompute inverse redundant discrete wavelet transform\n");
#endif
    MIRDWT(x, m, n, h, lh, L, yl, yh);

    /* Unprotect params */
    UNPROTECT(3);

#ifdef DEBUG_RWT
    REprintf("\tcreate list output object\n");
#endif
    PROTECT(vntOut = NEW_LIST(2));

#ifdef DEBUG_RWT
    REprintf("\tassigning value objects to list\n");
#endif
    SET_VECTOR_ELT(vntOut, 0, vntX);
    SET_VECTOR_ELT(vntOut, 1, vntLr);

    /* Unprotect value objects */
    UNPROTECT(2);

    {
        SEXP vntNames;

#ifdef DEBUG_RWT
        REprintf("\tassigning names to value objects in list\n");
#endif
        PROTECT(vntNames = NEW_CHARACTER(2));
        SET_STRING_ELT(vntNames, 0, CREATE_STRING_VECTOR("x"));
        SET_STRING_ELT(vntNames, 1, CREATE_STRING_VECTOR("L"));
        SET_NAMES(vntOut, vntNames);
        UNPROTECT(1);
    }

    /* Unprotect output object */
    UNPROTECT(1);

#ifdef DEBUG_RWT
    REprintf("\treturning output...\n");
#endif

    return vntOut;
}

