/*
 * $Id$
 *
 * .Call interface to redundant discrete wavelet transform method
 *
 * Copyright (c) 2004 MD Anderson Cancer Center. All rights reserved.
 * Created by Paul Roebuck, Department of Bioinformatics, MDACC.
 */

#include <R.h>
#include <Rdefines.h>
#include "do_mrdwt.h"
#include "do_util.h"


/*
 * Macros
 */
#define isint(x) ((x - floor(x)) > 0.0 ? 0 : 1)
#define min(a,b) (((a) < (b)) ? (a) : (b))


/*
 * Public
 */
SEXP do_mrdwt(SEXP vntX, SEXP vntH, SEXP vntL)
{
    SEXP vntOut;
    SEXP vntYl;
    SEXP vntYh;
    SEXP vntLr;
    double *x, *h, *yl, *yh;
    int m, n, lh, L;

#ifdef DEBUG_RWT
    REprintf("In do_mrdwt(x, h, L)...\n");
#endif

    /*
     * Handle first parameter (numeric matrix)
     */
#ifdef DEBUG_RWT
    REprintf("\tfirst param 'x'\n");
#endif
    if (GetMatrixDimen(vntX, &m, &n) != 2)
    {
        error("'x' is not a two dimensional matrix");
        /*NOTREACHED*/
    }
    PROTECT(vntX = AS_NUMERIC(vntX));
    x = NUMERIC_POINTER(vntX);
#ifdef DEBUG_RWT
    REprintf("x[%d][%d] = 0x%p\n", m, n, x);
#endif

    /*
     * Handle second parameter (numeric vector)
     */
#ifdef DEBUG_RWT
    REprintf("\tsecond param 'h'\n");
#endif
    PROTECT(vntH = AS_NUMERIC(vntH));
    h = NUMERIC_POINTER(vntH);
    lh = GET_LENGTH(vntH);
#ifdef DEBUG_RWT
    REprintf("h[%d] = 0x%p\n", GET_LENGTH(vntH), h);
#endif

    /*
     * Handle third parameter (integer scalar)
     */
#ifdef DEBUG_RWT
    REprintf("\tthird param 'L'\n");
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
    REprintf("\tcreating value objects\n");
#endif

    /* Create yl value object */
    {
#ifdef DEBUG_RWT
        REprintf("\tcreating 'yl' value object\n");
#endif
        PROTECT(vntYl = NEW_NUMERIC(m*n));
        yl = NUMERIC_POINTER(vntYl);

        /* Add dimension attribute to value object */
#ifdef DEBUG_RWT
        REprintf("\tconvert 'yl' value object to matrix\n");
#endif
        {
            SEXP vntDim;

            PROTECT(vntDim = NEW_INTEGER(2));
            INTEGER(vntDim)[0] = m;
            INTEGER(vntDim)[1] = n;
            SET_DIM(vntYl, vntDim);
            UNPROTECT(1);
        }
    }

    /* Create yh value object */
    {
        int cols = (min(m,n) == 1) ? (L * n) : (3 * L * n);

#ifdef DEBUG_RWT
        REprintf("\tcreating 'yh' value object\n");
#endif
        PROTECT(vntYh = NEW_NUMERIC(m*cols));
        yh = NUMERIC_POINTER(vntYh);

        /* Add dimension attribute to value object */
#ifdef DEBUG_RWT
        REprintf("\tconvert 'yh' value object to matrix\n");
#endif
        {
            SEXP vntDim;

            PROTECT(vntDim = NEW_INTEGER(2));
            INTEGER(vntDim)[0] = m;
            INTEGER(vntDim)[1] = cols;
            SET_DIM(vntYh, vntDim);
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
    REprintf("\tcompute redundant discrete wavelet transform\n");
#endif
    MRDWT(x, m, n, h, lh, L, yl, yh);

    /* Unprotect params */
    UNPROTECT(2);

#ifdef DEBUG_RWT
    REprintf("\tcreate list output object\n");
#endif
    PROTECT(vntOut = NEW_LIST(3));

#ifdef DEBUG_RWT
    REprintf("\tassigning value objects to list\n");
#endif
    SET_VECTOR_ELT(vntOut, 0, vntYl);
    SET_VECTOR_ELT(vntOut, 1, vntYh);
    SET_VECTOR_ELT(vntOut, 2, vntLr);

    /* Unprotect value objects */
    UNPROTECT(3);

    {
        SEXP vntNames;

#ifdef DEBUG_RWT
        REprintf("\tassigning names to value objects in list\n");
#endif
        PROTECT(vntNames = NEW_CHARACTER(3));
        SET_STRING_ELT(vntNames, 0, CREATE_STRING_VECTOR("yl"));
        SET_STRING_ELT(vntNames, 1, CREATE_STRING_VECTOR("yh"));
        SET_STRING_ELT(vntNames, 2, CREATE_STRING_VECTOR("L"));
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

