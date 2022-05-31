/*
 * $Id$
 *
 * .Call interface to inverse discrete wavelet transform method
 *
 * Copyright (c) 2004 MD Anderson Cancer Center. All rights reserved.
 * Created by Paul Roebuck, Department of Bioinformatics, MDACC.
 */

#include <R.h>
#include <Rdefines.h>
#include "do_midwt.h"
#include "do_util.h"


/*
 * Macros
 */
#define isint(x) ((x - floor(x)) > 0.0 ? 0 : 1)


/*
 * Public
 */
SEXP do_midwt(SEXP vntY, SEXP vntH, SEXP vntL)
{
    SEXP vntOut;
    SEXP vntX;
    SEXP vntLr;
    double *x, *h, *y;
    int m, n, lh, L;

#ifdef DEBUG_RWT
    REprintf("In do_midwt(y, h, L)...\n");
#endif

    /*
     * Handle first parameter (numeric matrix)
     */
#ifdef DEBUG_RWT
    REprintf("\tfirst param 'y'\n");
#endif
    if (GetMatrixDimen(vntY, &m, &n) != 2)
    {
        error("'y' is not a two dimensional matrix");
        /*NOTREACHED*/
    }
    PROTECT(vntY = AS_NUMERIC(vntY));
    y = NUMERIC_POINTER(vntY);
#ifdef DEBUG_RWT
    REprintf("y[%d][%d] = 0x%p\n", m, n, y);
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
    REprintf("\tcompute inverse discrete wavelet transform\n");
#endif
    MIDWT(x, m, n, h, lh, L, y);

    /* Unprotect params */
    UNPROTECT(2);

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

