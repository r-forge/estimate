/*
 * $Id$
 *
 * Utility routines associated with .Call interfaces
 *
 * Copyright (c) 2004 MD Anderson Cancer Center. All rights reserved.
 * Created by Paul Roebuck, Department of Bioinformatics, MDACC.
 */

#include <R.h>
#include <Rdefines.h>
#include "do_util.h"


/*
 * Public
 */
int GetMatrixDimen(SEXP vntX, int *nrow, int *ncol)
{
    SEXP vntXdim;
    int nX;

    PROTECT(vntXdim = GET_DIM(vntX));
    nX = GET_LENGTH(vntXdim);
    if (nX == 2)
    {
        int *piXdim = INTEGER_POINTER(vntXdim);

        *nrow = piXdim[0];
        *ncol = piXdim[1];
    }
    else
    {
        *nrow = -1;
        *ncol = -1;
    }
    UNPROTECT(1);

    return nX;
}

