/*
 * $Id$
 *
 * Package initialization
 *
 * Copyright (c) 2005 MD Anderson Cancer Center. All rights reserved.
 * Created by Paul Roebuck, Department of Bioinformatics, MDACC.
 */

#include <R.h>
#include <R_ext/Rdynload.h>
#include "do_mdwt.h"
#include "do_midwt.h"
#include "do_mirdwt.h"
#include "do_mrdwt.h"


/*
 * Module Variables
 */
static const R_CallMethodDef callMethods[] = {
    {"do_mdwt",   (DL_FUNC) &do_mdwt,   3},
    {"do_midwt",  (DL_FUNC) &do_midwt,  3},
    {"do_mirdwt", (DL_FUNC) &do_mirdwt, 4},
    {"do_mrdwt",  (DL_FUNC) &do_mrdwt,  3},
    {NULL, NULL, 0}
};


/*
 * Public
 */
void R_init_rwt(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

