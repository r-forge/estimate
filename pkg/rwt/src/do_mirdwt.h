/*
 * $Id$
 *
 * Public include for .Call interface to MIRDWT 
 *
 * Copyright (c) 2004 MD Anderson Cancer Center. All rights reserved.
 * Created by Paul Roebuck, Department of Bioinformatics, MDACC.
 */

#ifndef DO_MIRDWT_H
#define DO_MIRDWT_H	1

#include <Rdefines.h>
#include "mirdwt.h"


/*
 * Function Declarations
 */
extern SEXP do_mirdwt(
  SEXP vntYl,
  SEXP vntYh,
  SEXP vntH,
  SEXP vntL);

#endif /* DO_MIRDWT_H */

