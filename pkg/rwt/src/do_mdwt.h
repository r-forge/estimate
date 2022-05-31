/*
 * $Id$
 *
 * Public include for .Call interface to MDWT
 *
 * Copyright (c) 2004 MD Anderson Cancer Center. All rights reserved.
 * Created by Paul Roebuck, Department of Bioinformatics, MDACC.
 */

#ifndef DO_MDWT_H
#define DO_MDWT_H	1

#include <Rdefines.h>
#include "mdwt.h"


/*
 * Function Declarations
 */
extern SEXP do_mdwt(
  SEXP vntX,
  SEXP vntH,
  SEXP vntL);

#endif /* DO_MDWT_H */

