/*
 * $Id$
 *
 * Public include for .Call interface to MRDWT 
 *
 * Copyright (c) 2004 MD Anderson Cancer Center. All rights reserved.
 * Created by Paul Roebuck, Department of Bioinformatics, MDACC.
 */

#ifndef DO_MRDWT_H
#define DO_MRDWT_H	1

#include <Rdefines.h>
#include "mrdwt.h"


/*
 * Function Declarations
 */
extern SEXP do_mrdwt(
  SEXP vntX,
  SEXP vntH,
  SEXP vntL);

#endif /* DO_MRDWT_H */

