/*
 * $Id$
 *
 * Public include for .Call interface to MIDWT 
 *
 * Copyright (c) 2004 MD Anderson Cancer Center. All rights reserved.
 * Created by Paul Roebuck, Department of Bioinformatics, MDACC.
 */

#ifndef DO_MIDWT_H
#define DO_MIDWT_H	1

#include <Rdefines.h>
#include "midwt.h"


/*
 * Function Declarations
 */
extern SEXP do_midwt(
  SEXP vntX,
  SEXP vntH,
  SEXP vntL);

#endif /* DO_MIDWT_H */

