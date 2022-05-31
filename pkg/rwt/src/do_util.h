/*
 * $Id$
 *
 * Public include for utility routines associated with .Call interfaces 
 *
 * Copyright (c) 2004 MD Anderson Cancer Center. All rights reserved.
 * Created by Paul Roebuck, Department of Bioinformatics, MDACC.
 */

#ifndef DO_UTIL_H
#define DO_UTIL_H	1

#include <Rdefines.h>


/*
 * Function Declarations
 */
extern int GetMatrixDimen(
  SEXP vntX,
  int *nrow,
  int *ncol);

#endif /* DO_UTIL_H */

