/*
 * $Id$
 *
 * Public include for discrete wavelet transform method 
 *
 * Copyright (c) 2004 MD Anderson Cancer Center. All rights reserved.
 * Created by Paul Roebuck, Department of Bioinformatics, MDACC.
 */

#ifndef MDWT_H
#define MDWT_H	1


/*
 * Function Declarations
 */
extern void MDWT(
  double *  /* x  */,
  int       /* m  */,
  int       /* n  */,
  double *  /* h  */,
  int       /* lh */,
  int       /* L  */,
  double *  /* y  */
);

#endif /* MDWT_H */

