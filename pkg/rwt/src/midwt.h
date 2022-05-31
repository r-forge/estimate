/*
 * $Id$
 *
 * Public include for inverse discrete wavelet transform method 
 *
 * Copyright (c) 2004 MD Anderson Cancer Center. All rights reserved.
 * Created by Paul Roebuck, Department of Bioinformatics, MDACC.
 */

#ifndef MIDWT_H
#define MIDWT_H	1


/*
 * Function Declarations
 */
extern void MIDWT(
  double *  /* x  */,
  int       /* m  */,
  int       /* n  */,
  double *  /* h  */,
  int       /* lh */,
  int       /* L  */,
  double *  /* y  */
);

#endif /* MIDWT_H */

