/*
 * $Id$
 *
 * Public include for inverse redundant discrete wavelet transform method 
 *
 * Copyright (c) 2004 MD Anderson Cancer Center. All rights reserved.
 * Created by Paul Roebuck, Department of Bioinformatics, MDACC.
 */

#ifndef MIRDWT_H
#define MIRDWT_H	1


/*
 * Function Declarations
 */
extern void MIRDWT(
  double *  /* x  */,
  int       /* m  */,
  int       /* n  */,
  double *  /* h  */,
  int       /* lh */,
  int       /* L  */,
  double *  /* yl */,
  double *  /* yh */
);

#endif /* MIRDWT_H */

