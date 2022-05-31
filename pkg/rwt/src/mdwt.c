/*
 * $Id$
 *
 * Discrete wavelet transform method
 *
 * Copyright (c) 2000 Rice University. All rights reserved.
 * Created by Markus Lang, Department of ECE, Rice University (1994).
 * Modifications by Jan Erik Odegard, Department of ECE, Rice University (1995).
 * Modifications by Paul Roebuck, Department of Bioinformatics, MDACC (2004).

This software is distributed and licensed to you on a non-exclusive
basis, free-of-charge. Redistribution and use in source and binary forms,
with or without modification, are permitted provided that the following
conditions are met:

1. Redistribution of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistribution in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of the University nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY WILLIAM MARSH RICE UNIVERSITY, HOUSTON, TEXAS,
AND CONTRIBUTORS AS IS AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL RICE UNIVERSITY
OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTIONS) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE), PRODUCT LIABILITY, OR OTHERWISE ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE,  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 *
 */

#include <R.h>
#define printf	REprintf
#define calloc	R_alloc
#include "mdwt.h"


/*
 * Macros
 */
#define max(A,B) ((A) > (B) ? (A) : (B))
#define mat(a, i, j) (*(a + (m*(j)+i)))  /* macro for matrix indices */


/*
 * Private Function Prototypes
 */
static void fpsconv(
  double *  /* x_in   */,
  int       /* lx     */,
  double *  /* h0     */,
  double *  /* h1     */,
  int       /* lhm1   */,
  double *  /* x_outl */,
  double *  /* x_outh */
);


/*
 * Public Routines
 */
void MDWT(
  double *x,
  int m,
  int n,
  double *h,
  int lh,
  int L,
  double *y)
{
    double *h0, *h1, *ydummyl, *ydummyh, *xdummy;
    long i;
    int actual_L, actual_m, actual_n, r_o_a, c_o_a, ir, ic, lhm1;

#ifdef DEBUG_RWT
    printf("In MDWT()...\n");
#endif

    xdummy = (double *)calloc(max(m,n)+lh-1, sizeof(double));
    ydummyl = (double *)calloc(max(m,n), sizeof(double));
    ydummyh = (double *)calloc(max(m,n), sizeof(double));
    h0 = (double *)calloc(lh, sizeof(double));
    h1 = (double *)calloc(lh, sizeof(double));

    /* transpose if row vector */
    if (n == 1)
    {
        n = m;
        m = 1;
    }

    /* analysis lowpass and highpass */
    for (i = 0; i < lh; i++)
    {
        h0[i] = h[lh-i-1];
        h1[i] = h[i];
    }

    for (i = 0; i < lh; i += 2)
    {
        h1[i] = -h1[i];
    }

    lhm1 = lh - 1;
    actual_m = 2*m;
    actual_n = 2*n;

    /* main loop */
    for (actual_L = 1; actual_L <= L; actual_L++)
    {
        if (m == 1)
        {
            actual_m = 1;
        }
        else
        {
            actual_m = actual_m/2;
            r_o_a = actual_m/2;
        }
        actual_n = actual_n/2;
        c_o_a = actual_n/2;

        /* go by rows */
        for (ir = 0; ir < actual_m; ir++)             /* loop over rows */
        {
            /* store in dummy variable */
            for (i = 0; i < actual_n; i++)
            {
                if (actual_L == 1)
                    xdummy[i] = mat(x, ir, i);
                else
                    xdummy[i] = mat(y, ir, i);
            }

            /* perform filtering lowpass and highpass */
            fpsconv(xdummy, actual_n, h0, h1, lhm1, ydummyl, ydummyh);

            /* restore dummy variables in matrices */
            ic = c_o_a;
            for (i = 0; i < c_o_a; i++)
            {
                mat(y, ir, i) = ydummyl[i];
                mat(y, ir, ic++) = ydummyh[i];
            }
        }

        /* go by columns in case of a 2D signal */
        if (m > 1)
        {
            for (ic = 0; ic < actual_n; ic++)         /* loop over columns */
            {
                /* store in dummy variables */
                for (i = 0; i < actual_m; i++)
                    xdummy[i] = mat(y, i, ic);

                /* perform filtering lowpass and highpass*/
                fpsconv(xdummy, actual_m, h0, h1, lhm1, ydummyl, ydummyh);

                /* restore dummy variables in matrix */
                ir = r_o_a;
                for (i = 0; i < r_o_a; i++)
                {
                    mat(y, i, ic) = ydummyl[i];
                    mat(y, ir++, ic) = ydummyh[i];
                }
            }
        }
    }
}


/*
 * Private Routines
 */
static void fpsconv(
  double *x_in,
  int lx,
  double *h0,
  double *h1,
  int lhm1,
  double *x_outl,
  double *x_outh)
{
    register int i;
    int ind;

    for (i = lx; i < lx+lhm1; i++)
    {
        x_in[i] = *(x_in+(i-lx));
    }

    ind = 0;
    for (i = 0; i < lx; i += 2)
    {
        double x0, x1;
        register int j;

        x0 = 0.0;
        x1 = 0.0;
        for (j = 0; j <= lhm1; j++)
        {
            x0 = x0 + x_in[i+j]*h0[lhm1-j];
            x1 = x1 + x_in[i+j]*h1[lhm1-j];
        }
        x_outl[ind] = x0;
        x_outh[ind++] = x1;
    }
}

