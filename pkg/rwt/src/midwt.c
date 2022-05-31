/*
 * $Id$
 *
 * Inverse discrete wavelet transform method
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
#include "midwt.h"


/*
 * Macros
 */
#define max(A,B) ((A) > (B) ? (A) : (B))
#define mat(a, i, j) (*(a + (m*(j)+i)))  /* macro for matrix indices */


/*
 * Private Function Prototypes
 */
static void bpsconv(
  double *  /* x_out */,
  int       /* lx    */,
  double *  /* g0    */,
  double *  /* g1    */,
  int       /* lhm1  */,
  int       /* lhhm1 */,
  double *  /* x_inl */,
  double *  /* x_inh */
);


/*
 * Public Routines
 */
void MIDWT(
  double *x,
  int m,
  int n,
  double *h,
  int lh,
  int L,
  double *y)
{
    double *g0, *g1, *ydummyl, *ydummyh, *xdummy;
    long i;
    int actual_L, actual_m, actual_n, r_o_a, c_o_a;
    int ir, ic, lhm1, lhhm1, sample_f;

#ifdef DEBUG_RWT
    printf("In MIDWT()...\n");
#endif

    xdummy = (double *)calloc(max(m,n), sizeof(double));
    ydummyl = (double *)calloc(max(m,n)+lh/2-1, sizeof(double));
    ydummyh = (double *)calloc(max(m,n)+lh/2-1, sizeof(double));
    g0 = (double *)calloc(lh, sizeof(double));
    g1 = (double *)calloc(lh, sizeof(double));

    /* transpose if row vector */
    if (n == 1)
    {
        n = m;
        m = 1;
    }

    /* synthesis lowpass and highpass */
    for (i = 0; i < lh; i++)
    {
        g0[i] = h[i];
        g1[i] = h[lh-i-1];
    }

    for (i = 1; i <= lh; i += 2)
    {
        g1[i] = -g1[i];
    }

    lhm1 = lh - 1;
    lhhm1 = lh/2 - 1;

    /* 2^L */
    sample_f = 1;
    for (i = 1; i < L; i++)
    {
        sample_f = sample_f*2;
    }

    if (m > 1)
    {
        actual_m = m / sample_f;
    }
    else
    {
        actual_m = 1;
    }
    actual_n = n / sample_f;

    for (i = 0; i < (m*n); i++)
    {
        x[i] = y[i];
    }

    /* main loop */
    for (actual_L = L; actual_L >= 1; actual_L--)
    {
        r_o_a = actual_m / 2;
        c_o_a = actual_n / 2;

        /* go by columns in case of a 2D signal*/
        if (m > 1)
        {
            for (ic = 0; ic < actual_n; ic++)         /* loop over column */
            {
                /* store in dummy variables */
                ir = r_o_a;
                for (i = 0; i < r_o_a; i++)
                {
                    ydummyl[i+lhhm1] = mat(x, i, ic);
                    ydummyh[i+lhhm1] = mat(x, ir++, ic);
                }

                /* perform filtering lowpass and highpass*/
                bpsconv(xdummy, r_o_a, g0, g1, lhm1, lhhm1, ydummyl, ydummyh);

                /* restore dummy variables in matrix */
                for (i = 0; i < actual_m; i++)
                {
                    mat(x, i, ic) = xdummy[i];
                }
            }
        }

        /* go by rows */
        for (ir = 0; ir < actual_m; ir++)             /* loop over rows */
        {
            /* store in dummy variable */
            ic = c_o_a;
            for (i = 0; i < c_o_a; i++)
            {
                ydummyl[i+lhhm1] = mat(x, ir, i);
                ydummyh[i+lhhm1] = mat(x, ir, ic++);
            }

            /* perform filtering lowpass and highpass*/
            bpsconv(xdummy, c_o_a, g0, g1, lhm1, lhhm1, ydummyl, ydummyh);

            /* restore dummy variables in matrices */
            for (i = 0; i < actual_n; i++)
            {
                mat(x, ir, i) = xdummy[i];
            }
        }

        if (m == 1)
        {
            actual_m = 1;
        }
        else
        {
            actual_m = actual_m * 2;
        }
        actual_n = actual_n * 2;
    }
}


/*
 * Private Routines
 */
static void bpsconv(
  double *x_out,
  int lx,
  double *g0,
  double *g1,
  int lhm1,
  int lhhm1,
  double *x_inl,
  double *x_inh)
{
    register int i;
    int ind;

    for (i = lhhm1-1; i > -1; i--)
    {
        x_inl[i] = x_inl[lx+i];
        x_inh[i] = x_inh[lx+i];
    }

    ind = 0;
    for (i = 0; i < lx; i++)
    {
        double x0, x1;
        register int j, tj;

        x0 = 0.0;
        x1 = 0.0;
        tj = -2;
        for (j = 0; j <= lhhm1; j++)
        {
            tj += 2;
            x0 = x0 + x_inl[i+j]*g0[lhm1-1-tj] +
                      x_inh[i+j]*g1[lhm1-1-tj];
            x1 = x1 + x_inl[i+j]*g0[lhm1-tj] +
                      x_inh[i+j]*g1[lhm1-tj];
        }
        x_out[ind++] = x0;
        x_out[ind++] = x1;
    }
}

