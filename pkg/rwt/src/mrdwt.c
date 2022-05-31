/*
 * $Id$
 *
 * Redundant discrete wavelet transform method
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
#include "mrdwt.h"


/*
 * Macros
 */
#define mat(a, i, j) (*(a + (m*(j)+i)))
#define max(a, b) ((a) > (b) ? (a) : (b))


/*
 * Private Function Prototypes
 */
static void fpconv(
  double *  /* x_in   */,
  int       /* lx     */,
  double *  /* h0     */,
  double *  /* h1     */,
  int       /* lh     */,
  double *  /* x_outl */,
  double *  /* x_outh */
);


/*
 * Public Routines
 */
void MRDWT(
  double *x,
  int m,
  int n,
  double *h,
  int lh,
  int L,
  double *yl,
  double *yh)
{
    double *h0, *h1, *ydummyll, *ydummylh, *ydummyhl;
    double *ydummyhh, *xdummyl , *xdummyh;
    long i;
    int actual_L, actual_m, actual_n, c_o_a, ir, n_c, n_cb;
    int ic, n_r, n_rb, c_o_a_p2n, sample_f;

#ifdef DEBUG_RWT
    printf("In MRDWT()...\n");
#endif

    xdummyl = (double *)calloc(max(m,n)+lh-1, sizeof(double));
    xdummyh = (double *)calloc(max(m,n)+lh-1, sizeof(double));
    ydummyll = (double *)calloc(max(m,n), sizeof(double));
    ydummylh = (double *)calloc(max(m,n), sizeof(double));
    ydummyhl = (double *)calloc(max(m,n), sizeof(double));
    ydummyhh = (double *)calloc(max(m,n), sizeof(double));
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

    actual_m = 2 * m;
    actual_n = 2 * n;
    for (i = 0; i < m*n; i++)
    {
        yl[i] = x[i];
    }

    /* main loop */
    sample_f = 1;
    for (actual_L = 1; actual_L <= L; actual_L++)
    {
        actual_m = actual_m / 2;
        actual_n = actual_n / 2;

        /* actual (level dependent) column offset */
        if (m == 1)
        {
            c_o_a = n*(actual_L-1);
        }
        else
        {
            c_o_a = 3*n*(actual_L-1);
        }
        c_o_a_p2n = c_o_a + 2*n;

        /* go by rows */
        n_cb = n / actual_n;                   /* # of column blocks per row */
        for (ir = 0; ir < m; ir++)             /* loop over rows */
        {
            for (n_c = 0; n_c < n_cb; n_c++)     /* loop within one row */
            {
                /* store in dummy variable */
                ic = -sample_f + n_c;
                for (i = 0; i < actual_n; i++)
                {
                    ic = ic + sample_f;
                    xdummyl[i] = mat(yl, ir, ic);
                }

                /* perform filtering lowpass/highpass */
                fpconv(xdummyl, actual_n, h0, h1, lh, ydummyll, ydummyhh);

                /* restore dummy variables in matrices */
                ic = -sample_f + n_c;
                for (i = 0; i < actual_n; i++)
                {
                    ic = ic + sample_f;
                    mat(yl, ir, ic) = ydummyll[i];
                    mat(yh, ir, c_o_a+ic) = ydummyhh[i];
                }
            }
        }

        /* go by columns in case of a 2D signal*/
        if (m > 1)
        {
            n_rb = m / actual_m;                /* # of row blocks per column */
            for (ic = 0; ic < n; ic++)           /* loop over column */
            { 
                for (n_r = 0; n_r < n_rb; n_r++) /* loop within one column */
                {
                    /* store in dummy variables */
                    ir = -sample_f + n_r;
                    for (i = 0; i < actual_m; i++)
                    {
                        ir = ir + sample_f;
                        xdummyl[i] = mat(yl, ir, ic);
                        xdummyh[i] = mat(yh, ir,c_o_a+ic);
                    }

                    /* perform filtering: first LL/LH, then HL/HH */
                    fpconv(xdummyl, actual_m, h0, h1, lh, ydummyll, ydummylh);
                    fpconv(xdummyh, actual_m, h0, h1, lh, ydummyhl, ydummyhh);

                    /* restore dummy variables in matrices */
                    ir = -sample_f + n_r;
                    for (i = 0; i < actual_m; i++)
                    {
                        ir = ir + sample_f;
                        mat(yl, ir, ic) = ydummyll[i];
                        mat(yh, ir, c_o_a+ic) = ydummylh[i];
                        mat(yh, ir,c_o_a+n+ic) = ydummyhl[i];
                        mat(yh, ir, c_o_a_p2n+ic) = ydummyhh[i];
                    }
                }
            }
        }
        sample_f = sample_f*2;
    }
}


/*
 * Private Routines
 */
static void fpconv(
  double *x_in,
  int lx,
  double *h0,
  double *h1,
  int lh,
  double *x_outl,
  double *x_outh)
{
    register int i;

    for (i = lx; i < lx+lh-1; i++)
    {
        x_in[i] = x_in[i-lx];
    }

    for (i = 0; i < lx; i++)
    {
        double x0, x1;
        register int j;

        x0 = 0.0;
        x1 = 0.0;
        for (j = 0; j < lh; j++)
        {
            x0 = x0 + x_in[j+i]*h0[lh-1-j];
            x1 = x1 + x_in[j+i]*h1[lh-1-j];
        }
        x_outl[i] = x0;
        x_outh[i] = x1;
    }
}

