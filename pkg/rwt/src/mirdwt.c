/*
 * $Id$
 *
 * Inverse redundant discrete wavelet transform method
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
#include "mirdwt.h"


/*
 * Macros
 */
#define max(a, b) ((a) > (b) ? (a) : (b))
#define mat(a, i, j) (*(a + (m*(j)+i)))  /* macro for matrix indices */
#ifdef DEBUG_RWT
#define DumpArray(arr,larr)\
    {\
        int i;\
        printf("{\n");\
        for (i = 0; i < larr; i++)\
            printf("\t\t[%d]\t%lg\n", i, arr[i]);\
        printf("\t}\n");\
    }\

#endif


/*
 * Private Function Prototypes
 */
static void bpconv(
  double *  /* x_out */,
  int       /* lx    */,
  double *  /* g0    */,
  double *  /* g1    */,
  int       /* lh    */,
  double *  /* x_inl */,
  double *  /* x_inh */
);


/*
 * Public Routines
 */
void MIRDWT(
  double *x,
  int m,
  int n,
  double *h,
  int lh,
  int L,
  double *yl,
  double *yh)
{
    double *g0, *g1, *ydummyll, *ydummylh, *ydummyhl;
    double *ydummyhh, *xdummyl , *xdummyh, *xh;
    long i;
    int actual_L, actual_m, actual_n, c_o_a, ir, n_c, n_cb, lhm1;
    int ic, n_r, n_rb, c_o_a_p2n, sample_f;

#ifdef DEBUG_RWT
    printf("In MIRDWT()...\n");
#endif

#if 0
#ifdef DEBUG_RWT
    /* Debug written specific for MIRDWT example case */
    printf("\t**** Input parameters ****\n");
    printf("\tyl[%d] = ", m*n);
    DumpArray(yl, m*n);

    printf("\tyh[%d] = ", m*n);
    DumpArray(yh, m*n);

    printf("\th[%d] = ", lh);
    DumpArray(h, lh);

    printf("\tL = %d\n", L);
#endif
#endif

    xh = (double *)calloc(m*n, sizeof(double));
    xdummyl = (double *)calloc(max(m,n), sizeof(double));
    xdummyh = (double *)calloc(max(m,n), sizeof(double));
    ydummyll = (double *)calloc(max(m,n)+lh-1, sizeof(double));
    ydummylh = (double *)calloc(max(m,n)+lh-1, sizeof(double));
    ydummyhl = (double *)calloc(max(m,n)+lh-1, sizeof(double));
    ydummyhh = (double *)calloc(max(m,n)+lh-1, sizeof(double));
    g0 = (double *)calloc(lh, sizeof(double));
    g1 = (double *)calloc(lh, sizeof(double));

    /* transpose if row vector */
    if (n == 1)
    {
        n = m;
        m = 1;
    }

    /* analysis lowpass and highpass */
    for (i = 0; i < lh; i++)
    {
        g0[i] = h[i]/2;
        g1[i] = h[lh-i-1]/2;
    }
    for (i = 1; i <= lh; i += 2)
    {
        g1[i] = -g1[i];
    }

    lhm1 = lh - 1;

    /* 2^L */
    sample_f = 1;
    for (i = 1; i < L; i++)
    {
        sample_f = sample_f * 2;
    }
    actual_m = m / sample_f;
    actual_n = n / sample_f;

    /* restore yl in x */
    for (i = 0; i < m*n; i++)
    {
        x[i] = yl[i];
    }
 
    /* main loop */
    for (actual_L = L; actual_L >= 1; actual_L--)
    {
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

        /* go by columns in case of a 2D signal*/
        if (m > 1)
        {
            n_rb = m / actual_m;               /* # of row blocks per column */
            for (ic = 0; ic < n; ic++)         /* loop over column */
            {
                for (n_r = 0; n_r < n_rb; n_r++)   /* loop within one column */
                {
                    /* store in dummy variables */
                    ir = -sample_f + n_r;
                    for (i = 0; i < actual_m; i++)
                    {
                        ir = ir + sample_f;
                        ydummyll[i+lhm1] = mat(x, ir, ic);
                        ydummylh[i+lhm1] = mat(yh, ir, c_o_a+ic);
                        ydummyhl[i+lhm1] = mat(yh, ir,c_o_a+n+ic);
                        ydummyhh[i+lhm1] = mat(yh, ir, c_o_a_p2n+ic);
                    }

                    /* perform filtering and adding: first LL/LH, then HL/HH */
                    bpconv(xdummyl, actual_m, g0, g1, lh, ydummyll, ydummylh);
                    bpconv(xdummyh, actual_m, g0, g1, lh, ydummyhl, ydummyhh);

                    /* store dummy variables in matrices */
                    ir = -sample_f + n_r;
                    for (i = 0; i < actual_m; i++)
                    {
                        ir = ir + sample_f;
                        mat(x, ir, ic) = xdummyl[i];
                        mat(xh, ir, ic) = xdummyh[i];
                    }
                }
            }
        }

        /* go by rows */
        n_cb = n / actual_n;                 /* # of column blocks per row */
        for (ir = 0; ir < m; ir++)           /* loop over rows */
        {
            for (n_c = 0; n_c < n_cb; n_c++)     /* loop within one row */
            {
                /* store in dummy variable */
                ic = -sample_f + n_c;
                for (i = 0; i < actual_n; i++)
                {
                    ic = ic + sample_f;
                    ydummyll[i+lhm1] = mat(x, ir, ic);
                    if (m > 1)
                    {
                        ydummyhh[i+lhm1] = mat(xh, ir, ic);
                    }
                    else
                    {
                        ydummyhh[i+lhm1] = mat(yh, ir, c_o_a+ic);
                    }
                }

                /* perform filtering lowpass/highpass */
                bpconv(xdummyl, actual_n, g0, g1, lh, ydummyll, ydummyhh);

                /* restore dummy variables in matrices */
                ic = -sample_f + n_c;
                for (i = 0; i < actual_n; i++)
                {
                    ic = ic + sample_f;
                    mat(x, ir, ic) = xdummyl[i];
                }
            }
        }
        sample_f = sample_f / 2;
        actual_m = actual_m * 2;
        actual_n = actual_n * 2;
    }

#if 0
#ifdef DEBUG_RWT
    printf("\t**** Output parameters ****\n");
    printf("\tx[%d] = ", m*n);
    DumpArray(x, m*n);
#endif
#endif
}


/*
 * Private Routines
 */
static void bpconv(
  double *x_out,
  int lx,
  double *g0,
  double *g1,
  int lh,
  double *x_inl,
  double *x_inh)
{
    register int i;

    for (i = lh-2; i > -1; i--)
    {
        x_inl[i] = x_inl[lx+i];
        x_inh[i] = x_inh[lx+i];
    }

    for (i = 0; i < lx; i++)
    {
        double x0;
        register int j;

        x0 = 0.0;
        for (j = 0; j < lh; j++)
        {
            x0 = x0 + x_inl[j+i]*g0[lh-1-j] +
                      x_inh[j+i]*g1[lh-1-j];
        }
        x_out[i] = x0;
    }
}

