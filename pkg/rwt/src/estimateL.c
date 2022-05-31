/*
 * $Id$
 *
 * Level estimation routine
 *
 * Copyright (c) 2000 Rice University. All rights reserved.
 * Created by Markus Lang, Department of ECE, Rice University (1994).
 * Refactored by Paul Roebuck, Department of Bioinformatics, MDACC (2004).

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

#include "estimateL.h"


/*
 * Macros
 */
#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#define even(x)  (((x) & 1) ? 0 : 1)


/*
 * Public Routines
 */
int estimateL(int n, int m)
{
    register int L = m;
    register int i = n;
    register int j;

    for (j = 0; even(i); j++)
    {
        i = (i >> 1);
    }

    for (i = 0; even(L); i++)
    {
        L = (L >> 1);
    }

    L = (min(m,n) == 1) ? max(i,j) : min(i,j);

    return L;
}

