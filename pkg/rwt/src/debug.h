/*
 * $Id$
 *
 * Public include for debugging 
 *
 * Copyright (c) 2004 MD Anderson Cancer Center. All rights reserved.
 * Created by Paul Roebuck, Department of Bioinformatics, MDACC.
 */

#ifndef DEBUG_H
#define DEBUG_H	1


/*
 * Macros
 */
#define DumpArray(arr,larr)\
    {\
        int i;\
        printf("{\n");\
        for (i = 0; i < larr; i++)\
            printf("\t\t[%d]\t%lg\n", i, arr[i]);\
        printf("\t}\n");\
    }\


#endif /* DEBUG_H */

