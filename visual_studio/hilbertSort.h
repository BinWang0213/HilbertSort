#pragma once
#include <iostream>

static unsigned const int IMAX = ~(0U);

static unsigned const int idata2d[] =  /* 2 dimension to nkey conversion */
{ 0, 3, 1, 2,
 0, 1, 3, 2,
 2, 3, 1, 0,
 2, 1, 3, 0 };

static unsigned const int istate2d[] = /* 2 dimension to nkey state transitions */
{ 1, 2, 0, 0,
 0, 1, 3, 1,
 2, 0, 2, 3,
 3, 3, 1, 2 };

static unsigned const int data2d[] = /* nkey to 2 dimension conversion */
{ 0, 2, 3, 1,
 0, 1, 3, 2,
 3, 2, 0, 1,
 3, 1, 0, 2 };

static unsigned const int state2d[] = /* nkey to 2 dimension state transitions */
{ 1, 0, 0, 2,
 0, 1, 1, 3,
 3, 2, 2, 0,
 2, 3, 3, 1 };

static unsigned const data3d[] = {  /* nkey to 3 dimension conversion */
 0,  4,  6,  2,  3,  7,  5,  1,
 0,  1,  3,  2,  6,  7,  5,  4,
 0,  4,  5,  1,  3,  7,  6,  2,
 5,  4,  0,  1,  3,  2,  6,  7,
 6,  7,  3,  2,  0,  1,  5,  4,
 3,  7,  6,  2,  0,  4,  5,  1,
 5,  4,  6,  7,  3,  2,  0,  1,
 0,  1,  5,  4,  6,  7,  3,  2,
 5,  1,  0,  4,  6,  2,  3,  7,
 5,  1,  3,  7,  6,  2,  0,  4,
 0,  2,  6,  4,  5,  7,  3,  1,
 3,  1,  0,  2,  6,  4,  5,  7,
 5,  7,  6,  4,  0,  2,  3,  1,
 6,  7,  5,  4,  0,  1,  3,  2,
 3,  1,  5,  7,  6,  4,  0,  2,
 0,  2,  3,  1,  5,  7,  6,  4,
 3,  2,  0,  1,  5,  4,  6,  7,
 3,  2,  6,  7,  5,  4,  0,  1,
 6,  2,  0,  4,  5,  1,  3,  7,
 3,  7,  5,  1,  0,  4,  6,  2,
 5,  7,  3,  1,  0,  2,  6,  4,
 6,  2,  3,  7,  5,  1,  0,  4,
 6,  4,  0,  2,  3,  1,  5,  7,
 6,  4,  5,  7,  3,  1,  0,  2 };

static unsigned const state3d[] = { /* nkey to 3 dimension state transitions */
    1,  2,  0,  3,  4,  0,  5,  6,
    0,  7,  1,  8,  5,  1,  4,  9,
   15,  0,  2, 22, 20,  2, 19, 23,
   20,  6,  3, 23, 15,  3, 16, 22,
   22, 13,  4, 12, 11,  4,  1, 20,
   11, 19,  5, 20, 22,  5,  0, 12,
    9,  3,  6,  2, 21,  6, 17,  0,
   10,  1,  7, 11, 12,  7, 13, 14,
   12,  9,  8, 14, 10,  8, 18, 11,
    6,  8,  9,  7, 17,  9, 21,  1,
    7, 15, 10, 16, 13, 10, 12, 17,
    5, 14, 11,  9,  0, 11, 22,  8,
    8, 20, 12, 19, 18, 12, 10,  5,
   18,  4, 13,  5,  8, 13,  7, 19,
   17, 11, 14,  1,  6, 14, 23,  7,
    2, 10, 15, 18, 19, 15, 20, 21,
   19, 17, 16, 21,  2, 16,  3, 18,
   14, 16, 17, 15, 23, 17,  6, 10,
   13, 21, 18, 17,  7, 18,  8, 16,
   16,  5, 19,  4,  3, 19,  2, 13,
    3, 12, 20, 13, 16, 20, 15,  4,
   23, 18, 21, 10, 14, 21,  9, 15,
    4, 23, 22,  6,  1, 22, 11,  3,
   21, 22, 23,  0,  9, 23, 14,  2 };

static unsigned const idata3d[] = {   /* 3 dimension to nkey conversion */
 0,  7,  3,  4,  1,  6,  2,  5,
 0,  1,  3,  2,  7,  6,  4,  5,
 0,  3,  7,  4,  1,  2,  6,  5,
 2,  3,  5,  4,  1,  0,  6,  7,
 4,  5,  3,  2,  7,  6,  0,  1,
 4,  7,  3,  0,  5,  6,  2,  1,
 6,  7,  5,  4,  1,  0,  2,  3,
 0,  1,  7,  6,  3,  2,  4,  5,
 2,  1,  5,  6,  3,  0,  4,  7,
 6,  1,  5,  2,  7,  0,  4,  3,
 0,  7,  1,  6,  3,  4,  2,  5,
 2,  1,  3,  0,  5,  6,  4,  7,
 4,  7,  5,  6,  3,  0,  2,  1,
 4,  5,  7,  6,  3,  2,  0,  1,
 6,  1,  7,  0,  5,  2,  4,  3,
 0,  3,  1,  2,  7,  4,  6,  5,
 2,  3,  1,  0,  5,  4,  6,  7,
 6,  7,  1,  0,  5,  4,  2,  3,
 2,  5,  1,  6,  3,  4,  0,  7,
 4,  3,  7,  0,  5,  2,  6,  1,
 4,  3,  5,  2,  7,  0,  6,  1,
 6,  5,  1,  2,  7,  4,  0,  3,
 2,  5,  3,  4,  1,  6,  0,  7,
 6,  5,  7,  4,  1,  2,  0,  3 };

static unsigned const istate3d[] = { /* 3 dimension to nkey state transitions */
 1,  6,  3,  4,  2,  5,  0,  0,
 0,  7,  8,  1,  9,  4,  5,  1,
15, 22, 23, 20,  0,  2, 19,  2,
 3, 23,  3, 15,  6, 20, 16, 22,
11,  4, 12,  4, 20,  1, 22, 13,
22, 12, 20, 11,  5,  0,  5, 19,
17,  0,  6, 21,  3,  9,  6,  2,
10,  1, 14, 13, 11,  7, 12,  7,
 8,  9,  8, 18, 14, 12, 10, 11,
21,  8,  9,  9,  1,  6, 17,  7,
 7, 17, 15, 12, 16, 13, 10, 10,
11, 14,  9,  5, 11, 22,  0,  8,
18,  5, 12, 10, 19,  8, 12, 20,
 8, 13, 19,  7,  5, 13, 18,  4,
23, 11,  7, 17, 14, 14,  6,  1,
 2, 18, 10, 15, 21, 19, 20, 15,
16, 21, 17, 19, 16,  2,  3, 18,
 6, 10, 16, 14, 17, 23, 17, 15,
18, 18, 21,  8, 17,  7, 13, 16,
 3,  4, 13, 16, 19, 19,  2,  5,
16, 13, 20, 20,  4,  3, 15, 12,
 9, 21, 18, 21, 15, 14, 23, 10,
22, 22,  6,  1, 23, 11,  4,  3,
14, 23,  2,  9, 22, 23, 21,  0 };

/* Given a 1-d coordinate in [0,1], returns it as the Hilbert key) */
inline double Zoltan_HSFC_InvHilbert1d(double* coord)
{
    /* sanity check for input arguments */
    if (coord[0] < 0.0)
        printf("%s\n", "Spatial Coordinates out of range.");

    return *coord;
}



/* Given x,y coordinates in [0,1]x[0,1], returns the Hilbert key [0,1] */
inline double Zoltan_HSFC_InvHilbert2d(double* coord)
{
    static const unsigned* d[] = { idata2d,  idata2d + 4, idata2d + 8, idata2d + 12 };
    static const unsigned* s[] = { istate2d, istate2d + 4, istate2d + 8, istate2d + 12 };

    int level;
    unsigned int key[2], c[2], temp, state;
    const int MAXLEVEL = 28; /* 56 bits of significance, 28 per dimension */

    /* sanity check for input arguments */
    if ((coord[0] < 0.0) || (coord[0] > 1.0) || (coord[1] < 0.0)
        || (coord[1] > 1.0))
        printf("%s\n", "Spatial Coordinates out of range.");

    /* convert x,y coordinates to integers in range [0, IMAX] */
    c[0] = (unsigned int)(coord[0] * (double)IMAX);               /* x */
    c[1] = (unsigned int)(coord[1] * (double)IMAX);               /* y */

    /* use state tables to convert nested quadrant's coordinates level by level */
    key[0] = key[1] = 0;
    state = 0;
    for (level = 0; level < MAXLEVEL; level++) {
        temp = ((c[0] >> (30 - level)) & 2)    /* extract 2 bits at current level */
            | ((c[1] >> (31 - level)) & 1);

        /* treat key[] as long shift register, shift in converted coordinate */
        key[0] = (key[0] << 2) | (key[1] >> 30);
        key[1] = (key[1] << 2) | *(d[state] + temp);

        state = *(s[state] + temp);
    }

    /* convert 2 part Hilbert key to double and return */
    return ldexp((double)key[0], -24) + ldexp((double)key[1], -56);
}

/* Given x,y,z coordinates in [0,1]x[0,1]x[0,1], returns Hilbert key in [0,1] */
inline double Zoltan_HSFC_InvHilbert3d(double* coord)
{
    static const unsigned int* d[] =
    { idata3d,      idata3d + 8,   idata3d + 16,  idata3d + 24,
     idata3d + 32,  idata3d + 40,  idata3d + 48,  idata3d + 56,
     idata3d + 64,  idata3d + 72,  idata3d + 80,  idata3d + 88,
     idata3d + 96,  idata3d + 104, idata3d + 112, idata3d + 120,
     idata3d + 128, idata3d + 136, idata3d + 144, idata3d + 152,
     idata3d + 160, idata3d + 168, idata3d + 176, idata3d + 184 };

    static const unsigned int* s[] =
    { istate3d,      istate3d + 8,   istate3d + 16,  istate3d + 24,
     istate3d + 32,  istate3d + 40,  istate3d + 48,  istate3d + 56,
     istate3d + 64,  istate3d + 72,  istate3d + 80,  istate3d + 88,
     istate3d + 96,  istate3d + 104, istate3d + 112, istate3d + 120,
     istate3d + 128, istate3d + 136, istate3d + 144, istate3d + 152,
     istate3d + 160, istate3d + 168, istate3d + 176, istate3d + 184 };

    int level;
    unsigned int key[2], c[3], temp, state;
    const int MAXLEVEL = 19; /* 56 bits of significance, 18+ per dimension */
    //char* yo = "Zoltan_HSFC_InvHilbert3d";

    /* sanity check for input arguments */
    if ((coord[0] < 0.0) || (coord[0] > 1.0) || (coord[1] < 0.0)
        || (coord[1] > 1.0) || (coord[2] < 0.0) || (coord[2] > 1.0))
       printf("%s\n", "Spatial Coordinates out of range.");

    /* convert x,y,z coordinates to integers in range [0,IMAX] */
    c[0] = (unsigned int)(coord[0] * (double)IMAX);     /* x */
    c[1] = (unsigned int)(coord[1] * (double)IMAX);     /* y */
    c[2] = (unsigned int)(coord[2] * (double)IMAX);     /* z */

    /* use state tables to convert nested quadrant's coordinates level by level */
    key[0] = key[1] = 0;
    state = 0;
    for (level = 0; level < MAXLEVEL; level++) {
        temp = ((c[0] >> (29 - level)) & 4)  /* extract 3 bits at current level */
            | ((c[1] >> (30 - level)) & 2)
            | ((c[2] >> (31 - level)) & 1);

        /* treat key[] as long shift register, shift in converted coordinate */
        key[0] = (key[0] << 3) | (key[1] >> 29);
        key[1] = (key[1] << 3) | *(d[state] + temp);

        state = *(s[state] + temp);
    }

    /* convert 2 part Hilbert key to double and return */
    return ldexp((double)key[0], -25) + ldexp((double)key[1], -57);
}

/* Sorting pointers in decreasing order. Criteria (key) is double. */
inline static void quickpart_pointer_dec_double(
    int* sorted, double* val, int start, int end, int* equal, int* smaller
)
{
    int   i, next;
    double key = (val ? val[sorted[(end + start) / 2]] : 1.0);

    *equal = *smaller = start;
    for (i = start; i <= end; i++) {
        next = sorted[i];
        if ((val ? val[next] : 1.0) > key) {
            sorted[i] = sorted[*smaller];
            sorted[(*smaller)++] = sorted[*equal];
            sorted[(*equal)++] = next;
        }
        else if ((val ? val[next] : 1.0) == key) {
            sorted[i] = sorted[*smaller];
            sorted[(*smaller)++] = next;
        }
    }
}

inline void Zoltan_quicksort_pointer_dec_double(
    int* sorted, double* val, int start, int end
)
{
    int  equal, smaller;

    if (start < end) {
        quickpart_pointer_dec_double(sorted, val, start, end, &equal, &smaller);
        Zoltan_quicksort_pointer_dec_double(sorted, val, start, equal - 1);
        Zoltan_quicksort_pointer_dec_double(sorted, val, smaller, end);
    }
}

/* Sort in increasing order by first calling the decreasing sort,
   then reverse the order in linear time. */
inline void Zoltan_quicksort_pointer_inc_double(
    int* sorted, double* val, int start, int end
)
{
    int i, j;
    double temp;

    /* sort in decreasing order */
    Zoltan_quicksort_pointer_dec_double(sorted, val, start, end);
    /* reverse order */
    for (i = start, j = end; i < j; i++, j--) {
        temp = sorted[i];
        sorted[i] = sorted[j];
        sorted[j] = temp;
    }
}



