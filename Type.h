#ifndef _TYPE_H_
#define _TYPE_H_

struct double2 {
    double x, y;
};

struct float2 {
    float x, y;
};

#ifndef PRECISION
#define PRECISION 2
#endif
#if PRECISION==1
typedef float HS_float;
typedef float2 HS_float2;
#else
typedef double HS_float;
typedef double2 HS_float2;
#endif

typedef int HS_int;
typedef int HS_bigint;

#define INFINITY 1e30
#define VERYSMALL 1e-30

#define MAX_ARRAY 100000
#define EXCHANGEMAX 1000
#define BUFFMAX 100000

#define PBC_ANY_FLAG 0x1
#define PBC_POS_X 0x2
#define PBC_NEG_X 0x4
#define PBC_POS_Y 0x8
#define PBC_NEG_Y 0x10
#define PBC_POS_Z 0x20
#define PBC_NEG_Z 0x40

#endif
