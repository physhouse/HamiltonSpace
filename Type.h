#ifndef _TYPE_H_
#define _TYPE_H_

struct double2 {
    double x, y;
}
struct float2 {
    float x, y;
}

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

#endif
