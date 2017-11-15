#ifndef _RANDOM_H_
#define _RANDOM_H_

#include <Type.h>
#include <cstdlib>

namespace HamiltonSpace {

inline HS_float randomGauss()
{
    static HS_float n2 = 0.0;
    static int n2_cached = 0;
    HS_float mean = 0.0;
    HS_float stddev = 1.0;
    if (!n2_cached)
    {
        HS_float x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
        }
        while (r > 1.0);
        {
            HS_float d = sqrt(-2.0*log(r)/r);
            HS_float n1 = x*d;
            n2 = y*d;
            HS_float result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    } 
}

}

#endif
