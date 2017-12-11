#ifndef _MEMORY_H_
#define _MEMORY_H_

#include "Type.h"
#include <cstdlib>

template <typename TYPE>
void allocateArray1D(TYPE* &array, int n1)
{
    array = new TYPE [n1];
}

template <typename TYPE>
void destroy(TYPE *array)
{
    if (array == NULL) return;
    delete[] array;
}

template <typename TYPE>
void allocateMatrix2D(TYPE **&array, int n1, int n2)
{
    TYPE* data = new TYPE[n1*n2];
    array = new TYPE* [n1];

    bigint n = 0;
    for (int i = 0; i < n1; i++)
    {
        array[i] = &data[n];
        n += n2;
    }
}

template <typename TYPE>
void destroy(TYPE **array)
{
    if (array == NULL) return;
    delete[] array[0];
    delete[] array;
}

#endif
