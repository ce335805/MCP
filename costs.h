// header file which declares cost functions
#ifndef COSTS_H //inlcude guard
#define COSTS_H

#include <cstddef>//just to suppress the error that size_t is not known

//number of flops performed by a matrix matrix multiplication
unsigned long int constexpr FLOPS(const size_t d1, const size_t d2, const size_t d3) noexcept
{
    return d1 * d2 * d3;
}


//negative number of flops to find worst solution
int constexpr NEG_FLOPS(const size_t d1, const size_t d2, const size_t d3) noexcept
{
    return (-1) * static_cast<int>(d1 * d2 * d3);
}


//cost function that computes the actual time it takes dgemm to perform the multiplication
double WALL_TIME(const size_t d1, const size_t d2, const size_t d3);

#endif // costs.h

