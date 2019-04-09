// header file which declares cost functions
#ifndef UTILS_H //inlcude guard
#define UTILS_H
#include <assert.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <limits>


//function to print the found dimensions vectors to a file
void print_to_file(const double diff, const std::vector<size_t> dimensions, std::ofstream& output);

// just a convinient function for printing the output matricies
template <typename T>
void print_matrix(std::vector<T> to_print)
{

    //check that the passed value is a number somehow - thus function can be noexcept
    assert(std::is_arithmetic<T>::value);

    size_t N = static_cast<size_t>(std::sqrt(to_print.size()));

    std::cout << '\n';

    for (auto i = 0ul; i < N; ++i) {
        for (auto j = 0ul; j < N; ++j) {
            std::cout << std::setw(20) << to_print[i * N + j] << '\t';
        }
        std::cout << '\n';
    }

    std::cout << '\n';
}


//dont know why there is no factorial in std library - anyways here is one
template <typename T>
T fac(T n)
{
    T ret = 1;
    std::numeric_limits<T> lims;
    while (n > 0) {
				//just since it is very easy to pass too big values to factorial
        if (lims.max() / n < ret)
            return 0;
        else {
            ret = n * ret;
            n--;
        }
    }
		return ret;
}

#endif // utils.h
