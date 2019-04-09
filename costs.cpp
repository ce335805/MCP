#include <algorithm>
#include <assert.h>
//#include "mkl.h"
#include <cblas.h>
#include <chrono>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

#include "costs.h"

//cost function that computes the actual time it takes dgemm to perform the multiplication
double WALL_TIME(const size_t d1, const size_t d2, const size_t d3)
{
    //the time should be computed as an average - to have a more reliable result

    double wall_time{}; //return value
    int reps{ 100 }; //repetitions for average

    //create three matrices
    std::vector<double> M1(d1 * d2, 0.0);
    std::vector<double> M2(d2 * d3, 0.0);
    std::vector<double> M3(d1 * d3, 0.0);

    //same seed to always use the same matricies
    std::mt19937_64 generator(42);
    std::uniform_real_distribution<double> uniform(0.0, 100.0);

    //create random matrices - lambda must be mutable since the state of the generator changes
    //one could move this into the for-loop to average over different matricies
    //this however does not change execution time
    std::generate(M1.begin(), M1.end(), [uniform, generator]() mutable -> double { return uniform(generator); });
    std::generate(M2.begin(), M2.end(), [uniform, generator]() mutable -> double { return uniform(generator); });

    for (auto i = 0; i < reps; ++i) {

        //actual computation + timing
        auto start = std::chrono::steady_clock::now();
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, d1, d3, d2, 1.0, &M1[0], d2, &M2[0], d3, 0.0, &M3[0], d3);
        auto end = std::chrono::steady_clock::now();

        std::chrono::duration<double, std::milli> dur = end - start;
        wall_time += dur.count();
    }

    return wall_time / static_cast<double>(reps);
}
