#include <algorithm>
#include <assert.h>
#include "mkl.h"
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "costs.h"
#include "utils.h"

//computes time to evaluate a certain pranthesization
double compute_MC(std::vector<size_t> products, std::vector<size_t> dimensions);

//function to generate a vector of products for a given splits matrix - shift needed for recursion
std::vector<size_t> gen_prods(const std::vector<size_t> splits, const size_t shift);

//function to calculate the best paranthesization actually evaluating the matrix chain
std::vector<size_t> find_best_MC(const std::vector<size_t> dimensions);

//function to compute cost of a given paranthesization and cost function
template <typename T>
T compute_cost(const std::vector<size_t> dimensions, std::vector<size_t> splits, size_t i, size_t j,
    std::function<T(size_t, size_t, size_t)> const_func);

//function to compare two cost metrics for a given dimensions vector
template <typename T, typename S>
double compare_metrics(const std::vector<size_t> dimensions,
    std::function<T(size_t, size_t, size_t)> cost_func_1,
    std::function<S(size_t, size_t, size_t)> cost_func_2);

//length holds sizes of matrices
//costs(i, j) will hold the cost to evaluate the subchains A_i ... A_j
//splits holds the position where to split the chain to also save the actual solution
template <typename T>
void matrix_chain(const std::vector<size_t> dimensions, std::vector<T>& costs, std::vector<size_t>& splits,
    std::function<T(size_t, size_t, size_t)> cost_func)
{
    //check that the deduced type is a number somehow
    assert(std::is_arithmetic<T>::value);

    //number of matricies in the chain
    auto n = dimensions.size() - 1;

    //one could save only a triangular part - but this is additional work and n will be negligably small
    assert(costs.size() == n * n);
    assert(splits.size() == n * n);

    std::numeric_limits<T> lims;

    //initialize subchains of size 0 with 0 - otherswise with the maximum T
    std::generate(costs.begin(), costs.end(),
        [n, lims, i = 0ul]() mutable -> T {
            if (i / n == i % n) {
                //default initialization for arithmetic types is 0
                T zero{};
                ++i;
                return zero;
            } else
                ++i;
            return lims.max();
        });

    //on diagonal nothing can be split - else initalize at 0
    std::generate(splits.begin(), splits.end(),
        [n, i = 0ul]() mutable -> unsigned long int {
            if (i / n == i % n) {
                ++i;
                return (i - 1) % n;
            } else {
                ++i;
                return 0ul;
            }
        });

    //loop over dimensions of subchain
    for (auto l = 2ul; l < n + 1; ++l) {
        //left dim of first element in subchain when l is indexed
        for (auto i = 0ul; i < n - l + 1; ++i) {
            //number of the last matrix in the subchain - corresponds to left dimension
            auto j = i + l - 1;
            //number of matrix left of the splitup - corresponds to left dimension
            for (auto k = i; k < j; ++k) {
                //add cost of multiplying subchains
                auto temp_cost = cost_func(dimensions[i], dimensions[k + 1], dimensions[j + 1]);
                //add cost of computing each individual subchain
                temp_cost += costs[i * n + k] + costs[n * (k + 1) + j];
                //if the cost is actually lower - set entry in costs array
                if (temp_cost < costs[i * n + j]) {
                    costs[i * n + j] = temp_cost;
                    splits[i * n + j] = k;
                }
            }
        }
    }
}

int main()
{

    //time the program as a whole
    auto start = std::chrono::steady_clock::now();

    //for a test - check matrix chains of length 4 for differences between matrix chain algo to real life

    int len = 4;
    std::vector<size_t> dimensions(5, {});
    std::vector<size_t> costs(len * len, {});
    std::vector<size_t> splits(len * len, {});
    std::function<size_t(size_t, size_t, size_t)> flops = FLOPS;
    std::vector<size_t> prods_1{};
    std::vector<size_t> prods_2{};

    prods_1 = std::vector<size_t>{ 0, 1, 2 };
    dimensions = std::vector<size_t>{ 3, 5, 3, 5, 3 };
    //compute_MC(prods_1, dimensions);

    auto end = std::chrono::steady_clock::now();

    std::chrono::duration<double, std::milli> dur = end - start;

    double time = dur.count();
    std::cout << '\n';
    std::cout << "Program took : " << time / 1000 << "s" << '\n';
}

//function to calculate the best paranthesization actually evaluating the matrix chain
std::vector<size_t> find_best_MC(const std::vector<size_t> dimensions)
{

    std::vector<size_t> prod_ret(dimensions.size() - 2, {});
    std::vector<size_t> perm(dimensions.size() - 2, {});
    size_t k = 0;
    std::generate(perm.begin(), perm.end(), [&k]() { return k++; });

    std::numeric_limits<double> lims;
    double time{ lims.max() };

    for (auto i = 0ul; i < fac(dimensions.size() - 2); ++i) {

        //        for (auto elem : perm)
        //            std::cout << elem << '\t';
        //        std::cout << '\n';
        double temp_time{};
        for (auto rep = 0; rep < 1; ++rep)
            temp_time += compute_MC(perm, dimensions);
        std::cout << "Done with matrix chain" << '\n';

        if (temp_time < time) {
            std::copy(perm.begin(), perm.end(), prod_ret.begin());
            time = temp_time;
        }

        std::next_permutation(perm.begin(), perm.end());
    }
    return prod_ret;
}

//computes time to evaluate a certain pranthesization
double compute_MC(std::vector<size_t> products, std::vector<size_t> dimensions)
{
    assert(products.size() == dimensions.size() - 2);
    double time{};

    //first generate random matricies
    std::mt19937_64 generator(42);
    std::uniform_real_distribution<double> uniform(0.0, 100.0);

    //convert products vector in something convenient to evaluate the chain
    for (auto i = 0ul; i < products.size(); ++i) {
        for (auto j = i + 1; j < products.size(); ++j) {
            if (products[j] >= products[i])
                products[j]--;
        }
    }

    assert(products.back() == 0);

    //each entry points to one of the matricies
    std::vector<std::shared_ptr<double>> matricies(dimensions.size() - 1, nullptr);

    for (auto i = 0ul; i < matricies.size(); ++i) {
        //allocate memory
        matricies[i] = std::shared_ptr<double>(new double[dimensions[i] * dimensions[i + 1]]);
        for (auto k = 0; k < dimensions[i] * dimensions[i + 1]; ++k) {
            //initialize random matricies
            matricies[i].get()[k] = uniform(generator);
            //matricies[i].get()[k] = 1.0;
        }
    }

    std::vector<std::shared_ptr<double>> temps(products.size(), nullptr);

    //allocate temporarys in advance
    for (auto i = 0ul; i < products.size(); ++i)
        temps[i] = std::shared_ptr<double>(new double[dimensions[products[i]] * dimensions[products[i] + 2]]);

    //entry in products always stands for the left to be multilied matrix
    for (auto i = 0ul; i < products.size(); ++i) {
        std::cout << "in loop iter : " << i << '\n';
        auto d1 = dimensions[products[i]];
        auto d2 = dimensions[products[i] + 1];
        auto d3 = dimensions[products[i] + 2];

        auto start = std::chrono::steady_clock::now();
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, d1, d3, d2, 1.0, matricies[products[i]].get(), d2, matricies[products[i] + 1].get(), d3, 0.0, temps[i].get(), d3);
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double, std::milli> dur = end - start;
        time += dur.count();

        //update list of matricies
        //release memory of not needed pointers
        //        delete[] matricies[products[i]];
        //        delete[] matricies[products[i] + 1];
        matricies[products[i]] = temps[i];
        matricies.erase(matricies.begin() + products[i] + 1);
        dimensions.erase(dimensions.begin() + products[i] + 1);
    }

    //release memory
    //    for(auto ptr: matricies)
    //        delete[] ptr;
    std::cout << "returning from compute MC" << '\n';
    return time;
}

//function to generate a vector of products for a given splits matrix - shift needed for recursion
std::vector<size_t> gen_prods(const std::vector<size_t> splits, const size_t shift)
{

    size_t n = static_cast<size_t>(std::sqrt(splits.size()));

    if (n == 1)
        return std::vector<size_t>{};
    if (n == 2)
        return std::vector<size_t>{ splits[1] };

    //place where to split the matrix chain
    size_t split_elem = splits[n - 1] - shift;

    //create split matricies for the two sub chains
    std::vector<size_t> sub1((split_elem + 1) * (split_elem + 1), {});
    std::vector<size_t> sub2((n - (split_elem + 1)) * (n - (split_elem + 1)), {});

    //copy sub matricies from original splits matrix
    auto k = 0ul;
    std::copy_if(splits.begin(), splits.end(), sub1.begin(),
        [&k, split_elem, n](const size_t& elem) {
            ++k;
            if ((k - 1) / n <= split_elem && (k - 1) % n <= split_elem)
                return true;
            else {
                return false;
            }
        });

    k = 0;
    std::copy_if(splits.begin(), splits.end(), sub2.begin(),
        [&k, split_elem, n](const size_t& elem) {
            ++k;
            if ((k - 1) / n > split_elem && (k - 1) % n > split_elem)
                return true;
            else {
                return false;
            }
        });

    //calculate splits matricies from sub_matricies
    std::vector<size_t> sub_prods_1 = gen_prods(sub1, shift);
    std::vector<size_t> sub_prods_2 = gen_prods(sub2, shift + split_elem + 1);
    std::vector<size_t> prods{};
    //insert this into return value
    prods.insert(prods.end(), sub_prods_1.begin(), sub_prods_1.end());
    prods.insert(prods.end(), sub_prods_2.begin(), sub_prods_2.end());
    //at the end push back the original split which corresponds to last to be evaluated product
    prods.push_back(split_elem + shift);

    return prods;
}

//function to compute cost of a given paranthesization and cost function
template <typename T>
T compute_cost(const std::vector<size_t> dimensions, std::vector<size_t> splits, size_t i, size_t j,
    std::function<T(size_t, size_t, size_t)> cost_func)
{
    //since i is a left dimension and j a right dimension they should have a minimum seperation of 1
    assert(i < j);
    assert(splits.size() == (dimensions.size() - 1) * (dimensions.size() - 1));
    //check that the deduced type is a number somehow
    assert(std::is_arithmetic<T>::value);

    if (i == j - 1)
        return 0;
    else {
        return cost_func(dimensions[i], dimensions[splits[(dimensions.size() - 1) * i + j - 1] + 1], dimensions[j]) + compute_cost(dimensions, splits, i, splits[(dimensions.size() - 1) * i + j - 1] + 1, cost_func) + compute_cost(dimensions, splits, splits[(dimensions.size() - 1) * i + j - 1] + 1, j, cost_func);
    }
}

//function to compare two cost metrics for a given dimensions vector
template <typename T, typename S>
double compare_metrics(const std::vector<size_t> dimensions,
    std::function<T(size_t, size_t, size_t)> cost_func_1,
    std::function<S(size_t, size_t, size_t)> cost_func_2)
{

    //check that the deduced type is a number somehow
    assert(std::is_arithmetic<T>::value);
    assert(std::is_arithmetic<S>::value);

    //store results for matrix chain with two different metrics
    std::vector<T> costs_1((dimensions.size() - 1) * (dimensions.size() - 1), {});
    std::vector<size_t> splits_1((dimensions.size() - 1) * (dimensions.size() - 1), {});
    std::vector<S> costs_2((dimensions.size() - 1) * (dimensions.size() - 1), {});
    std::vector<size_t> splits_2((dimensions.size() - 1) * (dimensions.size() - 1), {});

    //call matrix chain algorithm for the two cost metrics
    matrix_chain(dimensions, costs_1, splits_1, cost_func_1);
    matrix_chain(dimensions, costs_2, splits_2, cost_func_2);

    //if the splits are the same the case is uninteresting
    //note that they can be different and still the matrix chain problem has found the same solution
    //if (std::equal(splits_1.begin(), splits_1.end(), splits_2.begin(), splits_2.end()))
    //    return 0.0;

    // ----- note - comparison only for second passed cost function!!! - this should thus be ex-time for project

    //previously calculated costs - just to save typing this is actually unnecessary
    S cost = costs_2[dimensions.size() - 2];

    //compare how the solution would perform the other cost_metric
    S cost_comp = compute_cost(dimensions, splits_1, 0, dimensions.size() - 1, cost_func_2);

    //if this cast fails the -compiler- will complain
    return static_cast<double>((cost_comp - cost) / cost);
}
