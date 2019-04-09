#include <assert.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <limits>

#include "utils.h"

//function to print the found dimensions vectors to a file
void print_to_file(const double diff, const std::vector<size_t> dimensions, std::ofstream& output)
{

    for (auto dim : dimensions)
        output << dim << '\t';
    output << diff << '\t';
    output << '\n';
}
