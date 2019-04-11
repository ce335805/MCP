# MCP
A program to solve the matrix chain problem is provided.
As a cost metric the number of FLOPS needed to compute the chain is given.
Any other cost function to multiply two matrices can be included by just passing the corresponding functin pointer.

The result for the fastest matrix chain is then compared to the actual best strategy to compute the matrix chain via a call to DGEMM.
The result might differ significantly, even for a perfect cost metric of mulitplying two matrices - due to caching of intermediate results.

#build
The project can be built using standard cmake
One needs to link to a BLAS implementation (eg. Intel's MKL)
This needs some modificatin in CMakeLists.txt and MY_MCP.cpp
