#include<iostream>
#include<vector>

int main(){

	std::vector<double*> ptrs (4, nullptr);
	for(auto i = 0; i < 4; i++)
		ptrs[i] = new double [i + 10];

	double *ptr = new double[20'000'000];

	ptrs[0] = ptr;

	for(auto p: ptrs)
		delete[] p;
	

}

