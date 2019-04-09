#include<iostream>
#include<vector>
#include<algorithm>

int main(){

	std::vector<int> v {1,2,3,4,5,6,7,8,9,5,11,12,13,14,15};

	//int i = 0;

	auto lam = [&i=0](const int &elem)mutable{
			std::cout << "i = " << i << '\n';
			++i;
		  if(elem == 5) return true;
			else{
				return false;
				} 
			};


	std::remove_if(v.begin(), v.end(), lam);

	for(auto elem: v)
		std::cout << elem << '\t';
	std::cout << '\n';



}
