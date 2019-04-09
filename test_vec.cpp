#include<iostream>
#include<vector>

void alter_vec(std::vector<int> v){
	v[0]++;
}

int main(){

	std::vector<int> v1 (4, 1);

	std::vector<int> v2 {};

	v2 = v1;

	v2[0]++;

	for(auto elem: v1)
		std::cout << elem << '\t';
	std::cout << '\n';

	for(auto elem: v2)
		std::cout << elem << '\t';
	std::cout << '\n';

	alter_vec(v2);

	for(auto elem: v1)
		std::cout << elem << '\t';
	std::cout << '\n';

	for(auto elem: v2)
		std::cout << elem << '\t';
	std::cout << '\n';


}
