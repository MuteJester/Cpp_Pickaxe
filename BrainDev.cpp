#include "CPPB_Lib.h"
#include <iostream>


int main() {
	Brain bt;
	std::vector<Coordinate> s,result;
	s.push_back(Coordinate{ 1,2,3 });
	s.push_back(Coordinate{ 1,2,3 });
	s.push_back(Coordinate{ 1,2,5 });
	s.push_back(Coordinate{ 2,2,3 });

	result = bt.K_Means<Coordinate>(s, 2, 100,2);

	for (auto i : result) {
		std::cout << "x: " <<i.x<<" y: "<<i.y << " z: "<<i.z << std::endl;
	}

	std::cout << "\n\n";
	system("pause");
	return 0;
}