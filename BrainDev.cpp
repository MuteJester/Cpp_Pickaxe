#include "CPPB_Lib.h"
#include <vector>
using namespace std;

int main() {
	CSV_File data;
	data.Load_CSV("data.csv");

	Pickaxe ml(data);
	
	CSV_File sp = data.Split_Data(10);
	cout << sp;

	cout << "\n\n";
	system("Pause");
	return 0;
}