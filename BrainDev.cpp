#include "CPPB_Lib.h"
using namespace std;

int main() {
	CSV_File data;
	data.Load_CSV("G");
	
	Matrix<double> x(3,IDENTITY_MATRIX);
	x[0][1] = 9;
	x[2][1] = 9;
	x.Random_Fill();
	cout << x << "\n\n";

	vector<Matrix<double> > lu = x.LU_Decomposition();

	cout << lu[0];


	cout << "\n\n";
	system("Pause");
	return 0;
}