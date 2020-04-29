#include "CPPB_Lib.h"
using namespace std;

int main() {
	CSV_File data;
	data.Load_CSV("G");
	String_Vector tx{ "1","2" };

	data.Add_Row(tx);
	data.Add_Row();
	data.Set_Value_At(1, data.Number_Of_Rows,"NewADD");
	data.Write_CSV("G_1");




	cout << "\n\n";
	system("Pause");
	return 0;
}