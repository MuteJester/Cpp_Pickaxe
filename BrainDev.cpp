#define PRINT_EPOCH_STATUS

#include "CPPB_Lib.h"
#include "Cpp_SIPL.h"

#include <vector>
using namespace std;

int main() {
	CSV_File data;
	data.Load_CSV("data.csv");
	cout << data.columns;

	//vector<double> X =  data["A Data"].get_Numeric_Values();
	//vector<double> Y = data["B Data"].get_Numeric_Values();
	Image Av;
	Av.Load_Image("Images/fruits-360/Test/Avocado/4_100.jpg");

	Av.Set_Scale(28, 28);

	Av.Grayscale(1);
	Av.Update_Pixel_Matrix();

	Matrix<double> dd = Av.Flatten();

	Training_Data training;
	training.Load_Data(dd[0], dd[0]);

	std::cout << "Loaded\n";

	Neural_Net encoder;
	encoder.Add_Layer(784);
	encoder.Add_Layer(500,"Sigmoid");
	encoder.Add_Layer(256, "Sigmoid");
	encoder.Add_Layer(100, "Sigmoid");
	encoder.Add_Layer(256, "Sigmoid");
	encoder.Add_Layer(500, "Sigmoid");
	encoder.Add_Layer(784, "Sigmoid");




	encoder.Compile(0.3, 0, 1, training);


	encoder.Start_Training(5);

	encoder.Print_Output_Neurons_Values();
	std::vector<double> vec = encoder.Get_Output_Values();
	std::vector<std::vector<double> > res = To_2D(vec, 28);

	Image result = Matrix_To_Graysacle(res);

	result.Write_Image("x.jpg");
	system("x.jpg");
	

	
	cout << "\n\n";
	system("Pause");
	return 0;
}