#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <random>
#include <limits>
#include <regex>
#include <cmath>


#define Mmax(a,b) ((a) > (b) ? (a) : (b))
#define Mmin(a,b) ((a) < (b) ? (a) : (b))
#define Mabs(x)  ( (x<0) ? -(x) : x )
#define IDENTITY_MATRIX 'i'
#define DEG_TO_RAD(A) A * PI/180
#define RAD_TO_DEG(A) A * 180/PI
#define SQUARE_IT(NUM) NUM*NUM
#define CUBE_IT(NUM) NUM*NUM*NUM

//typenames
typedef std::vector<std::string> String_Vector;

//general methods
String_Vector Split_String(std::string source, char delim) {
	String_Vector result;
	std::size_t current, previous = 0;
	current = source.find(delim);
	while (current != std::string::npos) {
		result.push_back(source.substr(previous, current - previous));
		previous = current + 1;
		current = source.find(delim, previous);
	}
	result.push_back(source.substr(previous, current - previous));
	return result;
}
template<class DType> std::ostream &operator<<(std::ostream &out, std::vector<DType> const &source) {
	for (unsigned i = 0; i < source.size(); i++) {
		if (i < source.size() - 1) {
			out << source[i] << ", ";
		}

		else {
			out << source[i] << "\n";
		}

	}
	return out;
}



class Random_Utilitis {
public:
	int Random_INT(int minimum_value,int maximum_value) {
		static std::random_device seed;
		static std::mt19937 random_number_generator(seed());
		std::uniform_int_distribution<size_t> indices(minimum_value, maximum_value-1);
		return indices(random_number_generator);
	}
	double Random_DOUBLE(double minimum_value, double maximum_value) {
		static std::random_device seed;
		static std::mt19937 random_number_generator(seed());
		std::uniform_real_distribution<> indices(minimum_value, maximum_value);
		return indices(random_number_generator);
	}
	
};


template<class MType> class Matrix {

public:
	Matrix(int Rows, int Cols) {
		Matrix_Body = std::vector<std::vector<MType> >(Rows, std::vector<MType>(Cols));
		for (int i = 0; i < Rows; i++) {
			for (int j = 0; j < Cols; j++) {
				Matrix_Body[i][j] = 0;
			}
		}
		this->Rows = Rows;
		this->Cols = Cols;
	}
	Matrix(int N, char x = IDENTITY_MATRIX) {
		Matrix_Body = std::vector<std::vector<MType> >(N, std::vector<MType>(N));
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				Matrix_Body[i][j] = 0;
				if (i == j) {
					Matrix_Body[i][j] = 1;
				}
			}
		}
		this->Rows = N;
		this->Cols = N;
	}
	Matrix() {
		this->Rows = 0;
		this->Cols = 0;
	}
	template<class MType> Matrix(Matrix<MType> copy) {
		this->Rows = copy.Rows;
		this->Cols = copy.Cols;
		for (int i = 0; i < this->Rows; i++) {
			this->Matrix_Body.push_back(std::vector<MType>(this->Cols));
			for (int j = 0; j < this->Cols; j++) {
				this->Matrix_Body[i][j] = copy[i][j];
			}
		}
	}

	~Matrix() {};
	int getRows() { return Rows; }
	int getCols() { return Cols; }
	void setRows(int &Rows) { this->Rows = Rows; }
	void setCols(int &Cols) { this->Cols = Cols; }
	std::vector<MType> &operator[](int a) {
		return Matrix_Body[a];
	}
	void Random_Fill();
	void operator=(Matrix<MType> B) {
		this->Rows = B.getRows();
		this->Cols = B.getCols();
		this->Matrix_Body = std::vector<std::vector<MType> >(Rows, std::vector<MType>(Cols));
		for (int i = 0; i < Rows; i++) {
			for (int j = 0; j < Cols; j++) {
				Matrix_Body[i][j] = B[i][j];
			}
		}

	}
	Matrix<MType> operator+(Matrix<MType> &B);
	Matrix<MType> operator-(Matrix<MType> &B);
	void operator/(int const &divide_by);
	Matrix<MType> Dot_Product(Matrix<MType> &B);
	Matrix<MType> Matrix_Transpose();
	void Multiply_By_Scalar(double const &scalar);
	Matrix<MType> Hadamard_Product(Matrix<MType> &Mul_By);
	Matrix<MType> Kronecker_Product(Matrix<MType> &Mul_By);
	void Horizontal_Matrix_Concatenation(Matrix<MType> &To_HConcat);
	void Convolve(Matrix<int> &Mask, int mask_h, int mask_w);
	std::vector<MType> Diagonal();
	Matrix<MType> Get_Upper_Matrix();
	Matrix<MType> Get_Lower_Matrix();
	std::vector<Matrix<MType> > LU_Decomposition();
	std::vector<Matrix<MType> >  QR_Decomposition();
	template<class MType> friend std::ostream &operator<<(std::ostream &out, Matrix<MType> const &mat);

protected:
	std::vector<std::vector<MType> > Matrix_Body;
	int Rows, Cols;

};



template<class MType> std::ostream &operator<<(std::ostream &out, Matrix<MType> const &mat) {
	out << std::left;
	for (int i = 0; i < mat.Rows; i++) {
		for (int j = 0; j < mat.Cols; j++) {
			if (j < mat.Cols - 1) {
				out << std::setw(9) << mat.Matrix_Body[i][j] << "|";

			}
			else {
				out << std::setw(9) << mat.Matrix_Body[i][j];

			}
		}
		out << "\n";
	}
	return out;
}
template<class MType> Matrix<MType> Matrix<MType>::Matrix_Transpose() {


	Matrix<MType> temp(Cols, Rows);



	for (int i = 0; i < Cols; i++) {
		for (int j = 0; j < Rows; j++) {
			temp[i][j] = Matrix_Body[j][i];
		}
	}
	return  temp;
}
template<class MType> void Matrix<MType>::Multiply_By_Scalar(double const &scalar) {
	for (int i = 0; i < Rows; i++) {
		for (int j = 0; j < Cols; j++) {
			Matrix_Body[i][j] *= scalar;
		}
	}
}
template<class MType> Matrix<MType> Matrix<MType>::operator+(Matrix<MType> &B) {
	if (Rows != B.getRows() || Cols != B.getCols()) {
		return;
	}
	else {
		Matrix<MType> Res(this->Rows, this->Cols);
		for (int i = 0; i < Rows; i++) {
			for (int j = 0; j < Cols; j++) {
				Res[i][j] = Matrix_Body[i][j] + B[i][j];
			}
		}
	}
}
template<class MType> Matrix<MType> Matrix<MType>::operator-(Matrix<MType> &B) {
	if (Rows != B.getRows() || Cols != B.getCols()) {
		return;
	}
	else {
		Matrix<MType> Res(this->Rows, this->Cols);
		for (int i = 0; i < Rows; i++) {
			for (int j = 0; j < Cols; j++) {
				Res[i][j] = Matrix_Body[i][j] - B[i][j];
			}
		}
	}
}
template<class MType> Matrix<MType> Matrix<MType>::Dot_Product(Matrix<MType> &B) {

	if (Cols != B.getRows()) {
		return;

	}
	else {
		MType sum = 0;
		Matrix<MType> temp(Rows, B.getCols());
		for (int i = 0; i < Rows; i++) {
			for (int j = 0; j < B.getCols(); j++) {
				for (int k = 0; k < Cols; k++) {
					sum += Matrix_Body[i][k] * B[k][j];

				}
				temp[i][j] = sum;
				sum = 0;
			}
		}
		return temp;
	}

}
template<class MType> Matrix<MType> Matrix<MType>::Hadamard_Product(Matrix<MType> &Mul_By) {

	if (Cols != Mul_By.getRows() || Rows != Mul_By.getRows()) {
		return;

	}
	else {
		Matrix<MType> Result(Rows, Cols);
		MType sum = 0;

		for (int i = 0; i < Rows; i++) {
			for (int j = 0; j < Mul_By.getCols(); j++) {
				Result[i][j] = Matrix_Body[i][j] * Mul_By.mat[i][j];
			}
		}
		return Result;

	}
}
template<class MType> Matrix<MType> Matrix<MType>::Kronecker_Product(Matrix<MType> &Mul_By) {
	Matrix<MType> Kronecker(Rows*Mul_By.getRows(), Cols*Mul_By.getCols());
	int startRow, startCol;
	for (int i = 0; i < Rows; i++) {
		for (int j = 0; j < Cols; j++) {
			startRow = i * Mul_By.getRows();
			startCol = j * Mul_By.getCols();
			for (int k = 0; k < Mul_By.getRows(); k++) {
				for (int l = 0; l < Mul_By.getCols(); l++) {
					Kronecker[startRow + k][startCol + l] = Matrix_Body[i][j] * Mul_By[k][l];
				}
			}
		}
	}
	return Kronecker;

}
template<class MType> void Matrix<MType>::Horizontal_Matrix_Concatenation(Matrix<MType> &To_HConcat) {
	if (this->Rows != To_HConcat.getRows())
		return;

	int  i, j, k, l = 0;
	Matrix<MType> ConcatH(Rows, Cols + To_HConcat.getCols());

	for (i = 0; i < Rows; i++) {
		for (j = 0; j < Cols; j++) {
			ConcatH[i][l] = Matrix_Body[i][j];
			l++;
		}
		for (k = 0; k < To_HConcat.getCols(); k++) {
			ConcatH[i][l] = To_HConcat[i][k];
			l++;
		}
		l = 0;
	}
	*this = ConcatH;
}
template<class MType> std::vector<MType> Matrix<MType>::Diagonal() {
	if (this->Cols != this->Rows) {
		std::cout << "Error Matrix Isnt Square (R != C)\n";
		return std::vector<MType>();
	}
	else {
		std::vector<MType> res(this->Cols);
		for (int i = 0; i < this->Cols; i++) {
			res[i] = this->Matrix_Body[i][i];
		}
		return res;
	}
	
}
template<class MType> void Matrix<MType>::Convolve(Matrix<int> &Mask, int mask_h, int mask_w) {

	Matrix<MType> sample(mask_h, mask_w);
	int accumulator = 0;
	for (int i = 0; i < this->Rows; i++) {
		for (int j = 0; j < this->Cols; j++) {
			for (int k = -(mask_h / 2); k < (mask_h / 2); k++) {
				for (int m = -(mask_w / 2); m < (mask_w / 2); m++) {

					if (i + k < 0 && j + m < 0) {
						sample[k + (mask_h / 2)][m + (mask_h / 2)] = Matrix_Body[this->Rows + k][this->Cols + m];

					}
					else if (i + k < 0 && j + m > 0) {
						sample[k + (mask_h / 2)][m + (mask_h / 2)] = Matrix_Body[this->Rows + k][j + m];

					}
					else if (j + m < 0 && i + k > 0) {
						sample[k + (mask_h / 2)][m + (mask_h / 2)] = Matrix_Body[i + k][this->Cols + m];

					}//till here before
					else if (i + k >= Rows && j + m >= Cols) {
						sample[k + (mask_h / 2)][m + (mask_h / 2)] = Matrix_Body[(i + k) - Rows][(j + m) - Cols];
					}
					else if (i + k >= Rows && j + m < Cols) {
						sample[k + (mask_h / 2)][m + (mask_h / 2)] = Matrix_Body[(i + k) - Rows][j + m];

					}
					else if (i + k < Rows && j + m >= Cols) {
						sample[k + (mask_h / 2)][m + (mask_h / 2)] = Matrix_Body[i + k][(j + m) - Cols];

					}
					else {
						sample[k + (mask_h / 2)][m + (mask_h / 2)] = Matrix_Body[i + k][j + m];

					}

				}
			}
			Mask.Matrix_Transpose();
			for (int k = 0; k < mask_h; k++) {
				for (int m = 0; m < mask_w; m++) {
					accumulator = +sample[k][m] * Mask[k][m];
				}
			}

			Matrix_Body[i][j] = accumulator;
			accumulator = 0;

		}
	}

}
template<class MType> void Matrix<MType>::operator/(int const &divide_by) {
	this->Multiply_By_Scalar((1.0 / divide_by));
}
template<class MType> Matrix<MType> Matrix<MType>::Get_Upper_Matrix() {
	if (this->Rows != this->Cols) {
		std::cout<<"Matrix Needs To Be NxN\n";
		return Matrix<MType>(1,1);
	}
	Matrix<MType> Upper(*this);
	for (int i = 0; i < this->Rows; i++) {
		for (int j = 0; j < i; j++) {
			Upper[i][j] = 0;
		}
	}
	return Upper;
}
template<class MType> Matrix<MType>  Matrix<MType>::Get_Lower_Matrix() {
	if (this->Rows != this->Cols) {
		std::cout<<"Matrix Needs To Be NxN\n";
		return Matrix<MType>(1,1);
	}
	Matrix Lower(*this);
	for (int i = 0; i < this->Rows; i++) {
		for (int j = i + 1; j < Cols; j++) {
			Lower.Matrix_Body[i][j] = 0;
		}
	}
	return Lower;
}
template<class MType> void Matrix<MType>::Random_Fill() {
	Random_Utilitis rand;
	for (int i = 0; i < this->Rows; i++) {
		for (int j = 0; j < this->Cols; j++) {
			this->Matrix_Body[i][j] = rand.Random_DOUBLE(0, 50);
		}
	}
}
template<class MType> std::vector<Matrix<MType> > Matrix<MType>::LU_Decomposition() {
	if (this->Cols != this->Rows) {
		std::cout << " Error  - Matrix Needs To Be NxN\n";
		return  std::vector<Matrix<MType> >();
	}
	Matrix<double> LU(*this);
	double *pivot = new double[this->Cols];
	for (int i = 0; i < this->Cols; i++) {
		pivot[i] = i;
	}
	int pivsign = 1;
	
	std::vector<MType> LUrowi;
	double *LUcolj = new double[this->Cols];
	//std::array<double, this->Cols> LUcolj;

	// Outer loop.

	for (int j = 0; j < this->Cols; j++) {
		for (int i = 0; i < this->Cols; i++) {
			LUcolj[i] = LU.Matrix_Body[i][j];
		}

		for (int i = 0; i < this->Cols; i++) {
			LUrowi = LU.Matrix_Body[i];
			int kmax = Mmin(i, j);
			double s = 0.0;
			for (int k = 0; k < kmax; k++) {
				s += LUrowi[k] * LUcolj[k];
			}

			LUrowi[j] = LUcolj[i] -= s;
		}

		int p = j;
		for (int i = j + 1; i < this->Cols; i++) {

			if (Mabs(LUcolj[i]) > Mabs(LUcolj[p])) {
				p = i;
			}
		}
		if (p != j) {
			for (int k = 0; k < this->Cols; k++) {
				double t = LU.Matrix_Body[p][k];
				LU.Matrix_Body[p][k] = LU.Matrix_Body[j][k];
				LU.Matrix_Body[j][k] = t;
			}
			double k = pivot[p];
			pivot[p] = pivot[j];
			pivot[j] = k;
			pivsign = -pivsign;
		}
		if ((j < this->Cols) & (LU[j][j] != 0.0)) {
			for (int i = j + 1; i < this->Cols; i++) {
				LU[i][j] /= LU[j][j];
			}
		}
	}


	//Matrix[] LUP = new Matrix[3];
	std::vector< Matrix<double> > LUP(3);
	Matrix<MType> Lower = LU.Get_Lower_Matrix();
	for (int i = 0; i < this->Cols; i++) {
		Lower[i][i] = 1;
	}
	Matrix<MType> pivot_mat(1, this->Cols);
	for (int i = 0; i < this->Cols; i++) {
		pivot_mat[0][i] = pivot[i];
	}
	LUP[0] = Lower;
	LUP[1] = LU.Get_Upper_Matrix();
	LUP[2] = pivot_mat;
	delete[] pivot;
	return LUP;
}
template<class MType> std::vector<Matrix<MType> >  Matrix<MType>::QR_Decomposition() {
	Matrix QR(*this);
	std::vector<double> R_diagonal(this->Rows);

	for (int k = 0; k < this->Rows; k++) {
		double nrm = 0;
		for (int i = k; i < this->Cols; i++) {
			nrm = std::sqrt(std::pow(nrm, 2) + std::pow(QR[i][k], 2));
		}

		if (nrm != 0.0) {
			if (QR[k][k] < 0) {
				nrm = -nrm;
			}
			for (int i = k; i < this->Cols; i++) {
				QR[i][k] /= nrm;
			}
			QR[k][k] += 1.0;

			for (int j = k + 1; j < this->Cols; j++) {
				double s = 0.0;
				for (int i = k; i < this->Cols; i++) {
					s += QR[i][k] * QR[i][j];
				}
				s = -s / QR.Matrix_Body[k][k];
				for (int i = k; i < this->Cols; i++) {
					QR[i][j] += s * QR[i][k];
				}
			}
		}
		R_diagonal[k] = -nrm;
	}

	int m = this->Rows;
	int n = this->Cols;
	Matrix XQ(m, n);
	Matrix<double> Q(m, n);
	for (int k = n - 1; k >= 0; k--) {
		for (int i = 0; i < m; i++) {
			Q[i][k] = 0.0;
		}
		Q[k][k] = 1.0;
		for (int j = k; j < n; j++) {
			if (QR[k][k] != 0) {
				double s = 0.0;
				for (int i = k; i < m; i++) {
					s += QR[i][k] * Q[i][j];
				}
				s = -s / QR[k][k];
				for (int i = k; i < m; i++) {
					Q[i][j] += s * QR[i][k];
				}
			}
		}
	}
	XQ = Q;
	Matrix XR(n, n);
	Matrix<double> R(n,n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i < j) {
				R[i][j] = QR[i][j];
			}
			else if (i == j) {
				R[i][j] = R_diagonal[i];
			}
			else {
				R[i][j] = 0.0;
			}
		}
	}
	XR = R;

	Matrix Xh(m, n);
	Matrix<double> H(m,n);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (i >= j) {
				H[i][j] = QR[i][j];
			}
			else {
				H[i][j] = 0.0;
			}
		}
	}
	Xh = H;
	std::vector<Matrix<MType> > Result(3);

	Result[0] = XQ;
	Result[1] = XR;
	Result[2] = Xh;
	return Result;
}













class Column {
private:
	void Determinate_Column_Type();
public:
	std::string Column_Name;
	int Column_Number;
	int Column_Type;
	int Missing;
	String_Vector Values;
	String_Vector Categories;
	std::string Get_Column_Type();
	void Rescan_Info();
	friend std::ostream& operator<<(std::ostream &output,Column const &column);
	std::string &operator[](int index);
};
std::ostream& operator<<(std::ostream &output, Column const &column) {
	output << column.Column_Name << "\n";
	for (int i = 0; i < (int)column.Values.size(); i++) {
		output << std::setw(2) <<column.Values[i] << "\n";
	}
	return output;
}
std::string &Column::operator[](int index) {
		if (index < 1 || index >(int)this->Values.size()) {
			std::cout << "\nIllegal Index Exception , Index -[" << index << "]\n";
			std::cout << "Valid Index Range Is - [ 1 - " << this->Values.size() << " ]\n";
			std::cout << "Returning Value At Index 1\n\n";
			return this->Values[0];
		}
		else {
			return this->Values[index - 1];
		}
	

	

}
 std::string Column::Get_Column_Type() {
	switch (this->Column_Type) {
	case 0:
		return "Numeric";
	case 1:
		return "Categorical";
	case 2:
		return "Text";
	case 3:
		return "Date-Time";
	}
	return "";
}
 void Column::Determinate_Column_Type() {
	 bool numeric = true;
	 //numeric test
	 std::string tm;
	 static std::regex num_regex("(^|\\s)([\\+-]?([0-9]+\\.?[0-9]*|\\.?[0-9]+))(\\s|$)");
	 static std::regex date_regex("([0-9]{4})/([0-9]{2})/([0-9]{2}) ([0-9]{2}):([0-9]{2}):([0-9]{2})");
	 for (int i = 0; i < (int)this->Values.size(); i++) {
		 tm = this->Values[i];
		 if (!std::regex_match(tm, num_regex)) {
			 if (!(tm == "")) {
				 numeric = false;
			 }
		 }
	 }

	 if (numeric == true) {
		 this->Column_Type = 0;
		 return;
	 }


	 //date time test
	 int num_of_dates = 0;
	 for (int i = 0; i < (int)this->Values.size(); i++) {
		 tm = this->Values[i];
		 if (std::regex_match(tm, date_regex)) {
			 num_of_dates++;
		 }
	 }
	 if (num_of_dates > 0) {
		 this->Column_Type = 3;
		 return;
	 }

	 String_Vector Categories = this->Values;

	 //categorical test
	 for (int i = 0; i < (int)Categories.size(); i++) {
		 tm = Categories[i];
		 if (tm != "0" && tm != "") {
			 for (int j = i + 1; j < (int)Categories.size(); j++) {
				 if (Categories[j] == (tm)) {
					 Categories[j] = "0";
				 }
			 }
		 }
	 }
	 int num_of_cat = 0;
	 for (int i = 0; i < (int)Categories.size(); i++) {
		 if (Categories[i] != "0") {
			 num_of_cat++;
		 }
	 }

	 if (num_of_cat < (int)(Categories.size() / 2) + 1) {
		 Categories.erase(std::remove(Categories.begin(), Categories.end(), "0"),Categories.end());
		 this->Categories = Categories;
		 this->Column_Type = 1;
		 return;
	 }




	 //else set as text
	 this->Column_Type = 2;
	 return;
 }
 void Column::Rescan_Info() {
	 this->Determinate_Column_Type();
 }






class CSV_File {
private:
	std::fstream RW;
public:
	CSV_File();
	std::vector<Column> Data;
	int Number_Of_Rows;
	int Number_Of_Columns;



	friend std::ostream& operator<<(std::ostream &output, CSV_File const &csv);
	void Load_CSV(std::string File_Location);
	void Write_CSV(std::string File_Name);
	Column Get_Column(int Column_Number);
	Column &operator[](int index);
	String_Vector Get_Row(int Row_Number);
	void Add_Row();
	void Add_Row(String_Vector Values);
	void Add_Column(std::string Column_Name);
	void Add_Column(std::string Column_Name, String_Vector const &Values);
	void Remove_Row(int Row_Number);
	void Remove_Row(int From_Row_Number, int To_Row_Number);
	void Remove_Column(int Column_Number);
	void Print_Column_Info();
	CSV_File Split_Data(int Split_Percentage);

};




CSV_File::CSV_File() {
	this->Number_Of_Columns = 0;
	this->Number_Of_Rows = 0;
}

void CSV_File::Write_CSV(std::string File_Name){
	File_Name.append(".csv");
	std::ofstream writer;
	writer.open(File_Name);
	for (int i = 0; i < (int)this->Data.size(); i++) {
		writer << this->Data[i].Column_Name << ",";
	}
	writer << "\n";
	for (int i = 0; i < (int)this->Data[0].Values.size(); i++) {
		for (int j = 0; j < (int)this->Data.size(); j++) {
		//	writer << this->Get_Value_At(j+1, i+1) << ",";
			writer << this->Data[j][i+1] << ",";

		}
		writer << "\n";
	}

	writer.close();
}


Column &CSV_File::operator[](int index) {
		if (index < 1 || index >(int)this->Data.size()) {
			std::cout << "\nIllegal Index Exception , Index -[" << index << "]\n";
			std::cout << "Valid Index Range Is - [ 1 - " << this->Data.size() << " ]\n";
			std::cout << "Returning Column At Index 1\n\n";
			return this->Data[0];
		}
		else {
			return this->Data[index - 1];
		}
	

}


String_Vector CSV_File::Get_Row(int Row_Number) {
	if (Row_Number < 1 || Row_Number > this->Number_Of_Rows) {
		std::cout << "==== Error -> [ Invalid Row ] ====\n";
		return String_Vector();
	}
	String_Vector Row;
	for (int i = 1; i <= (int)this->Data.size(); i++) {
		//Row.push_back(this->Get_Value_At(i, Row_Number));
		Row.push_back(this->Data[i-1][Row_Number]);
	}
	return Row;
}

std::ostream &operator<<(std::ostream &out, String_Vector vec) {
	out << std::left;
	out << "[";
	for (int i = 0; i < (int)vec.size(); i++) {
		if (i < (int)vec.size() - 1) {
			out << std::setw(3) << vec[i] << " | ";

		}
		else {
			out << std::setw(3) << vec[i];
			
		}
	}
	out << "]\n";
	return out;
}

std::ostream& operator<<(std::ostream &output, CSV_File const &csv) {
	for (int i = 0; i < (int)csv.Data.size(); i++) {
		output << std::setw(csv.Data[i].Column_Name.size()+2) << csv.Data[i].Column_Name;
	}
	output << "\n";
	int col = 0;
	for (int i = 0; i < (int)csv.Data[0].Values.size(); i++) {
		for (int j = 0; j < (int)csv.Data.size(); j++) {
			output << std::setw(csv.Data[col].Column_Name.size()) <<csv.Data[j].Values[i]<<"|";
			col++;
		}
		col = 0;
		output << "\n";

	}
	return output;
}

void CSV_File::Load_CSV(std::string File_Location){
	File_Location.append(".csv");
	this->RW.open(File_Location);
	std::string sample;
	std::vector<String_Vector> Loaded_Data;
	std::vector<std::string> line;
	while (!RW.eof()) {
		RW >> sample;
		line = Split_String(sample,',');
		Loaded_Data.push_back(line);
	}

	for (int i = 0; i < (int)Loaded_Data[0].size();i++) {
		Column temp;
		String_Vector values(Loaded_Data.size()-1);
		int aux = 0;
		int missing = 0;
		for (int j = 1; j < (int)Loaded_Data.size(); j++) {
			values[aux++] = Loaded_Data[j][i];
			if (Loaded_Data[j][i] == "") {
				missing++;
			}
		}
		temp.Column_Number = i;
		temp.Values = values;
		temp.Missing = missing;
		temp.Column_Name = Loaded_Data[0][i];
		temp.Rescan_Info();
		this->Data.push_back(temp);
	}
	this->Number_Of_Columns = this->Data.size();
	this->Number_Of_Rows = this->Data[0].Values.size();
	RW.close();

}

Column CSV_File::Get_Column(int Column_Number) {
	if (Column_Number < 1 || Column_Number >(int)this->Data.size()) {
		std::cout << "==== Error! -> [Invaid Column] ====";
		return Column();
	}
	else {
		return this->Data[Column_Number - 1];
	}
}

void CSV_File::Add_Row() {

	for (int i = 0; i < this->Number_Of_Columns; i++) {
		this->Data[i].Values.push_back("0");
	}
	this->Number_Of_Rows++;

}

void CSV_File::Add_Row(String_Vector Values) {
	if (Values.size() > this->Data.size()) {
		std::cout << "==== Error -> [ There Are More Values In Row Then Columns In The Dataset ] ====\n";
		return;
	}
	else {
		if (Values.size() < this->Data.size()) {
			int i = 0;
			for (i = 0; i < (int)Values.size(); i++) {
				this->Data[i].Values.push_back(Values[i]);
			}
			for (i; i < (int)this->Data.size(); i++) {
				this->Data[i].Values.push_back("0");
			}
			this->Number_Of_Rows++;
		}
		else {
			for (int i = 0; i < (int)this->Data.size(); i++) {
				this->Data[i].Values.push_back(Values[i]);
			}
			this->Number_Of_Rows++;
		}


	}
}

void CSV_File::Add_Column(std::string Column_Name) {
	Column new_column;
	this->Number_Of_Columns++;
	new_column.Column_Name = Column_Name;
	new_column.Column_Number = this->Number_Of_Columns;
	new_column.Missing = 0;
	new_column.Values.reserve(this->Number_Of_Rows);
	for (int i = 0; i < this->Number_Of_Rows; i++) {
		new_column.Values.push_back("0");
	}
	this->Data.push_back(new_column);
}

void CSV_File::Add_Column(std::string Column_Name, String_Vector const &Values) {
	Column new_column;
	this->Number_Of_Columns++;
	new_column.Column_Name = Column_Name;
	new_column.Column_Number = this->Number_Of_Columns;
	new_column.Missing = 0;
	if ((int)Values.size() > this->Number_Of_Rows) {
		int gap =Values.size()- this->Number_Of_Rows;
		for (int i = 0; i < gap; i++) {
			this->Add_Row();
		}
		new_column.Values = Values;
		this->Data.push_back(new_column);

	}
	else {
		if ((int)Values.size() < this->Number_Of_Rows) {
			new_column.Values = Values;
			int gap = this->Number_Of_Rows - Values.size();
			new_column.Missing = gap;
			for (int i = 0; i < gap; i++) {
				new_column.Values.push_back("0");
			}
			this->Data.push_back(new_column);
		}
		else {
			new_column.Values = Values;
			this->Data.push_back(new_column);
		}

	}
}

void CSV_File::Remove_Row(int Row_Number) {
	if (Row_Number < 1 || Row_Number > this->Number_Of_Rows) {
		std::cout << "==== Error -> [ Invalid Row ] ====\n";
		return;
	}

	for (int i = 0; i < this->Number_Of_Columns; i++) {
		this->Data[i].Values.erase(this->Data[i].Values.begin() + (Row_Number-1));
	}
	this->Number_Of_Rows--;
}
void CSV_File::Remove_Row(int From_Row_Number, int To_Row_Number) {
	if (From_Row_Number < 1 || From_Row_Number > this->Number_Of_Rows || To_Row_Number < 1 || To_Row_Number > this->Number_Of_Rows) {
		std::cout << "==== Error -> [ Invalid Row ] ====\n";
		return;
	}
	
	int dem = 0;
	for (int i = From_Row_Number; i <= To_Row_Number; i++) {
		this->Remove_Row(i - dem);
		dem++;
	}
	

}


void CSV_File::Remove_Column(int Column_Number) {
	if (Column_Number < 1 || Column_Number >(int)this->Data.size()) {
		std::cout << "==== Error! -> [Invaid Column] ====";
		return;
	}
	this->Data.erase(this->Data.begin() + (Column_Number - 1));
	this->Number_Of_Columns--;

}

void CSV_File::Print_Column_Info() {
	std::cout << std::right;

	for (int i = 0; i < this->Number_Of_Columns; i++) {
		Column sampler = this->Data[i];
		std::cout << "\n==================================\n";
		std::cout << "Column Number : " << std::setw(5) << "[ " << sampler.Column_Number << " ]" << "\n";
		std::cout << "Column Name : " << std::setw(7) << "[ " << sampler.Column_Name << " ]" << "\n";
		if (sampler.Column_Type == 1) {
			std::cout << "Column Categories : ";
			for (int j = 0; j < (int)sampler.Categories.size(); j++) {
				if (j < (int)sampler.Categories.size() - 1) {
					std::cout << "{" << sampler.Categories[j] << "}, ";

				}
				else {
					std::cout << "{" << sampler.Categories[j] << "}";

				}
			}
			std::cout << "\n";
		}
		else {
			std::cout << "Column Type : " << std::setw(7) << "[ " << sampler.Get_Column_Type() << " ]" << "\n";
		}
		std::cout << "Column Missing : " << std::setw(4) << "[ " << sampler.Missing << " ]" << "\n";



		std::cout << "\n==================================\n";

	}
}

CSV_File CSV_File::Split_Data(int Split_Percentage) {
	if (Split_Percentage > 100) {
		std::cout << "Cannot Split Into: " <<(100 - Split_Percentage) << "%  /  " << Split_Percentage << "% ";
		return CSV_File();
	}

	CSV_File Split;
	int number_of_rows;
	number_of_rows = (int)(((double)this->Number_Of_Rows / 100)*Split_Percentage);
	Split.Number_Of_Rows = number_of_rows + 1;
	Split.Number_Of_Columns = this->Number_Of_Columns;
	for (int i = 0; i < this->Number_Of_Columns; i++) {
		Split.Data.push_back(Column());
		Split.Data[Split.Data.size() - 1].Column_Name = this->Data[i].Column_Name;
		Split.Number_Of_Columns++;
	}

	for (int i = this->Number_Of_Rows - number_of_rows; i < this->Number_Of_Rows; i++) {
		Split.Add_Row(this->Get_Row(i));

	}


	this->Remove_Row(Number_Of_Rows - number_of_rows, Number_Of_Rows);

	this->Number_Of_Rows -= number_of_rows;
	for (unsigned i = 0; i < Split.Data.size(); i++) {
		Split.Data[i].Rescan_Info();
	}
	return Split;
}
