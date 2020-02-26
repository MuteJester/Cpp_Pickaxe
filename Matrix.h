#pragma once
#include <iostream>
#include <vector>
#define IDENTITY_MATRIX 'i'


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
	Matrix(int N,char x = IDENTITY_MATRIX) {
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
	~Matrix() {};
	int getRows() { return Rows; }
	int getCols() { return Cols; }
	void setRows(int &Rows) { this->Rows = Rows; }
	void setCols(int &Cols) { this->Cols = Cols; }
	void Print_Matrix();
	std::vector<MType> &operator[](int a){
		return Matrix_Body[a];
	}
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
	void operator+(Matrix<MType> &B);
	void operator-(Matrix<MType> &B);
	void Dot_Product(Matrix<MType> &B);
	void Matrix_Transpose();
	void Multiply_By_Scalar(int const &scalar);
	Matrix<MType> Hadamard_Product(Matrix<MType> &Mul_By);
	Matrix<MType> Kronecker_Product(Matrix<MType> &Mul_By);
	void Horizontal_Matrix_Concatenation(Matrix<MType> &To_HConcat);

protected:
	std::vector<std::vector<MType> > Matrix_Body;
	int Rows, Cols;

};

template<class MType> void  Matrix<MType>::Print_Matrix() {
	for (int i = 0; i < Rows; i++) {
		for (int j = 0; j < Cols; j++) {
			std::cout << Matrix_Body[i][j] << " ";
		}
		std::cout << "\n";
	}
}

template<class MType> void Matrix<MType>::Matrix_Transpose() {


	Matrix<MType> temp(Cols, Rows);

	

	for (int i = 0; i < Cols; i++) {
		for (int j = 0; j < Rows; j++) {
			temp[i][j] = Matrix_Body[j][i];
		}
	}
	*this = temp;
}

template<class MType> void Matrix<MType>::Multiply_By_Scalar(int const &scalar) {
	for (int i = 0; i < Rows; i++) {
		for (int j = 0; j < Cols; j++) {
			Matrix_Body[i][j] *= scalar;
		}
	}
}

template<class MType> void Matrix<MType>::operator+(Matrix<MType> &B) {
	if (Rows != B.getRows() || Cols != B.getCols()) {
		return;
	}
	else {
		for (int i = 0; i < Rows; i++) {
			for (int j = 0; j < Cols; j++) {
				Matrix_Body[i][j] += B[i][j];
			}
		}
	}
}

template<class MType> void Matrix<MType>::operator-(Matrix<MType> &B) {
	if (Rows != B.getRows() || Cols != B.getCols()) {
		return;
	}
	else {
		for (int i = 0; i < Rows; i++) {
			for (int j = 0; j < Cols; j++) {
				Matrix_Body[i][j] -= B[i][j];
			}
		}
	}
}

template<class MType> void Matrix<MType>::Dot_Product(Matrix<MType> &B) {

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
		*this = temp;
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
