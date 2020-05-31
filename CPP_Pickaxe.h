#pragma once
#include <iostream>
#include <stdio.h>
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
#include <algorithm>
#include <stack>
#include <exception>

namespace Pickaxe{

#define Mmax(a,b) ((a) > (b) ? (a) : (b))
#define Mmin(a,b) ((a) < (b) ? (a) : (b))
#define Mabs(x)  ( (x<0) ? -(x) : x )
#define DEG_TO_RAD(A) A * PI/180
#define RAD_TO_DEG(A) A * 180/PI
#define SQUARE_IT(NUM) NUM*NUM
#define CUBE_IT(NUM) NUM*NUM*NUM

//tokens
#define IDENTITY_MATRIX 'i'
#define RANDOM_MATRIX 'r'
#define SPEARMAN 18723
#define PEARSON 18722
#define FORMULA 'F'
#define POINT_STRUCTRE_THOTH
#define MATRIX_STRUCTRE_THOTH

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
static double Sigmoid(double x) {
	return 1.0 / (1.0 + std::exp(-x));
}
static double Euclidean_Distance(std::vector<double> X1, std::vector<double> X2) {
	if (X1.size() != X2.size()) {
		std::cout<<"Both Points Must Contain The Same Amount Of Values\n";
		return std::numeric_limits<double>::infinity();
	}
	else {
		double dis = 0;
		for (int i = 0; i < X1.size(); i++) {
			dis += std::pow((X2[i] - X1[i]), 2);
		}
		return std::sqrt(dis);
	}
}
static double Manhattan_Distance(std::vector<double> X1, std::vector<double> X2) {
	if (X1.size() != X2.size()) {
		std::cout<<"Both Points Must Contain The Same Amount Of Values\n";
		return std::numeric_limits<double>::infinity();
	}
	else {
		double dis = 0;
		for (int i = 0; i < X1.size(); i++) {
			dis += std::abs(X2[i] - X1[i]);
		}
		return dis;
	}
}
static double Rectified(double value) {
	return std::fmax(0, value);
}
template < class T >
std::vector< std::vector<T> > To_2D(const std::vector<T>& flat_vec, std::size_t ncols)
{
	// sanity check
	if (ncols == 0 || flat_vec.size() % ncols != 0) throw std::domain_error("bad #cols");

	const auto nrows = flat_vec.size() / ncols;

	std::vector< std::vector<T> > mtx;
	const auto begin = std::begin(flat_vec);
	for (std::size_t row = 0; row < nrows; ++row) mtx.push_back({ begin + row * ncols, begin + (row + 1)*ncols });
	return mtx;
}
 
class Random_Utilitis {
public:
	static int Random_INT(int minimum_value,int maximum_value) {
		static std::random_device seed;
		static std::mt19937 random_number_generator(seed());
		std::uniform_int_distribution<size_t> indices(minimum_value, maximum_value-1);
		return (int)indices(random_number_generator);
	}
	static double Random_DOUBLE(double minimum_value, double maximum_value) {
		static std::random_device seed;
		static std::mt19937 random_number_generator(seed());
		std::uniform_real_distribution<> indices(minimum_value, maximum_value);
		return indices(random_number_generator);
	}
	
};



template<class p_type> class Point {
public:
	p_type x, y, z;

	Point() {
		x = y = z = 0;
	}
	Point(p_type x, p_type y) {
		this->x = x;
		this->y = y;
		this->z = 0;
	}
	Point(p_type x, p_type y, p_type z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}
	double get_Distance(Point<p_type> const &point_B) {
		return sqrt(SQUARE_IT((point_B.x - this->x)) + SQUARE_IT((point_B.y - this->y)) + SQUARE_IT((point_B.z - this->z)));
	}
	static double get_Squared_Distance(Point<p_type> const &point_A,Point<p_type> const &point_B) {
		return (SQUARE_IT((point_B.x - point_A.x)) + SQUARE_IT((point_B.y - point_A.y)) + SQUARE_IT((point_B.z - point_A.z)));
	}

	double get_Slope(Point<p_type> const &point_B) {
		return (point_B.y - this->y) / (point_B.x - this->x);
	}
	double get_Magnitude() {
		return sqrt((SQUARE_IT(this->x)) + (SQUARE_IT(this->y)) + (SQUARE_IT(this->z)));
	}

};
template<class p_type> double X_Linear_Interpolation(Point<p_type> const &A, Point<p_type> const &B, p_type y_value) {
	return (((y_value - A.y)*(B.x - A.x)) / (B.y - A.y) + A.x);
}
template<class p_type> double Y_Linear_Interpolation(Point<p_type> const &A, Point<p_type> const &B, p_type x_value) {
	return (double)(x_value - A.x)*(B.y - A.y) / (B.x - A.x) + A.y;
}
void Quadratic_Equation(double const &a, double const &b, double const &c, double &result_root_a, double &result_root_b) {
	result_root_a = (-b + std::sqrt(((SQUARE_IT(b)) - 4 * a*c))) / 2 * a;
	result_root_b = ((-b - std::sqrt(((SQUARE_IT(b)) - 4 * a*c)))) / 2 * a;
}
template<class p_type> double Point_Dot_Product(Point<p_type> const &A, Point<p_type> const &B) {
	return A.x*B.x + A.y*B.y + A.z*B.z;
}
template<class p_type> std::ostream &operator<<(std::ostream &out, Point< p_type> const &body) {
	out << "[ " << body.x << ", " << body.y << ", " << body.z << " ]";
	return out;
}
template<class p_type> std::ostream &operator<<(std::ostream &out, std::vector<Point< p_type> > const &body) {
	for (int i = 0; i < body.size(); i++) {
		out << body[i] << std::endl;
	}
	return out;
}
template<class p_type> bool operator<(Point<p_type> const &A, Point<p_type> const &B) {
	if (A.x < B.x || A.y < B.y || A.x < B.x) {
		return true;
	}
	else {
		return false;
	}
}
template<class p_type> Point<p_type> Get_Smallest_Point_In_Vector(std::vector<Point<p_type> > const & Points) {
	Point<p_type> cur_min = Points[0];
	for (int i = 1; i < Points.size(); i++) {
		if (Points[i] < cur_min) {
			cur_min = Points[i];
		}
	}
	return cur_min;
}
template<class p_type> Point<p_type> Get_Largest_Point_In_Vector(std::vector<Point<p_type> > const & Points) {
	Point<p_type> cur_max = Points[0];
	for (int i = 1; i < Points.size(); i++) {
		if (!(Points[i] < cur_max)) {
			cur_max = Points[i];
		}
	}
	return cur_max;
}

inline double Squared_Point_Distance(Point<double> first, Point<double> second) {
	return (first.x - second.x)*(first.x - second.x) + (first.y - second.y)*(first.y - second.y) + (first.z - second.z)*(first.z - second.z);
}


//symbolic


class Monomial {
public:
	double Degree;
	double Coefficient;
	double Pry;
	bool is_monom;
	bool is_cos, is_sin, is_tan, is_lan, is_E;

	Monomial();
	Monomial(int const &Coeff, int const &Deg);
	Monomial(int const &Coeff);
	Monomial(std::string const &form);
	friend std::ostream &operator<<(std::ostream &out, Monomial const &source);
	int operator+(Monomial const&b);
	void operator-(Monomial const&b);
	void operator*(Monomial const&b);
	void operator^(int const&b);
	void operator/(Monomial const &b);
	void Derive();
	void Derive(int const &mag);
	double operator[](double const &x_value);

};


void get_Cof_Deg(std::string const &to_parse, double &cof, double &deg) {
	std::stringstream ss, ss2;
	std::string via, via2;
	std::size_t pos;
	pos = to_parse.find("*");
	if (pos == std::string::npos) {
		pos = to_parse.find("x");
		if (pos == std::string::npos) {
			return;
		}
		else {
			if (to_parse[pos + 1] == '^') {
				pos += 2;
				for (std::size_t i = pos; i < to_parse.size(); i++) {
					ss << to_parse[i];
				}
				via = ss.str();
				ss.str(std::string());
				deg = std::stod(via.c_str());

				cof = 1;
				return;
			}
			else {
				cof = 1;
				deg = 1;
			}
		}
	}
	else {
		for (std::size_t i = 0; i < pos; i++) {
			ss << to_parse[i];
		}
		via = ss.str();
		ss.str(std::string());
		cof = atoi(via.c_str());
		pos = to_parse.find("x");
		if (to_parse[pos + 1] == '^') {
			pos += 2;
			for (std::size_t i = pos; i < to_parse.size(); i++) {
				ss << to_parse[i];
			}
			via = ss.str();
			ss.str(std::string());
			std::string::size_type sz;
			deg = std::stod(via, &sz);
			return;
		}
		else {
			deg = 1;
			return;
		}
	}
}
Monomial::Monomial() {
	this->Degree = 0;
	this->Coefficient = 0;
	is_cos = is_sin = is_tan = is_lan = is_E = false;
	is_monom = false;
}
Monomial::Monomial(int const &Coeff, int const &Deg) {
	this->Coefficient = Coeff;
	this->Degree = Deg;
	is_cos = is_sin = is_tan = is_lan = is_E = false;
	is_monom = true;
	this->Pry = 1;
}
Monomial::Monomial(int const &Coeff) {
	this->Degree = 1;
	this->Coefficient = Coeff;
	is_cos = is_sin = is_tan = is_lan = is_E = false;
	is_monom = true;
	this->Pry = 1;
}
Monomial::Monomial(std::string const &form) {
	is_cos = is_sin = is_tan = is_lan = is_E = false;
	is_monom = false;
	this->Pry = 1;
	std::size_t pos, endp;
	std::stringstream ss, ss2;
	std::string via, via2, subs;
	double DEG, i = 0;
	double COF;
	if (form.find("cos(") != std::string::npos) {
		this->is_cos = true;
		pos = form.find("cos(");
		if (pos > 0) {
			if (form[pos - 1] == '*') {
				for (std::size_t i = 0; i < pos - 1; i++) {
					ss << form[i];
				}
				via = ss.str();
				this->Pry = atoi(via.c_str());

				pos += 3;
				endp = form.find(")", pos);
				subs = form.substr(pos + 1, endp - 1);
				get_Cof_Deg(subs, COF, DEG);
				this->Degree = DEG;
				this->Coefficient = COF;
			}
			else {
				pos += 3;
				endp = form.find(")", pos);
				subs = form.substr(pos + 1, endp - 1);
				get_Cof_Deg(subs, COF, DEG);
				this->Degree = DEG;
				this->Coefficient = COF;
			}
		}
		else {
			pos += 3;
			endp = form.find(")", pos);
			subs = form.substr(pos + 1, endp - 1);
			get_Cof_Deg(subs, COF, DEG);
			this->Degree = DEG;
			this->Coefficient = COF;
		}

	}
	else if (form.find("sin(") != std::string::npos) {
		this->is_sin = true;
		pos = form.find("sin(");
		if (pos > 0) {
			if (form[pos - 1] == '*') {
				for (std::size_t i = 0; i < pos - 1; i++) {
					ss << form[i];
				}
				via = ss.str();
				this->Pry = atoi(via.c_str());

				pos += 3;
				endp = form.find(")", pos);
				subs = form.substr(pos + 1, endp - 1);
				get_Cof_Deg(subs, COF, DEG);
				this->Degree = DEG;
				this->Coefficient = COF;
			}
			else {
				pos += 3;
				endp = form.find(")", pos);
				subs = form.substr(pos + 1, endp - 1);
				get_Cof_Deg(subs, COF, DEG);
				this->Degree = DEG;
				this->Coefficient = COF;
			}
		}
		else {
			pos += 3;
			endp = form.find(")", pos);
			subs = form.substr(pos + 1, endp - 1);
			get_Cof_Deg(subs, COF, DEG);
			this->Degree = DEG;
			this->Coefficient = COF;
		}
	}
	else if (form.find("tan(") != std::string::npos) {
		this->is_tan = true;
		pos = form.find("tan(");
		if (pos > 0) {
			if (form[pos - 1] == '*') {
				for (std::size_t i = 0; i < pos - 1; i++) {
					ss << form[i];
				}
				via = ss.str();
				this->Pry = atoi(via.c_str());

				pos += 3;
				endp = form.find(")", pos);
				subs = form.substr(pos + 1, endp - 1);
				get_Cof_Deg(subs, COF, DEG);
				this->Degree = DEG;
				this->Coefficient = COF;
			}
			else {
				pos += 3;
				endp = form.find(")", pos);
				subs = form.substr(pos + 1, endp - 1);
				get_Cof_Deg(subs, COF, DEG);
				this->Degree = DEG;
				this->Coefficient = COF;
			}
		}
		else {
			pos += 3;
			endp = form.find(")", pos);
			subs = form.substr(pos + 1, endp - 1);
			get_Cof_Deg(subs, COF, DEG);
			this->Degree = DEG;
			this->Coefficient = COF;
		}
	}
	else if (form.find("ln(") != std::string::npos) {
		this->is_lan = true;
		pos = form.find("ln(");
		if (pos > 0) {
			if (form[pos - 1] == '*') {
				for (std::size_t i = 0; i < pos - 1; i++) {
					ss << form[i];
				}
				via = ss.str();
				this->Pry = atoi(via.c_str());

				pos += 2;
				endp = form.find(")", pos);
				subs = form.substr(pos + 1, endp - 1);
				get_Cof_Deg(subs, COF, DEG);
				this->Degree = DEG;
				this->Coefficient = COF;
			}
			else {
				pos += 3;
				endp = form.find(")", pos);
				subs = form.substr(pos + 1, endp - 1);
				get_Cof_Deg(subs, COF, DEG);
				this->Degree = DEG;
				this->Coefficient = COF;
			}
		}
		else {
			pos += 3;
			endp = form.find(")", pos);
			subs = form.substr(pos + 1, endp - 1);
			get_Cof_Deg(subs, COF, DEG);
			this->Degree = DEG;
			this->Coefficient = COF;
		}
	}
	else {
		this->is_monom = true;
		get_Cof_Deg(form, COF, DEG);
		this->Degree = DEG;
		this->Coefficient = COF;
	}


}

std::ostream &operator<<(std::ostream &out, Monomial const &source) {
	if (source.is_cos == true) {
		if (source.Pry > 1) {
			std::cout << source.Pry << "*cos(";

		}
		else {
			std::cout << "cos(";
		}
		if (source.Degree == 0) {
			out << source.Coefficient;
			std::cout << ")";

			return out;

		}
		else if (source.Degree == 1) {
			if (source.Coefficient == 1) {
				out << "X";
				std::cout << ")";

				return out;
			}
			else {
				out << source.Coefficient << "*X";
				std::cout << ")";

				return out;
			}


		}
		else {
			out << source.Coefficient << "*X^" << source.Degree;
			std::cout << ")";

			return out;
		}
	}
	else if (source.is_sin == true) {
		if (source.Pry == 1) {
			std::cout << "sin(";
		}
		else {
			std::cout << source.Pry << "*sin(";

		}
		if (source.Degree == 0) {
			out << source.Coefficient;
			std::cout << ")";

			return out;

		}
		else if (source.Degree == 1) {
			if (source.Coefficient == 1) {
				out << "X";
				std::cout << ")";

				return out;
			}
			else {
				out << source.Coefficient << "*X";
				std::cout << ")";

				return out;
			}


		}
		else {
			out << source.Coefficient << "*X^" << source.Degree;
			std::cout << ")";

			return out;
		}
	}
	else if (source.is_tan == true) {
		if (source.Pry == 1) {
			std::cout << "tan(";
		}
		else {
			std::cout << source.Pry << "*tan(";

		}
		if (source.Degree == 0) {
			out << source.Coefficient;
			std::cout << ")";

			return out;

		}
		else if (source.Degree == 1) {
			if (source.Coefficient == 1) {
				out << "X";
				std::cout << ")";

				return out;
			}
			else {
				out << source.Coefficient << "*X";
				std::cout << ")";

				return out;
			}


		}
		else {
			out << source.Coefficient << "*X^" << source.Degree;
			std::cout << ")";

			return out;
		}
	}
	else if (source.is_monom == true) {
		if (source.Degree == 0) {
			out << source.Coefficient;

			return out;

		}
		else if (source.Degree == 1) {
			if (source.Coefficient == 1) {
				out << "X";

				return out;
			}
			else {
				out << source.Coefficient << "*X";

				return out;
			}


		}
		else {
			out << source.Coefficient << "*X^" << source.Degree;

			return out;
		}
	}
	else if (source.is_lan == true) {
		if (source.Pry == 1) {
			std::cout << "ln(";
		}
		else {
			std::cout << source.Pry << "*ln(";

		}
		if (source.Degree == 0) {
			out << source.Coefficient;
			std::cout << ")";

			return out;

		}
		else if (source.Degree == 1) {
			if (source.Coefficient == 1) {
				out << "X";
				std::cout << ")";

				return out;
			}
			else {
				out << source.Coefficient << "*X";
				std::cout << ")";

				return out;
			}

		}
		else {
			out << source.Coefficient << "*X^" << source.Degree;
			std::cout << ")";

			return out;
		}
	}
	else {
		return out;
	}





}
int Monomial::operator+(Monomial const&b) {
	if (this->is_monom == true) {
		if (this->Degree == b.Degree) {
			this->Coefficient += b.Coefficient;
			return 1;
		}
		else {
			return 0;
		}
	}
	else {
		if (this->Coefficient == b.Coefficient && this->Degree == this->Degree) {
			this->Pry++;
			return 1;
		}
		else {
			return 0;
		}
	}
}
void Monomial::operator-(Monomial const&b) {
	if (this->Degree == b.Degree) {
		this->Coefficient -= b.Coefficient;
	}
	else {
		return;
	}
}
void Monomial::operator*(Monomial const&b) {
	if (b.Degree == 0) {
		this->Coefficient *= b.Coefficient;
	}
	else {
		this->Coefficient *= b.Coefficient;
		this->Degree += b.Degree;
	}
}
void  Monomial::operator^(int const&b) {
	if (b == 0) {
		this->Coefficient = 0;
		this->Degree = 1;
	}
	else {
		this->Degree *= b;
		double p = this->Coefficient;
		for (int i = 1; i < b; i++) {
			this->Coefficient *= p;
		}
	}
}
void Monomial::operator/(Monomial const &b) {
	this->Coefficient /= b.Coefficient;
	this->Degree -= b.Degree;
}
void Monomial::Derive() {
	if (is_cos == true) {
		is_cos = false;
		is_sin = true;

		if (Degree == 1) {
			Pry *= -Coefficient;
			return;
		}
		else if (Degree == 0) {
			Coefficient = 0;
		}
		else {

		}
	}
	else if (is_sin == true) {
		is_cos = true;
		is_sin = false;

		if (Degree == 1) {
			Pry *= Coefficient;
			return;
		}
		else if (Degree == 0) {
			Coefficient = 0;
		}
		else {

		}
	}
	else if (is_tan == true) {
		is_tan = false;
		is_cos = true;

		if (Degree == 1) {
			Pry *= Coefficient;
			Degree = -2;

			return;
		}
		else if (Degree == 0) {
			Coefficient = 0;
		}
		else {

		}
	}
	else if (is_lan == true) {


		if (Degree == 1) {
			Pry *= Coefficient;
			Degree = -1;
			return;
		}
		else if (Degree == 0) {
			Coefficient = 0;
		}
		else {

		}
	}
	else if (is_monom == true) {
		if (Degree == 1) {
			Degree = 0;
			return;
		}
		else if (Degree == 0) {
			Coefficient = 0;
		}
		else {
			this->Coefficient *= Degree;
			this->Degree--;
		}
	}
	else {
		return;
	}
}
void Monomial::Derive(int const &mag) {
	for (int i = 0; i < mag; i++) {
		Derive();
	}
}
double Monomial::operator[](double const &x_value) {
	if (this->is_cos == true) {
		if (Degree == 0) {
			if (Pry == 1) {
				return cos(Coefficient);
			}
			else {
				return Pry * cos(Coefficient);

			}
		}
		else if (Degree == 1) {
			if (Pry == 1) {
				return cos(Coefficient * x_value);
			}
			else {
				return Pry * cos(Coefficient * x_value);

			}
		}
		else {
			double t_valv = x_value;

			for (int i = 1; i < Degree; i++) {
				t_valv *= x_value;
			}
			t_valv *= Coefficient;
			if (Pry == 1) {
				return cos(t_valv);
			}
			else {
				return Pry * cos(t_valv);

			}
		}
	}
	else if (this->is_monom == true) {

		if (Degree == 0) {
			return Coefficient;
		}
		else if (Degree == 1) {
			return Coefficient * x_value;
		}
		else {
			double t_valv = x_value;

			for (int i = 1; i < Degree; i++) {
				t_valv *= x_value;
			}
			t_valv *= Coefficient;
			return t_valv;
		}
	}
	else if (this->is_tan == true) {
		if (Degree == 0) {
			if (Pry == 1) {
				return std::tan(Coefficient);
			}
			else {
				return Pry * std::tan(Coefficient);

			}

		}
		else if (Degree == 1) {
			if (Pry == 1) {
				return std::tan(Coefficient * x_value);
			}
			else {
				return Pry * std::tan(Coefficient * x_value);

			}
		}
		else {
			double t_valv = x_value;

			for (int i = 1; i < Degree; i++) {
				t_valv *= x_value;
			}
			t_valv *= Coefficient;
			if (Pry == 1) {
				return std::tan(t_valv);
			}
			else {
				return Pry * std::tan(t_valv);

			}
		}
	}
	else if (this->is_lan == true) {
		if (Degree == 0) {
			if (Pry == 1) {
				return std::log(Coefficient);
			}
			else {
				return Pry * std::log(Coefficient);

			}

		}
		else if (Degree == 1) {
			if (Pry == 1) {
				return std::log(Coefficient * x_value);
			}
			else {
				return Pry * std::log(Coefficient * x_value);

			}
		}
		else {
			double t_valv = x_value;

			for (int i = 1; i < Degree; i++) {
				t_valv *= x_value;
			}
			t_valv *= Coefficient;
			if (Pry == 1) {
				return std::log(t_valv);
			}
			else {
				return Pry * std::log(t_valv);

			}
		}
	}
	else if (this->is_sin == true) {
		if (Degree == 0) {
			if (Pry == 1) {
				return std::sin(Coefficient);
			}
			else {
				return Pry * std::sin(Coefficient);

			}

		}
		else if (Degree == 1) {
			if (Pry == 1) {
				return std::sin((Coefficient)* x_value);
			}
			else {
				return Pry * std::sin(Coefficient * x_value);

			}
		}
		else {
			double t_valv = x_value;

			for (int i = 1; i < Degree; i++) {
				t_valv *= x_value;
			}
			t_valv *= Coefficient;
			if (Pry == 1) {
				return std::sin(t_valv);
			}
			else {
				return Pry * std::sin(t_valv);

			}
		}
	}
	else {
		return 0;
	}

}







class Function {
protected:

public:
	std::string input;
	std::vector<Monomial> Body;
	std::vector<char> signs;

	std::stack<char> operations;
	std::stack<Monomial> operands;
	Function();
	Function(const char *function);
	void operator=(Function const &B);
	double operator[](double const &x_value);
	friend std::ostream &operator<<(std::ostream &out, Function const &func);
	void Derive();
	void Derive(int const &mag);
	double Get_Highest_Degree();
	double Newton_Raphson_Method(double const &guess, double const &deg_of_accuracy, double const &check_slope);
	Function Taylor_Polynomial(int const &derivative_n);
	void operator+(Function const &B);
	double Bisection(double const &a, double const &b, double const &EPSILON);
	std::vector<double> Get_Roots(double const &EPSILON);
	std::vector<Point<double> > Get_Polynomial_Maximum();
	std::vector<Point<double> > Get_Polynomial_Minimum();
	void Remove_Zero();

};

double commit_operation(char const &op, double const &a, double const &b) {
	double sum = 0;
	switch (op)
	{
	case '+':
		sum = a + b;
		return sum;
		break;
	case '-':
		sum = a - b;
		return sum;
		break;

	case '*':
		sum = a * b;
		return sum;
		break;

	case '/':

		sum = a / b;
		return sum;
		break;

	case '^':
		sum = pow(a, b);
		return sum;
		break;



	default:
		return 0;
		break;
	}
}
int sign_pt(char const &op) {
	switch (op)
	{
	case '+':
		return 1;
		break;
	case '-':
		return 1;
		break;

	case '*':
		return 2;
		break;

	case '/':

		return 2;
		break;

	case '^':
		return 3;
		break;



	default:
		return 0;
		break;
	}
}
Function::Function() {
}
Function::Function(const char *function) {
	this->input = std::string(function);
	std::stringstream ss;
	std::string via;
	std::size_t pos, endp;
	if (input.find(' ') == std::string::npos) {

		operands.push(Monomial(input));
		Body.push_back(Monomial(input));

	}
	else {
		endp = input.find(' ');
		pos = 0;
		while (true) {
			for (std::size_t i = pos; i < endp; i++) {
				ss << input[i];
			}
			via = ss.str();
			ss.str(std::string());
			operands.push(Monomial(via));
			Body.push_back(Monomial(via));
			operations.push(input[endp + 1]);
			signs.push_back(input[endp + 1]);
			pos = endp + 2;

			if (input.find(' ', pos + 1) == std::string::npos) {
				for (std::size_t i = pos + 1; i < input.size(); i++) {
					ss << input[i];
				}
				via = ss.str();
				ss.str(std::string());
				operands.push(Monomial(via));
				Body.push_back(Monomial(via));
				break;
			}
			endp = input.find(' ', pos + 1);

		}
	}

}
double Function::operator[](double const &x_value) {
	std::vector<double> infix;
	std::stack<char> opts;
	std::stack<double> oprnds;
	char temp;
	double a, b, sum;
	unsigned k = 0;
	for (unsigned i = 0; i < Body.size(); i++) {
		oprnds.push(Body[i][x_value]);
		if (k < signs.size()) {
			if (opts.empty()) {
				opts.push(signs[k]);
			}
			else {
				if (sign_pt(opts.top()) < sign_pt(signs[k])) {
					opts.push(signs[k]);
				}
				else {
					temp = opts.top();
					opts.pop();
					a = oprnds.top();
					oprnds.pop();
					b = oprnds.top();
					oprnds.pop();
					sum = commit_operation(temp, b, a);
					oprnds.push(sum);
					opts.push(signs[k]);
				}
			}
			k++;
		}

	}

	while (!opts.empty()) {
		a = oprnds.top();
		oprnds.pop();
		b = oprnds.top();
		oprnds.pop();
		temp = opts.top();
		opts.pop();
		sum = commit_operation(temp, b, a);
		oprnds.push(sum);
	}
	return oprnds.top();

}
void Function::operator=(Function const &B) {
	this->Body = B.Body;
	this->signs = B.signs;
	this->input = B.input;
}
std::ostream &operator<<(std::ostream &out, Function const &func) {
	unsigned k = 0;
	for (unsigned i = 0; i < func.Body.size(); i++) {

		out << func.Body[i] << " ";
		if (k < func.signs.size()) {
			out << func.signs[k] << " ";
			k++;
		}


	}
	for (unsigned i = 0; i < func.Body.size(); i++) {
	}
	return out;
}
void  Function::Derive() {
	for (unsigned i = 0; i < Body.size(); i++) {
		Body[i].Derive();
	}
	std::vector<Monomial>::iterator it;
	it = Body.end();
	it--;
	if ((*it).Coefficient == 0 && (*it).Degree == 0) {
		Remove_Zero();
	}

}
void Function::Derive(int const &mag) {
	for (int i = 0; i < mag; i++) {
		this->Derive();
	}
}
double Function::Get_Highest_Degree() {
	double max = 0;
	for (auto i : Body) {
		max < i.Degree ? max = i.Degree : max = max;
	}
	return max;
}
double Function::Newton_Raphson_Method(double const &guess, double const &deg_of_accuracy, double const &check_slope = 0.5) {
	Function f, f1;
	double fv, f1v, x = guess, x1 = 0;
	int critical = 100, counter = 0;
	while (critical > counter) {
		f = *this;
		f.Derive();
		f1 = f;
		f = *this;


		fv = f[x];
		f1v = f1[x];

		if (abs(f1[x]) < check_slope) {
			//	std::cout << "Slope is too small" << std::endl;
			break;
		}

		x1 = x - (fv / f1v);

		if (std::abs((x1 - x) / x1) < deg_of_accuracy) {
			return  x1;
		}
		x = x1;
		counter++;
	}
	//	std::cout << "Method does not converge due to oscillation" << std::endl;
	return 666;
}

void  Function::operator+(Function const &B) {
	unsigned i = 0;
	unsigned j = 0;
	while (j < B.Body.size())
	{
		if (Body.size() == 0) {
			Body.push_back(B.Body[j]);
			return;
		}
		else {
			for (i; i < this->Body.size(); i++) {
				if ((Body[i] + B.Body[j]) == 1) {
					j++;
					break;
				}

			}
			if (i == Body.size()) {
				Body.push_back(B.Body[j]);
				signs.push_back('+');
				j++;
			}
		}
	}
}
Function Function::Taylor_Polynomial(int const &derivative_n) {
	Function TP;
	Function Temp = *this;
	for (int i = 0; i < derivative_n; i++) {
		TP + Temp;
		Temp.Derive();
	}
	return TP;
}
double Function::Bisection(double const &a, double const &b, double const &EPSILON) {
	Function temp = *this;
	double aa = a;
	double bb = b;
	if (temp[aa] * temp[bb] >= 0)
	{
		//std::cout << "You have not assumed right a and b\n";
		return 666;
	}

	double c = aa;
	while ((bb - aa) >= EPSILON)
	{
		c = (aa + bb) / 2;

		if (temp[c] == 0.0)
			break;

		else if (temp[c] * temp[aa] < 0)
			bb = c;
		else
			aa = c;
	}
	return c;

}
std::vector<double> Function::Get_Roots(double const &EPSILON) {
	std::vector<double> result, f_res;
	Function temp = *this;
	Function Reg = *this;
	for (double i = -50; i < 30; i += 1) {
		result.push_back(this->Bisection(i, 20, EPSILON));
		if (temp[i] == 0) {
			result.push_back(i);
		}
	}
	for (double i = -30; i < 20; i += 0.05) {
		result.push_back(this->Newton_Raphson_Method(i, EPSILON));
	}

	sort(result.begin(), result.end());

	for (unsigned i = 0; i < result.size() - 1; i++) {
		if (std::abs(result[i] - result[i + 1]) < EPSILON)
		{
			result[i] = 666;
		}
	}
	sort(result.begin(), result.end());
	for (unsigned i = 0; i < result.size(); i++) {
		if (result[i] == 666) {
			break;
		}
		else {
			f_res.push_back(result[i]);
		}
	}


	return f_res;
}
std::vector<Point<double> >  Function::Get_Polynomial_Maximum() {
	Function div = *this;
	Function Reg = *this;
	std::vector<Point<double> > Result;
	div.Derive();
	std::vector<double> roots = div.Get_Roots(0.001);
	std::vector<double> max_points_x;
	div.Derive(); //seconed diverative
	for (unsigned i = 0; i < roots.size(); i++) {
		if (div[roots[i]] < 0)
			max_points_x.push_back(roots[i]);
	}

	for (unsigned i = 0; i < max_points_x.size(); i++) {
		Result.push_back(Point<double>(max_points_x[i], Reg[max_points_x[i]]));
	}

	return Result;
}
std::vector<Point<double> > Function::Get_Polynomial_Minimum() {
	Function div = *this;
	Function Reg = *this;
	std::vector<Point<double> > Result;
	div.Derive();
	std::vector<double> roots = div.Get_Roots(0.00001);
	std::vector<double> max_points_x;
	div.Derive(); //seconed diverative
	for (unsigned i = 0; i < roots.size(); i++) {
		if (div[roots[i]] > 0)
			max_points_x.push_back(roots[i]);
	}

	for (unsigned i = 0; i < max_points_x.size(); i++) {
		Result.push_back(Point<double>(max_points_x[i], Reg[max_points_x[i]]));
	}

	return Result;
}
void Function::Remove_Zero() {
	std::vector<Monomial>::iterator it;
	it = Body.end();
	it--;
	this->signs.pop_back();
	Body.erase(it);
}






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
		if (x == IDENTITY_MATRIX) {
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
		else if (x == RANDOM_MATRIX) {
			Matrix_Body = std::vector<std::vector<MType> >(N, std::vector<MType>(N));
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++) {
					Matrix_Body[i][j] = Random_Utilitis::Random_DOUBLE(0, 1);
				}
			}
			this->Rows = N;
			this->Cols = N;
		}
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
	int Rows, Cols;

	~Matrix() {};
	std::vector<MType> &operator[](int a) {
		return Matrix_Body[a];
	}
	void Random_Fill();
	void operator=(Matrix<MType> B) {
		this->Rows = B.Rows;
		this->Cols = B.Cols;
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
	double Get_Determinant();
	std::vector<double> Get_Eigen_Values();
	Matrix<MType> Get_Eigen_Vectors();
	std::vector<double>  Gaussian_Elimination(std::vector<double> equals_to);
	template<class MType> friend std::ostream &operator<<(std::ostream &out, Matrix<MType> const &mat);
	Matrix<MType> Reshape(int Rows, int Cols);
	std::vector<MType> Flatten();


protected:
	std::vector<std::vector<MType> > Matrix_Body;
};



template<class MType> std::ostream &operator<<(std::ostream &out, Matrix<MType> const &mat) {
	std::ios_base::sync_with_stdio(false);
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
	if (Rows != B.Rows || Cols != B.Cols) {
		return Matrix<double>();
	}
	else {
		Matrix<MType> Res(this->Rows, this->Cols);
		for (int i = 0; i < Rows; i++) {
			for (int j = 0; j < Cols; j++) {
				Res[i][j] = Matrix_Body[i][j] - B[i][j];
			}
		}
		return Res;
	}
}
template<class MType> Matrix<MType> Matrix<MType>::Dot_Product(Matrix<MType> &B) {

	if (Cols != B.Rows) {
		std::cout << "Rows Of Matrix A Are != Cols Of Matrix B";
		return Matrix<MType>();

	}
	else {
		MType sum = 0;
		Matrix<MType> temp(Rows, B.Cols);
		for (int i = 0; i < Rows; i++) {
			for (int j = 0; j < B.Cols; j++) {
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
template<class MType>  double Matrix<MType>::Get_Determinant() {
	if (this->Cols != this->Rows) {
		std::cout << "Matrix Needs To Be NxN\n";
		return 0;
	}
	std::vector<Matrix<double> > LU = this->LU_Decomposition();
	double det = LU[1][0][0];
	for (int i = 1; i < this->Cols; i++) {
		det *= LU[1][i][i];
	}
	int v = (LU[2][0][this->Cols - 1] - this->Cols - 1);
	if ((v % 2) == 0) {
		return det;
	}
	else {
		return -det;
	}
}
template<class MType>  std::vector<double> Matrix<MType>::Get_Eigen_Values() {
	Matrix< MType> CopyOfOriginal(*this);
	Matrix< MType> A(*this);
	std::vector<Matrix<MType> > QR = this->QR_Decomposition();
	for (int i = 0; i < 20; i++) {
		A = QR[1].Dot_Product(QR[0]);
		this->Matrix_Body = A.Matrix_Body;
		QR = this->QR_Decomposition();
	}

	std::vector<double> EigenValues = this->Diagonal();
	this->Matrix_Body = CopyOfOriginal.Matrix_Body;
	return EigenValues;
}
template<class MType> Matrix<MType> Matrix<MType>::Get_Eigen_Vectors() {
	Matrix< MType> CopyOfOriginal(*this);
	Matrix< MType> A(*this);
	std::vector<Matrix<MType> > QR = this->QR_Decomposition();
	Matrix<MType> Q(QR[0]);

	for (int i = 0; i < 30; i++) {
		A = QR[1].Dot_Product(QR[0]);
		this->Matrix_Body = A.Matrix_Body;
		QR = this->QR_Decomposition();
		Q = Q.Dot_Product(QR[0]);


	}

	this->Matrix_Body = CopyOfOriginal.Matrix_Body;
	return Q;
}
template<class MType> std::vector<double>  Matrix<MType>::Gaussian_Elimination(std::vector<double> equals_to) {
	int n = this->Matrix_Body[0].size();
	for (int p = 0; p < n; p++) {
		int max = p;
		for (int i = p + 1; i < n; i++) {
			if (Mabs(this->Matrix_Body[i][p]) >Mabs(this->Matrix_Body[max][p])) {
				max = i;
			}
		}
		std::vector<double> temp(this->Matrix_Body[p]);
		this->Matrix_Body[p] = this->Matrix_Body[max];
		this->Matrix_Body[max] = temp;
		double t = equals_to[p];
		equals_to[p] = equals_to[max];
		equals_to[max] = t;
		if (Mabs(this->Matrix_Body[p][p]) <= std::numeric_limits<double>::min()) {
			std::cout<<"\nMatrix is singular or nearly singular\n";
			return std::vector<double>();

		}
		for (int i = p + 1; i < n; i++) {
			double alpha = this->Matrix_Body[i][p] / this->Matrix_Body[p][p];
			equals_to[i] -= alpha * equals_to[p];
			for (int j = p; j < n; j++) {
				this->Matrix_Body[i][j] -= alpha * this->Matrix_Body[p][j];
			}
		}
	}
	std::vector<double> x(n);
	for (int i = n - 1; i >= 0; i--) {
		double sum = 0.0;
		for (int j = i + 1; j < n; j++) {
			sum += this->Matrix_Body[i][j] * x[j];
		}
		x[i] = (equals_to[i] - sum) / this->Matrix_Body[i][i];
	}

	return x;
}
template<class MType> Matrix<MType> Matrix<MType>::Reshape(int Rows, int Cols) {


	Matrix<MType> res(Rows,Cols);
	if (this->Rows*this->Cols == 0 || Rows * Cols != this->Rows * this->Cols) {
		std::cout << "Error Cannot Reshape Into " << Rows << "X" << Cols << "\n";
		throw(std::exception("Invalid_Shape"));
	}
	std::vector<MType>  F = this->Flatten();

	int z = 0;
	for (int i = 0; i < Rows; i++) {
		for (int j = 0; j < Cols; j++) {
			res[i][j] = F[z++];
		}
	}
	return res;

}
template<class MType> std::vector<MType> Matrix<MType>::Flatten() {
	std::vector<MType> Result;
	for (int i = 0; i < this->Rows; i++) {
		for (int j = 0; j < this->Cols; j++) {

			Result.push_back(this->Matrix_Body[i][j]);
		}
	}
	return Result;
}








class Column {
	//types 
	// 0 = Numeric
	// 1 = Categorical
	// 2 = Text
	// 3 = Date Time
private:
	void Determinate_Column_Type();
public:
	Column();
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
	double Column_Mean();
	double Column_Median();
	double Column_Standard_Deviation();
	double Column_Variance();
	double Max_Value();
	double Min_Value();
	void Fill_Missing_With_Column_Mean();
	void Fill_Missing_With_Column_Median();
	void Fill_Missing_With_Pattern(std::string Pattern);
	void operator=(Column const &source);
	std::vector<double> get_Numeric_Values() {
		std::vector<double> r(Values.size());
		for (int i = 0; i < this->Values.size(); i++) {
			r[i] = std::stod(Values[i]);
		}
		return r;
	}

};
std::ostream& operator<<(std::ostream &output, Column const &column) {
	output << column.Column_Name << "\n";
	for (int i = 0; i < (int)column.Values.size(); i++) {
		output << std::setw(2) <<column.Values[i] << "\n";
	}
	return output;
}
std::string &Column::operator[](int index) {
		if (index < 0 || index >=(int)this->Values.size()) {
			std::cout << "\nIllegal Index Exception , Index -[" << index << "]\n";
			std::cout << "Valid Index Range Is - [ 0 - " << this->Values.size()-1 << " ]\n";
			std::cout << "Returning Value At Index 1\n\n";
			return this->Values[0];
		}
		else {
			return this->Values[index];
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
 double Column::Column_Mean() {
	 if (this->Column_Type != 0) {
		 std::cout << "Error ! - Column Is Not Numeric!\n";
		 return-1;
	 }
	 double mean = 0;
	 for (unsigned i = 0; i < this->Values.size(); i++) {
		 mean += std::stod(this->Values[i]);
	 }
	 return mean / this->Values.size();
 }
 double Column::Column_Median() {
	 if (this->Column_Type != 0) {
		 std::cout << "Error ! - Column Is Not Numeric!\n";
		 return-1;
	 }
	 std::vector<double> dt(this->Values.size());
	 for (unsigned i = 0; i < this->Values.size(); i++) {
		 dt[i] = std::stod(this->Values[i]);
	 }
	 std::sort(dt.begin(), dt.end());
	 if (dt.size() % 2 == 0) {
		 return (double)(dt[dt.size() / 2] + (double)dt[(dt.size() / 2) - 1]) / 2;
	 }
	 else {

		 return dt[dt.size() / 2];
	 }
 }
 double Column::Column_Standard_Deviation() {
	 if (this->Column_Type != 0) {
		 std::cout << "Error ! - Column Is Not Numeric!\n";
		 return-1;
	 }
	 std::vector<double> dt(this->Values.size());

	 for (unsigned i = 0; i < this->Values.size(); i++) {
		 dt[i] = std::stod(this->Values[i]);
	 }
	 double mean = this->Column_Mean();
	 double sum = 0;
	 for (unsigned i = 0; i < dt.size(); i++) {
		 sum += ((dt[i] - mean)*(dt[i] - mean));
	 }
	 sum /= (double)dt.size();
	 sum = std::sqrt(sum);
	 return sum;
 }
 double Column::Column_Variance() {
	 if (this->Column_Type != 0) {
		 std::cout << "Error ! - Column Is Not Numeric!\n";
		 return-1;
	 }
	 double Deviation = this->Column_Standard_Deviation();
	 return SQUARE_IT(Deviation);
 }
 double Column::Max_Value() {
	 if (this->Column_Type != 0) {
		 std::cout << "Error ! - Column Is Not Numeric!\n";
		 return-1;
	 }
	 double max = std::numeric_limits<double>::min();
	 double t;
	 for (unsigned i = 0; i < this->Values.size(); i++) {
		 t = std::stod(this->Values[i]);
		 if (t > max) {
			 max = t;
		 }
	 }
	 return max;
 }
 double Column::Min_Value() {
	 if (this->Column_Type != 0) {
		 std::cout << "Error ! - Column Is Not Numeric!\n";
		 return-1;
	 }
	 double min = std::numeric_limits<double>::max();
	 double t;
	 for (unsigned i = 0; i < this->Values.size(); i++) {
		 t = std::stod(this->Values[i]);
		 if (t < min) {
			 min = t;
		 }
	 }
	 return min;
 }
 void Column::Fill_Missing_With_Column_Mean() {
	 if (this->Column_Type != 0) {
		 std::cout << "Error ! - Column Is Not Numeric!\n";
		 return;
	 }
	 if (this->Missing == 0) {
		 std::cout << "There Are No Missing Values At Column [" << this->Column_Name << "]\n";
		 return;
	 }
	 double mean = 0;
	 std::vector<int> missing;
	 for (unsigned i = 0; i < this->Values.size(); i++) {
		 if (Values[i] == "") {
			 missing.push_back(i);
		 }
		 else {
			 mean += std::stod(this->Values[i]);
		 }
	 }
	 mean /= this->Values.size();
	 for (unsigned i = 0; i < missing.size(); i++) {
		 Values[missing[i]] = std::to_string(mean);
	 }
	 this->Missing = 0;

 }
 void Column::Fill_Missing_With_Column_Median() {
	 if (this->Column_Type != 0) {
		 std::cout << "Error ! - Column Is Not Numeric!\n";
		 return;
	 }
	 if (this->Missing == 0) {
		 std::cout << "There Are No Missing Values At Column [" << this->Column_Name << "]\n";
		 return;
	 }
	 std::vector<double> dt(this->Values.size());
	 std::vector<int> missing;
	 double median;
	 for (unsigned i = 0; i < this->Values.size(); i++) {
		 if (Values[i] == "") {
			 missing.push_back(i);
		 }
		 else {
			 dt[i] = std::stod(this->Values[i]);
		 }
	 }
	 std::sort(dt.begin(), dt.end());
	 if (dt.size() % 2 == 0) {
		 median =(double)(dt[dt.size() / 2] + (double)dt[(dt.size() / 2) - 1]) / 2;
	 }
	 else {

		 median =  dt[dt.size() / 2];
	 }
	 for (unsigned i = 0; i < missing.size(); i++) {
		 Values[missing[i]] = std::to_string(median);
	 }
	 this->Missing = 0;
 
 
 }
 void Column::Fill_Missing_With_Pattern(std::string Pattern) {
	 if (this->Missing == 0) {
		 std::cout << "There Are No Missing Values At Column [" << this->Column_Name << "]\n";
		 return;
	 }

	 for (unsigned i = 0; i < this->Values.size(); i++) {
		 if (Values[i] == "") {
			 Values[i] = Pattern;
		 }
	 }
 }
 void  Column::operator=(Column const &source) {
	 this->Missing = source.Missing;
	 this->Column_Name = source.Column_Name;
	 this->Column_Number = source.Column_Number;
	 this->Values = source.Values;
	 this->Categories = source.Categories;
 }
 Column::Column() {
	 this->Column_Name = "None";
	 this->Column_Type = 0;
	 this->Column_Number = 0;
	 this->Missing = 0;
}




class CSV_File {
private:
	std::fstream RW;
public:
	CSV_File();
	CSV_File(const CSV_File &copy);
	template<class Tp> CSV_File(Matrix<Tp> Source);

	std::vector<Column> Data;
	int Number_Of_Rows;
	int Number_Of_Columns;
	std::vector<std::string> columns;



	friend std::ostream& operator<<(std::ostream &output, CSV_File const &csv);
	void Load_CSV(std::string File_Location);
	void Write_CSV(std::string File_Name);
	Column Get_Column(int Column_Number);
	Column &operator[](std::string);
	Column &operator[](int index);
	void operator=(const CSV_File &copy);
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
	double Column_Covariance(int Column_X, int Column_Y);
	void Remove_Rows_With_Missing_Values();
	

	int Number_Of_Missing();
	
};




CSV_File::CSV_File() {
	this->Number_Of_Columns = 0;
	this->Number_Of_Rows = 0;
	this->Data = std::vector<Column>();
}
CSV_File::CSV_File(const CSV_File &copy) {
	this->Data = copy.Data;
	this->Number_Of_Columns = copy.Number_Of_Columns;
	this->Number_Of_Rows = copy.Number_Of_Rows;
	this->columns = copy.columns;
}
void CSV_File::operator=(const CSV_File &copy) {
	this->Data = copy.Data;
	this->Number_Of_Columns = copy.Number_Of_Columns;
	this->Number_Of_Rows = copy.Number_Of_Rows;
	this->columns = copy.columns;
}

void CSV_File::Write_CSV(std::string File_Name){
	File_Name.append(".csv");
	std::ofstream writer;
	writer.open(File_Name);
	for (int i = 0; i < (int)this->Data.size(); i++) {
		if (i < (int)this->Data.size() - 1) {
			writer << this->Data[i].Column_Name << ",";

		}
		else {
			writer << this->Data[i].Column_Name;

		}
	}
	writer << "\n";
	for (int i = 0; i < (int)this->Data[0].Values.size(); i++) {
		for (int j = 0; j < (int)this->Data.size(); j++) {
		//	writer << this->Get_Value_At(j+1, i+1) << ",";
			if ( j < (int)this->Data.size() - 1) {
				writer << this->Data[j][i] << ",";
			}
			else {
				writer << this->Data[j][i];
			}

		}
		writer << "\n";
	}

	writer.close();
}


Column &CSV_File::operator[](int index) {
		if (index < 0 || index >=(int)this->Data.size()) {
			std::cout << "\nIllegal Index Exception , Index -[" << index << "]\n";
			std::cout << "Valid Index Range Is - [ 0 - " << this->Data.size()-1 << " ]\n";
			std::cout << "Returning Column At Index 1\n\n";
			return this->Data[0];
		}
		else {
			return this->Data[index];
		}
	

}
Column &CSV_File::operator[](std::string column) {
	for (int i = 0; i < this->Number_Of_Columns; i++) {
		if (this->Data[i].Column_Name == column) {
			return this->Data[i];
		}
	}
	Column Empty;
	std::cout << "Error! Column Name Not Found\n";
	std::exception bad_c_name;
	bad_c_name.what();
	throw(bad_c_name);
	
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
	std::ios_base::sync_with_stdio(false);
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
		output << std::setw(csv.Data[i].Column_Name.size()+3) << csv.Data[i].Column_Name;
	}
	output << "\n";
	output << std::left;
	int col = 0;
	for (int i = 0; i < (int)csv.Data[0].Values.size(); i++) {
		output << "[ ";
		for (int j = 0; j < (int)csv.Data.size(); j++) {
			if (csv.Data[j].Values[i] == "") {
				output << std::setw(csv.Data[col].Column_Name.size()) << "$Missing$" << "|";
				col++;
			}
			else {
				output << std::setw(csv.Data[col].Column_Name.size()) << csv.Data[j].Values[i] << " | ";
				col++;
			}
		}
		col = 0;
		output << "] ";
		output << "\n";


	}
	return output;
}

void CSV_File::Load_CSV(std::string File_Location){
	this->RW.open(File_Location);
	std::string sample;
	std::vector<String_Vector> Loaded_Data;
	std::vector<std::string> line;

	while (std::getline(RW, sample,'\n')) {
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
	this->Number_Of_Columns = (int)this->Data.size();
	this->Number_Of_Rows = (int)this->Data[0].Values.size();

	for (int i = 0; i < this->Number_Of_Columns; i++) {
		this->columns.push_back(this->Get_Column(i).Column_Name);
	}
	RW.close();

}

Column CSV_File::Get_Column(int Column_Number) {
	if (Column_Number < 0 || Column_Number >(int)this->Data.size()-1) {
		std::cout << "==== Error! -> [Invaid Column] ====";
		return Column();
	}
	else {
		return this->Data[Column_Number];
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
	new_column.Column_Name = Column_Name;
	new_column.Column_Number = this->Number_Of_Columns;
	new_column.Missing = 0;
	if ((int)Values.size() > this->Number_Of_Rows) {
		std::size_t gap =Values.size() - this->Number_Of_Rows;
		for (int i = 0; i < gap; i++) {
			this->Add_Row();
		}
		new_column.Values = Values;
		this->Data.push_back(new_column);
		this->Number_Of_Columns++;

	}
	else {
		if ((int)Values.size() < this->Number_Of_Rows) {
			new_column.Values = Values;
			std::size_t gap = this->Number_Of_Rows - Values.size();
			new_column.Missing = (int)gap;
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
double CSV_File::Column_Covariance(int Column_X, int Column_Y) {
	double covar = 0;
	double x_mean = this->Data[Column_X-1].Column_Mean();
	double y_mean = this->Data[Column_Y-1].Column_Mean();
	for (int i = 1; i <= this->Number_Of_Rows; i++) {
		covar += (std::stod(this->Data[Column_X-1][i]) - x_mean)*(std::stod(this->Data[Column_Y-1][i]) - y_mean);
	}
	return covar / (this->Number_Of_Rows);
}
void CSV_File::Remove_Rows_With_Missing_Values() {
	std::string ms;
	std::vector<int> Miss;

	for (int j = 0; j < this->Number_Of_Columns; j++) {
		for (int i = 1; i <= this->Number_Of_Rows; i++) {
			ms = this->Data[j][i];
			if (ms == "") {
				Miss.push_back(i);
			}
		}
	}
	if (Miss.size() == 0)
	{
		return;
	}
	else
	{
		std::sort(Miss.begin(), Miss.end());
		Miss.erase(std::unique(Miss.begin(), Miss.end()), Miss.end());
		for (unsigned i = 0; i < Miss.size(); i++) {
			this->Remove_Row(Miss[i] - (i));
		}


	}
}
int CSV_File::Number_Of_Missing() {
	int missing = 0;
	for (int i = 0; i < this->Number_Of_Columns; i++) {
		missing += this->Data[i].Missing;
	}
	return missing;
}




template<class Tp> CSV_File::CSV_File(Matrix<Tp> Source) {
	this->Data = std::vector<Column>(Source.Cols);
	for (int i = 0; i < Source.Cols; i++) {
		String_Vector x;
		for (int j = 0; j < Source.Rows; j++) {
			x.push_back(std::to_string(Source[j][i]));
		}
		Data[i].Values = x;
		Data[i].Column_Name = "Col - " + i;
		Data[i].Column_Type = 0;
	}
	this->Number_Of_Columns = Source.Cols;
	this->Number_Of_Rows = Source.Rows;
}




class Pickaxe {
private:
	Matrix<double> Step_Gradient(Matrix<double> Current_Weights, double Learning_Rate, std::vector<int> Columns_Of_Sampels, int Result_Column);
public:
	CSV_File Dataset;

	Pickaxe();
	Pickaxe(CSV_File Dataset);
	Pickaxe(std::string CSV_Filepath);
	Matrix<double> Linear_Regression_Static_Formula(int X_Values_Column_Number, int Y_Values_Column_Number);
	Function Linear_Regression_Static_Formula(int X_Values_Column_Number, int Y_Values_Column_Number,char MODE);
	Matrix<double> Linear_Regression_Gradient_Descent(std::vector<int> Sample_Columns, int Reuslt_Column, double Leaning_Rate, int Number_Of_Iterations);
	std::vector<Point<double> > K_Means(std::vector<int> Target_Columns, int k, int number_of_iterations);
	static double MSE(String_Vector Y, String_Vector Y_Hat);
	static double MSE(std::vector<double> Y, std::vector<double> Y_Hat);
	double MSE(int Column_Y, int Column_Y_Hat);
	static double MAE(String_Vector Y, String_Vector Y_Hat);
	static double MAE(std::vector<double> Y, std::vector<double> Y_Hat);
	double MAE(int Column_Y, int Column_Y_Hat);

	static double R_Squared(String_Vector Y, String_Vector Y_Hat);
	static double Adjusted_R_Squared(String_Vector Y, String_Vector Y_Hat, int Indpendent_Variables);

	Matrix<int> Confusion_Matrix(Matrix<double> Weights, int Result_Column, std::vector<int> Sampled_Rows, double Decision_Boundary);
	void Print_Confusion_Matrix_List_Of_Rates(Matrix<double> Confusion_Matrix, int Number_Sampled);
	std::string KNN(int K, std::vector<double> Test_Values, std::vector<int> Sample_Columns, int Result_Column);
	Matrix<double> KNN(int K, std::vector<double> Test_Values, std::vector<int> Sample_Columns);
	Matrix<double> KNN(int K, Pickaxe Test_Dataset, std::vector<int> Sample_Columns);
	Matrix<double> Get_Variance_Covariance_Matirx(std::vector<int> Sample_Columns);
	Matrix<double> Logistic_Regression(std::vector<int> Value_Column_Numbers, int Binary_Category_Number, int number_of_iterations, double learning_rate);
	void Validate_Linear_Regression(Matrix<double> LR, int True_Column, std::vector<int> Samples_Column_Numbers);
	void Validate_Logistic_Regression(Matrix<double> LR_Weights, int Binary_Column, std::vector<int> Sampled_Rows);
	static double Get_Pearson_Correlation_Coefficient(String_Vector Y, String_Vector Y_Hat);
	double Get_Pearson_Correlation_Coefficient(int Column_X, int Column_Y);
	static double Get_Spearman_Correlation_Coefficient(String_Vector Y, String_Vector Y_Hat);
	double Get_Spearman_Correlation_Coefficient(int Column_X, int Column_Y);
	static std::vector<double> Rankify(String_Vector Y);
	CSV_File Compute_Column_Correlations(unsigned Correlation_Type);

	Matrix<double> PCA(std::vector<int> Selected_Columns);

};

Pickaxe::Pickaxe() {
} 
Pickaxe::Pickaxe(CSV_File Dataset) {
	this->Dataset = Dataset;
}
Pickaxe::Pickaxe(std::string CSV_Filepath) {
	Dataset.Load_CSV(CSV_Filepath);
}
Matrix<double>  Pickaxe::Linear_Regression_Static_Formula(int X_Values_Column_Number, int Y_Values_Column_Number) {


	if (Dataset.Data[X_Values_Column_Number-1].Column_Type != 0 || Dataset.Data[Y_Values_Column_Number - 1].Column_Type != 0) {
		std::cout<< "Error Both Columns Must Be Numeric!\n";
		return Matrix<double>();
	}

	double sum_xy = 0, sum_xsquared = 0, sum_ysquared = 0, sum_x = 0, sum_y = 0;
	std::cout << "X:  " << Dataset[X_Values_Column_Number][this->Dataset.Number_Of_Rows]<<std::endl;
	for (int i = 0; i < Dataset.Number_Of_Rows; i++) {
		double tx =std::stod(Dataset[X_Values_Column_Number][i]), ty = std::stod(Dataset[Y_Values_Column_Number][i]);
		sum_x += tx;
		sum_y += ty;
		sum_xy += ty * tx;
		sum_xsquared += tx * tx;
		sum_ysquared += ty * ty;
	}
	int nor = Dataset.Number_Of_Rows;
	Matrix<double> res(2, 1);
	res[0][0] = ((sum_y*sum_xsquared - sum_x * sum_xy)) / ((nor*sum_xsquared - sum_x * sum_x));
	res[1][0] = ((nor* sum_xy - sum_x * sum_y)) / ((nor*sum_xsquared - sum_x * sum_x));
	return res;

}
Function  Pickaxe::Linear_Regression_Static_Formula(int X_Values_Column_Number, int Y_Values_Column_Number,char MODE) {
	if (MODE == FORMULA) {

		if (Dataset.Data[X_Values_Column_Number - 1].Column_Type != 0 || Dataset.Data[Y_Values_Column_Number - 1].Column_Type != 0) {
			std::cout << "Error Both Columns Must Be Numeric!\n";
			return Function();
		}

		double sum_xy = 0, sum_xsquared = 0, sum_ysquared = 0, sum_x = 0, sum_y = 0;
		for (int i = 1; i <= Dataset.Number_Of_Rows; i++) {
			double tx = std::stod(Dataset[X_Values_Column_Number][i]), ty = std::stod(Dataset[Y_Values_Column_Number][i]);
			sum_x += tx;
			sum_y += ty;
			sum_xy += ty * tx;
			sum_xsquared += tx * tx;
			sum_ysquared += ty * ty;
		}
		int nor = Dataset.Number_Of_Rows;
		double ax = ((sum_y*sum_xsquared - sum_x * sum_xy)) / ((nor*sum_xsquared - sum_x * sum_x));
		double b = ((nor* sum_xy - sum_x * sum_y)) / ((nor*sum_xsquared - sum_x * sum_x));
		std::string fr = (std::to_string((int)ax) + "*" + "*x^1 + " + std::to_string((int)b) + "*x^0");
		Function res(fr.c_str());
		return res;
	}
	else {
		return Function("0");
	}
}
Matrix<double> Pickaxe::Step_Gradient(Matrix<double> Current_Weights, double Learning_Rate, std::vector<int> Columns_Of_Sampels, int Result_Column) {
	std::size_t nof = Columns_Of_Sampels.size() + 1;
	Matrix<double> Gradients((int)nof, 1);
	Matrix<double> Teta(Current_Weights);
	Matrix<double> Values((int)nof, 1);
	Matrix<double> Prediction(Current_Weights);
	Matrix<double> h0;
	Teta = Teta.Matrix_Transpose();

	for (int i = 0; i < this->Dataset.Number_Of_Rows; i++) {
		//h0  guess
		Values[0][0] = 1;
		for (int j = 1; j < nof; j++) {
			Values[j][0] = std::stod(this->Dataset[Columns_Of_Sampels[j - 1]][i]);
		}
		h0 = Teta.Dot_Product(Values);

		for (int j = 0; j < nof; j++) {
			Gradients[j][0] += (h0[0][0] - std::stod(this->Dataset[Result_Column][i]))*Values[j][0];

		}

	}
	for (int j = 0; j < nof; j++) {
		Gradients[j][0] *= (Learning_Rate / this->Dataset.Number_Of_Rows);

		Prediction[j][0] = Prediction[j][0] - Gradients[j][0];
	}
	return Prediction;
}
Matrix<double> Pickaxe::Linear_Regression_Gradient_Descent(std::vector<int> Sample_Columns, int Result_Column, double Leaning_Rate, int Number_Of_Iterations) {
	//y = mx + b - for slope calculation
	Matrix<double> LE((int)Sample_Columns.size() + 1, 1);
	for (int i = 0; i < Number_Of_Iterations; i++) {
		LE = this->Step_Gradient(LE, Leaning_Rate, Sample_Columns, Result_Column);
	}
	return LE;
}
std::vector<Point<double> > Pickaxe::K_Means(std::vector<int> Target_Columns, int k, int number_of_iterations) {

	if (Target_Columns.size() > 3 || Dataset[Target_Columns[0]].Column_Type != 0, Dataset[Target_Columns[1]].Column_Type != 0) {
		std::cout << "Error In Column Specification\n";
		return std::vector<Point<double> >();
	}

	std::vector<Point<double> > data;
	switch (Target_Columns.size())
	{
	case 2:
		for (int i = 0; i < Dataset.Number_Of_Rows; i++) {
			data.push_back(Point<double>(std::stod(Dataset[Target_Columns[0]][i]), std::stod(Dataset[Target_Columns[1]][i]), 0));
		}
		break;
	case 3:
		for (int i = 0; i < Dataset.Number_Of_Rows; i++) {
			data.push_back(Point<double>(std::stod(Dataset[Target_Columns[0]][i]), std::stod(Dataset[Target_Columns[2]][i]), std::stod(Dataset[Target_Columns[2]][i])));
		}
		break;
	default:
		return std::vector<Point<double> >();
		break;
	}


	Random_Utilitis R;
	std::vector<Point<double> > means(k);

	for (auto& cluster : means) {
		cluster = data[R.Random_INT(0,(int)data.size()-1)];
	}

	std::vector<size_t> assignments(data.size());

	for (std::size_t iteration = 0; iteration < number_of_iterations; ++iteration) {
		// Find assignments.
		for (size_t point = 0; point < data.size(); ++point) {
			double best_distance = std::numeric_limits<double>::max();
			std::size_t best_cluster = 0;
			for (std::size_t cluster = 0; cluster < k; ++cluster) {
				const double distance = Point<double>::get_Squared_Distance(data[point], means[cluster]);
				if (distance < best_distance) {
					best_distance = distance;
					best_cluster = cluster;
				}
			}
			assignments[point] = best_cluster;
		}

		// Sum up and count points for each cluster.
		std::vector<Point<double> > new_means(k);
		std::vector<size_t> counts(k, 0);

		for (std::size_t point = 0; point < data.size(); ++point) {
			const auto cluster = assignments[point];
			new_means[cluster].x += data[point].x;
			new_means[cluster].y += data[point].y;
			new_means[cluster].z += data[point].z;
			counts[cluster] += 1;
		}

		// Divide sums by counts to get new centroids.
		for (std::size_t cluster = 0; cluster < k; ++cluster) {
			// Turn 0/0 into 0/1 to avoid zero division.
			const auto count = std::max<size_t>(1, counts[cluster]);
			means[cluster].x = new_means[cluster].x / count;
			means[cluster].y = new_means[cluster].y / count;
			means[cluster].z = new_means[cluster].z / count;

		}
	}

	return means;
}
double Pickaxe::Pickaxe::MSE(String_Vector Y, String_Vector Y_Hat) {
	if (Y.size() != Y_Hat.size()) {
		std::cout<<"Test Groups Have To Be The Same Size\n";
		return std::numeric_limits<double>::infinity();
	}
	double es = 0;
	for (int i = 0; i < Y.size(); i++) {
		es += std::pow(std::stod(Y[i]) - std::stod(Y_Hat[i]), 2);
	}
	es /= Y.size();
	return es;
}
double Pickaxe::MSE(std::vector<double> Y, std::vector<double> Y_Hat) {
	if (Y.size() != Y_Hat.size()) {
		std::cout << "Test Groups Have To Be The Same Size\n";
		return std::numeric_limits<double>::infinity();
	}
	double es = 0;
	for (int i = 0; i < Y.size(); i++) {
		es += std::pow((Y[i]) - (Y_Hat[i]), 2);
	}
	es /= Y.size();
	return es;
}
double Pickaxe::MSE(int Column_Y, int Column_Y_Hat) {
	double es = 0;
	for (int i = 0; i < Dataset.Number_Of_Rows; i++) {
		es += std::pow(std::stod(Dataset[Column_Y][i]) - std::stod(Dataset[Column_Y_Hat][i]), 2);
		std::cout << std::stold(Dataset[Column_Y][i]) << std::endl;
	}
	es /= (Dataset.Number_Of_Rows);
	return es;
}
double Pickaxe::MAE(String_Vector Y, String_Vector Y_Hat) {
	if (Y.size() != Y_Hat.size()) {
		std::cout << "Test Groups Have To Be The Same Size\n";
		return std::numeric_limits<double>::infinity();
	}
	double es = 0;
	for (int i = 0; i < Y.size(); i++) {
		es += std::abs(std::stod(Y[i]) - std::stod(Y_Hat[i]));
	}
	es /= Y.size();
	return es;
}
double Pickaxe::MAE(std::vector<double> Y, std::vector<double> Y_Hat) {
	if (Y.size() != Y_Hat.size()) {
		std::cout << "Test Groups Have To Be The Same Size\n";
		return std::numeric_limits<double>::infinity();
	}
	double es = 0;
	for (int i = 0; i < Y.size(); i++) {
		es += std::abs((Y[i]) - (Y_Hat[i]));
	}
	es /= Y.size();
	return es;
}
double Pickaxe::MAE(int Column_Y, int Column_Y_Hat) {
	double es = 0;
	for (int i = 1; i <= Dataset.Number_Of_Rows; i++) {
		es += std::abs(std::stod(Dataset[Column_Y][i]) - std::stod(Dataset[Column_Y_Hat][i]));
	}
	es /= Dataset.Number_Of_Rows;
	return es;
}
std::vector<double> Pickaxe::Rankify(String_Vector Y) {
	std::vector<double> Rank_Y(Y.size());
	for (unsigned i = 0; i < Rank_Y.size(); i++)
	{
		int r = 1, s = 1;
		for (unsigned j = 0; j < i; j++) {
			if (Y[i] == "") {
				std::cout << "There Are Missing Values In One Of The Columns\nPlease Take Care Of Missing Values\n";
				return std::vector<double>();

			}

			if (std::stod(Y[j]) < std::stod(Y[i])) r++;
			if (std::stod(Y[j]) == std::stod(Y[i])) s++;
		}
		for (unsigned j = i + 1; j < Rank_Y.size(); j++) {
			if (Y[j] == "") {
				std::cout << "There Are Missing Values In One Of The Columns\nPlease Take Care Of Missing Values\n";
				return std::vector<double>();

			}
			if (std::stod(Y[j]) < std::stod(Y[i])) r++;
			if (std::stod(Y[j]) == std::stod(Y[i])) s++;
		}
		Rank_Y[i] = r + (s - 1) * 0.5;
	}
	return Rank_Y;
}
double Pickaxe::Get_Spearman_Correlation_Coefficient(String_Vector Y, String_Vector Y_Hat) {
	double SCC = 0;
	std::vector<double> rankY = Pickaxe::Rankify(Y), rankYhat = Pickaxe::Rankify(Y_Hat);
	double sum_X = 0, sum_Y = 0, sum_XY = 0;
	double squareSum_X = 0, squareSum_Y = 0;
	for (unsigned i = 0; i < rankY.size(); i++)
	{
		sum_X += +rankY[i];
		sum_Y += rankYhat[i];
		sum_XY += rankY[i] * rankYhat[i];
		squareSum_X += rankY[i] * rankY[i];
		squareSum_Y += rankYhat[i] * rankYhat[i];
	}

	SCC = (rankY.size() * sum_XY - sum_X * sum_Y) / std::sqrt((rankY.size() * squareSum_X - sum_X * sum_X) *  (rankY.size() * squareSum_Y - sum_Y * sum_Y));

	return SCC;
}
double  Pickaxe::Get_Spearman_Correlation_Coefficient(int Column_X, int Column_Y) {
	double SCC = 0;
	std::vector<double> rankY = Pickaxe::Rankify(this->Dataset[Column_X].Values), rankYhat = Pickaxe::Rankify(this->Dataset[Column_Y].Values);
	double sum_X = 0, sum_Y = 0, sum_XY = 0;
	double squareSum_X = 0, squareSum_Y = 0;
	for (unsigned i = 0; i < rankY.size(); i++)
	{
		sum_X += +rankY[i];
		sum_Y += rankYhat[i];
		sum_XY += rankY[i] * rankYhat[i];
		squareSum_X += rankY[i] * rankY[i];
		squareSum_Y += rankYhat[i] * rankYhat[i];
	}

	SCC = (rankY.size() * sum_XY - sum_X * sum_Y) / std::sqrt((rankY.size() * squareSum_X - sum_X * sum_X) *  (rankY.size() * squareSum_Y - sum_Y * sum_Y));

	return SCC;
}
CSV_File Pickaxe::Compute_Column_Correlations(unsigned Correlation_Type) {
	CSV_File result;
	std::vector<int> numeric;
	for (int i = 1; i <= this->Dataset.Number_Of_Columns; i++) {
		if (this->Dataset[i - 1].Column_Type == 0) {
			numeric.push_back(i - 1);
		}
	}
	if (numeric.size() < 2) {
		std::cout << "You Need Minimum 2 Numeric Columns To Find Correlations\n";
		return result;
	}
	else {

		int number_of_correlations = (int)std::round(((numeric.size())*(numeric.size() - 1)) / 2);
		result.Add_Column("Between");
		result.Add_Column("Correlation");

		for (int i = 0; i < number_of_correlations; i++) {
			result.Add_Row();
		}



		int pos = 1;
		for (unsigned i = 0; i < numeric.size(); i++) {
			for (unsigned j = i + 1; j < numeric.size(); j++) {
				result[1][pos] = this->Dataset[numeric[i]].Column_Name + " - " + this->Dataset[numeric[j]].Column_Name;
				if (Correlation_Type = PEARSON) {
					result[2][pos] = std::to_string(this->Get_Pearson_Correlation_Coefficient(numeric[i], numeric[j]));

				}
				else if (Correlation_Type == SPEARMAN) {
					result[2][pos] = std::to_string(this->Get_Pearson_Correlation_Coefficient(numeric[i], numeric[j]));
				}
				else {
					std::cout << "Unknown Correlation Type Please Specify - Pearson/Spearman\n";
					return CSV_File();
				}
				pos++;
			}
		}

		return result;
	}

}
double Pickaxe::Get_Pearson_Correlation_Coefficient(String_Vector Y, String_Vector Y_Hat) {
	double r = 0;
	double sigma_xy = 0, sigma_x = 0, sigma_y = 0, sigma_xs = 0, sigma_ys = 0;
	std::vector<double> y(Y.size()), y_h(Y.size());
	for (unsigned i = 0; i < Y.size(); i++) {
		y[i] = std::stod(Y[i]);
		y_h[i] = std::stod(Y_Hat[i]);
	}

	for (unsigned i = 0; i < Y.size(); i++) {
		sigma_xy += y[i] * y_h[i];
		sigma_x += y[i];
		sigma_y += y_h[i];
		sigma_xs += y[i] * y[i];
		sigma_ys += y_h[i] * y_h[i];
	}
	r = y.size()*sigma_xy - sigma_x * sigma_y;
	r /= std::sqrt(((Y.size()*sigma_xs - sigma_x * sigma_x)) * ((Y.size()*sigma_ys - sigma_y * sigma_y)));
	return r;
}
double  Pickaxe::Get_Pearson_Correlation_Coefficient(int Column_X, int Column_Y) {
	double r = 0;
	double sigma_xy = 0, sigma_x = 0, sigma_y = 0, sigma_xs = 0, sigma_ys = 0;
	std::vector<double>y(this->Dataset.Number_Of_Rows), y_h(this->Dataset.Number_Of_Rows);
	for (int i = 0; i <  this->Dataset.Number_Of_Rows; i++) {
		y[i] = std::stod(this->Dataset[Column_X][i]);
		y_h[i] = std::stod(this->Dataset[Column_Y][i]);
	}

	for (int i = 0; i < this->Dataset.Number_Of_Rows; i++) {
		sigma_xy += y[i] * y_h[i];
		sigma_x += y[i];
		sigma_y += y_h[i];
		sigma_xs += y[i] * y[i];
		sigma_ys += y_h[i] * y_h[i];
	}
	r = y.size()*sigma_xy - sigma_x * sigma_y;
	r /= std::sqrt(((y.size()*sigma_xs - sigma_x * sigma_x)) * ((y.size()*sigma_ys - sigma_y * sigma_y)));
	return r;
}
double Pickaxe::R_Squared(String_Vector Y, String_Vector Y_Hat) {
	double R_S = Pickaxe::Get_Pearson_Correlation_Coefficient(Y, Y_Hat);
	R_S *= R_S;
	return R_S;
}
double Pickaxe::Adjusted_R_Squared(String_Vector Y, String_Vector Y_Hat, int Indpendent_Variables) {
	double R_S = Pickaxe::R_Squared(Y, Y_Hat);
	double ARS = (1.0 - R_S)*(Y.size() - 1);
	ARS /= Y.size() - 1 - Indpendent_Variables;
	ARS = 1.0 - ARS;
	return ARS;
}

Matrix<int> Pickaxe::Confusion_Matrix(Matrix<double> Weights, int Result_Column, std::vector<int> Sampled_Rows, double Decision_Boundary) {
	Matrix<int> CM(2, 2);
	int TruePositive = 0;
	int TrueNegative = 0;
	int FalsePositive = 0;
	int FalseNegative = 0;
	for (int i = 0; i < Dataset.Number_Of_Rows; i++) {
		double pred_s = Weights[0][0];
		for (int j = 1; j <= Sampled_Rows.size(); j++) {
			pred_s += Weights[j][0] * std::stod(Dataset[Sampled_Rows[j - 1]][i]);
		}
		pred_s = Sigmoid(pred_s);
		double actual = std::stod(Dataset[Result_Column][i]);

		if (pred_s > Decision_Boundary) {
			pred_s = 1;
		}
		else {
			pred_s = 0;
		}


		if (pred_s == 0 && actual == 0) {
			TrueNegative++;
		}
		else if (pred_s == 1 && actual == 0) {
			FalsePositive++;
		}
		else if (pred_s == 1 && actual == 1) {
			TruePositive++;
		}
		else if (pred_s == 0 && actual == 1) {
			FalseNegative++;
		}


	}

	CM[0][0] = TruePositive;
	CM[1][0] = FalsePositive;
	CM[0][1] = FalseNegative;
	CM[1][1] = TrueNegative;
	return CM;
}
void Pickaxe::Print_Confusion_Matrix_List_Of_Rates(Matrix<double> Confusion_Matrix, int Number_Sampled) {
	double TP = Confusion_Matrix[0][0];
	double FP = Confusion_Matrix[1][0];
	double TN = Confusion_Matrix[1][1];
	double FN = Confusion_Matrix[0][1];
	printf("Accuracy: [ %f ]\n", ((TP + TN) / Number_Sampled));
	printf("Misclassification Rate: [ %f ]\n", ((FP + FN) / Number_Sampled));
	printf("Sensitivity: [ %f ]\n", ((TP) / Number_Sampled));
	printf("Precision: [ %f ]\n", ((TP) / TP + FP));
	printf("Recall: [ %f ]\n", ((TP) / TP + FN));
}
std::string Pickaxe::KNN(int K, std::vector<double> Test_Values, std::vector<int> Sample_Columns, int Result_Column) {
	std::vector<double> distances(Dataset.Number_Of_Rows);
	std::vector<double> x;

	for (int i = 0; i < Dataset.Number_Of_Rows; i++) {
		x = std::vector<double>(Sample_Columns.size());
		for (int j = 0; j < Sample_Columns.size(); j++) {
			x[j] = std::stod(Dataset[Sample_Columns[j]][i]);
		}
		distances[i - 1] = Euclidean_Distance(x, Test_Values);
	}
	std::vector<int> knn(K);
	double min = std::numeric_limits<double>::max();
	for (int i = 0; i < K; i++) {
		for (int j = 0; j < distances.size(); j++) {
			if (distances[j] < min) {
				min = distances[j];
				knn[i] = j + 1;
			}
		}
		distances[knn[i] - 1] = std::numeric_limits<double>::infinity();
		min = std::numeric_limits<double>::max();
	}

	String_Vector categories = Dataset[Result_Column - 1].Categories;
	std::vector<int> cats(categories.size(),0);

	for (int i = 0; i < K; i++) {
		std::string val = Dataset[Result_Column][knn[i]];
		for (int j = 0; j < categories.size(); j++) {
			if (val == (categories[j])) {
				cats[j]++;
				break;
			}
		}
	}

	std::sort(cats.begin(),cats.end());
	return categories[cats[cats.size() - 1] - 1];
}
Matrix<double> Pickaxe::KNN(int K, std::vector<double> Test_Values, std::vector<int> Sample_Columns) {
	std::vector<double> distances(Dataset.Number_Of_Rows);
	std::vector<double> x;
	Matrix<double> result(K, 2);

	for (int i = 0; i < Dataset.Number_Of_Rows; i++) {
		x = std::vector<double>(Sample_Columns.size());
		for (int j = 0; j < Sample_Columns.size(); j++) {
			x[j] = std::stod(Dataset[Sample_Columns[j]][i]);
		}
		distances[i - 1] = Euclidean_Distance(x, Test_Values);
	}
	std::vector<int> knn(K);
	double min = std::numeric_limits<double>::max();
	for (int i = 0; i < K; i++) {
		for (int j = 0; j < distances.size(); j++) {
			if (distances[j] < min) {
				min = distances[j];
				knn[i] = j + 1;
			}
		}
		result[i][1] = min;
		distances[knn[i] - 1] = std::numeric_limits<double>::infinity();
		min = std::numeric_limits<double>::max();
	}

	for (int i = 0; i < K; i++) {
		result[i][0] = knn[i];

	}
	return result;
}
Matrix<double> Pickaxe::KNN(int K, Pickaxe Test_Dataset, std::vector<int> Sample_Columns) {
	Matrix<double> Final_Res(Test_Dataset.Dataset.Number_Of_Rows, (1 + K));
	std::vector<double> xv(Sample_Columns.size());
	for (int i = 0; i < Test_Dataset.Dataset.Number_Of_Rows; i++) {
		for (int j = 0; j < xv.size(); j++) {
			xv[j] = std::stod(Test_Dataset.Dataset[Sample_Columns[j]][i]);
		}
		Matrix<double> res = this->KNN(K, xv, Sample_Columns);
		for (int j = 1; j < K + 1; j++) {
			Final_Res[i - 1][j] = res[j - 1][0];
		}

		Final_Res[i - 1][0] = i;

	}

	return Final_Res;
}



Matrix<double> Pickaxe::Get_Variance_Covariance_Matirx(std::vector<int> Sample_Columns){
	
	Matrix<double> data(Dataset.Number_Of_Rows, (int)Sample_Columns.size());
	Matrix<double> dataOnes(Dataset.Number_Of_Rows, 1);

	for (int i = 0; i < Dataset.Number_Of_Rows; i++) {
		dataOnes[i][0] = 1;

	}
	Matrix<double> dataOnestag(dataOnes);
	dataOnestag = dataOnestag.Matrix_Transpose();

	for (int i = 0; i < Dataset.Number_Of_Rows; i++) {
		for (int j = 0; j < Sample_Columns.size(); j++) {
			data[i][j] = std::stod(Dataset[Sample_Columns[j]][i]);
		}
	}

	//transformation a = A - 11'x/n
	dataOnestag = dataOnes.Dot_Product(dataOnestag);
	dataOnestag = dataOnestag.Dot_Product(data);
	dataOnestag/(Dataset.Number_Of_Rows);
	data = data - (dataOnestag);
	dataOnes = Matrix<double>(data);
	dataOnes = dataOnes.Matrix_Transpose();
	//(a*a')/n 
	data = dataOnes.Dot_Product(data);
	data/(Dataset.Number_Of_Rows);

	return data;
}


Matrix<double> Pickaxe::Logistic_Regression(std::vector<int> Value_Column_Numbers, int Binary_Category_Number, int number_of_iterations, double learning_rate) {
	CSV_File holder;
	for (int i = 0; i < Value_Column_Numbers.size(); i++) {
		holder.Add_Column(Dataset[Value_Column_Numbers[i]].Column_Name, Dataset[Value_Column_Numbers[i]].Values);

	}
	String_Vector Binary_Category =Dataset[Binary_Category_Number].Values;

	int number_of_features = holder.Number_Of_Columns + 1;
	int number_of_sampels = holder.Number_Of_Rows;
	Matrix<double> Weights(number_of_features, 1);

	for (int i = 0; i < number_of_features; i++) {
		Weights[i][0] = 0;
	}

	std::vector<double> Predictions;
	double Cost = 0, Lowest_Cost = std::numeric_limits<double>::max();
	for (int i = 0; i < number_of_iterations; i++) {

		//predictions
		for (int j = 0; j < number_of_sampels; j++) {

			double y = std::stod(Binary_Category[j]);
			double s_pred = Weights[0][0];
			for (int k = 1; k < number_of_features; k++) {
				double instance = std::stod(holder[k][j+1]);

				s_pred += Weights[k][0] * instance;
			}
			double pred = Sigmoid(s_pred);
			Predictions.push_back(pred);


			//updating weights
			Weights[0][0] = Weights[0][0] + learning_rate * (y - pred)*pred*(1.0 - pred)*1.0;

			for (int k = 1; k < number_of_features; k++) {
				Weights[k][0] = Weights[k][0] + learning_rate * (y - pred)*pred*(1.0 - pred)*std::stod(holder[k][j+1]);
			}

			//B0 = B0 + learning_rate*(y - pred)*pred*(1.0-pred)*1.0;
			//B1 = B1 + learning_rate*(y - pred)*pred*(1.0-pred)*Double.parseDouble(Values.get(j));
		}

		//cost

		for (int j = 0; j < Predictions.size(); j++) {
			double y = std::stod(Binary_Category[j]);
			Cost += -y * std::log(Predictions[j]) - (1 - y)*std::log(Predictions[j]);
		}

		Cost *= 1.0 / Predictions.size();

		if (Cost < Lowest_Cost) {
			Lowest_Cost = Cost;
		}


		Predictions.clear();
		Predictions = std::vector<double>();
	}



	return Weights;
}

void Pickaxe::Validate_Linear_Regression(Matrix<double> LR, int True_Column, std::vector<int> Samples_Column_Numbers) {
	for (int i = 1; i <= this->Dataset.Number_Of_Rows; i++) {
		double pred = LR[0][0];
		for (int j = 1; j < Samples_Column_Numbers.size() + 1; j++) {
			pred += LR[j][0] * std::stod(Dataset[Samples_Column_Numbers[j - 1]][i]);
		}
		printf("Predictions: { %f } ---- Actual: { %s }\n", pred, this->Dataset[True_Column][i].c_str());
	}
}

void Pickaxe::Validate_Logistic_Regression(Matrix<double> LR_Weights, int Binary_Column, std::vector<int> Sampled_Rows) {
	for (int i = 1; i <= Dataset.Number_Of_Rows; i++) {
		double pred_s = LR_Weights[0][0];
		for (int j = 1; j <= Sampled_Rows.size(); j++) {
			pred_s += LR_Weights[j][0] * std::stod(Dataset[Sampled_Rows[j - 1]][i]);
		}

		pred_s = Sigmoid(pred_s);

		printf("\n====================\n");
		printf("Values: ");
		for (int z = 0; z < Sampled_Rows.size(); z++) {
			printf("{%s} ", Dataset[Sampled_Rows[z]][i].c_str());

		}
		printf("\n====================\n");

		printf("Prediction:[%f] ---- Actual: [%s] ", pred_s, Dataset[Binary_Column][i].c_str());
		printf("\n---------------------\n");


	}
}


Matrix<double> Pickaxe::PCA(std::vector<int> Selected_Columns) {
	std::vector<double> means(Selected_Columns.size());
	for (int i = 0; i < Selected_Columns.size(); i++) {
		means[i] = Dataset[Selected_Columns[i]].Column_Mean();
	}
	Matrix<double> data(Dataset.Number_Of_Rows, (int)Selected_Columns.size());

	for (int i = 0; i < Dataset.Number_Of_Rows; i++) {
		for (int j = 0; j < Selected_Columns.size(); j++) {
			data[i][j] = std::stod(Dataset[Selected_Columns[j]][i]) - means[j];
		}
	}
	CSV_File tocsv(data);
	Pickaxe body;
	body.Dataset = tocsv;

	Matrix<double> Var_CO_Var = body.Get_Variance_Covariance_Matirx(Selected_Columns);

	Matrix<double> EV = Var_CO_Var.Get_Eigen_Vectors();

	data = data.Dot_Product(EV);
	
	//cout << data;

	return data;
}



//NN

class Neuron {
public:
	std::vector<double> weights;
	std::vector<double> backprop_weights;
	double gradient;
	double bias;
	double value = 0;

	static double minWeightValue;
	static double maxWeightValue;

	// Constructor for the hidden / output neurons
	Neuron(std::vector<double> weights, double bias, double value);
	// Constructor for the input neurons
	Neuron(double value);
	// Static function to set min and max weight for all variables
	static void setRangeWeight(double min, double max);
	// Function used at the end of the backprop to switch the calculated value in the
	// cache weight in the weights
	void Update_Weights();

};
double Neuron::maxWeightValue;
double Neuron::minWeightValue;
Neuron::Neuron(std::vector<double> weights, double bias, double value) {
	this->weights = weights;
	this->bias = bias;
	this->backprop_weights = this->weights;
	this->gradient = 0;
	this->value = value;
}
Neuron::Neuron(double value) {
	this->bias = -1;
	this->backprop_weights = this->weights;
	this->gradient = -1;
	this->value = value;
}
void Neuron::setRangeWeight(double min, double max) {
	minWeightValue = min;
	maxWeightValue = max;
	
}
void Neuron::Update_Weights() {
	this->weights = this->backprop_weights;
}

class Neuron_Layer {
public:
	bool ASigmoid = true;
	bool ATanh = false;
	bool ARelu = false;
	 std::vector<Neuron> neurons;
	 Neuron_Layer() {};


	// Constructor for the hidden and output layer
	 Neuron_Layer(int Neurons_Connections, int Number_of_Neurons);


	// Constructor for the input layer
	 Neuron_Layer(std::vector<double> input);

	void Set_Activation_Function(std::string Activation_Function_Name);


};

class Data_Cartridge {
	public:
	std::vector<double> Data;
	std::vector<double> Expected_Output;

	Data_Cartridge(std::vector<double> row_values, std::vector<double> result_values) {
		this->Data = row_values;
		this->Expected_Output = result_values;
	}

};

class Training_Data {
public:
	std::vector<Data_Cartridge> Data;


	Training_Data();

	void Load_Data(CSV_File Data_Set, std::vector<int> Data_Column_Number, std::vector<int> Result_Column);
	void Load_Data(std::vector<std::vector<double> > Data, std::vector< std::vector<double> > Result);
	void Load_Data(std::vector<double> Data, std::vector<double > Result);
	void Load_Data(Matrix<double> Data, Matrix<double> Result);



	void Print_Data();

};

Training_Data::Training_Data() {

}

void Training_Data::Load_Data(CSV_File Data_Set, std::vector<int> Data_Column_Number, std::vector<int> Result_Column) {
	this->Data.reserve(Data_Set.Number_Of_Rows);

	for (int i = 0; i < Data_Set.Number_Of_Rows; i++) {
		std::vector<double> row_values;
		row_values.reserve(Data_Column_Number.size());
		std::vector<double> result_values;
		result_values.reserve(Result_Column.size());

		for (int j = 0; j < Data_Column_Number.size(); j++) {
			row_values.push_back(std::stod(Data_Set[Data_Column_Number[j]][i]));
		}

		for (int q = 0; q < Result_Column.size(); q++) {
			result_values.push_back(std::stod(Data_Set[Result_Column[q]][i]));
		}
		this->Data.push_back(Data_Cartridge(row_values, result_values));

	}
}
void Training_Data::Load_Data(std::vector<std::vector<double> > Data, std::vector< std::vector<double> > Result) {
	for (int i = 0; i < Data.size(); i++) {
		this->Data.push_back(Data_Cartridge(Data[i], Result[i]));

	}
}
void Training_Data::Load_Data(Matrix<double> Data, Matrix<double> Result) {
	for (int i = 0; i < Data.Rows; i++) {
		this->Data.push_back(Data_Cartridge(Data[i], Result[i]));

	}
}
void Training_Data::Load_Data(std::vector<double> Data, std::vector<double > Result) {
	this->Data.push_back(Data_Cartridge(Data, Result));

}


void Training_Data::Print_Data() {
	for (int i = 0; i < Data.size(); i++) {
		std::cout << "Data: ";
		for (int j = 0; j < this->Data[i].Data.size(); j++) {
			std::cout<<  " |" << this->Data[i].Data[j];
		}
		std::cout << " | Target: ";
		for (int j = 0; j < this->Data[i].Expected_Output.size(); j++) {
			std::cout << this->Data[i].Expected_Output[j] << " ";
		}
		std::cout << "\n";
	}
}




void Neuron_Layer::Set_Activation_Function(std::string Activation_Function_Name) {
	if (Activation_Function_Name == ("Relu")) {
		this->ASigmoid = false;
		this->ATanh = false;
		this->ARelu = true;
	}
	else if (Activation_Function_Name == ("Sigmoid")) {
		this->ASigmoid = true;
		this->ATanh = false;
		this->ARelu = false;
	}
	else if (Activation_Function_Name == ("Tanh")) {
		this->ASigmoid = false;
		this->ATanh = true;
		this->ARelu = false;
	}
}

Neuron_Layer::Neuron_Layer(int Neurons_Connections, int Number_of_Neurons) {
	this->neurons.reserve(Neurons_Connections);

	for (int i = 0; i < Number_of_Neurons; i++) {
		std::vector<double> weights;
		weights.reserve(Neurons_Connections);

		for (int j = 0; j < Neurons_Connections; j++) {
			weights.push_back(Random_Utilitis::Random_DOUBLE(Neuron::minWeightValue, Neuron::maxWeightValue));
			
		}
		neurons.push_back(Neuron(weights, 1, Random_Utilitis::Random_DOUBLE(Neuron::minWeightValue, Neuron::maxWeightValue)));
	}
}


Neuron_Layer::Neuron_Layer(std::vector<double> input) {
	this->neurons.reserve(input.size());
	for (int i = 0; i < input.size(); i++) {
		this->neurons.push_back(Neuron(input[i]));
	}
}


class Neural_Net {

	std::vector<Neuron_Layer> Layers;
	Training_Data Training;

	//for model constructing and reconstructing 
	std::vector<int> Topology_Of_Neurons;
	std::vector<std::string> activations;
	double Learning_Rate;
	double Min_Weight;
	double Max_Weight;
public:

	Neural_Net(std::vector<int> Topology_Of_Neurons, Training_Data Training_Data, double Learning_Rate, double Min_Weight, double Max_Weight);
	Neural_Net(std::string Brain_File_Saved_Model);
	Neural_Net();
	void Add_Layer(int Number_Of_Neurons,std::string Activation);
	void Add_Layer(int Number_Of_Neurons);
	void Compile(double learning_rate, double min_weight, double max_weight, Training_Data training);
	Matrix<double> Predict(Matrix<double> inputs);
	std::vector<double>  Predict(std::vector<double> inputs);


	void Forward_Data(std::vector<double> inputs);
	// This function sums up all the gradient connecting a given neuron in a given layer
	double sumGradient(int neuron_index, int layer_index);
	void Back_Propagate(Data_Cartridge training_Data);
	void Start_Training(int Training_Iterations);
	std::vector<double> Get_Output_Values();
	void Print_Output_Neurons_Values();
	void Save_Model(std::string Model_Name);
	void Set_Activation_Function(std::string Activation_Function_Name, int Layer_Number);



};


Neural_Net::Neural_Net(std::vector<int> Topology_Of_Neurons, Training_Data Training_Data, double Learning_Rate, double Min_Weight, double Max_Weight) {
	Neuron::setRangeWeight(Min_Weight, Max_Weight);
	Layers.reserve(Topology_Of_Neurons.size() / 2 + 1);
	for (int i = 0; i < Topology_Of_Neurons.size() / 2 + 1; i++) {
		Layers.push_back(Neuron_Layer());
	}
	int j = 1;
	for (int i = 0; i < Topology_Of_Neurons.size(); i += 2) {
		Layers[j] = Neuron_Layer(Topology_Of_Neurons[i], Topology_Of_Neurons[i + 1]);
		j++;
	}

	this->Training = Training_Data;
	this->Learning_Rate = Learning_Rate;
	this->Min_Weight = Min_Weight;
	this->Max_Weight = Max_Weight;
	this->Topology_Of_Neurons = Topology_Of_Neurons;
}

Neural_Net::Neural_Net(std::string Brain_File_Saved_Model) {
	std::fstream in;
	if (!Brain_File_Saved_Model.find(".Axe")) {
		std::cout<<"Unkown File Format Specifed, Please Use Only Java_Brain Supported Formats";
		return;
	}
	in.open((Brain_File_Saved_Model));
		int sampler;
		std::string SAMPLER;
		std::getline(in, SAMPLER, '\n');
		this->Learning_Rate = std::stod(SAMPLER);
		std::getline(in, SAMPLER, '\n');
		this->Max_Weight = std::stod(SAMPLER);
		std::getline(in, SAMPLER, '\n');
		this->Min_Weight = std::stod(SAMPLER);
		std::getline(in, SAMPLER, '\n');
		sampler = std::stoi(SAMPLER);
		this->Topology_Of_Neurons = std::vector<int>(sampler, 0);
		for (int i = 0; i < sampler; i++) {
			std::getline(in, SAMPLER, '\n');
			this->Topology_Of_Neurons[i] = std::stoi(SAMPLER);
		}

		Neuron::setRangeWeight(Min_Weight, Max_Weight);
		
		Layers = std::vector<Neuron_Layer>(Topology_Of_Neurons.size() / 2 + 1,Neuron_Layer());
		Layers[0] = Neuron_Layer();
		int j = 1;
		for (int i = 0; i < Topology_Of_Neurons.size(); i += 2) {
			Layers[j] = Neuron_Layer(Topology_Of_Neurons[i], Topology_Of_Neurons[i + 1]);
			j++;
		}

		//reading layer structure

		for (int layer = 1; layer < this->Layers.size(); layer++) {
			for (int ners = 0; ners < this->Layers[layer].neurons.size(); ners++) {

				std::getline(in, SAMPLER, '\n');
				this->Layers[layer].neurons[ners].bias = std::stod(SAMPLER);
				std::getline(in, SAMPLER, '\n');
				this->Layers[layer].neurons[ners].gradient = std::stod(SAMPLER);
				std::getline(in, SAMPLER, '\n');
				this->Layers[layer].neurons[ners].value = std::stod(SAMPLER);

				std::getline(in, SAMPLER, '\n');
				this->Layers[layer].neurons[ners].backprop_weights = std::vector<double>(std::stoi(SAMPLER),0);
				double backprop = (double)this->Layers[layer].neurons[ners].backprop_weights.size();
				for (int m = 0; m < backprop; m++) {
					std::getline(in, SAMPLER, '\n');
					this->Layers[layer].neurons[ners].backprop_weights[m] = std::stod(SAMPLER);
				}
				std::getline(in, SAMPLER, '\n');
				this->Layers[layer].neurons[ners].weights = std::vector<double>(std::stoi(SAMPLER), 0);
				double weight = (double)this->Layers[layer].neurons[ners].weights.size();
				for (int m = 0; m < weight; m++) {
					std::getline(in, SAMPLER, '\n');
					this->Layers[layer].neurons[ners].weights[m] = std::stod(SAMPLER);
				}




			}
		}








	

};

void Neural_Net::Forward_Data(std::vector<double> inputs) {
	// First bring the inputs into the input layer layers[0]
	Layers[0] = Neuron_Layer(inputs);
	for (int i = 1; i < Layers.size(); i++) {
		for (int j = 0; j < Layers[i].neurons.size(); j++) {
			//summing the weights
			double sum = 0;
			for (int k = 0; k < Layers[i - 1].neurons.size(); k++) {
				sum += Layers[i - 1].neurons[k].value*Layers[i].neurons[j].weights[k];
			}
			sum += Layers[i].neurons[j].bias; // add in the bias 

			if (Layers[i].ASigmoid == true) {
				Layers[i].neurons[j].value = Sigmoid(sum);
			}
			else if (Layers[i].ATanh == true) {
				Layers[i].neurons[j].value = std::tanh(sum);

			}
			else if (Layers[i].ARelu == true) {
				Layers[i].neurons[j].value = Rectified(sum);
			}
			else {
				std::cout << "No Activaition Function Selected";
			}
		}
	}
}
// This function sums up all the gradient connecting a given neuron in a given layer
double Neural_Net::sumGradient(int neuron_index, int layer_index) {
	double gradient_sum = 0;
	Neuron_Layer current_layer = Layers[layer_index];
	for (int i = 0; i < current_layer.neurons.size(); i++) {
		Neuron current_neuron = current_layer.neurons[i];
		gradient_sum += current_neuron.weights[neuron_index] * current_neuron.gradient;
	}
	return gradient_sum;
}
void Neural_Net::Back_Propagate(Data_Cartridge training_Data) {

	int number_layers = (int)Layers.size();
	int out_index = number_layers - 1;
	double derivative = 0;
	double delta = 0;

	// Update the output layers 
// For each output
	for (int i = 0; i < Layers[out_index].neurons.size(); i++) {
		// and for each of their weights
		double output = Layers[out_index].neurons[i].value;
		double target = training_Data.Expected_Output[i];



		if (Layers[out_index].ASigmoid == true) {
			derivative = (output - target);
			delta = derivative * (output*(1 - output));
		}
		else if (Layers[out_index].ATanh == true) {
			derivative = (output - target);
			delta = derivative * (1 - output * output);
		}
		else if (Layers[out_index].ARelu == true) {
			if (output >= 0) {
				delta = (output - target);
			}
			else {
				delta = 0.01;
			}
		}



		Layers[out_index].neurons[i].gradient = delta;
		for (int j = 0; j < Layers[out_index].neurons[i].weights.size(); j++) {
			double previous_output = Layers[out_index - 1].neurons[j].value;
			double error = delta * previous_output;
			Layers[out_index].neurons[i].backprop_weights[j] = Layers[out_index].neurons[i].weights[j] - Learning_Rate * error;
		}
	}


	//Update all the subsequent hidden layers
	for (int i = out_index - 1; i > 0; i--) {
		// For all neurons in that layers
		for (int j = 0; j < Layers[i].neurons.size(); j++) {
			double output = Layers[i].neurons[j].value;
			double gradient_sum = sumGradient(j, i + 1);




			if (Layers[i].ASigmoid == true) {
				delta = gradient_sum * (output*(1 - output));
			}
			else if (Layers[i].ATanh == true) {
				delta = gradient_sum * (1 - output * output);

			}
			else if (Layers[i].ARelu == true) {
				if (output >= 0) {
					delta = gradient_sum;
				}
				else {
					delta = gradient_sum * 0.01;
				}
			}


			Layers[i].neurons[j].gradient = delta;
			// And for all their weights
			for (int k = 0; k < Layers[i].neurons[j].weights.size(); k++) {
				double previous_output = Layers[i - 1].neurons[k].value;
				double error = delta * previous_output;



				Layers[i].neurons[j].backprop_weights[k] = Layers[i].neurons[j].weights[k] - Learning_Rate * error;
			}
		}
	}


	// Here we do another pass where we update all the weights
	for (int i = 0; i < Layers.size(); i++) {
		for (int j = 0; j < Layers[i].neurons.size(); j++) {
			Layers[i].neurons[j].Update_Weights();
		}
	}

}
void Neural_Net::Start_Training(int Training_Iterations) {
	for (int i = 0; i < Training_Iterations; i++) {
		for (int j = 0; j < Training.Data.size(); j++) {
			Forward_Data(Training.Data[j].Data);
			Back_Propagate(Training.Data[j]);
		}
#ifdef PRINT_EPOCH_STATUS
		std::cout << "== Epoch Status : [" << (i+1) << " / " << Training_Iterations << " ] ==\n\n";

#endif // PRINT_EPOCH_STATUS

	}
	std::cout << "\n";
}
std::vector<double> Neural_Net::Get_Output_Values() {
	std::vector<double> results;
	results.reserve(this->Layers[Layers.size() - 1].neurons.size());
	
	for (int m = 0; m < this->Layers[Layers.size() - 1].neurons.size(); m++) {
		results.push_back(this->Layers[Layers.size() - 1].neurons[m].value);
	}
	return results;
}
void Neural_Net::Print_Output_Neurons_Values() {
	std::cout << "============\n";
	std::cout << ("   Output") << "\n";
	std::cout << ("============\n");
	for (int j = 0; j < Layers[Layers.size() - 1].neurons.size(); j++) {
		std::cout << Layers[Layers.size() - 1].neurons[j].value;
	}
}

void Neural_Net::Save_Model(std::string Model_Name) {

	std::ofstream writer((Model_Name+".Axe"));
	writer << this->Learning_Rate << '\n';
	writer << Max_Weight << '\n';;
	writer << Min_Weight << '\n';;
	writer << this->Topology_Of_Neurons.size() << '\n';;
		for (int i = 0; i < Topology_Of_Neurons.size(); i++) {
			writer << this->Topology_Of_Neurons[i] << '\n';;
		}

		for (int layer = 1; layer < this->Layers.size(); layer++) {
			for (int ners = 0; ners < this->Layers[layer].neurons.size(); ners++) {

				writer << this->Layers[layer].neurons[ners].bias << '\n';;
				writer << this->Layers[layer].neurons[ners].gradient << '\n';;
				writer << this->Layers[layer].neurons[ners].value << '\n';;

				writer << this->Layers[layer].neurons[ners].backprop_weights.size() << '\n';;
				double backprop = (double)this->Layers[layer].neurons[ners].backprop_weights.size();
				for (int m = 0; m < backprop; m++) {
					writer << this->Layers[layer].neurons[ners].backprop_weights[m] << '\n';;
				}

				writer << this->Layers[layer].neurons[ners].weights.size() << '\n';;
				double weight = (double)this->Layers[layer].neurons[ners].weights.size();
				for (int m = 0; m < weight; m++) {
					writer << this->Layers[layer].neurons[ners].weights[m] << '\n';;
				}




			}
		}




		writer.close();
		std::cout << "\n\nYour Model Was Successfully Saved As: [" + Model_Name + ".Axe]\n";
}
void Neural_Net::Set_Activation_Function(std::string Activation_Function_Name,int Layer_Number) {
	Layers[Layer_Number].Set_Activation_Function(Activation_Function_Name);
}

Neural_Net::Neural_Net() {

}
void Neural_Net::Add_Layer(int Number_Of_Neurons, std::string Activation) {
	if (Topology_Of_Neurons.size() == 1) {
		this->Topology_Of_Neurons.push_back(Number_Of_Neurons);
		this->activations.push_back(Activation);
	}
	else {
		this->Topology_Of_Neurons.push_back(Topology_Of_Neurons[Topology_Of_Neurons.size() - 1]);
		this->Topology_Of_Neurons.push_back(Number_Of_Neurons);

		this->activations.push_back(Activation);
	}
}
void Neural_Net::Add_Layer(int Number_Of_Neurons) {
	this->Topology_Of_Neurons.push_back(Number_Of_Neurons);
}
void Neural_Net::Compile(double learning_rate, double min_weight, double max_weight, Training_Data training) {
	Neuron::setRangeWeight(min_weight, max_weight);
	Layers.reserve(Topology_Of_Neurons.size() / 2 + 1);
	for (int i = 0; i < Topology_Of_Neurons.size() / 2 + 1; i++) {
		Layers.push_back(Neuron_Layer());
	}
	int j = 1;
	for (int i = 0; i < Topology_Of_Neurons.size(); i += 2) {
		Layers[j] = Neuron_Layer(Topology_Of_Neurons[i], Topology_Of_Neurons[i + 1]);
		Layers[j].Set_Activation_Function(this->activations[j - 1]);
		j++;
	}

	this->Training = training;
	this->Learning_Rate = learning_rate;
	this->Min_Weight = min_weight;
	this->Max_Weight = max_weight;
}
Matrix<double> Neural_Net::Predict(Matrix<double> inputs) {
	Matrix<double> Result(inputs.Rows, inputs.Cols + this->Topology_Of_Neurons[this->Topology_Of_Neurons.size()-1]);
	for (int i = 0; i < inputs.Rows; i++) {
		std::vector<double> res = this->Predict(inputs[i]);
		int z = 0;
		for (int j = 0; j < inputs.Cols; j++) {
			Result[i][z++] = inputs[i][j];
		}

		for (int j = 0; j < res.size(); j++) {
			Result[i][z++] = res[j];
		}
	}
	return Result;
}
std::vector<double>  Neural_Net::Predict(std::vector<double> inputs) {
	this->Forward_Data(inputs);
	return Get_Output_Values();

}



class CNN{
public:
	Neural_Net FCN_Network;
	std::vector<int> Topology_Of_Layers;
	std::vector<std::string> activations;



};

}