#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <sstream>
#include <fstream>
#include <iomanip>




typedef std::vector<std::string> String_Vector;

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





class Column {
public:
	std::string Column_Name;
	int Column_Number;
	int Column_Type;
	int Missing;
	String_Vector Values;
	String_Vector Categories;
	
	friend std::ostream& operator<<(std::ostream &output,Column const &column);
};
std::ostream& operator<<(std::ostream &output, Column const &column) {
	output << column.Column_Name << "\n";
	for (int i = 0; i < column.Values.size(); i++) {
		output << std::setw(2) <<column.Values[i] << "\n";
	}
	return output;
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
	std::string Get_Value_At(int Column, int Row);
	void Set_Value_At(int Column, int Row,std::string Item);
	String_Vector Get_Row(int Row_Number);
	void Add_Row();
	void Add_Row(String_Vector Values);
};




CSV_File::CSV_File() {
	
}

void CSV_File::Write_CSV(std::string File_Name){
	File_Name.append(".csv");
	std::ofstream writer;
	writer.open(File_Name);
	for (int i = 0; i < this->Data.size(); i++) {
		writer << this->Data[i].Column_Name << ",";
	}
	writer << "\n";
	for (int i = 0; i < this->Data[0].Values.size(); i++) {
		for (int j = 0; j < this->Data.size(); j++) {
			writer << this->Get_Value_At(j+1, i+1) << ",";

		}
		writer << "\n";
	}

	writer.close();
}

std::string CSV_File::Get_Value_At(int Column, int Row) {
	if (Column < 1 || Row < 1 || Column > this->Data.size() || Row > this->Data[0].Values.size()) {
		std::cout << "==== Error -> [Invalid Row/Column] ====\n";
		return "0";
	}
	return this->Data[Column-1].Values[Row-1];
}

void CSV_File::Set_Value_At(int Column, int Row,std::string Item) {
	if (Column < 1 || Row < 1 || Column > this->Data.size() || Row > this->Data[0].Values.size()) {
		std::cout << "==== Error -> [Invalid Row/Column] ====\n";
		return;
	}
	this->Data[Column-1].Values[Row-1] = Item;
}

String_Vector CSV_File::Get_Row(int Row_Number) {
	if (Row_Number < 1 || Row_Number > this->Number_Of_Rows) {
		std::cout << "==== Error -> [ Invalid Row ] ====\n";
		return String_Vector();
	}
	String_Vector Row;
	for (int i = 1; i <= this->Data.size(); i++) {
		Row.push_back(this->Get_Value_At(i, Row_Number));
	}
	return Row;
}

std::ostream &operator<<(std::ostream &out, String_Vector vec) {
	out << std::left;
	out << "[";
	for (int i = 0; i < vec.size(); i++) {
		if (i < vec.size() - 1) {
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
	for (int i = 0; i < csv.Data.size(); i++) {
		output << std::setw(csv.Data[i].Column_Name.size()+2) << csv.Data[i].Column_Name;
	}
	output << "\n";
	int col = 0;
	for (int i = 0; i < csv.Data[0].Values.size(); i++) {
		for (int j = 0; j < csv.Data.size(); j++) {
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

	for (int i = 0; i < Loaded_Data[0].size();i++) {
		Column temp;
		String_Vector values(Loaded_Data.size()-1);
		int aux = 0;
		int missing = 0;
		for (int j = 1; j < Loaded_Data.size(); j++) {
			values[aux++] = Loaded_Data[j][i];
			if (Loaded_Data[j][i] == "") {
				missing++;
			}
		}
		temp.Column_Number = i;
		temp.Values = values;
		temp.Missing = missing;
		temp.Column_Name = Loaded_Data[0][i];
		this->Data.push_back(temp);
	}
	this->Number_Of_Columns = this->Data.size();
	this->Number_Of_Rows = this->Data[0].Values.size();
	RW.close();

}

Column CSV_File::Get_Column(int Column_Number) {
	if (Column_Number < 1 || Column_Number > this->Data.size()) {
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
			for (i = 0; i < Values.size(); i++) {
				this->Data[i].Values.push_back(Values[i]);
			}
			for (i; i < this->Data.size(); i++) {
				this->Data[i].Values.push_back("0");
			}
			this->Number_Of_Rows++;
		}
		else {
			for (int i = 0; i < this->Data.size(); i++) {
				this->Data[i].Values.push_back(Values[i]);
			}
			this->Number_Of_Rows++;
		}


	}
}