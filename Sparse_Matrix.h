#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <vector>
#include <stdexcept> //so I can use std::out_of_range, I won't use it because I choose not to implement try...catch constructs
#include <iostream> //so I can use std::cout and and std::cerr, this one is in a comment in the function
//where we write a value in the matrix, delite the "//" part and comment the allocation
//if there's a preferene of returning an error instead of allocate a new value



//ASSUMPTIONS: - the number of rows and columns of the matrix can't change
//			   - working with pointers because we want to make changes to the matrix not to its copy


//creating an abstract class SparseMatrix, then implementing the concrete functions for the classes SaprseMatrixCOO and SparseMatrixCSR 

class SparseMatrix {
public:
	virtual ~SparseMatrix() {} //destruction of the matrix when going out of scope,
	//using vitrual because of the pointers


	//getting the number of columns and rows
	virtual int getNumberRows() const = 0;
	virtual int getNumberCols() const = 0;


	//getting the number of nonzeros
	virtual int getNonZeros() const = 0;

	//read entry of the matrix
	virtual double operator() (int NumRow, int NumCol) const= 0;

	//write values in the matrix
	virtual void operator() (double value, int NumRow, int NumCol)= 0;

	//compute Matrix-Vector product (y=Ax)
	virtual std::vector<double> operator*(const std::vector<double>& x) const = 0;

	//Printing the matrix
	virtual void print() const = 0;

	//getting the i-th value in the values vector
	virtual double getValue(int i) const = 0;

	//getting the i-th value in the columns vector
	virtual int getCol(int i) const = 0;

	//getting the i-th value in rows/row_idx vector
	virtual int getRow(int i) const = 0;

	//sum of the i-th row
	virtual double sumOfRow(int i) const = 0;

};


template <typename T>
class SparseMatrixCOO : public SparseMatrix {
public:


	//define the Costructor
	SparseMatrixCOO(std::vector<T>* values, std::vector<int>* rows, std::vector<int>* cols)
		: values(values), rows(rows), cols(cols) {
		numRows = getNumberRows();
		numCols = getNumberCols();
	}

	//getting the numbers of Rows
	int getNumberRows() const override {
		int RowSize = rows->size();
		if (RowSize == 0) {
			return 0;
		}

		int maxRow = 0;
		for (int i = 0; i < RowSize; ++i) {
			if ((*rows)[i] > maxRow) {
				maxRow = (*rows)[i];
			}
		}
		return maxRow + 1;
	}

	//getting the number of Columns
	int getNumberCols() const override {

		int ColSize = cols->size();

		if (ColSize == 0) {
			return 0;
		}

		int maxCol = 0;
		for (int i = 0; i < ColSize; ++i) {
			if ((*cols)[i] > maxCol) {
				maxCol = (*cols)[i];
			}
		}
		return maxCol + 1;
	}

	//getting all non zeros values
	int getNonZeros() const override {
		return values->size();
	}



	//read the entry of the matrix

	double operator() (int Row, int Col) const override {

		if (Row < 0 || Row >= numRows || Col < 0 || Col >= numCols) {
			throw std::out_of_range("Index out of range"); //I need to create an heandler (try...catch)
			//std::cerr << "Index out of range" << std::endl;
			//return 0.0;
		}

		int ValueSize = values->size();
		for (int i = 0; i < ValueSize; ++i) {
			if ((*rows)[i] == Row && (*cols)[i] == Col) {
				return (*values)[i];
			}
		}

		//if there's something wrong with the insertion of the element of the matrix
		//return an error to have a deterministic outcome

		//throw std::out_of_range("Element not found"); again need a try...catch heandler
		//std::cerr << "Value not found" << std::endl;
		return 0.0;
	}



	//Write an entry of the mmatrix  

	void operator()(double value, int row, int col) override {

		//checking if out of bounds
		if (row < 0 || row >= numRows || col >= numCols || col < 0) {
			//throw std::out_of_range("Index out of range"); again need a try...catch
			std::cerr << "Index out of range" << std::endl;
			return; //if Index out of range exit from the function
		}

		//CASE 1: there's already a value at A(row, col)

		//for cicle to see if there's a value, if so replace it
		int ValueSize = values->size();
		for (int i = 0; i < ValueSize; ++i) {
			if ((*rows)[i] == row && (*cols)[i] == col) {
				(*values)[i] = value;
			}
		}

		//CASE 2: there's no value at A(row, col)
		int InsertPosition = 0;
		while ((*rows)[InsertPosition] < row && InsertPosition < ValueSize) {
			InsertPosition++;
		}
		values->insert(values->begin() + InsertPosition, value);
		cols->insert(cols->begin() + InsertPosition, col);
		rows->insert(rows->begin() + InsertPosition, row);


		//alternativly we could send an error instead of inserting a new value like:
		// std::cerr << "Entry (" << row << "," << col << ") has not been inserted)" << std::endl;
	}



	//Matrix-vector computation

	std::vector<double> operator*(const std::vector<double>& x) const override {

		//if the length of the vector makes it impossible to execute the product, return an error
		int x_size = x.size();
		if (x_size != numCols) {
			// throw std::out_of_range("Impossible, size of vector differnt then number of columns in the matrix");
			std::cerr << "Impossible, size of vector differnt then number of columns in the matrix" << std::endl;
			return std::vector<double>();
		}

		//initializing the result vector
		std::vector<double> result(numRows, 0.0);

		//the result vector will have only the non zeros values multiplications
		//we can consider only the values vector and multiply those values for the x-coordinate
		//corresponding to the column (rows*colums)

		//for cicle to insert values in the result vector
		int result_size = result.size();
		for (int j = 0; j < result_size; ++j) {

			//for cicle to consider only one row at a time 
			int ValueSize=values->size();
			for (int i = 0; i < ValueSize; ++i) {
				if ((*rows)[i] == j) {
					result[j] += (*values)[i] * x[(*cols)[i]];
				}
			}
		}
		return result;
	}



	//to print I decided to allocate a null matrix, then inserting the non Zeros values, so print the complete matrix

	void print() const override {
		std::vector<std::vector<double>> matrix(numRows, std::vector<double>(numCols, 0.0)); //allocate a null matrix

		//Filling with the non Zeros values
		int ValueSize = values->size();
		for (int i = 0; i < ValueSize; ++i) {
			matrix[(*rows)[i]][(*cols)[i]] = (*values)[i];
		}

		//Printing the matrix
		for (int i = 0; i < numRows; ++i) {
			for (int j = 0; j < numCols; ++j) {
				std::cout << matrix[i][j] << " ";
			}
			std::cout << std::endl;
		}
	}


	//getting the i-th value in the values vector
	double getValue(int i) const override {
		return (*values)[i];
	}

	//getting the i-th value in the vector columns
	int getCol(int i) const override {
		return (*cols)[i];
	}

	//getting the i-th value in the vector rows
	int getRow(int i) const override {
		return (*rows)[i];
	}


	//sum of the i-th row

	double sumOfRow(int i) const override {
		double sum = 0.0;
		int RowSize = rows->size();
		for (int j = 0; j <RowSize; ++j) {
			if ((*rows)[j] == i) {
				sum += (*values)[j];
			}
		}
		return sum;
	}


private:
	int numRows;
	int numCols;
	std::vector<double>* values; //non Zeros values
	std::vector<int>* rows; //rows indexes of the non Zeros values
	std::vector<int>* cols; //cols indexes of the non Zeros values

};









template <typename T>
class SparseMatrixCSR : public SparseMatrix {

public:

	//defining the constructor

	SparseMatrixCSR(std::vector<T>* values, std::vector<int>* cols, std::vector<int>* row_idx)
		: values(values), cols(cols), row_idx(row_idx) {
		numRows = getNumberRows();
		numCols = getNumberCols();
	}



	//getting the number of rows

	int getNumberRows() const override {
		int Row_idx_size = row_idx->size();
		if (Row_idx_size<= 1) {
			return 0;
		}

		return  Row_idx_size - 1;
	}



	//getting the number of columns

	int getNumberCols() const override {
		int ColSize = cols->size();
		if (ColSize == 0) {
			return 0;
		}
		int maxCol = 0;
		for (int i = 0; i < ColSize; ++i) {
			if ((*cols)[i] > maxCol) {
				maxCol = (*cols)[i];
			}
		}
		return maxCol + 1;

	}


	//getting the number of Non Zeors elements

	int getNonZeros() const override {
		int ValueSize = values->size();
		return ValueSize;
	}

	//reading values in a matrix  

	double operator()(int Row, int Col) const override {

		int p = Row + 1;

		//if indexes out of bounds return an error
		if (Row < 0 || Row >= numRows || Col < 0 || Col >= numCols) {
			throw std::out_of_range("Index out of range");//I need to create an heandler (try...catch)
			//std::cerr << "Index out of range" << std::endl;
			//return 0.0;
		}

		//if the row passed in input as all zeros, so row_idx[i+1]-row_idx[i]==0, return zero without any calculation
		if ((*row_idx)[p] - (*row_idx)[Row] == 0) {
			return 0.0;
		}

		//if all constraints are satisfied => there's a non zero value in the row, at the position of the column, so read the value
		for (int i = (*row_idx)[Row]; i < (*row_idx)[p]; ++i) {
			if ((*cols)[i] == Col) {
				return (*values)[i];
			}
		}



		//if something's wrong throw an error
		// throw std::out_of_range("Index out of range"); better but I need to create an heandler (try...catch)
		//std::cerr << "Value not found" << std::endl;
		return 0.0;
	}



	//writnig an entry in the matrix 

	void operator()(double value, int row, int col) override {
		if (row < 0 || row >= numRows || col >= numCols || col < 0) { //checking if out of bounds
			// throw std::out_of_range("Index out of range"); better but need an heandler try...catch
			std::cerr << "Index out of range" << std::endl;
		}

		int p = row + 1;

		//CASE: all zeros row:

		//if the row has all values zeros allocate the value in input at column col and adjust row_idx
		if ((*row_idx)[p] - (*row_idx)[row] == 0) {

			//the position to insert the value is row_idx[row]
			int insertPosition = (*row_idx)[row];

			values->insert(values->begin() + insertPosition, value);
			cols->insert(cols->begin() + insertPosition, col);

			//every value in row_idx, where index is bigger then the index of the inserted value,  must be incremented by one
			for (int i = row + 1; i <= numRows; ++i) {
				(*row_idx)[i]++;
			}

		}

		//CASE: some values in the row:

		//if there's a value at the corresponding position, update the values vector
		for (int i = (*row_idx)[row]; i < (*row_idx)[p]; ++i) {
			if ((*cols)[i] == col) {
				(*values)[i] = value;
				break; //if not the programm won't exit the cicle
			}

			//if there's no value at the correpsonding position => allocate a new value and update the matrix-definition vectors
			else {
				//we need to find the correct place to insert the value
				int start = (*row_idx)[row];
				int end = (*row_idx)[p];

				int insertPosition = start;

				while (insertPosition < end && (*cols)[insertPosition] < col) {
					insertPosition++;
				}

				//when found the correct position, adjust: values, cols vectors

				values->insert(values->begin() + insertPosition, value);
				cols->insert(cols->begin() + insertPosition, col);

				//every value in row_idx, where index is bigger or equal then the index of the inserted value,  must be incremented by one
				for (int i = p; i <= numRows; ++i) {
					(*row_idx)[i]++;
				}
			}
		}
	};



	//Matrix-vector computation 
	std::vector<double> operator*(const std::vector<double>& x) const override {

		//if the length of the vector makes it impossible to execute the product, return an error
		int x_size = x.size();
		if (x_size != numCols) {
			std::cerr << "Impossible, size of vector differnt then number of columns in the matrix" << std::endl;
		}

		//initiate the values vector
		std::vector<double> result(numRows, 0.0);

		int l = result.size();

		//the result vector will have only the non zeros values multiplications
		//we can consider only the values vector and multiply those values for the x-coordinate
		//corresponding to the column (rows*colums)
		for (int i = 0; i < l; ++i) {
			int p = i + 1;
			for (int k = (*row_idx)[i]; k < (*row_idx)[p]; ++k) {
				result[i] += (*values)[k] * x[(*cols)[k]];
			}
		}
		return result;
	}



	//to print I decided to allocate a null matrix, then inserting the non Zeros values, so print the complete matrix

	void print() const override {

		//allocate a null matrix
		std::vector<std::vector<double>> matrix(numRows, std::vector<double>(numCols, 0.0));

		//Filling with the non Zeros values
		int Row_idx_size = row_idx->size();
		for (int j = 0; j < Row_idx_size - 1; j++) {

			for (int i = (*row_idx)[j]; i < (*row_idx)[j + 1]; ++i) {
				matrix[j][(*cols)[i]] = (*values)[i];

			}
		}

		//Printing the matrix
		for (int i = 0; i < numRows; ++i) {
			for (int j = 0; j < numCols; ++j) {
				std::cout << matrix[i][j] << " ";
			}
			std::cout << std::endl;
		}
	}


	//getting the i-th value in the values vector
	double getValue(int i) const override {
		return (*values)[i];
	}

	//getting the i-th value in the vector columns
	int getCol(int i) const override {
		return (*cols)[i];
	}

	//getting the i-th value in the vector row:idx
	int getRow(int i) const override {
		return (*row_idx)[i];
	}


	//sum of the i-th row

	double sumOfRow(int i) const override {
		double sum = 0.0;
		int q = i + 1;
		for (int j = (*row_idx)[i]; j < (*row_idx)[q]; ++j) {
				sum += (*values)[j];
			}
		return sum;
	}



private:

	int numRows;
	int numCols;
	std::vector<double>* values; //non Zeros values
	std::vector<int>* cols; //cols indexes of the non Zeros values
	std::vector<int>* row_idx; //Cumulative number of non Zeros for rows

};



//function to convert a matrix from COO format in CSR format

SparseMatrixCSR<double>* COO_to_CSR(const SparseMatrixCOO<double>& obj) {

	//getting the dimentions for the following vectors that will go to create the instance of CSR
	int RowSize = obj.getNumberRows();
	int ValSize = obj.getNonZeros();

	//dinamically allocating the vector for the CSR Instance, if not allocated they will be lost when out of scope
	//and the instance would be damaged or will work with non-vectors

	//must be deleted when useless in the main function

	std::vector<int>* Index = new std::vector<int>(RowSize + 1, 0);
	std::vector<int>* columns = new std::vector<int>(ValSize, 0);
	std::vector<double>* entries = new std::vector<double>(ValSize, 0.0);

	//inserting values in the row_idx (in this case Index) 
	//just counting how many times a certain i-value compare in the row vector of COO instance

	for (int i = 0; i < RowSize; ++i) {
		int count = 0;
		for (int j = 0; j < ValSize; ++j) {
			if (i == obj.getRow(j)) {
				count++;
			}
		}
		int q = i + 1;
		(*Index)[q] = (*Index)[i] + count;
	}

	//columns and values vectors are the same in both formats

	for (int i = 0; i < ValSize; ++i) {
		(*columns)[i] = obj.getCol(i);
		(*entries)[i] = obj.getValue(i);
	}

	//creating a pointer to an object SparseMatrixCSR that must be deleted when useless in the main function
	SparseMatrixCSR<double>* CSR_Matrix = new SparseMatrixCSR<double>(entries, columns, Index);

	return CSR_Matrix;

}


//function to convert a matrix from CSR format in COO format
SparseMatrixCOO<double>* CSR_to_COO(const SparseMatrixCSR<double>& obj) {

	//getting the dimentions of the vector that will be used in the computation
	int RowSize = obj.getNumberRows();
	int ValSize = obj.getNonZeros();

	//dinamically allocating the vector for the CSR Instance, if not allocated they will be lost when out of scope
	//and the instance would be damaged or will work with non-vectors

	//must be deleted when useless in the main function

	std::vector<int>* Rows = new std::vector<int>(ValSize, 0);
	std::vector<int>* Columns = new std::vector<int>(ValSize, 0);
	std::vector<double>* entries = new std::vector<double>(ValSize, 0.0);


	//columns and values vectors are the same in both formats

	for (int i = 0; i < ValSize; ++i) {
		(*Columns)[i] = obj.getCol(i);
		(*entries)[i] = obj.getValue(i);
	}

	//inserting values in the rows vector of COO
	//just conting how many non Zeros values there are in every row, so inserting a number of row indexes equal to this number

	for (int i = 0; i < RowSize; ++i) {
		int start = obj.getRow(i);
		int end = obj.getRow(i + 1);
		for (int j = start; j < end; ++j) {
			(*Rows)[j] = i;
		}
	}


	//creating a pointer to an object SparseMatrixCOO that must be deleted when useless in the main function

	SparseMatrixCOO<double>* COO_Matrix = new SparseMatrixCOO<double>(entries, Rows, Columns);

	return COO_Matrix;

}


#endif






