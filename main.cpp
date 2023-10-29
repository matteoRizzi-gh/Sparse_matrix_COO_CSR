//To understand better some concept I used the official website of the C++ STandards Commitee https://isocpp.org 

#include <iostream>
#include <vector>
#include <random>
#include "Sparse_Matrix.h"


//Test the code



void test_unity(const SparseMatrix* obj) {

    //flag to see if all test give out a good outcome
    bool check = true;

    //dimensions needed for the computation
    int NumCol = obj->getNumberCols();
    int NumRow = obj->getNumberRows();

    //unity vector unity=<1,...,1>=1_n
    std::vector<double> unity(NumCol, 1);

    //matrix-unity computation
    std::vector<double> verify = obj->operator*(unity);

    //setting a tolerance of 0.001 because of approximations
    double tolerance = 1e-3;

    //checking if all the i-ths values of the matrix-unity computation are equal to the sum of the rows

    for (int i = 0; i < NumRow; ++i) {
        double sum = obj->sumOfRow(i);
        if (std::abs(verify[i] - sum) > tolerance) {
            check = false;
            break; //if there an error the computation stops
        }
    }
    if (check == true) {
        std::cout << "Test unity: succeded" << std::endl;
    }
    else {
        std::cout << "Test unity: failed" << std::endl;
    }
}




//checking if given the i-th vector of the canonical basis the product A*e_i gives out the i-th column of A

void test_canonical(const SparseMatrix* obj){
    
    //There could be problems with the break, so using a flag that will decide the check
    
    bool flag = true;
    
    //dimensions needed for the computation
    int NumCol = obj->getNumberCols();
    int NumRow = obj->getNumberRows();

    //setting a tolerance of 0.001 because of approximations
    double tolerance = 1e-3;
    
 //getting the it-h vector of the canonical basis
   for (int i = 0; i < NumCol; ++i) {
       std::vector<double> canon(NumCol, 0);
        canon[i] = 1;

      //computing A*e_i
       std::vector<double> comp = obj->operator*(canon);

      //checking if every coordinate of the A*e_i is equal of the corresponding value of the i-th column of A
        for (int t = 0; t < NumRow; ++t) {
           double k = obj->operator()(t, i);
           if (std::abs(comp[t] - k)>tolerance) {
                flag=false;
                break; //if there's an element that differes more that 10^-3 stops the computation
            }
      }
        
    }

   if (flag == true) {
       std::cout << "Test canonical: succeded" << std::endl;
   }
   else {
       std::cout << "Test canonical: failed" << std::endl;
   }
}




void test_insert_read(SparseMatrix* obj) {
   
    bool flag = true;

    //dimensions needed for the computation
    int NumCol = obj->getNumberCols();
    int NumRow = obj->getNumberRows();
    double tolerance = 1e-3;

    // Initialize random number generator to get a random element of the matrix
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_int_distribution<int> colDistribution(0, NumCol - 1);
    std::uniform_int_distribution<int> rowDistribution(0, NumRow - 1);

    // Generate random values for i and j
    int j = colDistribution(generator);
    int i = rowDistribution(generator);

    //I'll use pi to insert at the random place
    double pi = 3.1415;

    //saving the actuall value so that the matrix won't be modified
    double safe = obj->operator()(i, j);

    //inserting pi
    obj->operator() (pi, i, j);

    //reading the value in A(i,j)
    double verify = obj->operator()(i, j);
    
    //checking if the insertion works 
    if (std::abs(verify - pi) > tolerance) {
        flag = false;
    }

    //restore matrix
    obj->operator() (safe, i, j);

    if (flag == true) {
        std::cout << "Test insert-read: succeded" << std::endl;
    }
    else {
        std::cout << "Test insert-read: failed" << std::endl;
    }
}



//testing if the convertion function works

void test_convertion(const SparseMatrix* obj_1, const SparseMatrix* obj_2) {

    bool flag = true;

    //dimensions needed for the computation
    int NumCol = obj_1->getNumberCols();
    int NumRow = obj_1->getNumberRows();

    int NumCol_2 = obj_2->getNumberCols();
    int NumRow_2 = obj_2->getNumberRows();

    //checking if the dimentions are correct, else return an error
    if (NumCol != NumCol_2 || NumRow != NumRow_2) {
        flag = false;
        std::cerr << "Wrong Matrix Dimentions" << std::endl;
    }

    double tolerance = 1e-3;

    //checking using the read operator for COO and CSR if the values in A(i,j) and A'(i,j) are the same
    for (int i = 0; i < NumRow; ++i) {
        for (int j = 0; j < NumCol; ++j) {

            double el_1 = obj_1->operator()(i, j);
            double el_2 = obj_2->operator()(i, j);

            if (std::abs(el_1 - el_2) > tolerance) {
                flag = false;
                break; //exits the cycle, we found an error
            }
        }
    }

    if (flag == true) {
        std::cout << "Test convertion: succeded" << std::endl;
    }
    else {
        std::cout << "Test convertion: failed" << std::endl;

    }

}



//testing if the exception in the overloaded read operator would be managed correctly

void test_errors(const SparseMatrix* obj) {

    //dimensions needed for the computation
    int NumCol = obj->getNumberCols();
    int NumRow = obj->getNumberRows();

    //getting out of bounds
    int i = NumCol + 2;
    int j = NumRow + 2;
    
    //try...catch to manage the exceprion
    try {
        obj->operator()(j, i);
    }
    catch (const std::out_of_range& e) {
        std::cerr << "Error: " << e.what();
        std::cout << " => Exception managed " << std::endl;
    }
    
}




int main() {
    // Instance of SparseMatrixCOO
    std::vector<double> values = { 3.1, 4, 5, 7.4, 2, 6 };
    std::vector<int> row_idx = { 0, 2, 4, 4, 6 };
    std::vector<int> cols = { 2, 4, 2, 4, 1, 3 };
    std::vector<int> rows = { 0, 0, 1, 1, 3, 3 };
    SparseMatrixCSR<double> myMatrixCSR(&values, &cols, &row_idx );
    SparseMatrixCSR<double>* MatCSR = &myMatrixCSR;
    SparseMatrixCOO<double>* MatCOO = CSR_to_COO(myMatrixCSR);
    MatCOO->print();
    test_canonical(MatCOO);
    test_insert_read(MatCOO);
    test_unity(MatCOO);
    test_convertion(MatCSR, MatCOO);
    test_errors(MatCSR);


    return 0;
}

