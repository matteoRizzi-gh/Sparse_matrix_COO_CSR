#include "pch.h"
#include "gtest/gtest.h"
#include "Sparse_Matrix.h"
#include <random>
#include <vector>
#include <iostream>
#include <cmath>

TEST(SparseMatrixTest, SumOfRowsTest) {

    // Instance of SparseMatrixCSR
    std::vector<double> values = { 3.1, 4, 5, 7.4, 2, 6 };
    std::vector<int> row_idx = { 0, 2, 4, 4, 6 };
    std::vector<int> cols = { 2, 4, 2, 4, 1, 3 };

    SparseMatrixCSR<double> obj(&values, &cols, &row_idx);

    int NumCol = obj.getNumberCols();
    int NumRow = obj.getNumberRows();

    // Unity vector unity=<1,...,1>=1_n
    std::vector<double> unity(NumCol, 1);

    //matrix-unity computation
    std::vector<double> verify = obj.operator*(unity);

    //setting a tolerance of 0.001 because of approximations
    double tolerance = 1e-3;

    // Checking if all the i-th values of the matrix-unity computation are equal to the sum of the rows
    for (int i = 0; i < NumRow; ++i) {
        double sum = obj.sumOfRow(i);
        EXPECT_NEAR(verify[i], sum, tolerance); // Use EXPECT_NEAR to check with tolerance
    }

}

TEST(SparseMatrixTest, CanonicalBasisComputation) {

    // Instance of SparseMatrixCSR
    std::vector<double> values = { 3.1, 4, 5, 7.4, 2, 6 };
    std::vector<int> row_idx = { 0, 2, 4, 4, 6 };
    std::vector<int> cols = { 2, 4, 2, 4, 1, 3 };

    SparseMatrixCSR<double> obj(&values, &cols, &row_idx);

    //dimensions needed for the computation
    int NumCol = obj.getNumberCols();
    int NumRow = obj.getNumberRows();

    // Setting a tolerance of 0.001 because of approximations
    double tolerance = 1e-3;

    bool check = true;
    bool flag = false;

    // Checking if given the i-th vector of the canonical basis, the product A*e_i gives out the i-th column of A
    for (int i = 0; i < NumCol; ++i) {
        std::vector<double> canon(NumCol, 0);
        canon[i] = 1;

        // Computing A*e_i
        std::vector<double> comp = obj.operator*(canon);

        // Checking if every coordinate of the A*e_i is equal to the corresponding value of the i-th column of A
        for (int t = 0; t < NumRow; ++t) {
            double k = obj.operator()(t, i);
            EXPECT_NEAR(comp[t], k, tolerance);
            if (std::abs(comp[t] - k) > tolerance) {
                flag = true;
                break; // If there's an element that differs more than 10^-3, stop the computation
            }
        }
    }
    // If the computation stopped, hence if the test failed, set the flag in false
    if (flag) {
        check = false;
    }

    EXPECT_TRUE(check);
}


TEST(SparseMatrixTest, InsertReadTest) {
    // Instance of SparseMatrixCSR
    std::vector<double> values = { 3.1, 4, 5, 7.4, 2, 6 };
    std::vector<int> row_idx = { 0, 2, 4, 4, 6 };
    std::vector<int> cols = { 2, 4, 2, 4, 1, 3 };

    SparseMatrixCSR<double> obj(&values, &cols, &row_idx);


    //dimensions needed for the computation
    int NumCol = obj.getNumberCols();
    int NumRow = obj.getNumberRows();
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

    //inserting pi
    obj.operator() (pi, i, j);

    //reading the value in A(i,j)
    double verify = obj.operator()(i, j);

    //checking if the insertion works 
    EXPECT_NEAR(obj.operator()(i, j), pi, tolerance);

}

TEST(SparseMatrixTest, ConvertionCSRTest) {

    // Instance of SparseMatrixCSR
    std::vector<double> values = { 3.1, 4, 5, 7.4, 2, 6 };
    std::vector<int> row_idx = { 0, 2, 4, 4, 6 };
    std::vector<int> cols = { 2, 4, 2, 4, 1, 3 };

    SparseMatrixCSR<double> obj(&values, &cols, &row_idx);

    SparseMatrixCOO<double>* MatCOO = CSR_to_COO(obj);

    int NumCol = obj.getNumberCols();
    int NumRow = obj.getNumberRows();

    double tolerance = 1e-3;

    for (int i = 0; i < NumRow; ++i) {
        for (int j = 0; j < NumCol; ++j) {
            double el_1 = obj.operator()(i, j);
            double el_2 = MatCOO->operator()(i, j);
            EXPECT_NEAR(el_1, el_2, tolerance);
        }
    }
}


TEST(SparseMatrixTest, ConvertionCOOTest) {

    // Instance of SparseMatrixCOO
    std::vector<double> values = { 3.1, 4, 5, 7.4, 2, 6 };
    std::vector<int> rows = { 0, 0, 1, 1, 3, 3 };
    std::vector<int> cols = { 2, 4, 2, 4, 1, 3 };

    SparseMatrixCOO<double> obj(&values, &rows, &cols);

    SparseMatrixCSR<double>* MatCSR = COO_to_CSR(obj);

    int NumCol = obj.getNumberCols();
    int NumRow = obj.getNumberRows();

    double tolerance = 1e-3;

    for (int i = 0; i < NumRow; ++i) {
        for (int j = 0; j < NumCol; ++j) {
            double el_1 = obj.operator()(i, j);
            double el_2 = MatCSR->operator()(i, j);
            EXPECT_NEAR(el_1, el_2, tolerance);
        }
    }
}



TEST(SparseMatrixTest, ErrorHandlingTest) {

    // Instance of SparseMatrixCSR
    std::vector<double> values = { 3.1, 4, 5, 7.4, 2, 6 };
    std::vector<int> row_idx = { 0, 2, 4, 4, 6 };
    std::vector<int> cols = { 2, 4, 2, 4, 1, 3 };

    SparseMatrixCSR<double> obj(&values, &cols, &row_idx);

    // Dimensions needed for the computation
    int NumCol = obj.getNumberCols();
    int NumRow = obj.getNumberRows();

    // Attempt to access an out-of-bounds element
    int i = NumCol + 2;
    int j = NumRow + 2;

    // Use EXPECT_THROW to check if an exception is thrown
    EXPECT_THROW(obj.operator()(j, i), std::out_of_range);
}



int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}