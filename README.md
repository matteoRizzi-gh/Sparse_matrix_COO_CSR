# Sparse_matrix_COO_CSR
Basic implementation of Sparse Matrix COO and CSR 

The codes file in this folder contains headers, google test and main file to handle Sparse Matrix in COO and CSR formats.
 
The `main.cpp` file includes tests to evaluate the correctness of the functions defined in the header files. Specifically, there is a simple instance of a CSR matrix that is converted into a COO matrix, printed, and all tests are run to verify the correctness of the implemented functions.

To certify the correctness of the code, Google Test is used. You can run the tests with the following commands in shell:

	g++ -o my_test test.cpp -lgtest -lgtest_main -pthread

or equivalently with the compilation_test script given, with:

	chmod -x compilation_test.sh
	./compilation_test.sh

The code has been successfully compiled and tested with the following commands:

	g++ -std=c++17 -Wall -Wextra -Wpedantic main.cpp -o sparse_matrix

	clang++ -std=c++11 -Wall -Wextra -pedantic main.cpp -o sparse_matrix

so evalueted looking at the file cretated this way.

Equivalently with the compile script given, with:

	chmod -x compile.sh
	./compile.sh

Rizzi Matteo
