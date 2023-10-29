#!/bin/bash

# Compiler and options
COMPILER="g++"
CPPFLAGS="-std=c++17"  # Modifica secondo le tue esigenze
LDFLAGS="-lgtest -lgtest_main -pthread"  # Flags di collegamento per Google Test

# Source file name
SOURCE_FILE="test.cpp"

# Output executable name
OUTPUT_EXECUTABLE="my_test"

# Compilation
$COMPILER $CPPFLAGS -o $OUTPUT_EXECUTABLE $SOURCE_FILE $LDFLAGS

# Check if compilation was successful
if [ $? -eq 0 ]; then
  echo "Compilation of google test completed. Execute ./$OUTPUT_EXECUTABLE to run the tests."
else
  echo "Errors during compilation of the test file."
fi
