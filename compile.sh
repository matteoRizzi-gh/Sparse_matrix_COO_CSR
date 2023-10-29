#!/bin/bash

# Specifica il comando del compilatore e le opzioni
COMPILER="g++"
OPTIONS="-std=c++17 -Wall -Wpedantic"

# Nome del file sorgente e nome dell'eseguibile
SOURCE_FILE="main.cpp"
OUTPUT_EXECUTABLE="sparse_matrix"

# Compila il file sorgente in un eseguibile
$COMPILER $OPTIONS $SOURCE_FILE -o $OUTPUT_EXECUTABLE

# Verifica se la compilazione è riuscita
if [ $? -eq 0 ]; then
  echo "Compilation Success. Execute ./$OUTPUT_EXECUTABLE to look up the programm."
else
  echo "Errors during compilation."
fi
