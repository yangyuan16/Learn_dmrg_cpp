#include <iostream>
using namespace std;
#include "case2_3.h"
#include "./../array/anyarray2d.h"
#include "./../array/array2d.h"

void runcase2_3()
{
    // Assuming you have a 2D array, for example:
    int rows = 3;
    int cols = 4;
    int** array2D = new int*[rows];
    for (int i = 0; i < rows; i++)
    {
        array2D[i] = new int[cols]{i * cols + 1, i * cols + 2, i * cols + 3, i * cols + 4};
    }

    // Allocate memory for the 1D array
    int* array1D = new int[rows * cols];

    // Covert 2D array to 1D array using pointers
    anyarray2d().arr2dtoarr1d(array2D, rows, cols, array1D);
    cout << "-------------------------------------------" << endl;
    cout << "2d array convert to 1d array,print the elements of the 1D array " << endl;
    // Print the elements of the 1D array
    cout << "1D Array: ";
    for (int i = 0; i< rows * cols; i++){
        cout << array1D[i] << " ";
    }
    cout << endl;
    // Free allocate memory
    for (int i = 0; i < rows; i++){
        delete[] array2D[i];
    }
    delete[] array2D;
    delete[] array1D;
}

