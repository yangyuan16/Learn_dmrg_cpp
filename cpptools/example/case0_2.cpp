#include <iostream>
using namespace std;
#include "case0_2.h"
#include "./../array/anyarray1d.h"

void runcase0_2()
{
    const int numRows = 3;
    const int numCols = 4;

    double* inputArray = new double[numRows * numCols]{1.1,2.1,3.1,4.1,5.2,6.2,7.2,8.2,9.3,10.3,11.3,12.3};

    double** outputArray = nullptr;
    anyarray1d().arr1dtoarr2d(inputArray, outputArray, numRows, numCols);

    // Print the 2D array
    cout << "----[convert 1d array to 2d array, print 2d array] ----" << endl;
    for (int i = 0; i < numRows; i++){
        for (int j = 0; j<numCols; j++){
            cout << outputArray[i][j] << " ";
        }
        cout << endl;
    }

    // Deallocate memory for the 2D array
    anyarray1d().delete2Darray(outputArray, numRows);

    // Deallocate memory for the 1D array
    delete[] inputArray;
}