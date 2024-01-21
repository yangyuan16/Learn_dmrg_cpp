#include <iostream>
using namespace std;
#include "case0_4.h"
#include "./../array/anyarray1d.h"
#include<vector>

void runcase0_4()
{
    const int numRows = 4;
    const int numCols = 4;

    double* inputArray = new double[numRows * numCols]{ 1.1,2.1,3.1,4.1,
                                                        5.1,6.1,7.1,8.1,
                                                        9.1,10.1,11.1,12.1,
                                                        13.1,14.1,15.1,16.1};
    std::vector<std::vector<double>> outputVector = anyarray1d().arr1dtovec2d(inputArray, numRows, numCols);
        
    cout << "---[Convert 1d array to 2d vector, Print the elements]---"<< endl;
    for (const auto& row : outputVector){
        for (double element : row){
            cout << element << " ";
        }
        cout << endl;
    }

    // Deallocate memory for the array

    delete[] inputArray;

}