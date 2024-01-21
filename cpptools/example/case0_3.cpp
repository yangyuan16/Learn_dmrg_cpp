#include <iostream>
using namespace std;
#include "case0_3.h"
#include "./../array/anyarray1d.h"
#include<vector>

void runcase0_3()
{
    const int size = 10;
    double* intputArray = new double[size]{1.1,2.1,3.1,4.1,5.1,6.1,7.1,8.1,9.1,10.1};

    std::vector<double> outputVector = anyarray1d().arr1dtovec1d(intputArray, size);

    cout << "[---Convert 1d array to 1d vector, Print the element:]" <<endl;
    for (double element : outputVector){
        cout << element << " ";
    }
    cout << endl;

    // Deallocate memory for the array
    delete[] intputArray;
}