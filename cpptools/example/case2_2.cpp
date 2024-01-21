#include <iostream>
using namespace std;
#include "case2_2.h"
#include "./../array/anyarray2d.h"
#include "./../array/array2d.h"
#include <vector>

void runcase2_2()
{
    double array2D[5][4]={{1.1,2.1,3.1,4.1},
                         {5.1,6.1,7.1,8.1},
                         {9.1,10.1,11.1,12.1},
                         {13.1,14.1,15.1,16.1},
                         {17.1,18.1,19.1,20.1}};
    // Get the dimension of the 2D array
    int rows = sizeof(array2D) / sizeof(array2D[0]);
    int cols = sizeof(array2D[0]) / sizeof(array2D[0][0]);

    // Convert 2D array to 1D vector using array pointers
    double *p = array2D[0];
    vector<double> vector1D = anyarray2d().arr2dtovec1d(p, rows, cols);
    cout << "----- Print the elements of the 1D vector ----" << endl;
    std::cout <<"1D vector: ";
    for (double element : vector1D){
        cout << element << " ";
    } 
    cout << endl;
}