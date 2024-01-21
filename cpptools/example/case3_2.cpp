#include <iostream>
using namespace std;
#include "case3_2.h"
#include "./../vec/vec.h"
#include <vector>

void runcase3_2()
{   
    const int rows = 1000;
    const int cols = 1000;
    // Initialize a large 2D array
    std::vector<std::vector<int>> large2DArray(rows, std::vector<int>(cols, 0));
    // Initialize a 1D array with the same total size
    std::vector<int> oneDArray(rows * cols);
    cout << "[---Convert 2D array to 1D array efficiently---]" << endl;
    vec().convert2DTo1D(large2DArray, oneDArray);
    //
    // Access elements in the 1D array
    cout << "First element in 1D array: " << oneDArray[0] << endl;
}