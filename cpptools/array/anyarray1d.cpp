#include<iostream>
using namespace std;
#include "anyarray1d.h"
#include <vector>
#include <cstring>

anyarray1d::anyarray1d(/* args */)
{
}

anyarray1d::~anyarray1d()
{
}

// -------------- print 1d array --------------------
template <typename T>
void anyarray1d::printArray(T* arr, int size)
{
    for (int i = 0; i < size; i++)
    {
        cout << *(arr + i) << " ";
    }
    cout << endl;
}
template void anyarray1d::printArray<double>(double*, int);
template void anyarray1d::printArray<int>(int*, int);
//
// ------------ convert 1d array to 2d array ------------
template <typename T>
void anyarray1d::arr1dtoarr2d(const T* inputArray, T**& outputArray, int numRows, int numCols)
{
    // Allocate memory for the 2D array
    outputArray = new T*[numRows];
    for (int i = 0; i< numRows; i++){
        outputArray[i] = new T[numCols];
    }

    // Copy elements from the 1D array to the 2D array
    int index = 0;
    for (int i = 0; i<numRows; i++){
        for (int j = 0; j<numCols; j++){
            outputArray[i][j] = inputArray[index++];
        }
    }
}
template void anyarray1d::arr1dtoarr2d<double>(const double*, double**&, int, int);
template void anyarray1d::arr1dtoarr2d<int>(const int*, int**&, int, int);
//
template <typename T>
void anyarray1d::delete2Darray(T**& array, int numRows){
    // Deallocate memory for the 2D array
    for (int i = 0; i<numRows; i++){
        delete[] array[i];
    }
    delete[] array;
    array = nullptr;
}
template void anyarray1d::delete2Darray<double>(double**&, int);
template void anyarray1d::delete2Darray<int>(int**&, int);
//
//---- convert 1d array to 1d vec
template <typename T>
std::vector<T> anyarray1d::arr1dtovec1d(const T* inputArray, int size)
{
    return std::vector<T>(inputArray, inputArray + size);
}
template std::vector<double> anyarray1d::arr1dtovec1d<double>(const double*, int);
template std::vector<int> anyarray1d::arr1dtovec1d<int>(const int*, int);
//
//---- convert 1d array to 1d vec
template <typename T>
std::vector<std::vector<T>> anyarray1d::arr1dtovec2d(const T* inputArray, int numRows, int numCols)
{
    std::vector<std::vector<T>> outputVector(numRows, std::vector<T>(numCols));
    int index = 0;
    for (int i = 0; i < numRows; i++){
        for (int j = 0; j< numCols; j++){
            outputVector[i][j] = inputArray[index++];
        }
    }
    return outputVector;
}
template std::vector<std::vector<double>> anyarray1d::arr1dtovec2d<double>(
    const double* inputArray, int numRows, int numCols);
template std::vector<std::vector<int>> anyarray1d::arr1dtovec2d<int>(
    const int* inputArray, int numRows, int numCols);