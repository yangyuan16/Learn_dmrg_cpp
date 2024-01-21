#ifndef _ANYARRAY1D_
#define _ANYARRAY1D_

#include<iostream>
using namespace std;
#include <vector>

class anyarray1d
{
private:
    /* data */
public:
    anyarray1d(/* args */);
    ~anyarray1d();
    //
    template <typename T>
    void printArray(T*, int);
    //
    template <typename T>
    void arr1dtoarr2d(const T* inputArray, T**& outputArray, int numRows, int numCols);
    //
    template <typename T>
    void delete2Darray(T**& array, int numRows);
    // 
    // 1d array to 1d vector
    template <typename T>
    std::vector<T> arr1dtovec1d(const T* inputArray, int size);
    //
    // 1d array to 2d vector
    template <typename T>
    std::vector<std::vector<T>> arr1dtovec2d(const T* inputArray, int numRows, int numCols);
    
    

};


#endif