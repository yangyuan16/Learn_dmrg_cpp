#include "random_list_generator.h"
#include <cstdlib>
#include <ctime>

int* RandomListGenerator::generateRandomList(int listSize, int desiredZeros, int seedValue){
	int* randomList = new int[listSize];
	int zerosCount = 0;
	int onesCount = 0;

	std::srand(seedValue);

	for (int i = 0; i < listSize; ++i){
	    if (zerosCount < desiredZeros && (onesCount == (listSize - desiredZeros) || std::rand() % 2 == 0)){
	        randomList[i] = 0;
		++zerosCount;
	    } else{
	       randomList[i] = 1;
	       ++onesCount;
	    }
	}

	return randomList;
}
