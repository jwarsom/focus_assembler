/*
 * SortEdges.cpp
 *
 *      Author: Julia Warnke
 */

#include <iostream>
#include <stdint.h>
#include <cstdio>
#include <algorithm>
#include <dirent.h>
#include<map>
using namespace std;
#include "MappingValues.h"
#include "SortEdges.h"



//int sortByNodeLabels(const void * a, const void * b)
//Description: This function compares the labels of the
//nodes belonging to the edge
//Output:None
//Input: two void *, which are cast to int *
//Return: an int, -1 if a < b, 0 if a == b, 1 if a > b
int sortByNodeLabels(const void * a, const void * b) {

	if(((uint32_t *)a)[0] < ((uint32_t *) b)[0])
	{
		return -1;
	}

	if( ((uint32_t *)a)[0] > ((uint32_t *)b)[0])
	{
		return 1;
	}

	if(((uint32_t *)a)[0] ==  ((uint32_t *)b)[0])
	{

		if(((uint32_t *)a)[1] < ((uint32_t *) b)[1])
		{
			return -1;
		}

		if( ((uint32_t *)a)[1] > ((uint32_t *)b)[1])
		{
			return 1;
		}
	}

	return 0;
}


//int sortByEdgeWeights(const void * a, const void * b)
//Description: This function compares the weight of
//two edges
//Output:None
//Input: two void *, which are cast to int *
//Return: an int, -1 if a < b, 0 if a == b, 1 if a > b
int sortByEdgeWeights(const void * a, const void * b) {

	if(((uint32_t *)a)[0] < ((uint32_t *) b)[0])
	{
		return -1;
	}

	if( ((uint32_t *)a)[0] > ((uint32_t *)b)[0])
	{
		return 1;
	}

	if(((uint32_t *)a)[0] == ((uint32_t *)b)[0])
	{
		int ovlLenA = ((uint32_t *)a)[2] >> 22;
		int ovlLenB = ((uint32_t *)b)[2] >> 22;
		if(ovlLenA > ovlLenB)
		{
			return -1;
		}

		if(ovlLenA < ovlLenB)
		{
			return 1;
		}
	}
	return 0;
}

//int sortByEdgeWeightsHeavy(const void * a, const void * b)
//Description: This function compares the weight of
//two edges
//Output:None
//Input: two void *, which are cast to int *
//Return: an int, -1 if a < b, 0 if a == b, 1 if a > b
int sortByEdgeWeightsHeavy(const void * a, const void * b) {

	if(((uint32_t *)a)[0] < ((uint32_t *) b)[0])
	{
		return -1;
	}

	if( ((uint32_t *)a)[0] > ((uint32_t *)b)[0])
	{
		return 1;
	}

	if(((uint32_t *)a)[0] == ((uint32_t *)b)[0])
	{
		int ovlLenA = ((uint32_t *)a)[2] >> 16;
		int ovlLenB = ((uint32_t *)b)[2] >> 16;
		if(ovlLenA > ovlLenB)
		{
			return -1;
		}

		if(ovlLenA < ovlLenB)
		{
			return 1;
		}
	}
}


//SortEdges(void)
//Description:The constructor
//Input:None
//Output:None
//Return:None
SortEdges::SortEdges(const int st)
	:sortType(st){;;;}

	//void sort(const char [])
	//Description: This function sorts an edge input file and
	//outputs it into
	//Input:workDir, const char [], the path to the directory
	//Output:None
	//Return:None
	void SortEdges::sort_edges(const char inFile[], const char outFile[]) {

		if(fileExists(inFile))
		{
			FILE * iFile = fopen(inFile, "r");

			fseek(iFile, 0, SEEK_END);
			long long int numElements = ftell(iFile)/sizeof(uint32_t);
			fseek(iFile, 0, SEEK_SET);

			uint32_t * buff = new uint32_t[numElements];
			fread(buff, sizeof(uint32_t), numElements, iFile);
			fclose(iFile);

			int numEdges = numElements/3;

			if(sortType == 0) //Sort by Node Labels
			{
				qsort(buff, numEdges, sizeof(uint32_t) * 3, sortByNodeLabels);
			}

			if(sortType == 1) //Sort by Edge Weight
			{
				qsort(buff, numEdges, sizeof(uint32_t) * 3, sortByEdgeWeights);
			}

			if(sortType == 2) //Sort by Heavy Edge Weight
			{
				qsort(buff, numEdges, sizeof(uint32_t) * 3, sortByEdgeWeightsHeavy);
			}
			

			FILE * oFile = fopen(outFile, "w");
			fwrite(buff, sizeof(uint32_t), numElements, oFile);
			fclose(oFile);

			delete [] buff;
		}
	}


//bool fileExists( char fileName [] )
//Description: This function checks to see if a file exists
//Input: fileName (char []) the name of the file
//Output:None
//Return:None
bool SortEdges::fileExists(const char fileName []) {

	FILE * pFile = fopen(fileName, "r");

	if(pFile != '\0')
	{
		fclose(pFile);
		return true;
	}

	return false;
}
