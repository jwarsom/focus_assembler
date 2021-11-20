/*
 * SortFragments.cpp
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
#include "Dictionary.h"
#include "Fragment_Index.h"
#include "Labels.h"
#include "SortFragments.h"


//SortFragments(void)
//Description:The constructor
//Input:None
//Output:None
//Return:None
SortFragments::SortFragments(){}

//void sort(const char [], const char [])
//Description: This function reads two input files 
//containing a fragment index and its corresponding label array.
//The fragment and label indexes are sorted by an ordering obtained
//from the overlap graph and then output into ouput files 
//Input:iFragments, const char [], path to fragment index file
//ofragments, const char [], path to fragment ouput file 
//Output:None
//Return:None
void SortFragments::sort(const char iFragments [], const char iLabels [], const char oFragments [], const char oLabels [],  int * & nodeMap) {
	if(fileExists(iFragments))
	{
		FILE * iFileF = fopen(iFragments, "r");
		FILE * iFileL = fopen(iLabels, "r");

		//Read in the input Fragment_Index
		FILE * pFile  = fopen(iFragments, "r");
		unsigned long long int fragmentOffset;
		fread(&fragmentOffset, sizeof(unsigned long long int), 1, pFile);

		unsigned long long int indexOffset;
		fread(&indexOffset, sizeof(unsigned long long int), 1, pFile);

		unsigned long long int intervalLength;
		fread(&intervalLength, sizeof(unsigned long long int), 1, pFile);

		unsigned long long int numIntervalFragments;
		fread(&numIntervalFragments, sizeof(unsigned long long int), 1, pFile);
		fclose(pFile);


		Fragment_Index iIndex;
		iIndex.reserve(numIntervalFragments, intervalLength, indexOffset, fragmentOffset);
		iIndex.read(iFileF);
		fclose(iFileF);		
		
		iIndex.set();

		//Read in the input Label
		int tagLen;
        	pFile = fopen(iLabels, "r");
        	fread(&tagLen, sizeof(int), 1, pFile);
        	fclose(pFile);
	
		Labels iLabel;
		iLabel.reserve(tagLen, numIntervalFragments, indexOffset);
		iLabel.read(iFileL);
		fclose(iFileL);

		Fragment_Index oIndex;
		oIndex.reserve(iIndex.numFragments(), iIndex.length(), 0, 0);

		Labels oLabel;

		//Adding 11 to the tag length because we are going to append the 
		//new rank of the fragment to the tag
		oLabel.reserve(iLabel.getTagLen()+11, iLabel.getNumLabels(), 0);

		//Sort the fragment labels according to their rank in the ordered graph
		unsigned long long int * sortBuff = new unsigned long long int [iIndex.numFragments()];
		for(int i = 0; i < iIndex.numFragments(); i++)
		{
			sortBuff[i] = iIndex.getIndexOffset() + i;
		}

		std::sort(sortBuff, sortBuff+iIndex.numFragments(), fragCmp(nodeMap));
	
		char * labelBuff = new char  [iLabel.getTagLen()+12];
		char * fragBuff = new char [10000];
		char * qualBuff = new char [10000];

		for(int i = 0; i < iIndex.numFragments(); i++)
		{
			//Extract and append the DNA seq and Qual values
			pair<unsigned long long int, unsigned long long int> myBounds;
			myBounds = iIndex.indexBounds(sortBuff[i] - iIndex.getIndexOffset());

			for(int j = myBounds.first; j < myBounds.second; j++)
			{
				fragBuff[j - myBounds.first] = iIndex.at(j);
			}

			for(int j = myBounds.first; j < myBounds.second; j++)
			{
				qualBuff[j-myBounds.first] = (iIndex.qualAt(j) + 33);
			}

			fragBuff[myBounds.second-myBounds.first] = '\0';
			qualBuff[myBounds.second-myBounds.first] = '\0';

			oIndex.append(fragBuff, qualBuff);

			//Extract and append the DNA seq label
			sprintf(labelBuff, "%s.%d", iLabel.labelAt(sortBuff[i]-iLabel.getOffset()).c_str(), nodeMap[sortBuff[i]]);
			oLabel.append(labelBuff);
		}
		
		oIndex.set();

		delete [] labelBuff;
		delete [] fragBuff;
		delete [] qualBuff;

		delete [] sortBuff;			

		FILE * oFileF = fopen(oFragments, "w");
		FILE * oFileL = fopen(oLabels, "w");

		oIndex.write(oFileF);
		oLabel.write(oFileL);

		fclose(oFileF);
		fclose(oFileL);
	}
}


//bool fileExists( char fileName [] )
//Description: This function checks to see if a file exists
//Input: fileName (char []) the name of the file
//Output:None
//Return:None
bool SortFragments::fileExists(const char fileName []) {

	FILE * pFile = fopen(fileName, "r");

	if(pFile != '\0')
	{
		fclose(pFile);
		return true;
	}

	return false;
}
