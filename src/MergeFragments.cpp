/*
 * MergeFragments.cpp
 *
 *      Author: Julia
 */

#include<iostream>
#include<stdint.h>
#include<cstdio>
#include<cstring>
#include <queue>
#include<cstdlib>
using namespace std;
#include "Dictionary.h"
#include "Fragment_Index.h"
#include "Labels.h"
#include "MergeFragments.h"

//Size unit
#define megabyte 1024000

//MergeFragments(void)
//Description: The ructor
//Input: None
//Output:None
//Return:None
MergeFragments::MergeFragments() {

}

//Merger(void)
//Description: The destructor
//Input:None
//Output:None
//Return:None
MergeFragments::~MergeFragments(){

}

//void merge( int,  int,  int,  int,  char [],  char [])
//Description: This function merges two sorted files
//Input: startSet (A  int): The starting file number of the first set, startSet2 (a  int): the
//starting number of the second set, numSets1 (A  int): The number of files in  set1,
//numSets2 (A  int): The number of files in set2, inDir ( char []): the name of the input directory
//outDir ( char []): The output directory
//Output:None
//Return:None
void MergeFragments::merge(int startSet,  int numSets,  char inDirF [],  char outDirF [],  char inDirL [] ,  char outDirL [])
{
	char * fragmentBuff = new char [22 * megabyte]; //This stores the sorted labels, fragments, and qual vals
	int buffPos = 0; //current position in the buffer

	Fragment_Index index1; 
	Fragment_Index index2; 
	Labels label1; 
	Labels label2; 

	int currSet1 = startSet;
	int currSet2 = startSet + numSets;
	int currOutSet = startSet;

	long int currFragment1 = 0;
	long int currFragment2 = 0;

	while((currSet1 < (startSet + numSets) || currFragment1 < index1.numFragments())  && (currSet2 < (startSet + 2 * numSets) || currFragment2 < index2.numFragments()))
	{
		//Make into functions
		if(currFragment1 == index1.numFragments())
		{
			char * inFile = new char [200];
			sprintf(inFile, "%s/sequences.%d", inDirF, currSet1);

			if(fileExists(inFile))
			{
				readData(inDirF, inDirL, index1, label1, currSet1);
				currFragment1 = 0;
				delete [] inFile;
				currSet1++;

			
			}else{
				delete [] inFile;
				break;
			}
		}

		if(currFragment2 == index2.numFragments())
		{
			char * inFile = new char [200];
			sprintf(inFile, "%s/sequences.%d", inDirF, currSet2);

			if(fileExists(inFile))
			{
				readData(inDirF, inDirL, index2, label2, currSet2);
				currFragment2 = 0; 
				delete [] inFile;
				currSet2++;
			}else{
				delete [] inFile;
				break;
			}
		}

		while(currFragment1 < index1.numFragments() && currFragment2 < index2.numFragments())
		{
			//Obtain the fragments
			pair<unsigned long long int, unsigned long long int> bounds1; 		
			pair<unsigned long long int, unsigned long long int> bounds2; 

			bounds1 = index1.indexBounds(currFragment1);
			bounds2 = index2.indexBounds(currFragment2);

			//Obtain the Labels 
			string l1 = label1.labelAt(currFragment1);
			string l2 = label2.labelAt(currFragment2);

			int l1Pos = l1.find_last_of('.');
			int l2Pos = l2.find_last_of('.');

			int l1Rank = atoi(l1.substr(l1Pos+1).c_str());
			int l2Rank = atoi(l2.substr(l2Pos+1).c_str());
		
			//Check to see if there will be buffer overflow ... write to file
			int addSize = 0;
			addSize += 2*((bounds1.second - bounds1.first) + 2*(bounds2.second - bounds2.first)); 
			addSize += (l1.size() + l2.size());
			addSize += 6;

			if(buffPos+addSize > 22 * megabyte)
			{
				if(currOutSet < startSet + numSets * 2)
				{
					writeBuff(fragmentBuff, buffPos, outDirF, currOutSet++);
				}else{
					appendBuff(fragmentBuff, buffPos, outDirF, currOutSet-1);
				}			
				buffPos = 0;
			}
			

			if(l1Rank < l2Rank)
			{
				appendFragment(fragmentBuff, buffPos, index1, bounds1, l1);
				currFragment1++;
			}
			else
			{
				appendFragment(fragmentBuff, buffPos, index2,  bounds2, l2);
				currFragment2++;
			}

		}
	}
		

	int track = 0;
	//This may not work for all cases if datasets are uneven
	while(currSet1 < (startSet + numSets) || currFragment1 < index1.numFragments())
	{
		if(currFragment1 == index1.numFragments())
		{
			char * inFile = new char [200];
			sprintf(inFile, "%s/sequences.%d", inDirF, currSet1);
			if(fileExists(inFile))
			{
				readData(inDirF, inDirL, index1, label1, currSet1);
				currFragment1 = 0; 
				delete [] inFile;
				currSet1++;
			}else{
				delete [] inFile;
				break;
			}
			track++;
		}else{
			char * inFile = new char [200];
			sprintf(inFile, "%s/sequences.%d", inDirF, currSet1-1);
			track++;
			delete [] inFile;
		}

		pair<unsigned long long int, unsigned long long int> bounds1; 		
		bounds1 = index1.indexBounds(currFragment1);
		string l1 = label1.labelAt(currFragment1);

		int addSize = 0;
		addSize += 2*(bounds1.second - bounds1.first); 
		addSize += l1.size();
		addSize += 6;

		if(buffPos+addSize > 22 * megabyte)
		{
			if(currOutSet < startSet + numSets * 2)
			{
				writeBuff(fragmentBuff, buffPos, outDirF, currOutSet++);
			}else{
				appendBuff(fragmentBuff, buffPos, outDirF, currOutSet-1);
			}			
			buffPos = 0;
		}
		appendFragment(fragmentBuff, buffPos, index1,  bounds1, l1);
		currFragment1++;
	}

	track = 0;
	while(currSet2 < (startSet + 2*numSets) || currFragment2 < index2.numFragments())
	{
		if(currFragment2 == index2.numFragments())
		{
			char * inFile = new char [200];
			sprintf(inFile, "%s/sequences.%d", inDirF, currSet2);
			if(fileExists(inFile))
			{
				readData(inDirF, inDirL, index2, label2, currSet2);
				currFragment2 = 0;
				delete [] inFile;
				currSet2++;
			}else{
				delete [] inFile;
				break;
			}
			track++;
		}else{
			
			char * inFile = new char [200];
			sprintf(inFile, "%s/sequences.%d", inDirF, currSet2-1);
			track++;
			delete [] inFile;

		}

		pair<unsigned long long int, unsigned long long int> bounds2; 		
		bounds2 = index2.indexBounds(currFragment2);
		string l2 = label2.labelAt(currFragment2);

		int addSize = 0;
		addSize += 2*(bounds2.second - bounds2.first); 
		addSize += l2.size();
		addSize += 6;

		if(buffPos+addSize > 22 * megabyte)
		{
			if(currOutSet < startSet + numSets * 2)
			{
				writeBuff(fragmentBuff, buffPos, outDirF, currOutSet++);
			}else{
				appendBuff(fragmentBuff, buffPos, outDirF, currOutSet-1);
			}			
			buffPos = 0;
		}
		appendFragment(fragmentBuff, buffPos, index2,  bounds2, l2);
		currFragment2++;
	}	

	if(buffPos > 0)
	{
		if(currOutSet < startSet + numSets * 2)
		{
			writeBuff(fragmentBuff, buffPos, outDirF, currOutSet++);
		}else{
			appendBuff(fragmentBuff, buffPos, outDirF, currOutSet-1);
		}			
	}

	delete [] fragmentBuff;
	makeIndexes(outDirF, outDirF, outDirL, startSet, startSet + 2 * numSets);
}

//bool fileExists( char fileName [] )
//Description: This function checks to see if a file exists
//Input: fileName (char []) the name of the file
//Output:None
//Return:None
bool MergeFragments::fileExists(const char fileName []) {
	
	FILE * pFile = fopen(fileName, "r");
	if(pFile != '\0')
	{
		fclose(pFile);
		return true;
	}
	return false;
}

//makeIndexes(int, int)
//Description: This function writes a buffer to a fragment index
//Input: startSet (int) : the starting sequence set, endSet (int) 
//the ending sequence set
//Output: None
//Return: None 
void MergeFragments::makeIndexes(char inDir [], char outDirF [], char outDirL [],  int startSet,  int endSet)
{
	long long int runningTotalFragments = 0; long long int runningTotalnumNucl = 0; int labelLen = 0;
	for(int i = startSet; i < endSet; i++)
	{
		char * fileName = new char [200];
		sprintf(fileName, "%s/tmp%d", inDir, i);

		if(fileExists(fileName))
		{
			FILE * pFile = fopen(fileName, "r");
			fseek(pFile, 0, SEEK_END);
			int size = ftell(pFile);
			fseek(pFile, 0, SEEK_SET);

			char * buff = new char [size+1];
			fread(buff, sizeof(char), size, pFile);
			fclose(pFile);

			buff[size] = '\0'; 	

			char * begin;
			begin = buff;

			char * end = begin;
			int numFragments = 0; long long int numNucl = 0;

			while(end != (buff+(size-1)))
			{
				end = strchr(begin, '\0');
				if(end-begin > labelLen)
					labelLen = (end-begin);

				begin = end+1;			
				end = strchr(begin, '\0');
				numNucl+=(end-begin);

				begin = end+1;
				end = strchr(begin, '\0');
				numFragments++;

				begin = end+1;
			}		

			Fragment_Index index;
			index.reserve(numFragments, numNucl, runningTotalFragments, runningTotalnumNucl);

			Labels label;
			label.reserve(labelLen, numFragments, runningTotalFragments);


			runningTotalFragments+=numFragments;
			runningTotalnumNucl+=numNucl;

			char * labelB = nullptr; char * labelE = nullptr; char * seqB = nullptr; 
			char * seqE = nullptr; char * qualB = nullptr; char * qualE = nullptr;

			labelB = buff; 

			while(qualE != (buff+(size-1)))
			{
				labelE = strchr(labelB, '\0');

				labelLen = (labelE-labelB);
				seqB = labelE + 1;
				seqE = strchr(seqB, '\0');
				qualB = seqE + 1;
				qualE = strchr(qualB, '\0');	

				index.append(seqB, qualB);
				label.append(labelB, labelLen);

				labelB = qualE+1;
			}		

			delete [] buff;		

			index.set();

			char * fileNameF = new char [200];
			char * fileNameL = new char [200];
			sprintf(fileNameF, "%s/sequences.%d", outDirF, i);
			sprintf(fileNameL, "%s/labels.%d", outDirL, i);

			FILE * pFile1 = fopen(fileNameF, "w");
			FILE * pFile2 = fopen(fileNameL, "w");

			index.write(pFile1);
			label.write(pFile2);

			fclose(pFile1);
			fclose(pFile2);

			delete [] fileNameF;
			delete [] fileNameL;
		}

		char * sysCmd = new char [500];
		sprintf(sysCmd, "rm -r %s", fileName);

		if(fileExists(fileName))
		{
			system(sysCmd);
		}

		delete [] sysCmd;
		delete [] fileName;
	} 
}

//writeFinalFragments(Index &, int, int, char * &, int, int, int)
//Description: This function writes the final fragments to the file
//Input: index (Index &): the index to be written to file, currFragment (int)
//the currentFragment, fragQualBuff (char * &): the buffer holding the fragments
//buggPos (int) : the current position in the buffer, currOutSet (int): the current
//output file, maxSet (int) : the maximum output file 
//Output: 
void MergeFragments::writeFinalFragments(Fragment_Index & index, Labels & label, int currFragment, char * & fragmentBuff, int buffPos, int currOutSet, int maxSet, char outDirF [])
{
	for(int i = currFragment; i < index.numFragments(); i++)
	{
		pair<unsigned long long int, unsigned long long int> bounds;
		bounds = index.indexBounds(i);

		string l = label.labelAt(i);

		int appSize = 2*(bounds.second - bounds.first);
		appSize += (l.size());
		appSize += 3;

		if(buffPos + appSize >= 22 * megabyte)
		{
			if(currOutSet < maxSet)
			{
				writeBuff(fragmentBuff, buffPos, outDirF, currOutSet++);
			}else{
				appendBuff(fragmentBuff, buffPos, outDirF, currOutSet-1);
			}			

			buffPos = 0;
		} 

		appendFragment(fragmentBuff, buffPos, index, bounds, l);
	}	

	if(buffPos > 0)
	{
		if(currOutSet < maxSet)
		{
			writeBuff(fragmentBuff, buffPos, outDirF, currOutSet++);
		}else{
			appendBuff(fragmentBuff, buffPos, outDirF, currOutSet-1);
		}			
	} 
}

//void appendFragment(char * &, int, pair<long long int, long long int>,
//Description: This appends a single fragment to the buffer
//pair<long long int, long long int>, string, string)
//Input: fragQualBuff (char * &), the buffer that hold the fragments, quals,
//and labels, buffPos (int) : the current position in the buffer,
//bounds (pair<long long int, long long int>) : fragment bounds in the index
//Output:None
//Return:None 
void MergeFragments::appendFragment(char * & fragQualBuff, int & buffPos, Fragment_Index & index, pair<unsigned long long int, unsigned long long int> & bounds, string & l)
{
	for(int i = 0; i < l.size(); i++)
		fragQualBuff[buffPos++] = l.at(i);

	fragQualBuff[buffPos++] = '\0';

	for(long long int i = bounds.first; i != bounds.second; i++)
		fragQualBuff[buffPos++] = index.at(i);

	fragQualBuff[buffPos++] = '\0';

	for(long long int i = bounds.first; i != bounds.second; i++)
		fragQualBuff[buffPos++] = index.qualAt(i) + 33;

	fragQualBuff[buffPos++] = '\0';
}

//appendBuff(char * &, int, char * &, int)
//Description: This function appends a buffer of sequences to a file
//Input: fragQualBuff (char * &) : the buff into which the fragments, labels,
//and quals are written, buffPos (int) : the current position in the buffer,
//outDirF (char * &) : the output directory, currOutSet (int) : The current
//output directory
//Output: None
//Return: None
void MergeFragments::appendBuff(char * & fragQualBuff,  int buffPos,  char outDirF [],  int currOutSet)
{
	char * fileName = new char [200];
	sprintf(fileName, "%s/tmp%d", outDirF, currOutSet);

	FILE * pFile = fopen(fileName, "a");
	fwrite(fragQualBuff, sizeof(char), buffPos, pFile);
	fclose(pFile);

	delete [] fileName;
}

//writeBuff(char * &, int, char * &, int)
//Description: This function writes a buffer to a file
//Input: fragQualBuff (char * &) : the buff into which the fragments, labels,
//and quals are written, buffPos (int) : the current position in the buffer,
//outDirF (char * &) : the output directory, currOutSet (int) : The current
//output directory
//Output: None
//Return: None
void MergeFragments::writeBuff( char * & fragQualBuff,  int buffPos,  char outDirF [],  int currOutSet)
{
	char * fileName = new char [200];
	sprintf(fileName, "%s/tmp%d", outDirF, currOutSet);

	FILE * pFile = fopen(fileName, "w");
	fwrite(fragQualBuff, sizeof(char), buffPos, pFile);
	fclose(pFile);

	delete [] fileName;
}

//bool readData( char [],  char [], Fragment_Index &, Labels &, int)
//Description: This function reads an index from a file
//Input: fileName (char []) the name of the file
//Output:None
//Return:None
void MergeFragments::readData(char inDirF [],  char inDirL [], Fragment_Index & index, Labels & label, int currSet)
{
	char * inFile = new char [200];
	sprintf(inFile, "%s/sequences.%d", inDirF, currSet);
	FILE * pFile = fopen(inFile, "r");
	
	unsigned long long int fragmentOffset;
	fread(&fragmentOffset, sizeof(unsigned long long int), 1, pFile);

	unsigned long long int indexOffset;
	fread(&indexOffset, sizeof(unsigned long long int), 1, pFile);

	unsigned long long int intervalLength;
	fread(&intervalLength, sizeof(unsigned long long int), 1, pFile);

	unsigned long long int numIntervalFragments;
	fread(&numIntervalFragments, sizeof(unsigned long long int), 1, pFile);
	fclose(pFile);

	index.reserve(numIntervalFragments, intervalLength, indexOffset, fragmentOffset);

	pFile = fopen(inFile, "r");
	index.read(pFile);
	fclose(pFile);

	index.set();

	sprintf(inFile, "%s/labels.%d", inDirL, currSet);

	unsigned int tagLen;
	unsigned long long int intervalOffset;
	unsigned long long int intervalNumLabels;

	pFile = fopen(inFile, "r");
	fread(&tagLen, sizeof(unsigned int), 1, pFile);
	fread(&intervalOffset, sizeof(unsigned long long int), 1, pFile);
	fread(&intervalNumLabels, sizeof(unsigned long long int), 1, pFile);
	fclose(pFile);

	label.reserve(tagLen, intervalNumLabels, intervalOffset);	

	pFile = fopen(inFile, "r");
	label.read(pFile);
	fclose(pFile);

	delete [] inFile;
}  

//bool fileExists( char [])
//Description: This function checks to see if a file exists
//Input: fileName (char []) the name of the file
//Output:None
//Return:None
bool MergeFragments::fileExists(char fileName []) {

	FILE * pFile = fopen(fileName, "r");

	if(pFile != '\0')
	{
		fclose(pFile);
		return true;
	}

	return false;
}
