/*
 * MergeEdges.cpp
 *
 *      Author: Julia
 */

#include<iostream>
#include<stdint.h>
#include<cstdio>
#include<cstring>
#include<queue>
#include<dirent.h>
using namespace std;
#include "Dictionary.h"
#include "Fragment_Index.h"
#include "MappingValues.h"
#include "MergeEdges.h"



//MergeEdges(void)
//Description: The constructor
//Input: mt (const int): The type of merging
//Output:None
//Return:None
MergeEdges::MergeEdges(const int mt)
	:pFile(new FILE * [10]), fileInfo(new uint32_t [60]), mergeType(mt) {;;;}

	//MergeEdges(void)
	//Description: The destructor
	//Input:None
	//Output:None
	//Return:None
	MergeEdges::~MergeEdges(){
		delete [] pFile;
		delete [] fileInfo;
	}

//void merge(const int, const int, const int, const int, const char [], const char [])
//Description: This function merges two sorted files
//Input: startSet (A const int): The starting file number of the first set, startSet2 (a const int): the
//starting number of the second set, numSets1 (A const int): The number of files in  set1,
//numSets2 (A const int): The number of files in set2, inDir (const char []): the name of the input directory
//outDir (const char []): The output directory
//Output:None
//Return:None
void MergeEdges::merge(const int startSet, const int numSets, const char inDir [], const char outDir [])
{
	fileCmp fCmp(fileInfo);
	priority_queue<int, vector<int>, fileCmp > myQueue(fCmp);

	uint32_t * iBuff = new uint32_t [60 * megabyte];
	uint32_t * oBuff = new uint32_t [33 * megabyte];
	int oBuffCount = 0; int oFileCount = 0;

	for(int i = 0; i < 2/*10*/; i++)
		if(!initiateQueue(myQueue, iBuff, inDir, startSet, i, numSets))
			break;

	while(!myQueue.empty())
	{
		int next = myQueue.top();
		myQueue.pop();

		if(oBuffCount >= 33 * megabyte)
		{
			int rollBack = 3;
			while(oBuff[33*megabyte-rollBack] == oBuff[33*megabyte-3] && rollBack < 33*megabyte)
				rollBack+=3;

			rollBack-=3;

			if(oFileCount < 2/*10*/ * numSets || !fileExists(inDir, startSet + oFileCount))
			{
				writeFile(oBuff, 33*megabyte-rollBack, outDir, startSet + oFileCount++);
			}
			else
			{
				appendFile(oBuff, 33*megabyte-rollBack, outDir, startSet + oFileCount-1);
			}

			memmove(oBuff, oBuff+(33*megabyte - rollBack), rollBack * sizeof(uint32_t));

			oBuffCount = rollBack;
		}
		
		addEdge(oBuff, oBuffCount, fileInfo[30+next], fileInfo[40+next], fileInfo[50+next]);
		getNextEdge(myQueue, inDir, iBuff, startSet, next, numSets);
	}

	if(oBuffCount > 0)
	{

		if((oFileCount < 2/*10*/ * numSets && fileExists(inDir, startSet + oFileCount)) || (!fileExists(inDir, startSet + oFileCount) && oFileCount > 11 * megabyte))
		{
			writeFile(oBuff, oBuffCount, outDir, startSet + oFileCount++);
		}
		else
		{
			appendFile(oBuff, oBuffCount, outDir, startSet + oFileCount-1);
		}
	}

	delete [] iBuff;
	delete [] oBuff;
}

//bool fileExists( char inDir [], const int set )
//Description: This function checks to see if a file exists
//Input: fileName (char []) the name of the file
//Output:None
//Return:None
bool MergeEdges::fileExists(const char inDir [], const int set) {

	char * fileName = new char [1000];
	sprintf(fileName, "%s/set%d", inDir, set);
	FILE * pFile = fopen(fileName, "r");

	if(pFile != '\0')
	{
		fclose(pFile);
		return true;
	}

	return false;
}

//void getNextEdge(priority_queue, char [], uint32_t [], int, int)
//Description: This function gets the next edge from a buffer and
//adds it to the queue
//Input: myQueue (priority_queue &): The queue, inDir (char []): the
//input directory, iBuff (uint32_t [] ): The input buffer, next (int):
//the current file being read, numSets (int): The number of sorted files per
//subset
//Output:None
//Return:None
void MergeEdges::getNextEdge(priority_queue<int, vector<int>, fileCmp> & myQueue, const char inDir [], uint32_t iBuff [], const int startSet, const int next, const int numSets){

	if(fileInfo[10+next] >= fileInfo[20+next])
	{
		readFile(iBuff, inDir, startSet, next, numSets);
	}

	if(fileInfo[20+next] > 0)
	{
		fileInfo[30 + next] = iBuff[(next*6*megabyte) + fileInfo[10+next]];
		fileInfo[40 + next] = iBuff[(next*6*megabyte) + fileInfo[10+next]+1];
		fileInfo[50 + next] = iBuff[(next*6*megabyte) + fileInfo[10+next]+2];

		fileInfo[10+next]+=3;
		myQueue.push(next);
	}
}

//FILE * initiateQueue(char [], int)
//Description: This function initiates the queue
//Input:inDir (char []): The input directory, i (int):
//the current file buffer
//Output:None
//Return: FILE *: The pointer to the file that was read from
FILE * MergeEdges::initiateQueue(priority_queue<int, vector<int>, fileCmp> & myQueue, uint32_t iBuff [], const char inDir [], const int startSet, const int i, const int numSets)
{
	char * fileName = new char [1000];
	sprintf(fileName, "%s/set%d", inDir, startSet + i * numSets);
	pFile[i] = fopen(fileName, "r");
	delete [] fileName;
	

	if(pFile[i] == '\0')
		return pFile[i];

	fileInfo[i] = 1; //file set
	fileInfo[10+i] = 3; //buffer index //TRY THIS 
	fileInfo[20+i] = fread(iBuff+(i*6*megabyte), sizeof(uint32_t), 6*megabyte, pFile[i]);
	fileInfo[30+i] = iBuff[(i*6*megabyte)]; //source
	fileInfo[40+i] = iBuff[(i*6*megabyte)+1]; //dest
	fileInfo[50+i] = iBuff[(i*6*megabyte)+2]; //weight

	myQueue.push(i);

	return pFile[i];
}

//void addEdge(uint32_t [], int &, uint32_t, uint32_t, uint32_t)
//Description: This function checks for duplicate edges and
//writes an edge to a buffer
//Input:buff (uint32_t []): The output buffer, index (int &): the current position in the buffer,
//source (uint32_t): The edge source, dest (uint32_t): The edge destination, weight (uint32_t): The
//weight of the edge
//Output:None
//Return:None
void MergeEdges::addEdge(uint32_t oBuff [], int & oBuffCount, uint32_t source, uint32_t dest, uint32_t weight)
{
	if(oBuffCount < 3 || oBuff[oBuffCount-3] != source || oBuff[oBuffCount-2] != dest)
	{
		oBuff[oBuffCount++] = source;
		oBuff[oBuffCount++] = dest;
		oBuff[oBuffCount++] = weight;

	}else
	{
		if(mergeType == 1) //Merge Weights
		{
			uint32_t weight1 = weight >> 16;
			uint32_t weight2 = oBuff[oBuffCount-1] >> 16;
			uint32_t mergedWeight = weight1+weight2;

			uint32_t rDovetail1 = weight << 19;
			uint32_t rDovetail2 = oBuff[oBuffCount-1] << 19;

			rDovetail1 = rDovetail1 >> 31;
			rDovetail2 = rDovetail2 >> 31;

			uint32_t sContained1 = weight << 17;
			uint32_t dContained1 = weight << 18;
			sContained1 = sContained1 >> 31;
			dContained1 = dContained1 >> 31;

			uint32_t sContained2 = oBuff[oBuffCount-1] << 17;
			uint32_t dContained2 = oBuff[oBuffCount-1] << 18;
			sContained2 = sContained2 >> 31;
			dContained2 = dContained2 >> 31;

			rDovetail1 = (rDovetail1 & ~sContained1 & ~dContained1) | (rDovetail2 & ((sContained1 | dContained1) & (~sContained2 & ~dContained2))); 
			rDovetail2 = (rDovetail2 & ~sContained2 & ~dContained2) | (rDovetail1 & ((sContained2 | dContained2) & (~sContained1 & ~dContained1))); 

			uint32_t sContained = sContained1 & sContained2;			
			uint32_t dContained = dContained1 & dContained2;
			uint32_t rDovetail = 0;

			if(!(sContained | dContained))
			{
				rDovetail = (weight1 > weight2) ? rDovetail1 : rDovetail2; 
			}else{
				rDovetail = rDovetail1 | rDovetail2;
			}

			oBuff[oBuffCount-1] = (mergedWeight << 16) | (sContained << 14) | (dContained << 13) | (rDovetail << 12);

		}else
		{
			uint32_t reverse = weight << 16;
			uint32_t sContained1 = weight << 17;
			uint32_t dContained1 = weight << 18;
			uint32_t rDovetail = weight << 19;
			uint32_t trimmings = weight << 20;
			sContained1 = sContained1 >> 31;
			dContained1 = dContained1 >> 31;
			reverse = reverse >> 31;
			rDovetail = rDovetail >> 31;
			trimmings = trimmings >> 20;

			uint32_t sContained2 = oBuff[oBuffCount-1] << 17;
			uint32_t dContained2 = oBuff[oBuffCount-1] << 18;
			sContained2 = sContained2 >> 31;
			dContained2 = dContained2 >> 31;

			weight = weight >> 16;

			oBuff[oBuffCount-1] = 0; 

			uint32_t rDovetail1 = 1;
			uint32_t rDovetail2 = 1;

			rDovetail1 = rDovetail1 ^ sContained1;
			rDovetail2 = rDovetail2 ^ dContained1;

			rDovetail = (rDovetail & rDovetail1 & rDovetail2);

			oBuff[oBuffCount-1] = (weight << 16) | (reverse << 15) | (sContained1 << 14) | (dContained1 << 13); 
			oBuff[oBuffCount-1] = oBuff[oBuffCount-1] | (sContained2 << 14) | (dContained2 << 13) | (rDovetail << 12) | (trimmings); 

		}
	}
}

//void readFile(uint32_t [], long long int &, FILE * &, const char [], int set)
//Description: This function reads edges from a file
//Input: buff (uint32_t []): the buffer that holds edges, fileSize (long long int): the size of the file,
//pFile (FILE *): the pointer to file, dir (const char): The files directory, set (int): the file number
//Output:None
//Return:None
void MergeEdges::readFile(uint32_t iBuff [], const char inDir [], const int startSet, const int next, const int numSets)
{
	fileInfo[20+next] = fread(iBuff+(next*6*megabyte), sizeof(uint32_t), 6*megabyte, pFile[next]);
	if(fileInfo[20+next] == 0)
	{
		fclose(pFile[next]);
		if(fileInfo[next] < numSets)
		{
			char * file = new char [500];
			sprintf(file, "%s/set%d", inDir, startSet + numSets * next + fileInfo[next]++); //add start set??
			pFile[next] = fopen(file, "r");
			delete [] file;

			if(pFile[next] != '\0')
				fileInfo[20+next] = fread(iBuff+(next*6*megabyte), sizeof(uint32_t), 6*megabyte, pFile[next]);
		}
	}

	fileInfo[10+next] = 0;
}

//void writeFile(uint32_t [], int size, const char [], int)
//Description: This function writes a buffer to a file
//Input: buff (uint32_t []): the buffer that holds edges, size (int): the size of the buffer,
// set (int): the file number
//Output:None
//Return:None
void MergeEdges::writeFile(uint32_t buff [], const int size, const char dir [], const int set){

	char * file = new char [500];
	sprintf(file, "%s/set%d", dir, set);
	FILE * oFile = fopen(file, "w");
	delete [] file;
	fwrite(buff, sizeof(uint32_t), size, oFile);
	fclose(oFile);
}

//void appendFile(uint32_t [], int size, const char [], int)
//Description: This function appends a buffer to a file
//Input: buff (uint32_t []): the buffer that holds edges, size (int): the size of the buffer,
// set (int): the file number
//Output:None
//Return:None
void MergeEdges::appendFile(uint32_t buff [], const int size, const char dir [], const int set){

	char * file = new char [500];
	sprintf(file, "%s/set%d", dir, set);
	FILE * oFile = fopen(file, "a");
	delete [] file;
	fwrite(buff, sizeof(uint32_t), size,  oFile);
	fclose(oFile);
}
