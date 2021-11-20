#include <iostream>
#include <map>
#include <string>
#include "mpi.h"
#include <algorithm>
#include <cstdio>
#include <stdint.h>
#include <vector>
#include <ctime>
#include <list>
#include <set>
#include <climits>
#include <queue>
#include <dirent.h>
using namespace std;
#include "MappingValues.h"
#include "Dictionary.h"
#include "Graph.h"
#include "Fragment_Index.h"




int main(int argc, char * argv [])
{	
	//User Input Flags
	string workDirS = "--workDir";               //Current Working Directory
	string timeS = "--time"; //Report the time

	string workDir = "";
	bool reportTime = false;

	for(int i = 1; i < argc; i+=2)
	{
		//The working directory
		if(argv[i] == workDirS)
		{
			workDir = argv[i+1];
			unsigned int pos = workDir.find_last_of("/");
			if(pos == workDir.size()-1)
				workDir.erase(pos);
		}

		//Display time
		if(argv[i] == timeS && reportTime != true)
		{
			reportTime = true;
			i--;
		}
	}

	char * fileName = new char [1000];

	sprintf(fileName, "%s/HybSpectrum/Graph0", workDir.c_str());

	int numFiles = countFiles(fileName);
	int graphsLoaded = numFiles-2;

	Graph ** graphs = new Graph * [graphsLoaded];
	for(int i = 0; i < graphsLoaded; i++)
		graphs[i] = new Graph();

	for(int i = 0; i < graphsLoaded; i++)
	{
		FILE * pFile;

		char * fileName = new char [200];
		sprintf(fileName, "%s/HybSpectrum/Graph0/graph%d", workDir.c_str(), i);

		pFile = fopen(fileName, "r");

		graphs[i]->read(pFile);

		fclose(pFile);
		delete [] fileName;
	}

	sprintf(fileName, "%s/
	

	for(int i = 0; i < graphsLoaded; i++)
	{
		for(int j = 0; j < graphs[i].numNodes(); j++)
		{
			if(graphs[i].nodeDegree(j) > 2)
			{
				
			}
		}
	}

	if(reportTime == true)
	{
		time_t rawtime;
		struct tm * timeinfo;
		time(&rawtime);
		timeinfo = localtime(&rawtime);
		cout<<"start time: "<<asctime(timeinfo)<<endl;
	}

			
	return 0;
}

//void countFiles(const char [])
//Description: This function counts files in a directory
//Input: dir (cosnt char []) : The directory where we are counting the files
//Output:None
//Return:None
int countFiles(const char dir [])
{
	DIR *pDir = opendir(dir);
	struct dirent * pDirent;

	//Count Files
	int numFiles = 0;
	while ((pDirent = readdir(pDir)) != NULL)
		if(strcmp(pDirent->d_name, ".") != 0 && strcmp(pDirent->d_name, "..") != 0)
			numFiles++;

	closedir(pDir);
