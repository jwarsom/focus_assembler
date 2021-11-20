#include<iostream>
#include<cstdlib>
#include<fstream>
#include <stdint.h>
#include <cstring>
#include <cstdio>
#include <dirent.h>
#include<cmath>
using namespace std;
#include "../src/MappingValues.h"
#include "../src/Dictionary.h"
#include "../src/Graph.h"

int countFiles(const char []);
float calcHval(Graph &, int, int);

int main(int argc, char ** argv)
{
	int numFiles = countFiles(argv[1]);

	int graphsLoaded = numFiles-2;
	/*
	   for(int i = 0; i < graphsLoaded; i++)
	   {
	   FILE * pFile;

	   char * fileName = new char [200];
	   sprintf(fileName, "%s/sequences.job%d.info", argv[1], i);

	   pFile = fopen(fileName, "r");

	   fseek(pFile, 0, SEEK_END);
	   int size = ftell(pFile)/sizeof(int);
	   fseek(pFile, 0, SEEK_SET);		

	   int * map = new int[size];
	   fread(map, sizeof(int), size, pFile);

	   for(int j = 1; j < size; j+=3)
	   cout<<map[j]<<endl;

	   fclose(pFile);
	   }*/

	int * nMap;
	{
		char * fileName = new char [1000];
		sprintf(fileName, "%s/nDen", argv[1]);
		FILE * pFile = fopen(fileName, "r");
		fseek(pFile, 0, SEEK_END);
		int size = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);

		nMap = new int [size];
		fread(nMap, sizeof(int), size, pFile);
		fclose(pFile);
		delete [] fileName;
	}


	for(int i = 0; i < graphsLoaded; i++)
	{
		FILE * pFile;

		char * fileName = new char [200];
		sprintf(fileName, "%s/graph%d", argv[1], i);

		pFile = fopen(fileName, "r");

		Graph graph;
		graph.read(pFile);
		fclose(pFile);


		for(int j = 0; j < graph.getNumNodes(); j++)
		{

			float Hval = calcHval(graph, j, 0);
		
			bool coverage = false;
			for(int k = 0; k < graph.nodeDegree(j); k++)
			if(graph.edgeOvlLen(j, k) >= 2000)
				coverage = true; 

			if(nMap[j] && coverage)
			{
				if(Hval < -2.5)
				{
					cout<<j<<" "<<-1 * Hval<<endl;
				}else{
				//	cout<<j<<" "<<0<<endl;
				}
			}

		}

		delete [] fileName;
	}

	delete [] nMap;

	return 0;

}

float calcHval(Graph & graph, int node, int numRounds)
{
	int numInGraph = 0;
	float totalWeight = 0;
	for(int k = 0; k < graph.nodeDegree(node); k++)
	{
		if(graph.edgeInGraph(node, k))
		{
			numInGraph++;
			totalWeight+=graph.edgeOvlLen(node, k);
		}
	}
	
	float Hval = 0;
	for(int k = 0; k < graph.nodeDegree(node); k++)
	{
		if(graph.edgeInGraph(node, k))
		{

			float edgeWeight = graph.edgeOvlLen(node, k);
			Hval += (edgeWeight/totalWeight * log(edgeWeight/totalWeight));
			//cout<<j+graph.getOffset()<<"\t"<<graph.edgeDest(j, k)<<endl;
		}

	}

	if(numRounds == 0)
	{
		return Hval;
	}else{
		float minHval = 0;
		for(int k = 0; k < graph.nodeDegree(node); k++)
		{
			if(graph.edgeInGraph(node, k))
			{
				float tmpHval = calcHval(graph, graph.edgeDest(node, k), numRounds-1);

				if(tmpHval < minHval)
					minHval = tmpHval;
			}


		}
		return minHval;
	}
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

	return numFiles;
}
