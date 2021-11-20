#include <iostream>
#include <unistd.h>
#include <sys/stat.h>
#include <stdint.h>
#include <cstdlib>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <ctime>
#include<climits>
#include <set>
#include <cmath>
#include <stack>
#include <queue>
#include <dirent.h>
#include <map>
using namespace std;
#include "../src/Dictionary.h"
#include "../src/Graph.h"

#define bufOverflow 792

long long int getEdgeCut(Graph ** &, int,  map<int, int> &);
int countFiles(const char []);

int main(int argc, char * argv [])
{
	string workDir = "NO_WORK_DIR";
	string distribDir = "DistributedGraph";

	string workDirS = "--workDir"; //Current Working Directory
	string distribDirS = "--distributedDir"; //The Distributed Graph Dir

	for(int i = 1; i < argc; i+=2)
	{
		//The working directory
		if(argv[i] == workDirS)
		{
			workDir = argv[i+1];
			unsigned int pos = workDir.find_last_of("/");
			if(pos == workDir.size()-1)
				workDir.erase(pos);

			//Don't want to segfault
			if(workDir.size() > 200)
			{
				cout<<"Please keep the working directory ";
				cout<<"name less than 200 characters in ";
				cout<<"length"<<endl;
				exit(bufOverflow);
			}
		}

		//Getting the distributed directory
		if(argv[i] == distribDirS)
		{
			distribDir = argv[i+1];		
		}
	}
	cout<<"HERE 1"<<endl;

	char * graphDir = new char[1000];
	sprintf(graphDir, "%s/Spectrum/Graph0", workDir.c_str());
	int numFiles = countFiles(graphDir);
        int graphsLoaded = numFiles-2;

	delete [] graphDir;

        Graph ** graphs = new Graph * [graphsLoaded];
        for(int i = 0; i < graphsLoaded; i++)
                graphs[i] = new Graph();

	cout<<"HERE "<<endl;

        for(int i = 0; i < graphsLoaded; i++)
        {
                FILE * pFile;

                char * fileName = new char [200];
                sprintf(fileName, "%s/Spectrum/Graph0/graph%d", workDir.c_str(), i);

                pFile = fopen(fileName, "r");

                graphs[i]->read(pFile);

                fclose(pFile);
                delete [] fileName;
        }	

	char * distHybGraphDir = new char [1000];
        sprintf(distHybGraphDir, "%s/%s/Graph0", workDir.c_str(), distribDir.c_str());
	
	int numPartitions = countFiles(distHybGraphDir);

	map<int, int> sections;
	cout<<"HERE 2"<<endl;

	for(int i = 0; i < numPartitions; i++)
	{
		char * fileName = new char[1000];
		int part = log2(numPartitions); int section = i;
                sprintf(fileName, "%s/%s/Graph0/Partition%d_%d", workDir.c_str(), distribDir.c_str(), part, section);

		FILE * pFile = fopen(fileName, "r");
		fseek(pFile, 0, SEEK_END);
		int size = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);

		int * buff = new int [size];
		fread(buff, sizeof(int), size, pFile);

		fclose(pFile);

		for(int j = 0; j < size; j++)
		{
			sections.insert(pair<int, int>(buff[j], i));
		}
	
		delete [] fileName;
		delete [] buff;
	}
	cout<<"HERE 3"<<endl;

	getEdgeCut(graphs, graphsLoaded, sections);


	delete [] distHybGraphDir;
	return 0;
}

long long int getEdgeCut(Graph ** & graphs, int graphsLoaded,  map<int, int> & sections)
{
	long long int edgeCut = 0;
	long long int totalEdgeCut = 0;

	long long int numTheSame = 0;
	long long int total = 0;
	
	for(int i = 0; i < graphsLoaded; i++)
	{
		for(int j = 0; j < graphs[i]->getNumEdges(); j++)
		{
			if(graphs[i]->edgeInGraph(j))
			{
				totalEdgeCut+=graphs[i]->edgeOvlLen(j);

				int source = graphs[i]->edgeSource(j)+ graphs[i]->getOffset();
				int dest = graphs[i]->edgeDest(j);

				int part1 = (*sections.find(source)).second;
				int part2 = (*sections.find(dest)).second;
				
				total++;

				if(part1 != part2)
				{
					edgeCut+=graphs[i]->edgeOvlLen(j);
				}else{
					numTheSame++;
				}
			}
		}
	}

	cout<<edgeCut<<" "<<totalEdgeCut<<endl;

	return edgeCut;
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

