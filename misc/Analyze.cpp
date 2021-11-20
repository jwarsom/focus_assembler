/*
 * Analyze.cpp
 *
 *  Created on: Aug. 6th, 2012
 *      Author: Julia Warnke
 */


#include<iostream>
#include<vector>
#include<string>
#include <dirent.h>
using namespace std;
#include "MappingValues.h"
#include "Dictionary.h"
#include "Fragment_Index.h"
#include "Graph.h"

/* This program analyzes a series of overlap graphs according to various graph theoretic characteristics

   A list of features that I can implement:

   GC content in read in nodes
   Node Degree/Overlap Node Degree
   Node Density (ie. Node Size)
   Edge Density
   Gene Alignment
   Average Alignment Identity
   Cycles in the overlap graph
   Output cluster

   Zoom in and out

*/

//The help message
void help_message(){

	cerr<<"\n\n";
	cerr<<"\t\t\t Overlap Graph Analyzer"<<endl;
	cerr<<"\t\t\tWritten by Julia Warnke\n"<<endl;
	cerr<<"Command line arguments\n"<<endl;
	cerr<<"--workDir :working directory :required"<<endl;

	cerr<<"Options\n"<<endl;

	cerr<<"help : This message"<<endl;
	cerr<<"displayCurrentGraph     : Display the currently selected graph level."<<endl;
	cerr<<"displayCurrentNode      : Display the currently selected node."<<endl;
	cerr<<"selectGraph  (int)      : Select a graph level."<<endl;
	cerr<<"selectNode (int)        : Select a node in a graph."<<endl;
	cerr<<"avgGraphDegree          : Display the average graph degree."<<endl;
	cerr<<"avgGraphOvlLen          : Display the average graph overlap length."<<endl;
	cerr<<"avgGraphOvlIden         : Display the average graph overlap identity."<<endl;
	cerr<<"avgGraphGCcontent       : Display the average graph GC content."<<endl;
	cerr<<"avgGraphEdgeDen         : Display the average graph edge density."<<endl;
	cerr<<"avgGraphClusterSize     : Display the average graph cluster size."<<endl;
	cerr<<"selectClustersGCcontent (int) (min) (max)  : Select "<<endl;
	cerr<<"selectClustersEdgeDen (int) (min) (max)    : Select "<<endl;
	cerr<<"selectClustersSize (int) (min) (max)       : Select"<<endl;
	cerr<<"selectClustersNodeDegree (int) (min) (max) : Select"<<endl;
	cerr<<"selectClustersAvgChildDegree (int) (min) (max) : Select"<<endl;
	cerr<<"selectClustersAvgIden (int) (min) (max)        : Select"<<endl;
	cerr<<"selectClustersAvgOvlLen (int) (min) (max)      : Select"<<endl;
	cerr<<"displayClusterGCcontent (int) (int)            : Display"<<endl;
	cerr<<"displayClusterEdgeDen (int) (int)              : Display"<<endl;
	cerr<<"displayClusterSize (int) (int)                 : Display"<<endl;
	cerr<<"displayClusterNodeDegree (int) (int)           : Display"<<endl;
	cerr<<"displayClusterAvgChildDegree (int) (int)       : Display"<<endl;
	cerr<<"displayClusterAvgIden (int) (int)              : Display"<<endl;
	cerr<<"displayClusterAvgOvlLen (int) (int)            : Display"<<endl;
	cerr<<"displayClusterNeighbors (int) (int)            : Display"<<endl;
	cerr<<"displayClusterChildren (int) (int)             : Display"<<endl;
	cerr<<""<<endl;

	cerr<<"\n\n"<<endl;

	cerr<<"Exiting program"<<endl;
	exit(help);
}

int  countFiles(const char []);
bool fileExists(const char []);
void avgGraphDegree(string &, const int);
void avgGraphOvlLen(string &);
void avgGraphOvlIden(string &);
void avgGraphGCcontent(string &);
void avgGraphEdgeDen(string &, const int);
void avgGraphClusterSize(string &, const int);
void selectClustersGCcontent(string &, const int, const float );
//void selectClustersGenes(const int );
void selectClustersEdgeDen(string &, const int, const float, const bool);
void selectClustersSize(string &, const int, const float, const bool);
void selectClustersNodeDegree(string &, const int, const int, const bool);
void selectClustersAvgChildDegree(string &, const int, const int, const bool);
void selectClustersAvgIden(string &, const int, const float minIden, const bool);
void selectClustersAvgOvlLen(string &, const int, const int ovlLen, const bool);
//void selectGraphCycles(const int );
void displayClusterGCcontent(string &, const int, const int );
//void displayClusterGenes(string &, const int, const int );
void displayClusterEdgeDen(string &, const int, const int );
void displayClusterSize(string &, const int, const int );
void displayClusterNodeDegree(string &, const int, const int );
void displayClusterAvgChildDegree(string &, const int, const int );
void displayClusterAvgIden(string &, const int , const int);
void displayClusterAvgOvlLen(string &, const int, const int );
void displayClusterNeighbors(string &, const int, const int );
void displayClusterChildren(string &, const int, const int );
int recoverClusters(string &, int * &, int, int);

//The main function is used to communicate with
//the user and retrieve commands for interrogating
//and manipulating an overlap graph.
int main (int argc, char * argv [])
{
	string workDirS = "--workDir";
	string workDir = "NO_WORK_DIR";

	for(int i = 1; i < argc; i+=2)
	{
		//The working directory
		if(argv[i] == workDirS)
		{
			workDir = argv[i+1];
			unsigned int pos = workDir.find_last_of("/");
			if(pos == workDir.size()-1)
				workDir.erase(pos);

			//Check for overflows
			if(workDir.size() > 200)
			{
				cout<<"Please keep the working directory name less than 200 characters in length"<<endl;
				exit(bufOverflow);
			}
		}
	}

	if(workDir == "NO_WORK_DIR")
	{
		help_message();
	}

	bool iterate = true;
	int graphLevel = 0; int currentNode = 0;

	char * cmd = new char [100];

	//Collect Queries
	while(iterate)
	{
		cin>>cmd;
		char * ptr = strtok(cmd, " ");

		if(strcmp(ptr, "help")  == 0)
			help_message();

		if(strcmp(ptr, "AvgGraphDegree") == 0)
			avgGraphDegree(workDir, graphLevel);

		if(strcmp(ptr, "AvgGraphOvlLen") == 0)
			avgGraphOvlLen(workDir);

		if(strcmp(ptr, "AvgGraphOvlIden") == 0)
			avgGraphOvlIden(workDir);

		if(strcmp(ptr, "AvgGraphGCcontent") == 0)
			avgGraphGCcontent(workDir);

		if(strcmp(ptr, "avgGraphClusterSize") == 0)
			avgGraphClusterSize(workDir, graphLevel);

		if(strcmp(ptr, "avgGraphEdgeDen") == 0)
			avgGraphEdgeDen(workDir, graphLevel);

		if(strcmp(ptr, "selectGraph") == 0)
		{
			while(ptr != NULL)
			{
				ptr = strtok (NULL, " ");
				graphLevel = atoi(ptr);

				if(graphLevel != 0)
					break;
			}
		}

		if(strcmp(ptr, "selectNode") == 0)
		{
			while(ptr != NULL)
			{
				ptr = strtok (NULL, " ");
				currentNode = atoi(ptr);

				if(currentNode != 0)
					break;
			}
		}

	if(strcmp(ptr, "selectClustersGCcontent") == 0)
	{
		float minGC = 0; bool greaterThan = true;
		while(ptr != NULL)
		{
			ptr = strtok (NULL, " ");
			minGC = atof(ptr);

			if(strcmp(ptr, "max") == 0)
				greaterThan = false;

			if(strcmp(ptr, "min") == 0)
				greaterThan = true;
		}

		selectClustersGCcontent(workDir, graphLevel, minGC, greaterThan);
	}

	//if(strcmp(ptr, "selectClustersGenes") == 0)
		//selectClustersGenes(graphLevel);

	if(strcmp(ptr, "selectClustersEdgeDen") == 0)
	{
		float minDensity = 0;
		while(ptr != NULL)
		{
			ptr = strtok (NULL, " ");
			minDensity = atof(ptr);

			if(minDensity != 0)
				break;
		}
		selectClustersEdgeDen(workDir, graphLevel, minDensity);
	}

	if(strcmp(ptr, "selectClustersSize") == 0)
	{
		float minDensity = 0;
		while(ptr != NULL)
		{
			ptr = strtok (NULL, " ");
			minDensity = atof(ptr);

			if(minDensity != 0)
				break;
		}

		selectClustersSize(workDir, graphLevel, minDensity);
	}

	if(strcmp(ptr, "selectClustersNodeDegree") == 0)
	{
		int minDegree = 0;
		while(ptr != NULL)
		{
			ptr = strtok (NULL, " ");
			minDegree = atoi(ptr);

			if(minDegree != 0)
				break;
		}
		selectClustersNodeDegree(workDir, graphLevel, minDegree);
	}

	if(strcmp(ptr, "selectClustersAvgChildDegree") == 0)
	{
		int minDegree = 0;
		while(ptr != NULL)
		{
			ptr = strtok (NULL, " ");
			minDegree = atoi(ptr);

			if(minDegree != 0)
				break;
		}

		selectClustersAvgChildDegree(workDir, graphLevel, minDegree);
	}

	if(strcmp(ptr, "selectClustersAvgIden") == 0)
	{
		float minIden = 0;
		while(ptr != NULL)
		{
			ptr = strtok (NULL, " ");
			minIden = atof(ptr);

			if(minIden != 0)
				break;
		}

		selectClustersAvgIden(workDir, graphLevel, minIden);
	}

	if(strcmp(ptr, "selectClustersAvgOvlLen") == 0)
	{
		int minOvlLen = 0;
		while(ptr != NULL)
		{
			ptr = strtok (NULL, " ");
			minOvlLen = atoi(ptr);

			if(minOvlLen != 0)
				break;
		}
		selectClustersAvgOvlLen(workDir, graphLevel, minOvlLen);
	}

	if(strcmp(ptr, "displayClusterGCcontent") == 0)
		displayClusterGCcontent(workDir, graphLevel, currentNode);

	//if(strcmp(ptr, "displayClusterGenes") == 0)
		//displayClusterGenes(workDir, graphLevel, currentNode);

	if(strcmp(ptr, "displayClusterEdgeDen") == 0)
		displayClusterEdgeDen(workDir, graphLevel, currentNode);

	if(strcmp(ptr, "displayClusterSize") == 0)
		displayClusterSize(workDir, graphLevel, currentNode);

	if(strcmp(ptr, "displayClusterNodeDegree") == 0)
		displayClusterNodeDegree(workDir, graphLevel, currentNode);

	if(strcmp(ptr, "displayClusterNeighbors") == 0)
		displayClusterNeighbors(workDir, graphLevel, currentNode);

	if(strcmp(ptr, "displayClusterAvgIden") == 0)
		displayClusterAvgIden(workDir, graphLevel, currentNode);

	if(strcmp(ptr, "displayClusterAvgOvlLen") == 0)
		displayClusterAvgOvlLen(workDir, graphLevel, currentNode);

	if(strcmp(ptr, "displayClusterChildren") == 0)
		displayClusterChildren(workDir, graphLevel, currentNode);

	if(strcmp(ptr, "displayClusterAvgChildDegree") == 0)
		displayClusterAvgChildDegree(workDir, graphLevel, currentNode);

	if(strcmp(ptr, "displayCurrentGraph") == 0)
		cout<<graphLevel<<endl;

	if(strcmp(ptr, "displayCurrentNode") == 0)
		cout<<currentNode<<endl;

	}

	delete [] cmd;

	return 0;
}

//void countFiles(const char [])
//Description: This function counts files in a directory
//Input: dir (cosnt char []) : The number of tasks that need to be finalized
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

//bool fileExists( char fileName [] )
//Description: This function checks to see if a file exists
//Input: fileName (char []) the name of the file
//Output:None
//Return:None
bool fileExists(const char fileName []) {

	FILE * pFile = fopen(fileName, "r");
	if(pFile != '\0')
	{
		fclose(pFile);
		return true;
	}
	return false;
}

//avgGraphDegree(string &, const int)
//Description: This function displays the average node degree of the
//current graph selected.
//Input: workDir (string &): The current work directory, graphLevel: The current graph level
//Output:None
//Return: None
void avgGraphDegree(string & workDir, const int graphLevel)
{
	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Edges/Graph%d", workDir.c_str(), graphLevel);

	int numFiles = countFiles(graphDir);

	long long int numNodes = 0; double numEdges = 0;

	for(int i = 0; i < numFiles; i++)
	{
		char * fileName = new char [1000];
		sprintf(fileName, "%s/graph%d", graphDir, i);

		if(fileExists(fileName))
		{
			long long int intervalNodes = 0; long long int intervalEdges = 0;
			FILE * pFile = fopen(fileName, "r");
			fread(&intervalNodes, sizeof(long long int), 1, pFile);
			fread(&intervalEdges, sizeof(long long int), 1, pFile);
			fclose(pFile);

			numNodes+=intervalNodes; numEdges+=intervalEdges;
		}

		delete [] fileName;
	}

	delete [] graphDir;

	float avgDegree = (2*numEdges)/numNodes;
	cout<<avgDegree<<endl;
}

//avgGraphOvlLen(string &, const int)
//Description: This function displays the average ovl length of a graph
//Input: workDir (string &): The current work directory
//Output:None
//Return: None
void avgGraphOvlLen(string & workDir)
{
	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Edges/Graph%d", workDir.c_str(), 0);

	int numFiles = countFiles(graphDir);
	double numNodes = 0; double totalOvlLen = 0;

	for(int i = 0; i < numFiles; i++)
	{
		char * fileName = new char [1000];
		sprintf(fileName, "%s/graph%d", graphDir, i);

		if(fileExists(fileName))
		{
			FILE * pFile = fopen(fileName, "r");
			Graph graph;

			graph.read(pFile);
			fclose(pFile);

			numNodes+=graph.getNumNodes();

			for(int j = 0; j < graph.getNumEdges(); j++)
			{
				totalOvlLen+=graph.edgeOvlLen(j);
			}
		}
		delete [] fileName;
	}

	delete [] graphDir;

	float avgOvlLen = totalOvlLen/numNodes;
	cout<<avgOvlLen<<endl;
}

//avgGraphOvlLen(string &, const int)
//Description: This function displays the average identity of the edges in a graph
//Input: workDir (string &): The current work directory
//Output:None
//Return: None
void avgGraphOvlIden(string & workDir)
{
	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Edges/Graph%d", workDir.c_str(), 0);

	int numFiles = countFiles(graphDir);
	double numNodes = 0; double totalOvlIden = 0;

	for(int i = 0; i < numFiles; i++)
	{
		char * fileName = new char [1000];
		sprintf(fileName, "%s/graph%d", graphDir, i);

		if(fileExists(fileName))
		{
			FILE * pFile = fopen(fileName, "r");
			Graph graph;

			graph.read(pFile);
			fclose(pFile);

			numNodes+=graph.getNumNodes();

			for(int j = 0; j < graph.getNumEdges(); j++)
			{
				totalOvlIden+=graph.edgeOvlIden(j);
			}
		}

		delete [] fileName;
	}

	delete [] graphDir;

	float avgOvlIden = totalOvlIden/numNodes;

	cout<<avgOvlIden<<endl;
}

//avgGraphGCcontent(string &, const int)
//Description: This function displays the average GC content of the reads in a graph
//Input: workDir (string &): The current work directory
//Output:None
//Return: None
void avgGraphGCcontent(string & workDir)
{
	char * fragmentDir = new char [1000];
	sprintf(fragmentDir, "%s/Edges/Graph%d", workDir.c_str(), 0);

	int numFiles = countFiles(fragmentDir);
	double numNucleotides; double totalGCcounts = 0;

	int numSets = 0;
	char * fileName = new char [1000];
	sprintf(fileName, "%s/sequences.job0.set%d", fragmentDir, numSets++);

	while(fileExists(fileName))
	{
		sprintf(fileName, "%s/sequences.job0.set%d", fragmentDir, numSets++);
	}

	numSets-=1;
	int numJobs = numFiles/numSets;

	for(int i = 0; i < numJobs; i++)
	{
		for(int j = 0; j < numSets; j++)
		{
			sprintf(fileName, "%s/sequences.job%d.set%d", fragmentDir, i, j);
			FILE * pFile = fopen(fileName, "r");
			Fragment_Index index;
			index.read(pFile);
			fclose(pFile);

			for(unsigned int j = 0; j < index.length(); j++)
			{
				if(index.at(j) == 'G' || index.at(j) == 'C')
					totalGCcounts++;

				numNucleotides++;
			}
		}
	}

	float percentGC = totalGCcounts/numNucleotides;

	cout<<percentGC<<endl;

	delete [] fragmentDir; delete [] fileName;
}

//avgGraphEdgeDen (string &, const int)
//Description: This returns the average edge density of the clusters
//Input: workDir (string &): The current work directory, graphLevel (const int): the current graph
//Output:None
//Return: None
void avgGraphEdgeDen(string & workDir, const int graphLevel)
{
	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Edges/Graph%d", workDir.c_str(), 0);

	char * fileName = new char [1000];
	sprintf(fileName, "%s/eDen", graphDir);

	int * eDen = '\0'; int * nDen = '\0'; int numNodes = 0;

	if(fileExists(fileName))
	{
		FILE * pFile = fopen(fileName, "r");
		fseek(pFile, 0, SEEK_END);
		numNodes = ftell(pFile)/sizeof(int);
		fseek(pFile, 0,  SEEK_SET);
		eDen = new int [numNodes];
		fread(eDen, sizeof(int), numNodes, pFile);
		fclose(pFile);

	}else
	{
		cout<<"No density files found"<<endl;
		return;
	}

	sprintf(fileName, "%s/nDen", graphDir);
	if(fileExists(fileName))
	{
		FILE * pFile = fopen(fileName, "r");
		nDen = new int [numNodes];
		fread(nDen, sizeof(int), numNodes, pFile);
		fclose(pFile);

	}else
	{
		cout<<"No density files found"<<endl;
		delete [] eDen;
		return;
	}

	double eDenTotal = 0;
	for(int i = 0; i < numNodes; i++)
	{
		float density = (2*eDen[i])/(nDen[i] * (nDen[i]-1));
		eDenTotal+=density;
	}

	delete [] eDen; delete [] nDen; delete [] graphDir; delete [] fileName;
	float avgDensity = eDenTotal/numNodes;

	cout<<avgDensity<<endl;
}

//avgGraphClusterSize(string &, const int)
//Description: This function displays the average cluster size of a graph
//Input: workDir (string &): The current work directory, graphLevel (const int): The current graph level
//Output:None
//Return: None
void avgGraphClusterSize(string & workDir, const int graphLevel) {

	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Edges/Graph%d", workDir.c_str(), 0);

	char * fileName = new char [1000];
	sprintf(fileName, "%s/nDen", graphDir);

	int * nDen = '\0'; int numNodes = 0;
	if(fileExists(fileName))
	{
		FILE * pFile = fopen(fileName, "r");
		fseek(pFile, 0, SEEK_END);
		numNodes = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);
		nDen = new int [numNodes];
		fread(nDen, sizeof(int), numNodes, pFile);
		fclose(pFile);

	}else
	{
		cout<<"No density files found"<<endl;
		return;
	}

	double nWeightsTotal = 0;
	for(int i = 0; i < numNodes; i++)
	{
		nWeightsTotal+=nDen[i];
	}

	delete [] nDen; delete [] graphDir; delete [] fileName;
}

//selectClustersGCcontent(string &, const int)
//Description: This function selects clusters by their GC content
//Input: workDir (string &): The current work directory, graphLevel (const int):
//Output:None
//Return: None
void selectClustersGCcontent(string &  workDir, const int graphLevel, const float minGC, const bool greaterThan)
{
	char * fragmentDir = new char [1000];
	sprintf(fragmentDir, "%s/Fragments", workDir.c_str());

	int numFiles = countFiles(fragmentDir);

	int numSets = 0;
	char * fileName = new char [1000];
	sprintf(fileName, "%s/sequences.job0.set%d", fragmentDir, numSets++);

	while(fileExists(fileName))
	{
		sprintf(fileName, "%s/sequences.job0.set%d", fragmentDir, numSets++);
	}

	numSets-=1;
	int numJobs = numFiles/numSets;

	Fragment_Index index;
	for(int i = 0; i < numJobs; i++)
	{
		for(int j = 0; j < numSets; j++)
		{
			sprintf(fileName, "%s/sequences.job%d.set%d", fragmentDir, i, j);
			FILE * pFile = fopen(fileName, "r");
			index.read(pFile);
			fclose(pFile);
		}
	}

	int size = 2 * index.numFragments();

	int * contigs = new int [size];
	int count = recoverClusters(workDir, contigs, size, graphLevel);

	int currentCluster = 0;
	float GCcount = 0; float length = 0;

	cout<<"Cluster, GC content"<<endl;
	for(int i = 0; i < count; i++)
	{
		if(contigs[i] != -1)
		{
			pair<long long int, long long int> bounds;
			bounds = index.indexBounds(contigs[i]);
			length+=bounds.second-bounds.first;
			for(int j = bounds.first; j < bounds.second; j++)
			{
				if(index.at(j) == 'G' || index.at(j) == 'C')
					GCcount++;
			}

		}else{

			float GCcontent = GCcount/length;
			if((GCcontent >= minGC && greaterThan) || (GCcontent < minGC && !greaterThan))
			{
				cout<<currentCluster<<","<<GCcontent<<endl;
			}
			currentCluster++; length = 0; GCcount = 0;
		}
	}

	if(GCcount != 0 && length != 0)
	{
		float GCcontent = GCcount/length;
		if((GCcontent >= minGC && greaterThan) || (GCcontent < minGC && !greaterThan))
		{
			cout<<currentCluster<<","<<GCcontent<<endl;
		}
	}
	delete [] contigs; delete [] fileName; delete [] fragmentDir;
}

//void selectClustersGenes(const int graphLevel);

//selectClustersEdgeDen(string &, const int, const float)
//Description: This function selects nodes based on their edge density
//Input: workDir (string &): The working directory, graphLevel (const int): The current graph we are on,
//minDensity (const float): The minimum density for nodes to be returned.
//Output:None
//Return: None
void selectClustersEdgeDen(string & workDir, const int graphLevel, const float minDensity, const bool greaterThan)
{
	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Edges/Graph%d", workDir.c_str(), graphLevel);

	char * fileName = new char [1000];
	sprintf(fileName, "%s/eDen", graphDir);

	int * eDen = '\0'; int * nDen = '\0'; int numNodes = 0;

	if(fileExists(fileName))
	{
		FILE * pFile = fopen(fileName, "r");
		fseek(pFile, 0, SEEK_END);
		numNodes = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);
		eDen = new int [numNodes];
		fread(eDen, sizeof(int), numNodes, pFile);
		fclose(pFile);

	}else
	{
		cout<<"No density files found"<<endl;
		return;
	}

	sprintf(fileName, "%s/nDen", graphDir);
	if(fileExists(fileName))
	{
		FILE * pFile = fopen(fileName, "r");
		nDen = new int [numNodes];
		fread(nDen, sizeof(int), numNodes, pFile);
		fclose(pFile);

	}else
	{
		cout<<"No density files found"<<endl;
		delete [] eDen;
		return;
	}

	cout<<"Cluster, Edge Density"<<endl;
	for(int i = 0; i < numNodes; i++)
	{
		float density = (2*eDen[i])/(nDen[i] * (nDen[i]-1));

		if((density >= minDensity && greaterThan) || (density < minDensity && !greaterThan))
			cout<<i<<","<<density<<endl;
	}

	delete [] eDen; delete [] nDen; delete graphDir; delete [] fileName;
}

//selectGraphClusterSize(string &, const int, const float)
//Description: This function selects nodes based on their edge density
//Input: workDir (string &): The working directory, graphLevel (const int): The current graph we are on,
//minWeight (const float): The minimum weight (number of child nodes) for a parent node
//Output:None
//Return: None
void selectClustersSize(string & workDir, const int graphLevel, const float minWeight, const bool greaterThan)
{
	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Edges/Graph%d", workDir.c_str(), 0);

	char * fileName = new char [1000];
	sprintf(fileName, "%s/nDen", graphDir);

	int * nDen = '\0'; int numNodes = 0;
	if(fileExists(fileName))
	{
		FILE * pFile = fopen(fileName, "r");
		fseek(pFile, 0, SEEK_END);
		numNodes = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);
		nDen = new int [numNodes];
		fread(nDen, sizeof(int), numNodes, pFile);
		fclose(pFile);

	}else
	{
		cout<<"No density files found"<<endl;
		return;
	}

	cout<<"Cluster, Weight"<<endl;
	for(int i = 0; i < numNodes; i++)
	{
		if((nDen[i] >= minWeight && greaterThan) || (nDen[i] < minWeight & !greaterThan))
			cout<<i<<","<<nDen[i]<<endl;
	}

	delete [] nDen; delete [] graphDir; delete [] fileName;
}

//selectGraphAvgChildDegree(string &, const int, const float)
//Description: This function selects nodes their children nodes average node degree in the graph.
//Input: workDir (string &): The working directory, graphLevel (const int): The current graph we are on,
//minDegree (const float): The minimum average child degree that a parent node can have
//Output:None
//Return: None
void selectClustersAvgChildDegree(string & workDir, const int graphLevel, const int minDegree, const bool greaterThan)
{
	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Edges/Graph0", workDir.c_str());

	int numFiles = countFiles(graphDir) - 2;
	long int numNodes = 0;
	{
		Graph graph;
		char * fileName = new char [1000];
		sprintf(fileName, "%s/graph%d", graphDir, numFiles-1);
		FILE * pFile = fopen(fileName, "r");
		delete [] fileName;

		graph.read(pFile);
		numNodes = graph.getOffset() + graph.getNumNodes();
	}

	int size = 2 * numNodes;

	int * contigs = new int [size];
	int count = recoverClusters(workDir, contigs, size, graphLevel);

	int * nodeMap = new int [numNodes];
	int numClusters = 0;

	for(int i = 0; i < count; i++)
	{
		if(contigs[i] != -1)
		{
			nodeMap[contigs[i]] = numClusters;
		}else
		{
			numClusters++;
		}
	}

	delete [] contigs;

	float * degreeTotals = new float [numClusters];
	for(int i = 0; i < numClusters; i++)
		degreeTotals[i] = 0;

	for(int i = 0; i < numFiles; i++)
	{
		char * fileName = new char [1000];
		sprintf(fileName, "%s/graph%d", graphDir, i);
		FILE * pFile = fopen(fileName, "r");

		Graph graph;
		graph.read(pFile);
		fclose(pFile);

		graph.set();

		for(int j = 0; j < graph.getNumNodes(); j++)
		{
			degreeTotals[nodeMap[j+graph.getOffset()]]+=graph.nodeDegree(j);
		}

		delete [] fileName;
	}

	delete [] nodeMap;

	sprintf(graphDir, "%s/Edges/Graph%d", workDir.c_str(), graphLevel);

	char * fileName = new char [1000];
	sprintf(fileName, "%s/nDen", graphDir);

	int * nDen = '\0';
	if(fileExists(fileName))
	{
		FILE * pFile = fopen(fileName, "r");
		fseek(pFile, 0, SEEK_END);
		numNodes = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);
		nDen = new int [numNodes];
		fread(nDen, sizeof(int), numNodes, pFile);
		fclose(pFile);

	}else
	{
		cout<<"No density files found"<<endl;
		return;
	}

	cout<<"Cluster, Average Child Node Degree"<<endl;
	for(int i = 0; i < numNodes; i++)
	{
		float avgDegree = degreeTotals[i]/nDen[i];
		if((avgDegree >= minDegree && greaterThan) || (avgDegree < minDegree && !greaterThan))
		{
			cout<<i<<","<<avgDegree<<endl;
		}
	}

	delete [] degreeTotals; delete [] nDen; delete [] fileName; delete [] graphDir;
}

//selectClustersNodeDegree(string &, const int, const ing)
//Description: This function selects nodes that have a node degree greater than a minimum
//Input: workDir (string &): The working directory, graphLevel (const int): The current graph we are on,
//minDegree (const int): The minimum degree that a node can have
//Output:None
//Return: None
void selectClustersNodeDegree(string & workDir, const int graphLevel, const int minDegree, const bool greaterThan)
{
	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Edges/Graph%d", workDir.c_str(), graphLevel);

	int numFiles = countFiles(graphDir) - 2;

	for(int i = 0; i < numFiles; i++)
	{
		char * fileName = new char [1000];
		sprintf(fileName, "%s/graph%d", graphDir, i);
		FILE * pFile = fopen(fileName, "r");

		Graph graph;
		graph.read(pFile);
		fclose(pFile);

		graph.set();

		cout<<"Cluster, Degree"<<endl;
		for(int j = 0; j < graph.getNumNodes(); j++)
		{
			if((graph.nodeDegree(j) >= minDegree && greaterThan) || (graph.nodeDegree(j) < minDegree && !greaterThan))
				cout<<graph.getOffset() + j<<","<<graph.nodeDegree(j)<<endl;
		}

		delete [] fileName;
	}

	delete [] graphDir;
}

//selectClustersAvgIden(string &, const int, const float)
//Description: This function selects nodes that have a node degree greater than a minimum
//Input: workDir (string &): The working directory, graphLevel (const int): The current graph we are on,
//minDegree (const float): The minimum degree that a node can have
//Output:None
//Return: None
void selectClustersAvgIden(string & workDir, const int graphLevel, const float minIden, const bool greaterThan)
{
	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Edges/Graph0", workDir.c_str());

	int numFiles = countFiles(graphDir) - 2;

	long int numNodes = 0;
	{
		Graph graph;
		char * fileName = new char [1000];
		sprintf(fileName, "%s/graph%d", graphDir, numFiles-1);
		FILE * pFile = fopen(fileName, "r");
		delete [] fileName;

		graph.read(pFile);
		numNodes = graph.getOffset() + graph.getNumNodes();
	}

	int size = 2 * numNodes;

	int * contigs = new int [size];
	int count = recoverClusters(workDir, contigs, size, graphLevel);

	int * nodeMap = new int [numNodes];
	int numClusters = 0;

	for(int i = 0; i < count; i++)
	{
		if(contigs[i] != -1)
		{
			nodeMap[contigs[i]] = numClusters;
		}else
		{
			numClusters++;
		}
	}

	delete [] contigs;

	float * idenTotals = new float [numClusters];
	for(int i = 0; i < numClusters; i++)
		idenTotals[i] = 0;

	for(int i = 0; i < numFiles; i++)
	{
		char * fileName = new char [1000];
		sprintf(fileName, "%s/graph%d", graphDir, i);
		FILE * pFile = fopen(fileName, "r");

		Graph graph;
		graph.read(pFile);
		fclose(pFile);

		graph.set();

		for(int j = 0; j < graph.getNumNodes(); j++)
		{
			for(int k = 0; k < graph.nodeDegree(j); k++)
				idenTotals[nodeMap[j+graph.getOffset()]] += graph.edgeOvlIden(j, k);
		}

		delete [] fileName;
	}

	delete [] nodeMap;

	sprintf(graphDir, "%s/Edges/Graph%d", workDir.c_str(), graphLevel);

	char * fileName = new char [1000];
	sprintf(fileName, "%s/eDen", graphDir);

	int * eDen = '\0';
	if(fileExists(fileName))
	{
		FILE * pFile = fopen(fileName, "r");
		fseek(pFile, 0, SEEK_END);
		numNodes = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);
		eDen = new int [numNodes];
		fread(eDen, sizeof(int), numNodes, pFile);
		fclose(pFile);

	}else
	{
		cout<<"No density files found"<<endl;
		return;
	}

	cout<<"Cluster, Average Overlap Identity"<<endl;
	for(int i = 0; i < numNodes; i++)
	{
		float avgIden = idenTotals[i]/eDen[i];
		if((avgIden >= minIden && greaterThan) || (avgIden < minIden && !greaterThan))
		{
			cout<<i<<","<<avgIden<<endl;
		}
	}

	delete [] idenTotals; delete [] eDen; delete [] fileName; delete [] graphDir;
}

//selectClustersAvgOvlLen(string &, const int, const ing)
//Description: This function selects nodes whose child nodes have an average overlap length greater
//than a provided minimum
//Input: workDir (string &): The working directory, graphLevel (const int): The current graph we are on,
//minOvlLen (const float): The minimum average overlap length of the child node edges
//Output:None
//Return: None
void selectClustersAvgOvlLen(string & workDir, const int graphLevel, const int minOvlLen, const bool greaterThan)
{
	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Edges/Graph0", workDir.c_str());

	int numFiles = countFiles(graphDir) - 2;

	long int numNodes = 0;
	{
		Graph graph;
		char * fileName = new char [1000];
		sprintf(fileName, "%s/graph%d", graphDir, numFiles-1);
		FILE * pFile = fopen(fileName, "r");
		delete [] fileName;

		graph.read(pFile);
		numNodes = graph.getOffset() + graph.getNumNodes();
	}

	int size = 2 * numNodes;

	int * contigs = new int [size];
	int count = recoverClusters(workDir, contigs, size, graphLevel);

	int * nodeMap = new int [numNodes];
	int numClusters = 0;

	for(int i = 0; i < count; i++)
	{
		if(contigs[i] != -1)
		{
			nodeMap[contigs[i]] = numClusters;
		}else
		{
			numClusters++;
		}
	}

	delete [] contigs;

	float * ovlLenTotals = new float [numClusters];
	for(int i = 0; i < numClusters; i++)
		ovlLenTotals[i] = 0;

	for(int i = 0; i < numFiles; i++)
	{
		char * fileName = new char [1000];
		sprintf(fileName, "%s/graph%d", graphDir, i);
		FILE * pFile = fopen(fileName, "r");

		Graph graph;
		graph.read(pFile);
		fclose(pFile);

		graph.set();

		for(int j = 0; j < graph.getNumNodes(); j++)
		{
			for(int k = 0; k < graph.nodeDegree(j); k++)
				ovlLenTotals[nodeMap[j+graph.getOffset()]] += graph.edgeOvlLen(j, k);
		}

		delete [] fileName;
	}

	delete [] nodeMap;
	sprintf(graphDir, "%s/Edges/Graph%d", workDir.c_str(), graphLevel);

	char * fileName = new char [1000];
	sprintf(fileName, "%s/eDen", graphDir);

	int * eDen = '\0';
	if(fileExists(fileName))
	{
		FILE * pFile = fopen(fileName, "r");
		fseek(pFile, 0, SEEK_END);
		numNodes = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);
		eDen = new int [numNodes];
		fread(eDen, sizeof(int), numNodes, pFile);
		fclose(pFile);

	}else
	{
		cout<<"No density files found"<<endl;
		return;
	}

	cout<<"Cluster, Average Overlap Identity"<<endl;
	for(int i = 0; i < numNodes; i++)
	{
		float avgOvlLen = ovlLenTotals[i]/eDen[i];
		if((avgOvlLen >= minOvlLen && greaterThan) || (avgOvlLen < minOvlLen && !greaterThan))
		{
			cout<<i<<","<<avgOvlLen<<endl;
		}
	}

	delete [] ovlLenTotals; delete [] eDen; delete [] fileName; delete [] graphDir;
}

//void selectGraphCycles(const int graphLevel);

//selectClustersGCcontent(string &, const int, const float)
//Description: This function selects nodes that have a node degree greater than a minimum
//Input: workDir (string &): The working directory, graphLevel (const int): The current graph we are on,
//minDegree (const float): The minimum degree that a node can have
//Output:None
//Return: None
void displayClusterGCcontent(string & workDir, const int graphLevel, const int currentNode)
{
	char * fragmentDir = new char [1000];
	sprintf(fragmentDir, "%s/Fragments", workDir.c_str());

	int numFiles = countFiles(fragmentDir);

	int numSets = 0;
	char * fileName = new char [1000];
	sprintf(fileName, "%s/sequences.job0.set%d", fragmentDir, numSets++);

	while(fileExists(fileName))
	{
		sprintf(fileName, "%s/sequences.job0.set%d", fragmentDir, numSets++);
	}

	numSets-=1;
	int numJobs = numFiles/numSets;

	Fragment_Index index;
	for(int i = 0; i < numJobs; i++)
	{
		for(int j = 0; j < numSets; j++)
		{
			sprintf(fileName, "%s/sequences.job%d.set%d", fragmentDir, i, j);
			FILE * pFile = fopen(fileName, "r");
			index.read(pFile);
			fclose(pFile);
		}
	}

		int size = 2 * index.numFragments();

		int * contigs = new int [size];
		int count = recoverClusters(workDir, contigs, size, graphLevel);

		int currentCluster = 0;
		float GCcount = 0; float length = 0;

		for(int i = 0; i < count; i++)
		{
			if(contigs[i] != -1)
			{
				if(currentCluster == currentNode)
				{
					pair<long long int, long long int> bounds;
					bounds = index.indexBounds(contigs[i]);
					length+=bounds.second-bounds.first;
					for(int j = bounds.first; j < bounds.second; j++)
					{
						if(index.at(j) == 'G' || index.at(j) == 'C')
							GCcount++;
					}
				}
			}else{

				if(currentCluster == currentNode)
				{
					float GCcontent = GCcount/length;
					cout<<GCcontent<<endl;

					break;
				}
				currentCluster++;
			}
		}

		delete [] contigs; delete [] fileName; delete [] fragmentDir;
}

//void displayClusterEdgeDen(string &, const int, const int)
//Description: This displays the edge density of a cluster
//Input: workDir (string &): the working directory, graphLevel (const int): the
//current graph, currentNode (const int): The current node we are on
//Output:None
//Return:None
void displayClusterEdgeDen(string & workDir, const int graphLevel, const int currentNode)
{
	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Edges/Graph%d", workDir.c_str(), graphLevel);

	char * fileName = new char [1000];
	sprintf(fileName, "%s/eDen", graphDir);

	int * eDen = '\0'; int * nDen = '\0'; int numNodes = 0;
	if(fileExists(fileName))
	{
		FILE * pFile = fopen(fileName, "r");
		fseek(pFile, 0, SEEK_END);
		numNodes = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);
		eDen = new int [numNodes];
		fread(eDen, sizeof(int), numNodes, pFile);
		fclose(pFile);

	}else
	{
		cout<<"No density files found"<<endl;
		return;
	}

	sprintf(fileName, "%s/nDen", graphDir);
	if(fileExists(fileName))
	{
		FILE * pFile = fopen(fileName, "r");
		nDen = new int [numNodes];
		fread(nDen, sizeof(int), numNodes, pFile);
		fclose(pFile);
	}else
	{
		cout<<"No density files found"<<endl;
		delete [] eDen;
		return;
	}

	cout<<"Cluster, Edge Density"<<endl;
	float density = (2*eDen[currentNode])/(nDen[currentNode] * (nDen[currentNode]-1));

	cout<<density<<endl;

	delete [] eDen; delete [] nDen; delete graphDir; delete [] fileName;
}

//void displayClusterSize(string &, const int, const int)
//Description: This displays the number of child nodes in a cluster
//Input: workDir (string &): the working directory, graphLevel (const int): the
//current graph, currentNode (const int): The current node we are on
//Output:None
//Return:None
void displayClusterSize(string & workDir, const int graphLevel, const int currentNode)
{
	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Edges/Graph%d", workDir.c_str(), 0);

	char * fileName = new char [1000];
	sprintf(fileName, "%s/nDen", graphDir);

	int * nDen = '\0'; int numNodes = 0;
	if(fileExists(fileName))
	{
		FILE * pFile = fopen(fileName, "r");
		fseek(pFile, 0, SEEK_END);
		numNodes = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);
		nDen = new int [numNodes];
		fread(nDen, sizeof(int), numNodes, pFile);
		fclose(pFile);

	}else
	{
		cout<<"No density files found"<<endl;
		return;
	}

	cout<<nDen[currentNode]<<endl;

	delete [] nDen; delete [] graphDir; delete [] fileName;
}

//void displayClusterNodeDegree(string &, const int, const int)
//Description: This displays a cluster's child node's average degree
//Input: workDir (string &): the working directory, graphLevel (const int): the
//current graph, currentNode (const int): The current node we are on
//Output:None
//Return:None
void displayClusterNodeDegree(string & workDir, const int graphLevel, const int currentNode)
{
	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Edges/Graph%d", workDir.c_str(), graphLevel);

	int numFiles = countFiles(graphDir) - 2;

	for(int i = 0; i < numFiles; i++)
	{
		char * fileName = new char [1000];
		sprintf(fileName, "%s/graph%d", graphDir, i);
		FILE * pFile = fopen(fileName, "r");

		long long int numNodes;
		for(int i = 0; i < 2; i++)
		fread(&numNodes, sizeof(long long int), 1, pFile);

		long long int offset;
		for(int i = 0; i < 4; i++)
			fread(&offset, sizeof(long long int), 1, pFile);

		if(currentNode >= offset && currentNode < offset+numNodes)
		{
			fseek(pFile, 0, SEEK_SET);

			Graph graph;
			graph.read(pFile);
			fclose(pFile);

			graph.set();
			cout<<graph.nodeDegree(currentNode-offset)<<endl;
		}

		delete [] fileName;
	}

	delete [] graphDir;
}

//void displayClusterAvgIden(string &, const int, const int)
//Description: This displays a node's average identity in a graph
//Input: workDir (string &): the working directory, graphLevel (const int): the
//current graph, currentNode (const int): The current node we are on
//Output:None
//Return:None
void displayClusterAvgChildDegree(string & workDir, const int graphLevel, const int currentNode)
{
	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Edges/Graph0", workDir.c_str());

	int numFiles = countFiles(graphDir) - 2;

	long int numNodes = 0;
	{
		Graph graph;
		char * fileName = new char [1000];
		sprintf(fileName, "%s/graph%d", graphDir, numFiles-1);
		FILE * pFile = fopen(fileName, "r");
		delete [] fileName;

		graph.read(pFile);
		numNodes = graph.getOffset() + graph.getNumNodes();
	}

	int size = 2 * numNodes;

	int * contigs = new int [size];
	int count = recoverClusters(workDir, contigs, size, graphLevel);

	int * nodeMap = new int [numNodes];
	int numClusters = 0;

	for(int i = 0; i < count; i++)
	{
		if(contigs[i] != -1)
		{
			nodeMap[contigs[i]] = numClusters;
		}else
		{
			numClusters++;
		}
	}

	delete [] contigs;

	float degreeTotal = 0;

	for(int i = 0; i < numFiles; i++)
	{
		char * fileName = new char [1000];
		sprintf(fileName, "%s/graph%d", graphDir, i);
		FILE * pFile = fopen(fileName, "r");

		Graph graph;
		graph.read(pFile);
		fclose(pFile);

		graph.set();

		for(int j = 0; j < graph.getNumNodes(); j++)
		{
			if(nodeMap[j+graph.getOffset()] == currentNode)
				degreeTotal+=graph.nodeDegree(j);
		}

		delete [] fileName;
	}

	delete [] nodeMap;

	sprintf(graphDir, "%s/Edges/Graph%d", workDir.c_str(), graphLevel);

	char * fileName = new char [1000];
	sprintf(fileName, "%s/nDen", graphDir);

	int * nDen = '\0';
	if(fileExists(fileName))
	{
		FILE * pFile = fopen(fileName, "r");
		fseek(pFile, 0, SEEK_END);
		numNodes = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);
		nDen = new int [numNodes];
		fread(nDen, sizeof(int), numNodes, pFile);
		fclose(pFile);

	}else
	{
		cout<<"No density files found"<<endl;
		return;
	}


	float avgDegree = degreeTotal/nDen[currentNode];
	cout<<avgDegree<<endl;

	delete [] nDen; delete [] fileName; delete [] graphDir;
}

//void displayClusterAvgIden(string &, const int, const int)
//Description: This displays a node's average identity in a graph
//Input: workDir (string &): the working directory, graphLevel (const int): the
//current graph, currentNode (const int): The current node we are on
//Output:None
//Return:None
void displayClusterAvgIden(string & workDir, const int graphLevel, const int currentNode)
{
	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Edges/Graph0", workDir.c_str());

	int numFiles = countFiles(graphDir) - 2;

	long int numNodes = 0;
	{
		Graph graph;
		char * fileName = new char [1000];
		sprintf(fileName, "%s/graph%d", graphDir, numFiles-1);
		FILE * pFile = fopen(fileName, "r");
		delete [] fileName;

		graph.read(pFile);
		numNodes = graph.getOffset() + graph.getNumNodes();
	}

	int size = 2 * numNodes;

	int * contigs = new int [size];
	int count = recoverClusters(workDir, contigs, size, graphLevel);

	int * nodeMap = new int [numNodes];
	int numClusters = 0;

	for(int i = 0; i < count; i++)
	{
		if(contigs[i] != -1)
		{
			nodeMap[contigs[i]] = numClusters;
		}else
		{
			numClusters++;
		}
	}

	delete [] contigs;

	float idenTotal = 0;

	for(int i = 0; i < numFiles; i++)
	{
		char * fileName = new char [1000];
		sprintf(fileName, "%s/graph%d", graphDir, i);
		FILE * pFile = fopen(fileName, "r");

		Graph graph;
		graph.read(pFile);
		fclose(pFile);

		graph.set();

		for(int j = 0; j < graph.getNumNodes(); j++)
		{
			if(j + graph.getOffset() == currentNode)
				for(int k = 0; k < graph.nodeDegree(j); k++)
					idenTotal+= graph.edgeOvlIden(j, k);
		}

		delete [] fileName;
	}

	delete [] nodeMap;

	sprintf(graphDir, "%s/Edges/Graph%d", workDir.c_str(), graphLevel);

	char * fileName = new char [1000];
	sprintf(fileName, "%s/eDen", graphDir);

	int * eDen = '\0';
	if(fileExists(fileName))
	{

		FILE * pFile = fopen(fileName, "r");
		fseek(pFile, 0, SEEK_END);
		numNodes = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);
		eDen = new int [numNodes];
		fread(eDen, sizeof(int), numNodes, pFile);
		fclose(pFile);

	}else
	{
		cout<<"No density files found"<<endl;
		return;
	}

	cout<<"Cluster, Average Overlap Identity"<<endl;

	float avgIden = idenTotal/eDen[currentNode];
	cout<<avgIden<<endl;


	delete [] eDen; delete [] fileName; delete [] graphDir;
}

//void displayClusterAvgOvlLen(string &, const int, const int)
//Description: This displays a clusters average overlap lengths
//Input: workDir (string &): the working directory, graphLevel (const int): the
//current graph, currentNode (const int): The current node we are on
//Output:None
//Return:None
void displayClusterAvgOvlLen(string & workDir, const int graphLevel, const int currentNode)
{
	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Edges/Graph0", workDir.c_str());

	int numFiles = countFiles(graphDir) - 2;

	long int numNodes = 0;
	{
		Graph graph;
		char * fileName = new char [1000];
		sprintf(fileName, "%s/graph%d", graphDir, numFiles-1);
		FILE * pFile = fopen(fileName, "r");
		delete [] fileName;

		graph.read(pFile);
		numNodes = graph.getOffset() + graph.getNumNodes();
	}

	int size = 2 * numNodes;

	int * contigs = new int [size];
	int count = recoverClusters(workDir, contigs, size, graphLevel);

	int * nodeMap = new int [numNodes];
	int numClusters = 0;

	for(int i = 0; i < count; i++)
	{
		if(contigs[i] != -1)
		{
			nodeMap[contigs[i]] = numClusters;
		}else
		{
			numClusters++;
		}
	}

	delete [] contigs;

	float ovlLenTotal = 0;

	for(int i = 0; i < numFiles; i++)
	{
		char * fileName = new char [1000];
		sprintf(fileName, "%s/graph%d", graphDir, i);
		FILE * pFile = fopen(fileName, "r");

		Graph graph;
		graph.read(pFile);
		fclose(pFile);

		graph.set();

		for(int j = 0; j < graph.getNumNodes(); j++)
		{
			if(j+graph.getOffset() == currentNode)
				for(int k = 0; k < graph.nodeDegree(j); k++)
					ovlLenTotal += graph.edgeOvlLen(j, k);
		}

		delete [] fileName;
	}

	delete [] nodeMap;
	sprintf(graphDir, "%s/Edges/Graph%d", workDir.c_str(), graphLevel);

	char * fileName = new char [1000];
	sprintf(fileName, "%s/eDen", graphDir);

	int * eDen = '\0';
	if(fileExists(fileName))
	{
		FILE * pFile = fopen(fileName, "r");
		fseek(pFile, 0, SEEK_END);
		numNodes = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);
		eDen = new int [numNodes];
		fread(eDen, sizeof(int), numNodes, pFile);
		fclose(pFile);

	}else
	{
		cout<<"No density files found"<<endl;
		return;
	}

	float avgOvlLen = ovlLenTotal/eDen[currentNode];
	cout<<avgOvlLen<<endl;

	delete [] eDen; delete [] fileName; delete [] graphDir;
}

//void displayClusterNeighbors(string &, const int, const int)
//Description: This displays a node's neighbors in the current graph
//Input: workDir (string &): the working directory, graphLevel (const int): the
//current graph, currentNode (const int): The current node we are on
//Output:None
//Return:None
void displayClusterNeighbors(string & workDir, const int graphLevel, const int currentNode)
{
	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Edges/Graph%d", workDir.c_str(), graphLevel);

	int numFiles = countFiles(graphDir) - 2;

	for(int i = 0; i < numFiles; i++)
	{
		char * fileName = new char [1000];
		sprintf(fileName, "%s/graph%d", graphDir, i);
		FILE * pFile = fopen(fileName, "r");

		long long int numNodes;
		for(int i = 0; i < 2; i++)
		fread(&numNodes, sizeof(long long int), 1, pFile);

		long long int offset;
		for(int i = 0; i < 4; i++)
			fread(&offset, sizeof(long long int), 1, pFile);

		if(currentNode >= offset && currentNode < offset+numNodes)
		{
			fseek(pFile, 0, SEEK_SET);

			Graph graph;
			graph.read(pFile);
			fclose(pFile);

			graph.set();

			int j = graph.nodeDegree(currentNode-offset);
			for(int k = 0; k < j; k++)
				cout<<graph.edgeDest(currentNode-offset, k);
		}
		delete [] fileName;
	}
}

//void displayClusterChildren(string &, const int, const int)
//Description: This displays a node's immediate children
//Input: workDir (string &): the working directory, graphLevel (const int): the
//current graph, currentNode (const int): The current node we are on
//Output:None
//Return:None
void displayClusterChildren(string & workDir, const int graphLevel, const int currentNode)
{

	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Edges/Graph%d", workDir.c_str(), graphLevel);
	char * fileName = new char [1000];
	sprintf(fileName, "%s/map", graphDir);

	int size = 0; int * map = '\0';
	if(fileExists(fileName))
	{
		FILE * pFile = fopen(fileName, "r");
		fseek(pFile, 0, SEEK_END);
		size = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);
		map = new int [size];
		fread(map, sizeof(int), size, pFile);
		fclose(pFile);
	}else
	{
		cout<<"No density files found"<<endl;
		return;
	}

	cout<<map[currentNode * 2];
	if(map[currentNode * 2 + 1] != -1)
		cout<<" "<<map[currentNode * 2 + 1];
	cout<<endl;

	delete [] fileName; delete [] graphDir; delete [] map;
}

//void recoverClusters()
//Description: This runs the the parent thread in a parallel program
//Input: workDir (string &): the working directory, minMerge (int):
//the minimum overlap length for merging, minDensity (double): The minimum inter-node
//density
//Output:None
//Return:None
int recoverClusters(string & workDir, int * & contigs, int size, int graphLevel) {

	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Edges/Graph%d", workDir.c_str(), graphLevel);

	int count = 0;
	for(int i = graphLevel; i > 0; i--)
	{
		char * fileName = new char [1000];
		sprintf(fileName, "%s/map", graphDir);

		FILE * pFile = fopen(fileName, "r");
		fseek(pFile, 0, SEEK_END);

		int buffSize = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);

		int * buff = new int [buffSize];
		fread(buff, sizeof(int), buffSize, pFile);
		fclose(pFile);

		if(i == graphLevel)
		{
			for(int j = 0; j < buffSize; j+=2)
			{
				contigs[j] = j/2;
				contigs[j+1] = -1;
			}
			count = buffSize;
		}

		int * tmp = new int [size];

		int cHold = count; count = 0;
		for(int j = 0; j < cHold; j++)
		{
			if(contigs[j]  != -1)
			{
				tmp[count++] = buff[contigs[j]*2];

				if(buff[contigs[j]*2+1] != -1)
					tmp[count++] = buff[contigs[j]*2+1];
			}else{
				tmp[count++] = -1;
			}
		}

		delete [] contigs;

		contigs = tmp;

		delete [] buff;
		delete [] fileName;
	}

	delete [] graphDir;

	return count;
}
