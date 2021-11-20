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
using namespace std;
#include <set>
#include <stack>
#include <queue>
#include <dirent.h>
#include <map>
#include "mpi.h"
#include "MappingValues.h"
#include "Dictionary.h"
#include "Graph.h"
#include "Fragment_Index.h"
#include "Labels.h"
#include "SortEdges.h"
#include "MergeEdges.h"
#include "SortFragments.h"
#include "MergeFragments.h"

//MPI Tags
#define sortTagE 0
#define mergeTagE 1
#define sortTagF 2
#define mergeTagF 3
#define quitTag 4
#define graphTag 5

//The help message
void help_message(){

	cout<<"\n\n";
	cout<<"\t\t\t Merge and Traverse Focus Assembler - Coarsen Graph"<<endl;
	cout<<"\t\t\tWritten by Julia Warnke-Sommer\n"<<endl;
	cout<<"Command line arguments\n"<<endl;
	cout<<"--workDir :working directory :required"<<endl;
	cout<<"--minMerge : The minimum edge ovl length for merging"<<endl;
	cout<<"--minDensity : The minimum inter-cluster density for merging "<<endl;
	cout<<"--percentMerged : The minimum percentage of nodes that have been merged "<<endl;
	cout<<"for the algorithm to continue merging"<<endl;
	cout<<"--verbose : verbose mode"<<endl;
	cout<<"--oClusters : Output the clusters of fragments"<<endl;
	cout<<"--time : Report the runtime"<<endl;
	cout<<"\n\n";

	cout<<"\n\n"<<endl;

	cout<<"Exiting program"<<endl;
	MPI_Abort(MPI_COMM_WORLD, help);
}

//Function Declarations
void createGraphDir(string &, const int);
void moveEdges(string &);
void runParentTask(string &, int, double, int, float, bool);
void runChildTask(string &);
void sortEdges(string &, int, int, int, int);
void mergeEdges(string &, int, int, int);
void sortFragments(string &, int);
void mergeFragments(string &, int);
void quitTasks(const int);
void makeGraph(string &, int, int, int);
bool fileExists(const char []);
void finalizeTasks(const int);
int  countFiles(const char []);
unsigned long long int countFragments(string &);
void cleanTmps(const char[] , const char []);
void renameFiles(const char []);
void makeDir(const char []);
void moveDir(const char [], const char []);
void removeDir(const char []);
void calculateIntervals(string &, const int, int * &, const int, const int);
float findMatching(string &, int, int, float, int * &, int [], int, int);
void makeEdges(string &, int, int, int, int * &);
void makeFinalEdges(string &, int, int * &);
void initiateDensity(string &);
int recoverClusters(string &, int * &);
void scrambleNodes(int  * & , int * &, int  * &, int * &, int);
void createFinalMapping(string &, int * &, int, int * &);
void outputClusters(string &, int * &, int, int);

//Main Function
int main(int argc, char * argv [])
{
	//Set up the MPI
	int rank;
	int nTasks;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nTasks);

	//User Input Flags
	string workDirS = "--workDir"; //Current Working Directory
	string minMergeS = "--minMerge"; //Minimum Overlap for Merging
	string minDensityS = "--minDensity"; //Minimum Node Density for Merging
	string verboseS = "--verbose";
	string percentMergedS = "--percentMerged"; //Minimum percentage of nodes merged
	string oClustersS = "--oClusters"; //Output the clusters at each level
	string timeS = "--time"; //Report the time

	//Input values
	string workDir = "NO_WORK_DIR";
	int minMerge = 50;
	double minDensity = 50;
	bool verbose = false;
	float percentMerged = .01;
	bool oClusters = false;
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

			//Don't want to segfault
			if(workDir.size() > 200)
			{
				cout<<"Please keep the working directory ";
				cout<<"name less than 200 characters in ";
				cout<<"length"<<endl;
				exit(bufOverflow);
			}
		}

		//The minimum merging value
		if(argv[i] == minMergeS){
			minMerge = atoi(argv[i+1]);
		}

		//The minimum merging value
		if(argv[i] == minDensityS){
			minDensity = atof(argv[i+1]);
		}

		//The minimum percent of node merged to continue
		if(argv[i] == percentMergedS){
			percentMerged = atof(argv[i+1]);
		}

		//Verbose or not
		if(argv[i] == verboseS && verbose != true)
		{
			verbose = true;
			i--;
		}

		//Verbose or not
		if(argv[i] == oClustersS && oClusters != true)
		{
			oClusters = true;
			i--;
		}

		//Display time
		if(argv[i] == timeS && reportTime != true)
		{
			reportTime = true;
			i--;
		}
	}

	//User didn't provided the working directory
	if(workDir == "NO_WORK_DIR")
	{
		help_message();
		exit(help);
	}

	//Log file stuff
	if(rank == 0)
	{
		cerr<<"Graph Coarsening Algorithm"<<endl;
		cerr<<"Written by Julia Warnke"<<endl;
		cerr<<endl;
		cerr<<"Graph Coarsening Parameters"<<endl;
		cerr<<endl;
		cerr<<"Working directory : "<<workDir<<endl;
		cerr<<"Minimum merging value : "<<minMerge<<endl;
		cerr<<"Minimum node density : "<<minDensity<<endl;
		cerr<<"Verbose mode : ";

		if(verbose) 
		{ 
			cerr<<"yes"<<endl; 
		}else{
			cerr<<"no"<<endl;			
		}

		cerr<<"Minimum node percent merged : "<<percentMerged<<endl;
		cerr<<"Output clusters : ";

		if(oClusters) 
		{ 
			cerr<<"yes"<<endl; 
		}else{
			cerr<<"no"<<endl;			
		}

		cerr<<"Reporting time : ";
		if(reportTime) 
		{ 
			cerr<<"yes"<<endl; 
		}else{
			cerr<<"no"<<endl;			
		}
	}

	if(reportTime == true && rank == 0)
	{
		time_t rawtime;
		struct tm * timeinfo;
		time(&rawtime);
		timeinfo = localtime(&rawtime);
		cout<<"start time: "<<asctime(timeinfo)<<endl;
	}

	//Run the program
	if(rank == 0)
		runParentTask(workDir, minMerge, minDensity, nTasks, percentMerged, oClusters);
	else
		runChildTask(workDir);

	if(reportTime == true && rank == 0)
	{
		time_t rawtime;
		struct tm * timeinfo;
		time(&rawtime);
		timeinfo = localtime(&rawtime);
		cout<<"end time: "<<asctime(timeinfo)<<endl;
	}

	//Log file stuff
	if(rank == 0)
	{
		cerr<<"Job is complete "<<endl;
		cerr<<"Thank you for running the Focus algorithm"<<endl;
	}

	MPI_Finalize();

	return 0;
}

//void runParentTask()
//Description: This runs the the parent thread in a parallel program
//Input: workDir (string &): the working directory, minMerge (int):
//the minimum overlap length for merging, minDensity (double): The minimum inter-node
//outputClusters (bool): Output the clusters at each level if true
//density
//Output:None
//Return:None
void runParentTask(string & workDir, int minMerge, double minDensity, int nTasks, float minPercent, bool oClusters){

	int graphIter = 0;
	//Create first graph directory
	//Organize/move edges to it

	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Spectrum/Graph%d", workDir.c_str(), graphIter);
	mkdir(graphDir, S_IRWXU);
	delete [] graphDir;

	moveEdges(workDir); 

	//Create the sub-node and inter-edge density of each node
	//Collapse Graph/Merge Nodes
	float percentMerged = 0;

	do{
		int auxFiles = 3;

		if(graphIter == 0)
			auxFiles = 0;

		cout<<"SORTING "<<endl;
		sortEdges(workDir, 0, graphIter, auxFiles, nTasks);
		cout<<"DONE SORTING "<<endl;
		mergeEdges(workDir, graphIter, auxFiles, nTasks);
		cout<<"DONE MERGING "<<endl;

		if(graphIter == 0)
			sortEdges(workDir, 1, graphIter, auxFiles, nTasks);
		else
			sortEdges(workDir, 2, graphIter, auxFiles, nTasks);

		cout<<"MAKING GRAPH "<<endl;
		makeGraph(workDir, graphIter, auxFiles, nTasks);
		cout<<"DONE MAKING "<<endl;

		if(graphIter > 0 && oClusters)
		{
			int * contigs  = '\0';
			int count = recoverClusters(workDir, contigs);

			outputClusters(workDir, contigs, count, graphIter);
			delete [] contigs;
		}

		if(graphIter == 0)
		{

			initiateDensity(workDir); auxFiles = 2;
		}

	
		char * graphDir = new char [1000];
        	sprintf(graphDir, "%s/Spectrum/Graph%d", workDir.c_str(), graphIter+1);
        	mkdir(graphDir, S_IRWXU);
        	delete [] graphDir;

		int * intervals = new int [5]; int numIntervals = 4;
		calculateIntervals(workDir, graphIter, intervals, auxFiles, numIntervals);

		int * nodeMap = '\0';
		percentMerged = findMatching(workDir, graphIter, minMerge, minDensity, nodeMap, intervals, auxFiles, numIntervals);
		makeEdges(workDir, graphIter, graphIter+1, auxFiles, nodeMap);

		delete [] nodeMap;
		delete [] intervals;

		graphIter++;

	}while(percentMerged >= minPercent);

	int auxFiles = 3;
	sortEdges(workDir, 0, graphIter, auxFiles, nTasks);
	mergeEdges(workDir, graphIter, auxFiles, nTasks);

	if(graphIter == 0)
		sortEdges(workDir, 1, graphIter, auxFiles, nTasks);
	else
		sortEdges(workDir, 2, graphIter, auxFiles, nTasks);

	makeGraph(workDir, graphIter, auxFiles, nTasks);

	int * contigs  = '\0';
	int count = recoverClusters(workDir, contigs);

	if(oClusters)
		outputClusters(workDir, contigs, count, graphIter++);

	int * nodeMap = '\0';
	createFinalMapping(workDir, contigs, count, nodeMap);
	delete [] contigs;

	graphIter = -1; //Making the final OvlGraph

	auxFiles = 2;
	makeFinalEdges(workDir, auxFiles, nodeMap);
	delete [] nodeMap;

	auxFiles = 1;
	sortEdges(workDir, 0, graphIter, auxFiles, nTasks);
	mergeEdges(workDir, graphIter, auxFiles, nTasks);
	sortEdges(workDir, 1, graphIter, auxFiles, nTasks);

	makeGraph(workDir, graphIter, auxFiles, nTasks);

	sortFragments(workDir, nTasks);
	mergeFragments(workDir, nTasks);

	//Clean up
	for(int i = 1; i < nTasks; i++)
		quitTasks(i);

}

//void outputClusters(string &, int * &, int, int)
//Description: This outputs the fragments in a clustering to a file
//Input: workDir (string &):the working directory, contigs (int * &),
//count (int),
//Output:None
//Return:None
void outputClusters(string & workDir, int * & contigs, int count, int graphIter)
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

	sprintf(fileName, "%s/sequences.job%d.set%d", fragmentDir, numJobs-1, numSets-1);

	FILE * pFile  = fopen(fileName, "r");
	unsigned long long int fragmentOffset;
	fread(&fragmentOffset, sizeof(unsigned long long int), 1, pFile);

	unsigned long long int indexOffset;
	fread(&indexOffset, sizeof(unsigned long long int), 1, pFile);

	unsigned long long int intervalLength;
	fread(&intervalLength, sizeof(unsigned long long int), 1, pFile);

	unsigned long long int numIntervalFragments;
	fread(&numIntervalFragments, sizeof(unsigned long long int), 1, pFile);
	fclose(pFile);

	unsigned long long int numFragments = numIntervalFragments + indexOffset;
	unsigned long long int length = intervalLength + fragmentOffset;

	Fragment_Index index;
	index.reserve(numFragments, length, 0, 0);

	Labels label; int tagLen;
	sprintf(fileName, "%s/Labels/labels.job%d.set%d", workDir.c_str(), numJobs-1, numSets-1);
	pFile = fopen(fileName, "r");
	fread(&tagLen, sizeof(int), 1, pFile);
	fclose(pFile);

	label.reserve(tagLen, numFragments, 0);

	for(int i = 0; i < numJobs; i++)
	{
		for(int j = 0; j < numSets; j++)
		{
			sprintf(fileName, "%s/sequences.job%d.set%d", fragmentDir, i, j);
			pFile = fopen(fileName, "r");
			index.read(pFile);
			fclose(pFile);

			sprintf(fileName, "%s/Labels/labels.job%d.set%d", workDir.c_str(), i, j);
			pFile = fopen(fileName, "r");
			label.read(pFile);
			fclose(pFile);
		}
	}

	index.set();

	sprintf(fileName, "%s/Clusters/clusters%d", workDir.c_str(), graphIter);
	pFile = fopen(fileName, "a");

	char * buff = new char [10 * megabyte];
	int curr = 0;
	for(int i = 0; i < count; i++)
	{
		if(contigs[i] != -1)
		{
			string fLabel = label.labelAt(contigs[i]);
			pair<unsigned long long int, unsigned long long int> bounds = index.indexBounds(contigs[i]);
			int appLen = fLabel.length() + 2;

			if(curr+appLen >= 10 * megabyte-1)
			{
				fwrite(buff, sizeof(char), curr, pFile);
				curr = 0;
			}

			buff[curr++] = '>';
			for(unsigned int j = 0; j < fLabel.length(); j++)
			{
				buff[curr++] = fLabel.at(j);
			}
			buff[curr++] = '\n';


			for(unsigned int j = bounds.first; j < bounds.second; j++)
			{
				buff[curr++] = index.at(j);
			}

			buff[curr++] = '\n';

		}else{

			if(curr >= 10 * megabyte-1)
			{
				fwrite(buff, sizeof(char), curr, pFile);
				curr = 0;
			}

			buff[curr++] = '*'; buff[curr++] = '\n';
		}
	}

	if(curr > 0)
	{
		fwrite(buff, sizeof(char), curr, pFile);
	}

	fclose(pFile);

	delete [] buff;

	delete [] fragmentDir;
	delete [] fileName;
}


//void createFinalMapping(string &, int * &, int, int * &, int)
//Description: This remaps the nodes according to the clustering
//Input: workDir (string &): the working directory, minMerge (int):
//contigs (int * &), count (int), nodeMap (int * &), graphIter (int)
//Output:None
//Return:None
void createFinalMapping(string & workDir, int * & contigs, int count, int * & nodeMap) {

	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Spectrum/Graph0", workDir.c_str());

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

	nodeMap = new int [numNodes];

	int numBounds = 0;
	for(int i = 0; i < count; i++)
	{
		if(contigs[i] != -1)
		{
			nodeMap[contigs[i]] = i - numBounds;
		}else
		{
			numBounds++;
		}
	}

	char * mapFile = new char [200];
	sprintf(mapFile, "%s/OvlGraph/finalMap", workDir.c_str());

	FILE * pFile = fopen(mapFile, "w");
	fwrite(nodeMap, sizeof(int), numNodes, pFile);
	fclose(pFile);
}

//void recoverClusters()
//Description: This runs the the parent thread in a parallel program
//Input: workDir (string &): the working directory, minMerge (int):
//the minimum overlap length for merging, minDensity (double): The minimum inter-node
//density
//Output:None
//Return:None
int recoverClusters(string & workDir, int * & contigs) {

	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Spectrum/Graph0", workDir.c_str());

	int numFiles = countFiles(graphDir) - 2;
	long int numNodes = 0;

	{
		Graph graph;
		char * fileName = new char [1000];
		sprintf(fileName, "%s/graph%d", graphDir, numFiles-1);
		FILE * pFile = fopen(fileName, "r");
		delete [] fileName;

		graph.read(pFile);
		fclose(pFile);
		numNodes = graph.getOffset() + graph.getNumNodes();
	}

	sprintf(graphDir, "%s/Spectrum", workDir.c_str());
	numFiles = countFiles(graphDir);

	contigs = new int [2 * numNodes];

	int count = 0;
	for(int i = numFiles-1; i > 0; i--)
	{
		char * fileName = new char [1000];
		sprintf(fileName, "%s/Graph%d/map", graphDir, i);

		FILE * pFile = fopen(fileName, "r");
		fseek(pFile, 0, SEEK_END);

		int buffSize = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);

		int * buff = new int [buffSize];
		fread(buff, sizeof(int), buffSize, pFile);
		fclose(pFile);

		if(i == numFiles - 1)
		{
			for(int j = 0; j < buffSize; j+=2)
			{
				contigs[j] = j/2;
				contigs[j+1] = -1;
			}
			count = buffSize;
		}

		int * tmp = new int [2 * numNodes];

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

//void makeEdges(string &, int, int * &)
//Description: This creates new edges for the collapsed graph
//Input: workDir (string &): the working directory, graphIter (int): The current
//graph iteration, nodeMap (int *): the mapping from the previous graph
//to the new one
//Output:None
//Return:None
void makeFinalEdges(string & workDir, int auxFiles, int * & nodeMap) {

	char * iGraphDir = new char [1000];
	char * oGraphDir = new char [1000];
	char * fileName = new char [1000];

	sprintf(iGraphDir, "%s/Spectrum/Graph0", workDir.c_str());
	sprintf(oGraphDir, "%s/OvlGraph", workDir.c_str());

	int numFiles = countFiles(iGraphDir);
	numFiles-=auxFiles;

	uint32_t * buff = new uint32_t[33 * megabyte];
	int count = 0; int numSets = 0;

	long long int numEdges = 0;
	for(int i = 0; i < numFiles; i++)
	{
		sprintf(fileName, "%s/graph%d", iGraphDir, i);
		FILE * pFile = fopen(fileName, "r");

		Graph graph;
		graph.read(pFile);
		fclose(pFile);

		for(long int j = 0; j < graph.getNumEdges(); j++)
		{
			uint32_t source = nodeMap[graph.edgeSource(j) + graph.getOffset()];
			uint32_t dest = nodeMap[graph.edgeDest(j)];

			uint32_t ovlLen = graph.edgeOvlLen(j)/10;

			uint32_t ovlIden = (graph.edgeOvlIden(j) - 70)/.5;
			uint32_t reverse = 0; uint32_t sContained = 0; uint32_t dContained = 0;
			uint32_t rDovetail = 0; uint32_t matePair = 0;

			if(graph.isReverse(j)) 
				reverse = 1;

			if(graph.sContained(j))
			{
				sContained = 1;
			}

			if(graph.dContained(j))
			{
				dContained = 1;
			}

			if(graph.rDovetail(j))
				rDovetail = 1;

			if(graph.matePair(j))
				matePair = 1;	


			if(count >= 33 * megabyte)
			{
				sprintf(fileName, "%s/set%d", oGraphDir, numSets++);
				pFile = fopen(fileName, "w");
				fwrite(buff, sizeof(uint32_t), count, pFile);
				fclose(pFile); count = 0;
			}

			if(source != dest)
			{
				numEdges++;
				buff[count++] = source;
				buff[count++] = dest;
				buff[count++] = ((ovlLen << 22) | (ovlIden << 16) | (reverse << 15) | (sContained << 14) | (dContained << 13) | (rDovetail << 12));
			}
		}
	}

	if(count > 11 * megabyte || numSets == 0)
	{
		sprintf(fileName, "%s/set%d", oGraphDir, numSets++);
		FILE * pFile = fopen(fileName, "w");
		fwrite(buff, sizeof(uint32_t), count, pFile);
		fclose(pFile);

	}else{

		if(count > 0)
		{
			sprintf(fileName, "%s/set%d", oGraphDir, numSets-1);
			FILE * pFile = fopen(fileName, "a");
			fwrite(buff, sizeof(uint32_t), count, pFile);
			fclose(pFile);
		}
	}

	delete [] iGraphDir;
	delete [] oGraphDir;
	delete [] fileName;
}


//void makeEdges(string &, int, int * &)
//Description: This creates new edges for the collapsed graph
//Input: workDir (string &): the working directory, graphIter (int): The current
//graph iteration, nodeMap (int *): the mapping from the previous graph
//to the new one
//Output:None
//Return:None
void makeEdges(string & workDir, int iGraph, int oGraph, int auxFiles, int * & nodeMap) {

	char * iGraphDir = new char [1000];
	char * oGraphDir = new char [1000];
	char * fileName = new char [1000];

	sprintf(iGraphDir, "%s/Spectrum/Graph%d", workDir.c_str(), iGraph);
	sprintf(oGraphDir, "%s/Spectrum/Graph%d", workDir.c_str(), oGraph);

	int numFiles = countFiles(iGraphDir);
	numFiles-=auxFiles;

	uint32_t * buff = new uint32_t[33 * megabyte];
	int count = 0; int numSets = 0;

	long long int numEdges = 0;
	for(int i = 0; i < numFiles; i++)
	{
		sprintf(fileName, "%s/graph%d", iGraphDir, i);
		FILE * pFile = fopen(fileName, "r");

		Graph graph;
		graph.read(pFile);
		fclose(pFile);

		for(long int j = 0; j < graph.getNumEdges(); j++)
		{
			uint32_t source = nodeMap[graph.edgeSource(j) + graph.getOffset()];
			uint32_t dest = nodeMap[graph.edgeDest(j)];
			uint32_t weight = graph.edgeOvlLen(j)/10;
			weight = weight << 16;

			uint32_t sContained = 0; uint32_t dContained = 0;
			uint32_t rDovetail = 0; uint32_t matePair = 0;

			if(graph.sContained(j))
			{
				sContained = 1;
			}

			if(graph.dContained(j))
			{
				dContained = 1;
			}

			if(graph.rDovetail(j))
				rDovetail = 1;

			if(count >= 33 * megabyte)
			{
				sprintf(fileName, "%s/set%d", oGraphDir, numSets++);
				pFile = fopen(fileName, "w");
				fwrite(buff, sizeof(uint32_t), count, pFile);
				fclose(pFile); count = 0;
			}

			if(source != dest)
			{
				buff[count++] = source;
				buff[count++] = dest;
				buff[count++] = weight | ((sContained << 14) | (dContained << 13) | (rDovetail << 12));
			}
		}
	}

	if(count > 11 * megabyte || numSets == 0)
	{
		sprintf(fileName, "%s/set%d", oGraphDir, numSets++);
		FILE * pFile = fopen(fileName, "w");
		fwrite(buff, sizeof(uint32_t), count, pFile);
		fclose(pFile);

	}else{

		if(count > 0)
		{
			sprintf(fileName, "%s/set%d", oGraphDir, numSets-1);
			FILE * pFile = fopen(fileName, "a");
			fwrite(buff, sizeof(uint32_t), count, pFile);
			fclose(pFile);
		}
	}

	delete [] iGraphDir;
	delete [] oGraphDir;
	delete [] fileName;
}

//void initiateDensity(string)
//Description: This initiates the nodes and edges densities
//Input: workDir (string &): the working directory, graphIter
//Output:None
//Return:None
void initiateDensity(string & workDir) {

	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Spectrum/Graph0", workDir.c_str());

	int numFiles = countFiles(graphDir);

	Graph graph;

	char * fileName = new char [1000];
	sprintf(fileName, "%s/graph%d", graphDir, numFiles-1);
	FILE * pFile = fopen(fileName, "r");

	graph.read(pFile);

	long int numNodes = graph.getOffset() + graph.getNumNodes();

	fclose(pFile);

	int * eDen = new int [numNodes];
	memset(eDen, 0, sizeof(int) * numNodes);

	sprintf(fileName, "%s/eDen", graphDir);
	pFile = fopen(fileName, "w");
	fwrite(eDen, sizeof(int), numNodes, pFile);
	fclose(pFile);

	delete [] eDen;

	int * nDen = new int [numNodes];

	for(long int i = 0; i < numNodes; i++)
	{
		nDen[i] = 1;
	}

	char * ovlInfoDir = new char [1000];
	sprintf(ovlInfoDir, "%s/OvlInfo", workDir.c_str());

	numFiles = countFiles(ovlInfoDir);

	for(int i = 0; i < numFiles; i++)
	{
		sprintf(fileName, "%s/sequences.job%d.info", ovlInfoDir, i);	
		pFile = fopen(fileName, "r");
		fseek(pFile, 0, SEEK_END);
		int size = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);

		int * cOvl = new int [size];
		fread(cOvl, sizeof(int), size, pFile);
		fclose(pFile);

		for(int j = 1; j < size; j+=3)
		{
			if(((j-1)/3)+cOvl[0] != cOvl[j])
				nDen[((j-1)/3)+cOvl[0]] = 0; 	
		}

		delete [] cOvl;
	}

	sprintf(fileName, "%s/nDen", graphDir);
	pFile = fopen(fileName, "w");

	fwrite(nDen, sizeof(int), numNodes, pFile);
	fclose(pFile);

	delete [] nDen;

	delete [] graphDir;
	delete [] fileName;
}

//void findMatching(string &, int, int, float, uint32_t * &, int [], int)
//Description: This finds a matching in the graph
//Input: workDir (string &): the working directory, graphIter (int): The current
//graph iteration, minMerge (int): The minimum overlapping for merging, minDensity (int)
//min Node density for merging, nodeMap (int *): the mapping from the previous graph,
//intervals(int []) : holds quartile information, numIntervals (int): the number
//of intervals
//Output:None
//Return:None
float findMatching(string & workDir, int graphIter, int minMerge, float minDensity, int * & nodeMap, int intervals [], int auxFiles, int numIntervals)
{
	char * iGraphDir = new char [1000];
	char * oGraphDir = new char [1000];

	sprintf(iGraphDir, "%s/Spectrum/Graph%d", workDir.c_str(), graphIter);
	sprintf(oGraphDir, "%s/Spectrum/Graph%d", workDir.c_str(), graphIter+1);

	int numFiles = countFiles(iGraphDir);
	numFiles-=auxFiles;

	char * fileName = new char [1000];
	sprintf(fileName, "%s/eDen", iGraphDir);

	FILE * pFile = fopen(fileName, "r");
	fseek(pFile, 0, SEEK_END);

	int numNodes = ftell(pFile)/sizeof(int);
	fseek(pFile, 0, SEEK_SET);

	int * eDen = new int [numNodes];
	fread(eDen, sizeof(int), numNodes, pFile);
	fclose(pFile);

	sprintf(fileName, "%s/nDen", iGraphDir);
	pFile = fopen(fileName, "r");

	int * nDen = new int [numNodes];
	fread(nDen, sizeof(int), numNodes, pFile);
	fclose(pFile);

	int * eDen_new = new int [numNodes];
	int * nDen_new = new int [numNodes];

	nodeMap = new int [numNodes];
	int * map_inv = new int [2*numNodes];

	for(int i = 0; i < numNodes; i++)
		nodeMap[i] = -1;

	for(int i = 0; i < 2*numNodes; i++)
		map_inv[i] = -1;

	int numMerged = 0; int numMatched = 0;

	long long int numEdges = 0;
	for(int i = numIntervals; i > 0; i--)
	{
		for(int j = 0; j < numFiles; j++)
		{
			sprintf(fileName, "%s/graph%d", iGraphDir, j);
			FILE * pFile = fopen(fileName, "r");

			Graph graph;
			graph.read(pFile);
			fclose(pFile);

			if(i == numIntervals)
				numEdges+=graph.getNumEdges();

			for(int k = 0; k < graph.getNumNodes(); k++)
			{
				uint32_t source = k + graph.getOffset();

				if(nodeMap[source] == -1)
				{
					for(int h = 0; h < graph.nodeDegree(k); h++)
					{
						if(graph.edgeOvlLen(k, h) > intervals[i - 1] && graph.edgeOvlLen(k, h) <= intervals[i])
						{
							uint32_t dest = graph.edgeDest(k, h);

							float nWeight = nDen[source] + nDen[dest];
							float eWeight = eDen[source] + eDen[dest] + graph.edgeOvlLen(k, h);

							float density = (2*eWeight)/(nWeight * (nWeight-1));

							if(graph.edgeOvlLen(k, h) < minMerge)
							{
								break;
							}

							if(density >= minDensity && nodeMap[dest] == -1)
							{
								nodeMap[source] = numMerged;
								nodeMap[dest] = numMerged;

								map_inv[2* numMerged] = source;
								map_inv[2* numMerged+1] = dest;

								eDen_new [numMerged] = eWeight;
								nDen_new[numMerged++] = nWeight;

								numMatched++;
								break;
							}

						}
					}

					if(i == 1 && nodeMap[source] == -1)
					{
						nodeMap[source] = numMerged;

						map_inv[2* numMerged] = source;
						eDen_new [numMerged] = eDen[source];
						nDen_new[numMerged++] = nDen[source];
					}
				}
			}
		}
	}

	scrambleNodes(nodeMap, map_inv, eDen_new, nDen_new, numMerged);

	sprintf(fileName, "%s/map", oGraphDir);
	pFile = fopen(fileName, "w");
	fwrite(map_inv, sizeof(int), numMerged*2, pFile);
	fclose(pFile);

	sprintf(fileName, "%s/eDen", oGraphDir);
	pFile = fopen(fileName, "w");
	fwrite(eDen_new, sizeof(int), numMerged, pFile);
	fclose(pFile);

	sprintf(fileName, "%s/nDen", oGraphDir);
	pFile = fopen(fileName, "w");
	fwrite(nDen_new, sizeof(int), numMerged, pFile);
	fclose(pFile);

	delete [] eDen;
	delete [] nDen;

	delete [] eDen_new;
	delete [] nDen_new;

	delete [] map_inv;

	delete [] iGraphDir;
	delete [] oGraphDir;
	delete [] fileName;

	return (float) (2*numMatched)/(float) numNodes;
}

//void scrambleNodes(int [], int [], int [], int [], int)
//Description: This randomizes the nodes for the next iteration
//Input: nodeMap (int): the node mappings, map_inv (int): the node map inverse,
//Input: eDen_new (int): the edge densities (int), nDen_new (int): the node densities
//Input: numMerges (int): number of nodes merged
//Output:None
//Return:None
void scrambleNodes(int * & nodeMap, int * & map_inv, int * & eDen_new, int * & nDen_new, int numMerged) {

	int * randomNodes = new int [numMerged];
	for(int i = 0; i < numMerged; i++)
		randomNodes[i] = i;

	random_shuffle(randomNodes, randomNodes+numMerged);
	int * tmp = new int [2*numMerged];
	for(int i = 0; i < numMerged; i++)
	{
		tmp[2*i] = map_inv[2*randomNodes[i]];
		tmp[2*i+1] = map_inv[2*randomNodes[i]+1];
	}

	delete [] map_inv;
	map_inv = tmp;

	tmp = new int [numMerged];
	for(int i = 0; i < numMerged; i++)
	{
		tmp[i] = eDen_new[randomNodes[i]];
	}

	delete [] eDen_new;
	eDen_new = tmp;

	tmp = new int [numMerged];
	for(int i = 0; i < numMerged; i++)
	{
		tmp[i] = nDen_new[randomNodes[i]];
	}

	delete [] nDen_new;
	nDen_new = tmp;

	for(int i = 0; i < 2 * numMerged; i++)
	{
		if(map_inv[i] != -1)
		{
			nodeMap[map_inv[i]] = i/2;
		}
	}
	delete [] randomNodes;
}

//void calculateIntervals(string &, const int, int * &, const int)
//Description: This calculates ovl intervals for iterative matching
//Input: workDir (string &): the working directory, graphIter (int): The current
//graph iteration, intervals(int []) : holds quartile information,
//numIntervals (int): the number of intervals
//Output:None
//Return:None
void calculateIntervals(string & workDir, const int graphIter, int * & intervals, const int auxFiles, const int numIntervals) {

	map<int, int> weightCounts;
	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Spectrum/Graph%d", workDir.c_str(), graphIter);

	int numFiles = countFiles(graphDir);
	numFiles-=auxFiles;

	long long int numEdges = 0;

	for(int i = 0; i < numFiles && i < 30; i++)
	{
		char * fileName = new char [1000];
		sprintf(fileName, "%s/graph%d",  graphDir, i);

		Graph graph;

		FILE * pFile = fopen(fileName, "r");
		graph.read(pFile);
		fclose(pFile);

		for (long int j = 0; j < graph.getNumEdges(); j++)
		{
			pair<map<int, int>::iterator, bool > wasInserted;
			wasInserted = weightCounts.insert(pair<int, int>(graph.edgeOvlLen(j), 1 ));
			if(!wasInserted.second)
			{
				(*wasInserted.first).second++;
			}

			numEdges++;
		}
		delete [] fileName;
	}

	int chunkSize = numEdges/numIntervals; int numChunks = 1;

	if(chunkSize = 0) chunkSize = 1;

	intervals[0] = 0; intervals[numIntervals] = INT_MAX;
	long long int countedEdges = 0;

	for(map<int, int>::iterator it = weightCounts.begin(); it != weightCounts.end(); it++)
	{
		countedEdges+=(*it).second;
		if(countedEdges >= numChunks * chunkSize)
		{
			intervals[numChunks++] = (*it).first;
		}

		if(numChunks >= numIntervals) //changed from greater than to greater than equal
			break;
	}

	for(int i = numChunks; i < numIntervals; i++)
	{
		intervals[i] = intervals[numChunks-1];
	} 

	delete [] graphDir;
}

//void quitTasks (const int)
//Description: This function quits all of the worker processors
//Input: i (const int) This is the index of the worker processor that we want to quit
//Output:None
//Return:None
void quitTasks(const int i) {
	int message []  = {0,0,0,0};
	MPI_Send(message, 4, MPI_INT, i, quitTag, MPI_COMM_WORLD);
}

//void sortFragments()
//Description: This sorts a subset of fragments  
//Input: workDir (string &): the working directory, graphIter (int), the current graph iteration
//nTasks (int) the number of available processors
//density
//Output:None
//Return:None
void sortFragments(string & workDir, int nTasks) {

	char * fragmentDir = new char [200];
	sprintf(fragmentDir, "%s/Fragments", workDir.c_str());

	//Count Files
	int numFiles = countFiles(fragmentDir);

	int numSets = 0;
	char * fileName = new char [1000];
	sprintf(fileName, "%s/sequences.job0.set%d", fragmentDir, numSets++);

	while(fileExists(fileName))
	{
		sprintf(fileName, "%s/sequences.job0.set%d", fragmentDir, numSets++);
	}
	delete [] fragmentDir;
	numSets-=1;

	int numJobs = numFiles/numSets;

	//Determine task set size
	int setSize = numFiles/(nTasks-1);
	if(setSize == 0)  setSize++;

	//Seed Tasks
	int currFile = 0; int taskNum = 1; int message [] = {0, 0, 0, 0};
	for(; currFile < numFiles && taskNum < nTasks; currFile+=setSize)
	{
		message[0] = currFile; message[1] = setSize;
		message[2] = numJobs; message[3] = numSets;

		MPI_Send(message, 4, MPI_INT, taskNum++, sortTagF, MPI_COMM_WORLD);
	}

	//Receive/Send Tasks
	MPI_Status status; int result;
	while(currFile < numFiles)
	{
		message[0] = currFile; message[1] = setSize;
		message[2] = numJobs; message[3] = numSets;

		MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		MPI_Send(message, 4, MPI_INT, status.MPI_SOURCE, sortTagF, MPI_COMM_WORLD);

		currFile+=setSize;
	}

	//Clean up
	finalizeTasks(taskNum);
}

//void sortEdges()
//Description: This sorts a subset of files
//Input: workDir (string &): the working directory, graphIter (int), the current graph iteration
//nTasks (int) the number of available processors
//density
//Output:None
//Return:None
void sortEdges(string & workDir, int sortType, int graphIter, int auxFiles, int nTasks) {

	char * graphDir = new char [1000];
	if(graphIter != -1)
		sprintf(graphDir, "%s/Spectrum/Graph%d", workDir.c_str(), graphIter);
	else
		sprintf(graphDir, "%s/OvlGraph", workDir.c_str());

	//Count Files
	int numFiles = countFiles(graphDir);
	numFiles-=auxFiles;

	delete [] graphDir;

	//Determine task set size
	int setSize = numFiles/(nTasks-1);
	if(setSize == 0)  setSize++;

	//Seed Tasks
	int currFile = 0; int taskNum = 1; int message [] = {0, 0, 0, 0};
	for(; currFile < numFiles && taskNum < nTasks; currFile+=setSize)
	{
		message[0] = currFile; message[1] = setSize;
		message[2] = graphIter; message[3] = sortType;

		MPI_Send(message, 4, MPI_INT, taskNum++, sortTagE, MPI_COMM_WORLD);
	}

	//Receive/Send Tasks
	MPI_Status status; int result;
	while(currFile < numFiles)
	{
		message[0] = currFile; message[1] = setSize;
		message[2] = graphIter; message[3] = sortType;

		MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		MPI_Send(message, 4, MPI_INT, status.MPI_SOURCE, sortTagE, MPI_COMM_WORLD);

		currFile+=setSize;
	}

	//Clean up
	finalizeTasks(taskNum);
}

//void makeGraph()
//Description: This function builds graphs from the edges
//Input: workDir (string &): the working directory, graphIter (int), the current graph iteration
//nTasks (int) the number of available processors
//Output:None
//Return:None
void makeGraph(string & workDir, int graphIter, int auxFiles, int nTasks) {

	char * graphDir = new char [1000];
	if(graphIter != -1)
		sprintf(graphDir, "%s/Spectrum/Graph%d", workDir.c_str(), graphIter);
	else
		sprintf(graphDir, "%s/OvlGraph", workDir.c_str());


	//Count Files
	int numFiles =  countFiles(graphDir);
	numFiles-=auxFiles;

	delete [] graphDir;

	//Determine task set size
	int setSize = numFiles/(nTasks-1);
	if(setSize == 0)  setSize++;

	//Seed Tasks
	int currFile = 0; int taskNum = 1; int message [] = {0, 0, 0, 0};
	for(; currFile < numFiles && taskNum < nTasks; currFile+=setSize)
	{
		message[0] = currFile; message[1] = setSize;
		message[2] = graphIter;

		MPI_Send(message, 4, MPI_INT, taskNum++, graphTag, MPI_COMM_WORLD);
	}

	//Receive/Send Tasks
	MPI_Status status; int result;
	while(currFile < numFiles)
	{
		message[0] = currFile; message[1] = setSize;
		message[2] = graphIter;

		MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		MPI_Send(message, 4, MPI_INT, status.MPI_SOURCE, graphTag, MPI_COMM_WORLD);

		currFile+=setSize;
	}

	//Clean up
	finalizeTasks(taskNum);

	//Remove all the set files
	char * dirName = new char [1000]; char * fileName = new char [1000]; int eSet = 0; 
	if(graphIter != -1)
		sprintf(dirName,"%s/Spectrum/Graph%d", workDir.c_str(), graphIter);
	else
		sprintf(dirName,"%s/OvlGraph/", workDir.c_str());

	sprintf(fileName, "%s/set%d", dirName, eSet++);

	while(remove(fileName) == 0)
	{
		sprintf(fileName,"%s/set%d", dirName, eSet++);
	}


	delete [] dirName; delete [] fileName;

}

//void mergeEdges()
//Description: This sorts a subset of files
//Input: workDir (string &): the working directory, graphIter (int), the current graph iteration
//nTasks (int) the number of available processors
//density
//Output:None
//Return:None
void mergeEdges(string & workDir, int graphIter, int auxFiles, int nTasks){

	char * graphDir = new char[1000];
	if(graphIter != -1)
		sprintf(graphDir, "%s/Spectrum/Graph%d", workDir.c_str(), graphIter);
	else
		sprintf(graphDir, "%s/OvlGraph", workDir.c_str());

	char * tmpDir = new char [1000];
	sprintf(tmpDir, "%s/Spectrum/Tmp", workDir.c_str());

	mkdir(tmpDir, S_IRWXU);

	int numFiles =  countFiles(graphDir);
	numFiles-=auxFiles;

	//Determine Merge Type
	int mergeType = 0;
	if(graphIter > 0) mergeType = 1;

	//Seed Tasks
	int currFile = 0; int taskNum = 1; int setSize = 1; int message [] = {0, 0, 0, 0};
	for(; currFile < numFiles && taskNum < nTasks; currFile+=(setSize * 2/*10*/))
	{
		message[0] = currFile; message[1] = setSize;
		message[2] = graphIter; message[3] = mergeType;
		MPI_Send(message, 4, MPI_INT, taskNum++, mergeTagE, MPI_COMM_WORLD);
	}

	//Receive/Send Tasks
	MPI_Status status; int result; bool multipleIterations = false;
	while(setSize < numFiles)
	{
		multipleIterations = true;
		while(currFile < numFiles)
		{
			message[0] = currFile; message[1] = setSize;
			message[2] = graphIter; message[3] = mergeType;

			MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			MPI_Send(message, 4, MPI_INT, status.MPI_SOURCE, mergeTagE, MPI_COMM_WORLD);

			currFile+=(setSize * 2/*10*/);
		}

		//Clean Up
		finalizeTasks(taskNum); cleanTmps(graphDir, tmpDir);
		setSize*=2/*10*/; currFile = 0; taskNum = 1;
		for(; setSize < numFiles && currFile < numFiles && taskNum < nTasks; currFile+=(setSize * 2/*10*/))
		{
			message[0] = currFile; message[1] = setSize;
			message[2] = graphIter; message[3] = mergeType;

			MPI_Send(message, 4, MPI_INT, taskNum++, mergeTagE, MPI_COMM_WORLD);
		}
	}

	//Clean up
	finalizeTasks(taskNum);

	if(!multipleIterations)
		cleanTmps(graphDir, tmpDir);

	removeDir(tmpDir);
	renameFiles(graphDir);

	delete [] tmpDir; delete [] graphDir;
}

//void mergeFragments()
//Description: This merges 
//Input: workDir (string &): the working directory, graphIter (int), the current graph iteration
//nTasks (int) the number of available processors
//density
//Output:None
//Return:None
void mergeFragments(string & workDir, int nTasks){

	char * fragmentDir = new char[200];
	sprintf(fragmentDir, "%s/Fragments", workDir.c_str());

	char * tmpDirF = new char [200];
	sprintf(tmpDirF, "%s/TmpF", workDir.c_str());
	mkdir(tmpDirF, S_IRWXU);

	char * labelDir = new char[200];
	sprintf(labelDir, "%s/Labels", workDir.c_str());

	char * tmpDirL = new char [200];
	sprintf(tmpDirL, "%s/TmpL", workDir.c_str());
	mkdir(tmpDirL, S_IRWXU);

	int numFiles =  countFiles(fragmentDir);
	//Seed Tasks
	int currFile = 0; int taskNum = 1; int setSize = 1; int message [] = {0, 0, 0, 0};
	for(; currFile < numFiles && taskNum < nTasks; currFile+=(setSize * 2/*10*/))
	{
		message[0] = currFile; message[1] = setSize;
		MPI_Send(message, 4, MPI_INT, taskNum++, mergeTagF, MPI_COMM_WORLD);
	}
	//Receive/Send Tasks
	MPI_Status status; int result; bool multipleIterations = false;
	while(setSize < numFiles)
	{
		multipleIterations = true;
		while(currFile < numFiles)
		{
			message[0] = currFile; message[1] = setSize;

			MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			MPI_Send(message, 4, MPI_INT, status.MPI_SOURCE, mergeTagF, MPI_COMM_WORLD);

			currFile+=(setSize * 2/*10*/);
		}

		//Clean Up
		finalizeTasks(taskNum);
		removeDir(fragmentDir); 
		rename(tmpDirF, fragmentDir); mkdir(tmpDirF, S_IRWXU);
		removeDir(labelDir); 
		rename(tmpDirL, labelDir); mkdir(tmpDirL, S_IRWXU);
		setSize*=2/*10*/; currFile = 0; taskNum = 1;
		for(; setSize < numFiles && currFile < numFiles && taskNum < nTasks; currFile+=(setSize * 2/*10*/))
		{
			message[0] = currFile; message[1] = setSize;
			MPI_Send(message, 4, MPI_INT, taskNum++, mergeTagF, MPI_COMM_WORLD);
		}
	}

	//Clean up
	finalizeTasks(taskNum);
	if(!multipleIterations)
	{
		removeDir(fragmentDir); removeDir(labelDir);
		rename(tmpDirF, fragmentDir);
		rename(tmpDirL, labelDir);
	}
	else
	{
		removeDir(tmpDirF);
		removeDir(tmpDirL);
	}

	delete [] tmpDirF; delete [] tmpDirL; delete [] fragmentDir; delete [] labelDir;
}


//void cleanEdgeTmps(const char [], const char [])
//Description: This function cleans up the TMP and Graph directories
//Input: graphDir (const char []), tmpDir (const char [])
//Output:None
//Return:None
void cleanTmps(const char graphDir [], const char tmpDir [] ) {

	//Deleting/Moving directories
	char * fileName = new char[1000];
	sprintf(fileName, "%s/map", graphDir);
	char * dest = new char [1000];

	if(fileExists(fileName))
	{
		sprintf(dest, "%s/map", tmpDir);
		rename(fileName, dest);
	}

	sprintf(fileName, "%s/eDen", graphDir);
	if(fileExists(fileName))
	{
		sprintf(dest, "%s/eDen", tmpDir);
		rename(fileName, dest);
	}

	sprintf(fileName, "%s/nDen", graphDir);
	if(fileExists(fileName))
	{
		sprintf(dest, "%s/nDen", tmpDir);
		rename(fileName, dest);
	}

	sprintf(fileName, "%s/finalMap", graphDir);
	if(fileExists(fileName))
	{
		sprintf(dest, "%s/finalMap", tmpDir);
		rename(fileName, dest);
	}

	delete [] dest;
	delete [] fileName;

	removeDir(graphDir); 
	rename(tmpDir, graphDir);
	mkdir(tmpDir, S_IRWXU);
}

//void cleanTmps(const char [], const char [])
//Description: This function cleans up the TMP and Graph directories
//Input: graphDir (const char []), tmpDir (const char [])
//Output:None
//Return:None
void renameFiles(const char graphDir []) {
	//Renaming files
	char * fileName = new char[1000];
	int numFiles = countFiles(graphDir);
	int fileIndex = 0;

	for (int i = 0; i < numFiles; i++)
	{
		sprintf(fileName, "%s/set%d", graphDir, i);

		if(fileExists(fileName))
		{
			sprintf(fileName, "mv %s/set%d %s/set%d", graphDir, i, graphDir, fileIndex);
			char * name1 = new char [1000]; char * name2 = new char[1000];

			sprintf(name1, "%s/set%d", graphDir, i);
			sprintf(name2, "%s/set%d", graphDir, fileIndex);

			if(fileIndex != i)
			{
				rename(name1, name2);
			}

			fileIndex++;
		}
	}

	delete [] fileName;
}

//void countFiles(const char [])
//Description: This function counts the number of fragments in the fragment set
//Input: workDir (string &) : The working directory
//Output:None
//Return:None
unsigned long long int countFragments(string & workDir)
{
	char * fileName = new char [1000];
	sprintf(fileName, "%s/Fragments", workDir.c_str());

	int numFiles = countFiles(fileName);
	int numJobs = 0; int numSets = 0;

	for(int i = 0; i < numFiles; i++)
	{
		sprintf(fileName, "%s/Fragments/sequences.job0.set%d", workDir.c_str(), i);
		if(!fileExists(fileName))
			break;

		numSets++;
	}

	numJobs = numFiles/numSets;

	sprintf(fileName, "%s/Fragments/sequences.job%d.set%d", workDir.c_str(), numJobs-1, numSets-1);

	unsigned long long int indexOffset = 0; unsigned long long int totalFragments;

	FILE * pFile = fopen(fileName, "r");
	fread(&indexOffset, sizeof(unsigned long long int), 1, pFile);
	fread(&indexOffset, sizeof(unsigned long long int), 1, pFile);

	fread(&totalFragments, sizeof(unsigned long long int), 1, pFile);
	fread(&totalFragments, sizeof(unsigned long long int), 1, pFile);
	fclose(pFile); //Added this 

	delete [] fileName;

	return (indexOffset + totalFragments);
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

//void makeDir(const char [])
//Description: This function calls the linux mkdir command
//Input: dirName (const char [])
//Output:None
//Return:None
void makeDir(const char dirName []) {

	char * sysCmd = new char [1000];
	sprintf(sysCmd, "mkdir %s", dirName);
	system(sysCmd);
	delete [] sysCmd;
}

//void moveDir(const char [], const char [])
//Description: This function calls the linux mv command
//Input: source (const char []), dest (const char [])
//Output:None
//Return:None
void moveDir(const char source [], const char dest []) {

	char * sysCmd = new char [1000];
	sprintf(sysCmd, "mv %s %s", source, dest);
	system(sysCmd);
	delete [] sysCmd;
}

//void removeDir(const char [])
//Description: This function calls the linux rm -r command
//Input: dirName (const char [])
//Output:None
//Return:None
void removeDir(const char dirName []) {

	DIR *pDir = opendir(dirName);
	struct dirent * pDirent;
	char * fileName = new char [1000];

	while ((pDirent = readdir(pDir)) != NULL)
		if(strcmp(pDirent->d_name, ".") != 0 && strcmp(pDirent->d_name, "..") != 0)
		{
			sprintf(fileName, "%s/%s", dirName, pDirent->d_name);
			remove(fileName);
		}

	rmdir(dirName);
	delete [] fileName;
}

//void finalizeTasks(const int )
//Description: This function finalizes child tasks that are running
//Input: taskNum (cosnt int) : The number of tasks that need to be finalized
//Output:None
//Return:None
void finalizeTasks(const int taskNum)
{
	MPI_Status status; int result;
	for(int i = 1; i < taskNum; i++)
	{
		MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	}
}

//void runChildTask()
//Description: This runs the the child thread in a parallel program
//Input: workDir (string &): the working directory
//density
//Output:None
//Return:None
void runChildTask(string & workDir){

	while(true)
	{
		MPI_Status status;
		int message [] = {0,0,0,0};
		MPI_Recv(message, 4, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		if(status.MPI_TAG == sortTagE)
		{
			int currFile = message[0]; 
			int setSize = message[1]; 
			int graphIter = message[2]; 
			int sortType = message[3];

			for(int i = 0; i < setSize; i++)
			{
				char * fileName = new char [1000];
				if(graphIter != -1)
					sprintf(fileName, "%s/Spectrum/Graph%d/set%d", workDir.c_str(), graphIter, currFile+i);
				else
					sprintf(fileName, "%s/OvlGraph/set%d", workDir.c_str(),  currFile+i);

				SortEdges mySort(sortType);
				mySort.sort_edges(fileName, fileName);

				delete [] fileName;
			}
		}

		if(status.MPI_TAG == sortTagF)
		{
			int currFile = message[0]; 
			int setSize = message[1]; 
			int numJobs = message[2]; 
			int numSets = message[3];

			char * mapFile = new char [200];

			sprintf(mapFile, "%s/OvlGraph/finalMap", workDir.c_str());
			FILE * pFile = fopen(mapFile, "r");
			fseek(pFile, 0, SEEK_END);
			int size = ftell(pFile)/sizeof(int);
			fseek(pFile, 0, SEEK_SET);

			int * nodeMap = new int[size];
			fread(nodeMap, sizeof(int), size, pFile);

			fclose(pFile);
			delete [] mapFile;

			int i = currFile/numSets; int j = currFile%numSets;
			for(; i < numJobs && (i * numSets + j) < currFile+setSize; i++)
			{
				for(; j < numSets && (i * numSets + j) < currFile+setSize; j++)
				{
					char * iFragments = new char [200];
					sprintf(iFragments, "%s/Fragments/sequences.job%d.set%d", workDir.c_str(), i, j);
					char * oFragments = new char [200];
					sprintf(oFragments, "%s/Fragments/sequences.%d", workDir.c_str(), (i*numSets)+j);
					char * iLabels = new char [200];
					sprintf(iLabels, "%s/Labels/labels.job%d.set%d", workDir.c_str(), i, j);
					char * oLabels = new char [200];
					sprintf(oLabels, "%s/Labels/labels.%d", workDir.c_str(), (i*numSets)+j);

					SortFragments mySort;
					mySort.sort(iFragments, iLabels, oFragments, oLabels, nodeMap);

					remove(iFragments);
					remove(iLabels);

					delete [] iFragments;
					delete [] oFragments;

					delete [] iLabels;
					delete [] oLabels;
				}

				j = 0;
			}

			delete [] nodeMap;
		}

		if(status.MPI_TAG == mergeTagE)
		{
			int currFile = message[0]; 
			int setSize = message[1]; 
			int graphIter = message[2]; 
			int sortType = message[3];

			char * inDir = new char [1000];
			char * outDir = new char [1000];
			if(graphIter != -1)
				sprintf(inDir, "%s/Spectrum/Graph%d", workDir.c_str(), graphIter);
			else
				sprintf(inDir, "%s/OvlGraph", workDir.c_str());

			sprintf(outDir, "%s/Spectrum/Tmp", workDir.c_str());

			MergeEdges myMerge(sortType);
			myMerge.merge(currFile, setSize, inDir, outDir);

			delete [] inDir;
			delete [] outDir;
		}

		if(status.MPI_TAG == mergeTagF)
		{
			int currFile = message[0]; 
			int setSize = message[1]; 

			char * inDirF = new char [1000];
			char * outDirF = new char [1000];
			sprintf(inDirF, "%s/Fragments", workDir.c_str());
			sprintf(outDirF, "%s/TmpF", workDir.c_str());

			char * inDirL = new char [1000];
			char * outDirL = new char [1000];
			sprintf(inDirL, "%s/Labels", workDir.c_str());
			sprintf(outDirL, "%s/TmpL", workDir.c_str());

			MergeFragments myMerge;
			myMerge.merge(currFile, setSize, inDirF, outDirF, inDirL, outDirL);

			delete [] inDirF;
			delete [] outDirF;
			delete [] inDirL;
			delete [] outDirL;
		}

		if(status.MPI_TAG == graphTag) //Fix this for Parallel Computing it will not erase files right
		{
			int currFile = message[0]; 
			int setSize = message[1]; 
			int graphIter = message[2]; 
			int sortType = message[3];

			for(int i = currFile; i < currFile + setSize; i++)
			{
				uint32_t offset = 0; uint32_t currNode = 0;
				if(i != 0)
				{
					char * fileName = new char [1000];
					if(graphIter != -1)
						sprintf(fileName,"%s/Spectrum/Graph%d/set%d", workDir.c_str(), graphIter,  i-1);
					else
						sprintf(fileName,"%s/OvlGraph/set%d", workDir.c_str(), i-1);

					FILE * pFile = fopen(fileName, "r");

					fseek(pFile, 0, SEEK_END);
					int fileSize = ftell(pFile)/sizeof(uint32_t);
					fseek(pFile, 0, SEEK_SET);

					uint32_t * arr = new uint32_t [fileSize];
					fread(arr, sizeof(uint32_t), fileSize, pFile);
					fclose(pFile);

					offset = arr[fileSize - 3];
					currNode = ++offset;

					delete [] arr;
					delete [] fileName;
				}

				char * fileName = new char [1000];
				if(graphIter != -1)
					sprintf(fileName,"%s/Spectrum/Graph%d/set%d", workDir.c_str(), graphIter,  i);
				else
					sprintf(fileName,"%s/OvlGraph/set%d", workDir.c_str(), i);


				if(!fileExists(fileName))
				{
					delete [] fileName;
					break;
				}	

				FILE * pFile = fopen(fileName, "r");

				fseek(pFile, 0, SEEK_END);
				int fileSize = ftell(pFile)/sizeof(uint32_t);
				fseek(pFile, 0, SEEK_SET);

				uint32_t * arr = new uint32_t [fileSize];
				fread(arr, sizeof(uint32_t), fileSize, pFile);
				fclose(pFile);

				int numEdges = fileSize/3;
				int numNodes = arr[fileSize-3] - currNode + 1;

				if(graphIter != -1)
					sprintf(fileName,"%s/Spectrum/Graph%d/set%d", workDir.c_str(), graphIter,  i+1);
				else
					sprintf(fileName,"%s/OvlGraph/set%d", workDir.c_str(),  i+1);

				if(!fileExists(fileName))
				{
					if(graphIter != -1)
						sprintf(fileName,"%s/Spectrum/Graph%d/map", workDir.c_str(), graphIter);
					else
						sprintf(fileName,"%s/OvlGraph/map", workDir.c_str());

					if(fileExists(fileName))
					{
						FILE * pFile = fopen(fileName, "r");
						fseek(pFile, 0, SEEK_END);
						numNodes = (ftell(pFile)/sizeof(int))/2 - currNode;
						fclose(pFile);
					}else
					{
						numNodes = countFragments(workDir) - currNode;
					}
				}

				Graph graph;
				graph.reserveNumEdges(numEdges);
				graph.reserveNumNodes(numNodes);
				graph.setOffset(offset);
				if(graphIter != 0 && graphIter != -1)
					graph.setHeavyGraph();

				if(arr[0] == currNode)
				{
					graph.newEdgeSet();
				}
				else
					currNode--;


				for(int j = 0; j < fileSize; j+=3)
				{

					if(arr[j] != currNode)
					{
						graph.addSilentNodes(arr[j] - currNode - 1);
						graph.newEdgeSet();
					}

					graph.addEdge(arr[j+1], arr[j+2]);
					currNode = arr[j];
				}

				if(graphIter != -1)
					sprintf(fileName,"%s/Spectrum/Graph%d/set%d", workDir.c_str(), graphIter,  i+1);
				else
					sprintf(fileName,"%s/OvlGraph/set%d", workDir.c_str(), i+1);

				if(!fileExists(fileName))
				{
					if(graphIter != -1)
						sprintf(fileName,"%s/Spectrum/Graph%d/map", workDir.c_str(), graphIter);
					else
						sprintf(fileName,"%s/OvlGraph/map", workDir.c_str());

					if(fileExists(fileName))
					{
						FILE * pFile = fopen(fileName, "r");
						fseek(pFile, 0, SEEK_END);
						numNodes = (ftell(pFile)/sizeof(int))/2;
						fclose(pFile);

						graph.addSilentNodes(numNodes - currNode - 1);

					}else
					{
						numNodes = countFragments(workDir);
						graph.addSilentNodes(numNodes - currNode - 1);
					}

				}

				graph.set();

				if(graphIter != -1)
					sprintf(fileName,"%s/Spectrum/Graph%d/graph%d", workDir.c_str(), graphIter,  i);
				else
					sprintf(fileName,"%s/OvlGraph/graph%d", workDir.c_str(), i);

				pFile = fopen(fileName, "w");
				graph.write(pFile);
				fclose(pFile);

				delete [] arr;
				delete [] fileName;
			}

		}

		if(status.MPI_TAG == quitTag)
			break;

		MPI_Send(0, 0, MPI_INT, 0, status.MPI_TAG, MPI_COMM_WORLD);
	}
}

//void createGraphDir(string &, int)
//Description: This function creates the graph directory
//Input: workDir(string): This is the name of the working directory,
//Input: graph (int), the current graph
//Output:None
//Return:None
void createGraphDir(string & workDir, int graphIter) {

	char * graphDir = new char[1000];
	sprintf(graphDir, " mkdir %s/Spectrum/Graph%d", workDir.c_str(), graphIter);
	system(graphDir);

	delete [] graphDir;
}

//void moveEdges(string &)
//Description: This function moves edges from the overlap directory
//Input: workDir(string): This is the name of the working directory,
//Output:None
//Return:None
void moveEdges(string & workDir)
{
	char * ovlDir = new char[1000];
	char * graphDir = new char[1000];
	sprintf(ovlDir, "%s/Overlaps", workDir.c_str());
	sprintf(graphDir, "%s/Spectrum/Graph0", workDir.c_str());

	DIR *pDir = opendir(ovlDir);
	struct dirent * pDirent;

	uint32_t * buff = new uint32_t[33 * megabyte];
	int index = 0; int set = 0;

	if(pDir != NULL)
	{
		while ((pDirent = readdir(pDir)) != NULL)
		{
			if(strcmp(pDirent->d_name, ".") != 0 && strcmp(pDirent->d_name, "..") != 0)
			{
				char * file = new char[1000];
				sprintf(file, "%s/%s", ovlDir, pDirent->d_name);
				FILE * pFile = fopen(file, "r");
				int count = 0;

				while((count = fread(buff+index, sizeof(uint32_t), 33*megabyte-index, pFile)) == 33*megabyte-index)
				{
					char * file2 = new char [1000];
					sprintf(file2, "%s/set%d", graphDir, set++);
					FILE * pFile2 = fopen(file2, "w");
					fwrite(buff, sizeof(uint32_t), index+count, pFile2);
					fclose(pFile2);
					delete [] file2;
					index = 0; count = 0;
				}
				index+=count;

				fclose(pFile);
				delete [] file;
			}
		}

		if(index > 11 * megabyte || set == 0)
		{
			char * file = new char [1000];
			sprintf(file, "%s/set%d", graphDir, set++);
			FILE * pFile = fopen(file, "w");
			fwrite(buff, sizeof(uint32_t), index, pFile);
			fclose(pFile);
			delete [] file;

		}else{

			if(index > 0)
			{
				char * file = new char [1000];
				sprintf(file, "%s/set%d", graphDir, set-1);
				FILE * pFile = fopen(file, "a");
				fwrite(buff, sizeof(uint32_t), index, pFile);
				fclose(pFile);
				delete [] file;
			}
		}

		closedir(pDir);
	}

	delete [] buff;
	delete [] ovlDir;
	delete [] graphDir;
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
