#include <iostream>
#include <map>
#include <string>
#include "mpi.h"
#include <algorithm>
#include <cstdio>
#include <stdint.h>
#include <vector>
#include <list>
#include <set>
#include <climits>
#include <queue>
#include <dirent.h>
using namespace std;
#include "Dictionary.h"
#include "Fragment_Index.h"
#include "Labels.h"
#include "Graph.h"
#include "BandedAlignment.h"
#include "SortEdges.h"
#include "MergeEdges.h"

//Size unit
#define megabyte 1024000

char arr []  = {'A', 'C', 'T', 'G', '-'};


#define assembleTag 1
#define sortTag 2
#define mergeTag 3
#define graphTag 4
#define workTag 5
#define endTag 0
#define help 911

void help_message(){

	cout<<"\n\n";
	cout<<"\t\t\t Merge and Traverse Assembler"<<endl;
	cout<<"\t\t\tWritten by Julia Warnke\n"<<endl;
	cout<<"Command line arguments\n"<<endl;
	cout<<"--workDir :working directory :required"<<endl;
	cout<<"--graphIter : The current graph iteration"<<endl;
	cout<<"--minDensity : The minimum inter-cluster d"<<endl;
	cout<<"--percentMerged : The minimum percentage of nodes that have been merged "<<endl;
	cout<<"for the algorithm to continue merging"<<endl;
	cout<<"--verbose : verbose mode"<<endl;
	cout<<"\n\n";

	cout<<"\n\n"<<endl;

	cout<<"Exiting program"<<endl;
	MPI_Abort(MPI_COMM_WORLD, help);
}

int sortByClusterSize(const void * a, const void * b) {

	if(((int *)a)[1] < ((int *) b)[1])
	{
		return 1;
	}else{

		if(((int *)a)[1] > ((int *) b)[1])
		{
			return -1;
		}
	}

	return 0;
}


void initHybSpec(string & workDir);
void makeDir(const char dirName []);
int countFiles(const char []);
int recoverClusters(string &, int * &, int * &, int);
int recoverClusterIndex(string &, int * &, int);
int recoverHybClusters(string &, int * &, int);
void outputClusters(string &, int * &, int, int);
void loadFragments(string &, Fragment_Index ** &, int &, int &, int &, int * &, int, int, int * &);
void loadGraphs(string &, Graph ** &, int &, int &, int &, int * &, int, int, int * &);
void loadLabels(string &, Labels ** &, int &, int &, int &, int * &, int, int, int * &);
void removeContainments(Graph ** &, int, int * &, int, int, int * &, int, int);
void transitiveReduce(Graph ** &, int, uint32_t * &, uint32_t * &, uint32_t * &, int * &, int, int, int * &, int, int);
int trimGraph(uint32_t * &, uint32_t * &, uint32_t * &, uint32_t * &, int, int, int);
int popBubbles(uint32_t * &, uint32_t * &, uint32_t * &, uint32_t * &, int, int, int);
int traverseGraph(uint32_t * &, uint32_t * &, uint32_t * &, uint32_t * &, uint32_t * &, int, int);
void makeEdges(string &, int, int, int * &);
void sortEdges(string &, int, int, int, int);
void mergeEdges(string &, int, int, int);
void makeGraph(string &, int, int, int, int);
void finalizeTasks(const int);
void cleanTmps(const char [], const char []);
void removeDir(const char []);
void moveDir(const char [], const char []);
bool fileExists(const char []);
void renameFiles(const char []);

int main(int argc, char * argv [])
{	
	//Set up the MPI
	int rank;
	int nTasks;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nTasks);

	//User Input Flags
	string workDirS = "--workDir";               //Current Working Directory
	string minPathLenS = "--minPathLen";         //Minimum path length
	string minBubbleLenS = "--minBubbleLen";     //Min length of the bubble
	string samS = "--sam";                       //output sam format 
	string verboseS = "--verbose";		     //verbose mode

	string workDir = "";
	int minPathLen = 2;
	int minBubbleLen = 2;
	bool sam = false;
	bool verbose = false;	

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

		//The minimum path length
		if(argv[i] == minPathLenS){
			minPathLen = atoi(argv[i+1]);
		}

		//The minimum length of the bubble
		if(argv[i] == minBubbleLenS){
			minBubbleLen = 2;
		}	

		//The minimum percent of node merged to continue
		if(argv[i] == samS && sam != true){
			sam = true;
			i--;
		}

		//Verbose or not
		if(argv[i] == verboseS && verbose != true)
		{
			verbose = true;
			i--;
		}
	}	

	if(rank == 0)
	{
		initHybSpec(workDir);
		char * graphDir = new char [1000];
		sprintf(graphDir, "%s/Spectrum", workDir.c_str());
		int numGraphs = countFiles(graphDir);

		int clusterCount = 0;
		for(int i = numGraphs-1; i > 0; i--)
		{
			int * cIndex;
			long int numNodes = recoverClusterIndex(workDir, cIndex, i);

			int * contigs;
			int count = recoverClusters(workDir, cIndex, contigs, i);

			int * cSizes = new int [numNodes * 2];

			//This hold currPos in binArr, binArr max Size, current total content
			int * binStats = new int [(nTasks-1) * 3];
			for(int j = 0; j < (nTasks-1)*3; j+=3)
			{
				binStats[j] = 0; binStats[j+1] = 0; 
				binStats[j+2] = numNodes/(nTasks-1) + numNodes/((nTasks-1)*(nTasks-1));
			}

			int ** bins = new int * [nTasks-1];
			for(int j = 0; j < nTasks-1; j++)
				bins[j] = new int[numNodes/(nTasks-1)+numNodes/((nTasks-1)*(nTasks-1))];

			char * fileName = new char [200];
			sprintf(fileName, "%s/HybSpectrum/Graph%d/tmpMap", workDir.c_str(), i);
			FILE * pFile = fopen(fileName, "r");
			fseek(pFile, 0, SEEK_END);
			int mCount1 = ftell(pFile)/sizeof(int);
			fseek(pFile, 0, SEEK_SET); 

			int * map1 = new int[mCount1];						
			fread(map1, sizeof(int), mCount1, pFile);
			fclose(pFile);

			int currNode = 0;
			int currSize = 0; int numRecorded = 0;
			for(int j = 0; j < count; j++)
			{
				if(contigs[j] == -1)
				{	
					if(map1[cIndex[currNode]] == -1)
					{
						cSizes[2* numRecorded] = cIndex[currNode];
						cSizes[2* numRecorded+1] = currSize;
						numRecorded++;
					}
					currNode++; currSize = 0;
				}else{
					currSize++;
				}		
			}										

			qsort(cSizes, numRecorded, sizeof(int) * 2, sortByClusterSize);

			for(int j = 0; j < numRecorded * 2; j+=2)
			{
				int min = INT_MAX; int minPos = 0;
				for(int k = 0; k < 3*(nTasks-1); k+=3)
				{			
					if(binStats[k] < min)
					{
						min = binStats[k]; minPos = k;
					}  
				}	

				int currSize = binStats[minPos+1]; 
				int maxSize = binStats[minPos+2];
				int bin = minPos/3;

				if(currSize + 1 <= maxSize)
				{
					bins[bin][currSize] = cSizes[j]; 
					binStats[minPos]+=cSizes[j+1];
					binStats[minPos+1]++; 
				}else{
					int * tmpBuff = new int[maxSize + numNodes/(nTasks-1)];
					memcpy(tmpBuff, bins[bin], sizeof(int) * currSize);

					delete [] bins[bin];

					bins[bin] = tmpBuff;

					binStats[minPos+2] = maxSize + numNodes/(nTasks-1);

					bins[bin][currSize] = cSizes[j]; 
					binStats[minPos]+=cSizes[j+1];
					binStats[minPos+1]++; 
				} 
			}			

			for(int j = 0; j < nTasks-1; j++)
			{
				MPI_Send(&i, 1, MPI_INT, j+1, assembleTag, MPI_COMM_WORLD);
				MPI_Send(&i, 1, MPI_INT, j+1, assembleTag, MPI_COMM_WORLD);
				MPI_Send(&binStats[j*3+1], 1, MPI_INT, j+1, assembleTag, MPI_COMM_WORLD);
				MPI_Send(bins[j], binStats[j*3+1], MPI_INT, j+1, assembleTag, MPI_COMM_WORLD);	
			}

			delete [] cSizes;
			delete [] binStats;
			delete [] contigs;
			delete [] cIndex;

			for(int j = 0; j < nTasks-1; j++)
				delete [] bins[j];

			delete [] bins;

			sprintf(fileName, "%s/Spectrum/Graph%d/map", workDir.c_str(), i);
			pFile = fopen(fileName, "r");

			fseek(pFile, 0, SEEK_END);
			int ncCount = ftell(pFile)/sizeof(int);
			fseek(pFile, 0, SEEK_SET);

			int * nodeMap = new int[ncCount];						
			fread(nodeMap, sizeof(int), ncCount, pFile);
			fclose(pFile);

			sprintf(fileName, "%s/HybSpectrum/Graph%d/tmpMap", workDir.c_str(), i-1);
			pFile = fopen(fileName, "r");
			fseek(pFile, 0, SEEK_END);
			int mCount2 = ftell(pFile)/sizeof(int);
			fseek(pFile, 0, SEEK_SET); 

			int * map2 = new int[mCount2];						
			fread(map2, sizeof(int), mCount2, pFile);
			fclose(pFile);

			for(int j = 0; j < mCount1; j++)
			{	
				if(map1[j] != -1)
				{
					int child1 = nodeMap[2*j];
					map2[child1] = map1[j];
					if(nodeMap[2*j+1] != -1)
					{
						int child2 = nodeMap[2*j+1];
						map2[child2] = map1[j];
					}
				}
			}

			int numReceived = 0; int handshake = 0; MPI_Status status;
			while(numReceived < (nTasks-1))
			{
				int buffSize = 0;
				MPI_Recv(&buffSize, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

				if(status.MPI_TAG == workTag)
				{
					int * buff = new int [buffSize]; 
					MPI_Send(&handshake, 1, MPI_INT, status.MPI_SOURCE, workTag, MPI_COMM_WORLD);
					MPI_Recv(buff, buffSize, MPI_INT, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);			
					MPI_Send(&handshake, 1, MPI_INT, status.MPI_SOURCE, workTag, MPI_COMM_WORLD);

					for(int j = 0; j < buffSize; j++)
					{
						map1[buff[j]] = clusterCount;

						int child1 = nodeMap[2*buff[j]];
						map2[child1] = clusterCount;
						if(nodeMap[2*buff[j]+1] != -1)
						{
							int child2 = nodeMap[2*buff[j]+1];
							map2[child2] = clusterCount;
						}

						clusterCount++;
					}

					MPI_Recv(&handshake, 1, MPI_INT, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
					delete [] buff;

				}else{
					if(status.MPI_TAG == endTag)
						numReceived++;
				}
			} 

			int relable = clusterCount;
			for(int j = 0; j < mCount1; j++)
			{
				if(map1[j] == -1)
				{
					map1[j] = relable++;
				}	
			}
			sprintf(fileName, "%s/HybSpectrum/Graph%d/tmpMap", workDir.c_str(), i);
			pFile = fopen(fileName, "w");
			fwrite(map1, sizeof(int), mCount1, pFile);
			fclose(pFile);

			sprintf(fileName, "%s/HybSpectrum/Graph%d/tmpMap", workDir.c_str(), i-1);
			pFile = fopen(fileName, "w");
			fwrite(map2, sizeof(int), mCount2, pFile);
			fclose(pFile);

			delete [] map2;

			int auxFiles = 3;
			makeEdges(workDir, i, 3, map1);				
			sortEdges(workDir, 0, i, 1, nTasks);
			mergeEdges(workDir, i, 1, nTasks);
			sortEdges(workDir, 1, i, 1, nTasks);

			makeGraph(workDir, i, 1, relable, nTasks); 

			delete [] map1;
			delete [] nodeMap;

		}

		char * fileName = new char [1000];
		sprintf(fileName, "%s/HybSpectrum/Graph0/tmpMap", workDir.c_str());

		FILE * pFile = fopen(fileName, "r");
		fseek(pFile, 0, SEEK_END);
		int mapSize = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);

		int * map = new int [mapSize];
		fread(map, sizeof(int), mapSize, pFile);		
		fclose(pFile);

		int relabel = clusterCount;
		for(int i = 0; i < mapSize; i++)
		{
			if(map[i] == -1)
			{
				map[i] = relabel++;
			}
		}

		pFile = fopen(fileName, "w");
		fwrite(map, sizeof(int), mapSize, pFile);
		fclose(pFile);

		delete [] fileName;

		makeEdges(workDir, 0, 2, map);
		sortEdges(workDir, 0, 0, 1, nTasks);
		mergeEdges(workDir, 0, 1, nTasks);
		sortEdges(workDir, 1, 0, 1, nTasks);

		makeGraph(workDir, 0, 1, relabel, nTasks);

		for(int i = numGraphs -1; i >= 0; i--)
		{ 
			int * contigs;
			int count = recoverHybClusters(workDir, contigs, i);
			outputClusters(workDir, contigs, count, i);
			delete [] contigs;
		}

		delete [] graphDir; int handshake = 0;
		for(int i = 1; i < nTasks; i++)
		{
			MPI_Send(&handshake, 1, MPI_INT, i, endTag, MPI_COMM_WORLD);
		}
	} 
	else{
		while(true)
		{
			MPI_Status status;
			int buffSize = 0; int handshake = 0; int graphIter = 0; int workType = 0;
			MPI_Recv(&workType, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			if(status.MPI_TAG == assembleTag)
			{
				MPI_Recv(&graphIter, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&buffSize, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				int * buff = new int [buffSize];
				MPI_Recv(buff, buffSize, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

				set<int> assembleSet;
				for(int i = 0; i < buffSize; i++)
					assembleSet.insert(buff[i]);

				buffSize = 0;

				int * cIndex;
				long int numNodes = recoverClusterIndex(workDir, cIndex, graphIter);

				int * contigs;
				int count = recoverClusters(workDir, cIndex, contigs, graphIter);

				char * fileName = new char [200];
				sprintf(fileName, "%s/OvlGraph/finalMap", workDir.c_str());
				FILE * pFile = fopen(fileName, "r");
				fseek(pFile, 0, SEEK_END);
				int mSize = ftell(pFile)/sizeof(int);	
				fseek(pFile, 0, SEEK_SET);				

				int * nodeMap = new int [mSize];
				fread(nodeMap, sizeof(int), mSize, pFile);
				fclose(pFile);	

				delete [] fileName;
				int currNode = 0;

				int gSize = 2; int graphsLoaded = 0; int currGraphSet = 0;
				Graph ** graphs = new Graph *[gSize];

				int iSize = 2; int indexesLoaded = 0; int currIndexSet = 0;
				Fragment_Index ** indexes = new Fragment_Index * [iSize];

				uint32_t * orderings = new uint32_t [2000];
				uint32_t * sandBox = new uint32_t [1100];

				uint32_t * nodes = new uint32_t [1000];
				uint32_t * forwardEdges = new uint32_t [5000];
				uint32_t * backEdges = new uint32_t [5000];

				int tally = 0;
				for(int i = 0; i < count; i++)
				{
					if(i == 0 || contigs[i] == -1)
					{	
						int nIndex = cIndex[currNode];
						if(assembleSet.find(nIndex) != assembleSet.end())
						{
							int startPos = i+1;
							if(i == 0)
								startPos--;

							int maxNode = 0; int minNode = INT_MAX; int numNodes = 0;
							for(int j = startPos; j < count && contigs[j] != -1; j++)
							{
								if(nodeMap[contigs[j]] <  minNode) minNode = nodeMap[contigs[j]];	
								if(nodeMap[contigs[j]] >  maxNode) maxNode = nodeMap[contigs[j]];	

								numNodes++;									
							} 

							if(numNodes > 1000)
							{
								delete [] orderings; orderings = new uint32_t [2*numNodes]; 
								delete [] sandBox; sandBox = new uint32_t [100+numNodes];
								delete [] nodes; nodes = new uint32_t [numNodes];
								delete [] forwardEdges; forwardEdges = new uint32_t [5*numNodes];
								delete [] backEdges; backEdges = new uint32_t [5*numNodes];
							}


							loadGraphs(workDir, graphs, gSize, graphsLoaded, currGraphSet, contigs, startPos, count, nodeMap);
							loadFragments(workDir, indexes, iSize, indexesLoaded, currIndexSet, contigs, startPos, count, nodeMap);


							removeContainments(graphs, graphsLoaded, contigs, startPos, count, nodeMap, minNode, maxNode);

							for(int j = startPos; j < count && contigs[j] != -1; j++)
							{
								int currNode = nodeMap[contigs[j]];
								Graph *  graph;
								for(int k = 0; k < graphsLoaded; k++)
								{
									graph = graphs[k];
									if(graph->getOffset() <= currNode && currNode < graph->getOffset() + graph->getNumNodes())
										break;
								}

								if(graph->nodeInGraph(currNode-graph->getOffset()))
								{
									nodes[currNode-minNode] = 1;
								}else{
									nodes[currNode-minNode] = 0;
								}

								forwardEdges[5*(currNode-minNode)] = 0;
								backEdges[5*(currNode-minNode)] = 0;
								sandBox[100+(currNode-minNode)] = 0;
							}

							transitiveReduce(graphs, graphsLoaded, nodes, forwardEdges, backEdges, contigs, startPos, count, nodeMap, minNode, maxNode);

							trimGraph(nodes, forwardEdges, backEdges, sandBox, minNode, maxNode, minPathLen);
							popBubbles(nodes, forwardEdges, backEdges, sandBox, minNode, maxNode,minPathLen);
							float numContigs = traverseGraph(nodes, forwardEdges, backEdges, orderings, sandBox, minNode, maxNode);
							if(numContigs < 4)
							{
								buff[buffSize++] = nIndex;
							}							

							if(numNodes > 1000)
							{
								delete [] orderings; orderings = new uint32_t [2000]; 
								delete [] sandBox; sandBox = new uint32_t [1100];
								delete [] nodes; nodes = new uint32_t [1000];
								delete [] forwardEdges; forwardEdges = new uint32_t [5000];
								delete [] backEdges; backEdges = new uint32_t [5000];
							}
						}

						currNode++;
					}
				}

				delete [] sandBox;
				delete [] orderings;

				delete [] nodes;
				delete [] forwardEdges;
				delete [] backEdges;

				delete [] contigs;
				delete [] nodeMap;
				delete [] cIndex;

				for(int i = 0; i < graphsLoaded; i++)
					delete graphs[i];

				delete [] graphs;

				for(int i = 0; i < indexesLoaded; i++)
					delete indexes[i];

				delete [] indexes;

				if(buffSize > 0)
				{
					MPI_Send(&buffSize, 1, MPI_INT, 0, workTag, MPI_COMM_WORLD);
					MPI_Recv(&handshake, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
					MPI_Send(buff, buffSize, MPI_INT, 0, workTag, MPI_COMM_WORLD);

					MPI_Recv(&handshake, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

					MPI_Send(&handshake, 1, MPI_INT, 0, workTag, MPI_COMM_WORLD);
				}

				delete [] buff;

				MPI_Send(&handshake, 1, MPI_INT, 0, endTag, MPI_COMM_WORLD);
			}

			if(status.MPI_TAG == sortTag)
			{
				int message [] = {0,0,0,0};
				MPI_Recv(message, 4, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

				int currFile = message[0]; 
				int setSize = message[1]; 
				int graphIter = message[2]; 
				int sortType = message[3];

				for(int i = 0; i < setSize; i++)
				{
					char * fileName = new char [1000];
					sprintf(fileName, "%s/HybSpectrum/Graph%d/set%d", workDir.c_str(), graphIter, currFile+i);

					SortEdges mySort(sortType);
					mySort.sort_edges(fileName, fileName);

					delete [] fileName;
				}

				MPI_Send(0, 0, MPI_INT, 0, status.MPI_TAG, MPI_COMM_WORLD);
			}

			if(status.MPI_TAG == mergeTag)
			{
				int message [] = {0,0,0,0};
				MPI_Recv(message, 4, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

				int currFile = message[0]; 
				int setSize = message[1]; 
				int graphIter = message[2]; 
				int sortType = message[3];

				char * inDir = new char [1000];
				char * outDir = new char [1000];

				sprintf(inDir, "%s/HybSpectrum/Graph%d", workDir.c_str(), graphIter);
				sprintf(outDir, "%s/HybSpectrum/Tmp", workDir.c_str());

				MergeEdges myMerge(sortType);
				myMerge.merge(currFile, setSize, inDir, outDir);

				delete [] inDir;
				delete [] outDir;

				MPI_Send(0, 0, MPI_INT, 0, status.MPI_TAG, MPI_COMM_WORLD);
			}

			if(status.MPI_TAG == graphTag)
			{	
				int message [] = {0,0,0,0};
				MPI_Recv(message, 4, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

				int currFile = message[0]; 
				int setSize = message[1]; 
				int graphIter = message[2]; 
				int totalNodes = message[3];

				for(int i = currFile; i < currFile + setSize; i++)
				{
					uint32_t offset = 0; uint32_t currNode = 0;
					if(i != 0)
					{
						char * fileName = new char [1000];
						sprintf(fileName,"%s/HybSpectrum/Graph%d/set%d", workDir.c_str(), graphIter,  i-1);
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
					sprintf(fileName,"%s/HybSpectrum/Graph%d/set%d", workDir.c_str(), graphIter,  i);

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

					sprintf(fileName,"%s/HybSpectrum/Graph%d/set%d", workDir.c_str(), graphIter,  i+1);

					if(!fileExists(fileName))
					{
						numNodes = totalNodes - currNode;
					}

					Graph graph;
					graph.reserveNumEdges(numEdges);
					graph.reserveNumNodes(numNodes);
					graph.setOffset(offset);
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

					sprintf(fileName,"%s/HybSpectrum/Graph%d/set%d", workDir.c_str(), graphIter,  i+1);

					if(!fileExists(fileName))
					{
						numNodes = totalNodes;

						graph.addSilentNodes(numNodes - currNode - 1);
					}

					graph.set();

					sprintf(fileName,"%s/HybSpectrum/Graph%d/graph%d", workDir.c_str(), graphIter,  i);

					pFile = fopen(fileName, "w");
					graph.write(pFile);
					fclose(pFile);
				}

				MPI_Send(0, 0, MPI_INT, 0, status.MPI_TAG, MPI_COMM_WORLD);
			}

			if(status.MPI_TAG == endTag)
				break;
		}
	}

	MPI_Finalize();

	return 0;
}

//void outputClusters(string &, int * &, int, int)
//Description: This outputs the fragments in a clustering to a file
//Input: workDir (string &):the working directory, contigs (int * &),
//count (int),
//Output:None
//Return:None
void outputClusters(string & workDir, int * & contigs, int count, int graphIter)
{
	int gSize = 2; int graphsLoaded = 0; int currGraphSet = 0;
	Graph ** graphs = new Graph *[gSize];

	int iSize = 2; int indexesLoaded = 0; int currIndexSet = 0;
	Fragment_Index ** indexes = new Fragment_Index * [iSize];

	int lSize = 2; int labelsLoaded = 0; int currLabelSet = 0;
	Labels ** labels = new Labels * [iSize];

	char * graphDir = new char [500]; 
	sprintf(graphDir, "%s/OvlGraph", workDir.c_str());

	char * fileName = new char [1000];
	sprintf(fileName, "%s/finalMap", graphDir);

	FILE * pFile = fopen(fileName, "r");
	fseek(pFile, 0, SEEK_END);

	int mapSize = ftell(pFile)/sizeof(int);
	fseek(pFile, 0, SEEK_SET);

	int * nodeMap = new int[mapSize];
	fread(nodeMap, sizeof(int), mapSize, pFile);
	fclose(pFile);
	char * contigDir = new char [500];
	sprintf(contigDir, "%s/Contigs", workDir.c_str());

	sprintf(fileName, "%s/contigs.%d", contigDir, graphIter);
	pFile = fopen(fileName, "w");

	char * buff = new char [10 * megabyte];
	int curr = 0;

	loadFragments(workDir, indexes, iSize, indexesLoaded, currIndexSet, contigs, 0, count, nodeMap);
	loadLabels(workDir, labels, lSize, labelsLoaded, currLabelSet, contigs, 0, count, nodeMap);

	for(int i = 0; i < count; i++)
	{
		if(contigs[i] != -1)
		{	
			int currNode = nodeMap[contigs[i]];
			Fragment_Index * index; Labels * label;
			for(int j = 0; j < indexesLoaded; j++)
			{
				index = indexes[j];
				if(currNode >= index->getIndexOffset() && currNode < index->getIndexOffset() + index->numFragments())
					break; 
			}			

			for(int j = 0; j < labelsLoaded; j++)
			{
				label = labels[j];
				if(currNode >= label->getOffset() && currNode < label->getOffset() + label->getNumLabels())
					break; 
			}			

			string fLabel = label->labelAt(currNode-label->getOffset());
			pair<unsigned long long int, unsigned long long int> bounds = index->indexBounds(currNode-index->getIndexOffset());
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
				buff[curr++] = index->at(j);
			}

			buff[curr++] = '\n';
		}else{

			if(curr >= 10 * megabyte-1)
			{
				fwrite(buff, sizeof(char), curr, pFile);
				curr = 0;
			}
			buff[curr++] = '*'; buff[curr++] = '\n';

			loadFragments(workDir, indexes, iSize, indexesLoaded, currIndexSet, contigs, i+1, count, nodeMap);
			loadLabels(workDir, labels, lSize, labelsLoaded, currLabelSet, contigs, i+1, count, nodeMap);
		}
	}

	if(curr > 0)
	{
		fwrite(buff, sizeof(char), curr, pFile);
	}

	fclose(pFile);

	delete [] buff;

	delete [] contigDir;
	delete [] graphDir;
	delete [] fileName;
}


//void cleanEdgeTmps(const char [], const char [])
//Description: This function cleans up the TMP and Graph directories
//Input: graphDir (const char []), tmpDir (const char [])
//Output:None
//Return:None
void cleanTmps(const char graphDir [], const char tmpDir [] ) {

	//Deleting/Moving directories
	char * fileName = new char[1000];
	sprintf(fileName, "%s/tmpMap", graphDir);
	char * dest = new char [1000];

	if(fileExists(fileName))
	{
		sprintf(dest, "%s/tmpMap", tmpDir);
		moveDir(fileName, dest);
	}

	delete [] dest;
	delete [] fileName;

	removeDir(graphDir); 
	moveDir(tmpDir, graphDir);
	makeDir(tmpDir);
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

			if(fileIndex != i)
			{
				system(fileName);
			}

			fileIndex++;
		}
	}

	delete [] fileName;
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

//void makeGraph()
//Description: This function builds graphs from the edges
//Input: workDir (string &): the working directory, graphIter (int), the current graph iteration
//nTasks (int) the number of available processors
//Output:None
//Return:None
void makeGraph(string & workDir, int graphIter, int auxFiles, int relabel, int nTasks) {

	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/HybSpectrum/Graph%d", workDir.c_str(), graphIter);

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
		message[2] = graphIter; message[3] = relabel;

		MPI_Send(&graphIter, 1, MPI_INT, taskNum, graphTag, MPI_COMM_WORLD);
		MPI_Send(message, 4, MPI_INT, taskNum++, workTag, MPI_COMM_WORLD);
	}

	//Receive/Send Tasks
	MPI_Status status; int result;
	while(currFile < numFiles)
	{
		message[0] = currFile; message[1] = setSize;
		message[2] = graphIter; message[3] = relabel;

		MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		MPI_Send(&graphIter, 1, MPI_INT, status.MPI_SOURCE, graphTag, MPI_COMM_WORLD);
		MPI_Send(message, 4, MPI_INT, status.MPI_SOURCE, workTag, MPI_COMM_WORLD);

		currFile+=setSize;
	}

	//Clean up
	finalizeTasks(taskNum);

	//Remove all the set files
	char * fileName = new char [1000];
	sprintf(fileName,"%s/HybSpectrum/Graph%d/set*", workDir.c_str(), graphIter);

	removeDir(fileName);
	delete [] fileName;
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
	sprintf(graphDir, "%s/HybSpectrum/Graph%d", workDir.c_str(), graphIter);

	char * tmpDir = new char [1000];
	sprintf(tmpDir, "%s/HybSpectrum/Tmp", workDir.c_str());

	makeDir(tmpDir);

	int numFiles =  countFiles(graphDir);
	numFiles-=auxFiles;

	//Determine Merge Type
	int mergeType = 1;

	//Seed Tasks
	int currFile = 0; int taskNum = 1; int setSize = 1; int message [] = {0, 0, 0, 0};
	for(; currFile < numFiles && taskNum < nTasks; currFile+=(setSize * 2/*10*/))
	{
		message[0] = currFile; message[1] = setSize;
		message[2] = graphIter; message[3] = mergeType;

		MPI_Send(&graphIter, 1, MPI_INT, taskNum, mergeTag, MPI_COMM_WORLD);
		MPI_Send(message, 4, MPI_INT, taskNum++, workTag, MPI_COMM_WORLD);
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
			MPI_Send(&graphIter, 1, MPI_INT, status.MPI_SOURCE, mergeTag, MPI_COMM_WORLD);
			MPI_Send(message, 4, MPI_INT, status.MPI_SOURCE, workTag, MPI_COMM_WORLD);

			currFile+=(setSize * 2/*10*/);

		}

		//Clean Up
		finalizeTasks(taskNum); cleanTmps(graphDir, tmpDir);
		setSize*=2/*10*/; currFile = 0; taskNum = 1;
		for(; setSize < numFiles && currFile < numFiles && taskNum < nTasks; currFile+=(setSize * 2/*10*/))
		{
			message[0] = currFile; message[1] = setSize;
			message[2] = graphIter; message[3] = mergeType;

			MPI_Send(&graphIter, 1, MPI_INT, taskNum, mergeTag, MPI_COMM_WORLD);
			MPI_Send(message, 4, MPI_INT, taskNum++, workTag, MPI_COMM_WORLD);
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


//void sortEdges()
//Description: This sorts a subset of files
//Input: workDir (string &): the working directory, graphIter (int), the current graph iteration
//nTasks (int) the number of available processors
//density
//Output:None
//Return:None
void sortEdges(string & workDir, int sortType, int graphIter, int auxFiles, int nTasks) {

	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/HybSpectrum/Graph%d", workDir.c_str(), graphIter);

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

		MPI_Send(&graphIter, 1, MPI_INT, taskNum, sortTag, MPI_COMM_WORLD);
		MPI_Send(message, 4, MPI_INT, taskNum++, workTag, MPI_COMM_WORLD);
	}

	//Receive/Send Tasks
	MPI_Status status; int result;
	while(currFile < numFiles)
	{
		message[0] = currFile; message[1] = setSize;
		message[2] = graphIter; message[3] = sortType;

		MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		MPI_Send(&graphIter, 1, MPI_INT, status.MPI_SOURCE, sortTag, MPI_COMM_WORLD);
		MPI_Send(message, 4, MPI_INT, status.MPI_SOURCE, workTag, MPI_COMM_WORLD);

		currFile+=setSize;
	}

	//Clean up
	finalizeTasks(taskNum);
}

//void makeEdges(string &, int, int * &)
//Description: This creates new edges for the collapsed graph
//Input: workDir (string &): the working directory, graphIter (int): The current
//graph iteration, nodeMap (int *): the mapping from the previous graph
//to the new one
//Output:None
//Return:None
void makeEdges(string & workDir, int graphIter, int auxFiles, int * & nodeMap) {

	char * iGraphDir = new char [1000];
	char * oGraphDir = new char [1000];
	char * fileName = new char [1000];

	sprintf(iGraphDir, "%s/Spectrum/Graph%d", workDir.c_str(), graphIter);
	sprintf(oGraphDir, "%s/HybSpectrum/Graph%d", workDir.c_str(), graphIter);

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
				buff[count++] = weight;
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

//void loadFragments()
//Description: This function loads the fragment indexes
//Input: graph (graph &): the overlap graph, pathLen (int): The length of the ends
//that can be trimmed.
//Output:None
//Return:None
void loadFragments(string & workDir, Fragment_Index ** & indexes, int & iSize, int & indexesLoaded, int & currIndexSet, int * & contigs, int startPos, int count, int * & nodeMap)
{
	for(int i = startPos; i < count && contigs[i] != -1; i++)
	{
		int currNode = nodeMap[contigs[i]];
		bool foundIndex = false;

		while(!foundIndex)
		{
			for(int j = 0; j < indexesLoaded; j++)
			{
				if(indexes[j]->getIndexOffset() <= currNode && currNode < indexes[j]->getIndexOffset() + indexes[j]->numFragments())
				{
					foundIndex = true; 
					if(i == startPos && j != 0)
					{
						for(int k = 0; k < j; k++)
							delete indexes[k];

						for(int k = 0; k+j < indexesLoaded; k++)
						{
							indexes[k] = indexes[k+j]; 
						}

						indexesLoaded = indexesLoaded - j;
					}
					break;		
				}
			}

			if(!foundIndex)
			{
				if(i == startPos)
				{
					for(int j = 0; j < indexesLoaded; j++)
						delete indexes[j];

					indexesLoaded = 0;
				}else{

					if(indexesLoaded == iSize)
					{
						Fragment_Index ** tmp = new Fragment_Index * [iSize++];
						for(int j = 0; j < indexesLoaded; j++)
							tmp[j] = indexes[j];

						delete [] indexes;

						indexes = tmp;
					}
				}

				char * fileName = new char [200];
				sprintf(fileName, "%s/Fragments/sequences.%d", workDir.c_str(), currIndexSet++);
				FILE * pFile = fopen(fileName, "r");

				unsigned long long int fragmentOffset;
				fread(&fragmentOffset, sizeof(unsigned long long int), 1, pFile);

				unsigned long long int indexOffset;
				fread(&indexOffset, sizeof(unsigned long long int), 1, pFile);

				unsigned long long int length;
				fread(&length, sizeof(unsigned long long int), 1, pFile);

				unsigned long long int numFragments;
				fread(&numFragments, sizeof(unsigned long long int), 1, pFile);
				fclose(pFile);

				pFile = fopen(fileName, "r");

				indexes[indexesLoaded] = new Fragment_Index();
				indexes[indexesLoaded]->reserve(numFragments, length, indexOffset, fragmentOffset);
				indexes[indexesLoaded]->read(pFile);
				indexes[indexesLoaded++]->set();

				fclose(pFile);
				delete [] fileName;
			}
		}
	}		
}

//void laodGraphs()
//Description: This function loads the overlap graphs
//Input: graph (graph &): the overlap graph, pathLen (int): The length of the ends
//that can be trimmed.
//Output:None
//Return:None
void loadLabels(string & workDir, Labels ** & labels, int & lSize, int & labelsLoaded, int & currLabelSet, int * & contigs, int startPos, int count, int * & nodeMap)
{
	for(int i = startPos; i < count && contigs[i] != -1; i++)
	{
		int currNode = nodeMap[contigs[i]];
		bool foundLabel = false;

		while(!foundLabel)
		{
			for(int j = 0; j < labelsLoaded; j++)
			{
				if(labels[j]->getOffset() <= currNode && currNode < labels[j]->getOffset() + labels[j]->getNumLabels())
				{
					foundLabel = true; 
					if(i == startPos && j != 0)
					{
						for(int k = 0; k < j; k++)
							delete labels[k];

						for(int k = 0; k+j < labelsLoaded; k++)
						{
							labels[k] = labels[k+j]; 
						}

						labelsLoaded = labelsLoaded - j;
					}

					break;		
				}
			}

			if(!foundLabel)
			{
				if(i == startPos)
				{
					for(int j = 0; j < labelsLoaded; j++)
						delete labels[j];

					labelsLoaded = 0;
				}else{

					if(labelsLoaded == lSize)
					{
						Labels ** tmp = new Labels * [lSize++];
						for(int j = 0; j < labelsLoaded; j++)
							tmp[j] = labels[j];

						delete [] labels;

						labels = tmp;
					}
				}

				char * fileName = new char [200];
				sprintf(fileName, "%s/Labels/labels.%d", workDir.c_str(), currLabelSet++);
				FILE * pFile = fopen(fileName, "r");

				unsigned int tagLen; unsigned long long int offset;
				unsigned long long int numLabels;

				fread(&tagLen, sizeof(unsigned int), 1, pFile);
				fread(&offset, sizeof(unsigned long long int), 1, pFile);
				fread(&numLabels, sizeof(unsigned long long int), 1, pFile);
				fseek(pFile, 0, SEEK_SET);

				labels[labelsLoaded] = new Labels();
				labels[labelsLoaded]->reserve(tagLen, numLabels, offset);
				labels[labelsLoaded++]->read(pFile);

				fclose(pFile);
				delete [] fileName;
			}
		}
	}		
}


//void laodGraphs()
//Description: This function loads the overlap graphs
//Input: graph (graph &): the overlap graph, pathLen (int): The length of the ends
//that can be trimmed.
//Output:None
//Return:None
void loadGraphs(string & workDir, Graph ** & graphs, int & gSize, int & graphsLoaded, int & currGraphSet, int * & contigs, int startPos, int count, int * & nodeMap)
{
	for(int i = startPos; i < count && contigs[i] != -1; i++)
	{
		int currNode = nodeMap[contigs[i]];
		bool foundGraph = false;

		while(!foundGraph)
		{
			for(int j = 0; j < graphsLoaded; j++)
			{
				if(graphs[j]->getOffset() <= currNode && currNode < graphs[j]->getOffset() + graphs[j]->getNumNodes())
				{
					foundGraph = true; 
					if(i == startPos && j != 0)
					{
						for(int k = 0; k < j; k++)
							delete graphs[k];

						for(int k = 0; k+j < graphsLoaded; k++)
						{
							graphs[k] = graphs[k+j]; 
						}

						graphsLoaded = graphsLoaded - j;
					}

					break;		
				}
			}

			if(!foundGraph)
			{
				if(i == startPos)
				{
					for(int j = 0; j < graphsLoaded; j++)
						delete graphs[j];

					graphsLoaded = 0;
				}else{

					if(graphsLoaded == gSize)
					{
						Graph ** tmp = new Graph * [gSize++];
						for(int j = 0; j < graphsLoaded; j++)
							tmp[j] = graphs[j];

						delete [] graphs;

						graphs = tmp;
					}
				}

				char * fileName = new char [200];
				sprintf(fileName, "%s/OvlGraph/graph%d", workDir.c_str(), currGraphSet++);
				FILE * pFile = fopen(fileName, "r");

				graphs[graphsLoaded] = new Graph();
				graphs[graphsLoaded++]->read(pFile);

				fclose(pFile);
				delete [] fileName;
			}
		}
	}		
}

//void trimGraph()
//Description: This function trims the overlap graph
//Input: graph (graph &): the overlap graph, pathLen (int): The length of the ends
//that can be trimmed.
//Output:None
//Return:None
int traverseGraph(uint32_t * & nodes, uint32_t * & forwardEdges, uint32_t * & backEdges, uint32_t * & orderings, uint32_t * & sandBox, int minNode, int maxNode)
{	
	int numContigs = 0;
	int startPos = 0; orderings[startPos] = 0;
	for(int i = 0; i < maxNode-minNode+1; i++)
	{
		if(nodes[i] && !sandBox[i+100])
		{
			orderings[startPos+1] = i+minNode;			
			orderings[startPos]++;	
			sandBox[100+i] = 1;

			int fIndex = 0; int bIndex = 0;
			int fTally = 0; int bTally = 0;
			for(int j = 0; j < forwardEdges[5*i]; j++)
			{
				int n = forwardEdges[5*i+j+1];
				if(nodes[n-minNode])
				{
					fIndex = j; fTally++;
				}
			}

			for(int j = 0; j < backEdges[5*i]; j++)
			{
				int n = backEdges[5*i+j+1];
				if(nodes[n-minNode])
				{
					bIndex = j; bTally++;
				}
			}

			if(fTally == 1)
			{		
				int nextNode = forwardEdges[5*i+1+fIndex];

				fTally = 0; bTally = 0;
				for(int j = 0; j < forwardEdges[5*(nextNode-minNode)]; j++)
				{
					int n = forwardEdges[5*(nextNode-minNode)+j+1];
					if(nodes[n-minNode])
					{
						fIndex = j; fTally++;
					}
				}

				for(int j = 0; j < backEdges[5*(nextNode-minNode)]; j++)
				{
					int n = backEdges[5*(nextNode-minNode)+j+1];
					if(nodes[n-minNode])
					{
						bIndex = j; bTally++;
					}
				}

				while(!sandBox[(nextNode-minNode)+100] && fTally == 1 && bTally == 1)
				{
					int currPos = orderings[startPos];
					orderings[startPos+currPos+1] = nextNode;
					sandBox[(nextNode-minNode)+100] = 1;
					nextNode = forwardEdges[5*(nextNode-minNode)+1+fIndex];

					fTally = 0; bTally = 0;
					for(int j = 0; j < forwardEdges[5*(nextNode-minNode)]; j++)
					{
						int n = forwardEdges[5*(nextNode-minNode)+j+1];
						if(nodes[n-minNode])
						{
							fIndex = j; fTally++;
						}
					}

					for(int j = 0; j < backEdges[5*(nextNode-minNode)]; j++)
					{
						int n = backEdges[5*(nextNode-minNode)+j+1];
						if(nodes[n-minNode])
						{
							bIndex = j; bTally++;
						}
					}

					orderings[startPos]++;	
				}

				if(!sandBox[(nextNode-minNode)+100] && bTally == 1)
				{
					int currPos = orderings[startPos];
					orderings[startPos+currPos+1] = nextNode;
					sandBox[(nextNode-minNode)+100] = 1;
					orderings[startPos]++;	
				}
			}

			startPos = orderings[startPos] + startPos + 1;		
			orderings[startPos] = 0;			

			fTally = 0; bTally = 0;
			for(int j = 0; j < forwardEdges[5*i]; j++)
			{
				int n = forwardEdges[5*i+j+1];
				if(nodes[n-minNode])
				{
					fIndex = j; fTally++;
				}
			}

			for(int j = 0; j < backEdges[5*i]; j++)
			{
				int n = backEdges[5*i+j+1];
				if(nodes[n-minNode])
				{
					bIndex = j; bTally++;
				}
			}

			if(bTally == 1)
			{
				int nextNode = backEdges[5*i+1+bIndex];
				fTally = 0; bTally = 0;
				for(int j = 0; j < forwardEdges[5*(nextNode-minNode)]; j++)
				{
					int n = forwardEdges[5*(nextNode-minNode)+j+1];
					if(nodes[n-minNode])
					{
						fIndex = j; fTally++;
					}
				}

				for(int j = 0; j < backEdges[5*(nextNode-minNode)]; j++)
				{
					int n = backEdges[5*(nextNode-minNode)+j+1];
					if(nodes[n-minNode])
					{
						bIndex = j; bTally++;
					}
				}
				while(!sandBox[(nextNode-minNode)+100] && fTally == 1 && bTally == 1)
				{
					int currPos = orderings[startPos];
					orderings[startPos+currPos+1] = nextNode;
					sandBox[(nextNode-minNode)+100] = 1;
					nextNode = backEdges[5*(nextNode-minNode)+1+bIndex];

					fTally = 0; bTally = 0;
					for(int j = 0; j < forwardEdges[5*(nextNode-minNode)]; j++)
					{
						int n = forwardEdges[5*(nextNode-minNode)+j+1];
						if(nodes[n-minNode])
						{
							fIndex = j; fTally++;
						}
					}

					for(int j = 0; j < backEdges[5*(nextNode-minNode)]; j++)
					{
						int n = backEdges[5*(nextNode-minNode)+j+1];
						if(nodes[n-minNode])
						{
							bIndex = j; bTally++;
						}
					}
					orderings[startPos]++;	
				}

				if(!sandBox[(nextNode-minNode)+100] && fTally == 1)
				{
					int currPos = orderings[startPos];
					orderings[startPos+currPos+1] = nextNode;
					sandBox[(nextNode-minNode)+100] = 1;
					orderings[startPos]++;	
				}
			}

			startPos = orderings[startPos] + startPos + 1;		
			orderings[startPos] = 0;			
			numContigs++;
		}		
	}

	return numContigs;
}

//void trimGraph()
//Description: This function trims the overlap graph
//Input: graph (graph &): the overlap graph, pathLen (int): The length of the ends
//that can be trimmed.
//Output:None
//Return:None
int trimGraph(uint32_t * & nodes, uint32_t * & forwardEdges, uint32_t * & backEdges, uint32_t * & sandBox, int minNode, int maxNode, int pathLen)
{
	int numTrimmed = 0; int numNodesInGraph = 0;
	for(int i = 0; i < maxNode-minNode+1; i++)
	{
		if(nodes[i])
		{
			int numExplored = 0; numNodesInGraph++;
			int fIndex = 0; int bIndex = 0;
			int fTally = 0; int bTally = 0;
			for(int j = 0; j < forwardEdges[5*i]; j++)
			{
				int n = forwardEdges[5*i+j+1];
				if(nodes[n-minNode])
				{
					fIndex = j; fTally++;
				}
			}

			for(int j = 0; j < backEdges[5*i]; j++)
			{
				int n = backEdges[5*i+j+1];
				if(nodes[n-minNode])
				{
					bIndex = j; bTally++;
				}
			}

			if(fTally == 0 && bTally == 1)
			{
				sandBox[numExplored++] = i + minNode;
				int nextNode = backEdges[5*i+bIndex+1];
				fTally = 0; bTally = 0;

				for(int j = 0; j < forwardEdges[5*(nextNode-minNode)]; j++)
				{
					int n = forwardEdges[5*(nextNode-minNode)+j+1];
					if(nodes[n-minNode])
					{
						fIndex = j; fTally++;
					}
				}

				for(int j = 0; j < backEdges[5*(nextNode-minNode)]; j++)
				{
					int n = backEdges[5*(nextNode-minNode)+j+1];
					if(nodes[n-minNode])
					{
						bIndex = j; bTally++;
					}
				}


				while(fTally == 1 && bTally == 1)
				{
					sandBox[numExplored++] = nextNode;
					nextNode = backEdges[5*(nextNode-minNode)+ bIndex + 1];

					fTally = 0; bTally = 0;
					for(int j = 0; j < forwardEdges[5*(nextNode-minNode)]; j++)
					{
						int n = forwardEdges[5*(nextNode-minNode)+j+1];
						if(nodes[n-minNode])
						{
							fIndex = j; fTally++;
						}
					}

					for(int j = 0; j < backEdges[5*(nextNode-minNode)]; j++)
					{
						int n = backEdges[5*(nextNode-minNode)+j+1];
						if(nodes[n-minNode])
						{
							bIndex = j; bTally++;
						}
					}


					if(numExplored > pathLen)
						break;
				}

				fTally = 0; bTally = 0;
			}

			if(fTally == 1 && bTally == 0)
			{
				sandBox[numExplored++] = i + minNode;
				int nextNode = forwardEdges[5*i+1];
				fTally = 0; bTally = 0;

				for(int j = 0; j < forwardEdges[5*(nextNode-minNode)]; j++)
				{
					int n = forwardEdges[5*(nextNode-minNode)+j+1];
					if(nodes[n-minNode])
					{
						fIndex = j; fTally++;
					}
				}

				for(int j = 0; j < backEdges[5*(nextNode-minNode)]; j++)
				{
					int n = backEdges[5*(nextNode-minNode)+j+1];
					if(nodes[n-minNode])
					{
						bIndex = j; bTally++;
					}
				}

				while(fTally == 1 && bTally  == 1)
				{
					sandBox[numExplored++] = nextNode;
					nextNode = forwardEdges[5*(nextNode-minNode)+fIndex + 1];

					fTally = 0; bTally = 0;
					for(int j = 0; j < forwardEdges[5*(nextNode-minNode)]; j++)
					{
						int n = forwardEdges[5*(nextNode-minNode)+j+1];
						if(nodes[n-minNode])
						{
							fIndex = j; fTally++;
						}
					}

					for(int j = 0; j < backEdges[5*(nextNode-minNode)]; j++)
					{
						int n = backEdges[5*(nextNode-minNode)+j+1];
						if(nodes[n-minNode])
						{
							bIndex = j; bTally++;
						}
					}

					if(numExplored > pathLen)
						break;
				}
			}

			float p1 = numExplored; float g1 = numNodesInGraph;
			if(numExplored < pathLen && p1/g1 < .5)
			{
				for(int j = 0; j < numExplored; j++)
				{
					nodes[sandBox[j]-minNode] = 0;
					numTrimmed++;
				}					
			}
		}
	}

	return numTrimmed;
}


//void popBubbles()
//Description: This function makes the graph
//Input: cMap (map<uint32_t, uint32_t>): the map from containments, graphArr (Graph **), pointers to 
//the graphs, numGraphs (int): the number of graphs that are loaded, graph (Graph &) the graph we are 
//creating.  
//Output:None
//Return:None
int popBubbles(uint32_t * & nodes, uint32_t * & forwardEdges, uint32_t * & backEdges, uint32_t * & sandBox, int minNode, int maxNode, int pathLen)
{
	int numPopped = 0;
	for(int i = 0; i < maxNode - minNode + 1; i++)
	{
		if(nodes[i])
		{
			if(forwardEdges[i*5] == 2)
			{
				int node1 = forwardEdges[i*5+1];
				int node2 = forwardEdges[i*5+2];

				if(nodes[node1-minNode] && nodes[node2-minNode])
				{
					int path1 = 0; int path2 = 2 * (pathLen+1);
					sandBox[path1++] = node1;
					sandBox[path2++] = node2;
					sandBox[100+node1-minNode] = 1;
					sandBox[100+node2-minNode] = 1;

					int fIndex1 = 0; int bIndex1 = 0;
					int fTally1 = 0; int bTally1 = 0;
					for(int j = 0; j < forwardEdges[5*(node1-minNode)]; j++)
					{
						int n = forwardEdges[5*(node1-minNode)+j+1];
						if(nodes[n-minNode])
						{
							fIndex1 = j; fTally1++;
						}
					}

					for(int j = 0; j < backEdges[5*(node1-minNode)]; j++)
					{
						int n = backEdges[5*(node1-minNode)+j+1];
						if(nodes[n-minNode])
						{
							bIndex1 = j; bTally1++;
						}
					}

					int fIndex2 = 0; int bIndex2 = 0;
					int fTally2 = 0; int bTally2 = 0;
					for(int j = 0; j < forwardEdges[5*(node2-minNode)]; j++)
					{
						int n = forwardEdges[5*(node2-minNode)+j+1];
						if(nodes[n-minNode])
						{
							fIndex2 = j; fTally2++;
						}
					}

					for(int j = 0; j < backEdges[5*(node2-minNode)]; j++)
					{
						int n = backEdges[5*(node2-minNode)+j+1];
						if(nodes[n-minNode])
						{
							bIndex2 = j; bTally2++;
						}
					}

					while(fTally1 == 1 && fTally2 == 1 && bTally1 == 1 && bTally2 == 1)
					{
						node1 = forwardEdges[5*(node1-minNode)+1+fIndex1];
						node2 = forwardEdges[5*(node2-minNode)+1+fIndex2];

						fTally1 = 0; bTally1 = 0;
						for(int j = 0; j < forwardEdges[5*(node1-minNode)]; j++)
						{
							int n = forwardEdges[5*(node1-minNode)+j+1];
							if(nodes[n-minNode])
							{
								fIndex1 = j; fTally1++;
							}
						}

						for(int j = 0; j < backEdges[5*(node1-minNode)]; j++)
						{
							int n = backEdges[5*(node1-minNode)+j+1];
							if(nodes[n-minNode])
							{
								bIndex1 = j; bTally1++;
							}
						}

						fTally2 = 0; bTally2 = 0;
						for(int j = 0; j < forwardEdges[5*(node2-minNode)]; j++)
						{
							int n = forwardEdges[5*(node2-minNode)+j+1];
							if(nodes[n-minNode])
							{
								fIndex2 = j; fTally2++;
							}
						}

						for(int j = 0; j < backEdges[5*(node2-minNode)]; j++)
						{
							int n = backEdges[5*(node2-minNode)+j+1];
							if(nodes[n-minNode])
							{
								bIndex2 = j; bTally2++;
							}
						}

						if(sandBox[100+(node1-minNode)] == 1)
						{
							for(int j = path1-1; j >= 0; j--)
							{
								nodes[sandBox[j]-minNode] = 0;
								numPopped++;
							}

							break;
						}

						if(sandBox[100+(node2-minNode)] == 1)
						{
							for(int j = path2-1; j >= 2*(pathLen+1); j--)
							{
								nodes[sandBox[j]-minNode] = 0;
								numPopped++;
							}

							break;
						}

						sandBox[path1++] = node1;
						sandBox[path2++] = node2;

						sandBox[100+(node1-minNode)] = 1;
						sandBox[100+(node2-minNode)] = 1;

						if(path1 > pathLen || path2 > pathLen)
							break;
					}

					for(int j = 0; j < path1; j++)
					{
						nodes[sandBox[j]-minNode] = 0;
					}

					for(int j = 2*(pathLen+1); j < path2; j++)
					{
						nodes[sandBox[j]-minNode] = 0;
					}
				}
			}
		}
	}	

	return numPopped;
}

//void transitiveReduce()
//Description: This function transitively reduces a graph
//Input: 
//Output:None
//Return:None
void transitiveReduce(Graph ** & graphs, int graphsLoaded, uint32_t * & nodes, uint32_t * & forwardEdges, uint32_t * & backEdges, int * & contigs, int startPos, int count, int * & nodeMap, int minNode, int maxNode)
{
	for(int i = startPos; i < count && contigs[i] != -1; i++)
	{
		uint32_t node1 = nodeMap[contigs[i]];
		if(nodes[node1-minNode])
		{
			Graph * graph1;
			for(int l = 0; l < graphsLoaded; l++)
			{
				graph1 = graphs[l];
				if(graph1->getOffset() <= node1 && node1 < graph1->getOffset() + graph1->getNumNodes())
					break;
			}

			for(int j = graph1->nodeDegree(node1-graph1->getOffset())-1; j >= 0; j--)
			{
				uint32_t node3 = graph1->edgeDest(node1-graph1->getOffset(), j);

				if(!graph1->rDovetail(node1-graph1->getOffset(), j) && node3 >= minNode && node3 <= maxNode && nodes[node3-minNode] )
				{
					bool reduced = false;
					for(int k = 0; k < graph1->nodeDegree(node1-graph1->getOffset()); k++)
					{
						uint32_t node2 = graph1->edgeDest(node1-graph1->getOffset(), k);
						if(!graph1->rDovetail(node1-graph1->getOffset(), k) && node2 >= minNode && node2 <= maxNode && nodes[node2-minNode])
						{
							Graph * graph2;
							for(int l = 0; l < graphsLoaded; l++)
							{
								graph2 = graphs[l];
								if(graph2->getOffset() <= node2 && node2 < graph2->getOffset() + graph2->getNumNodes())
									break;
							}

							for(int h = 0; h < graph2->nodeDegree(node2-graph2->getOffset()); h++)
							{
								if(!graph2->rDovetail(node2-graph2->getOffset(), h) && graph2->edgeDest(node2-graph2->getOffset(), h) == node3)
								{
									reduced = true; break;
								}
							}

						}

						if(reduced)
							break;
					}

					if(!reduced && node3 >= minNode && node3 <= maxNode && nodes[node3-minNode])
					{
						int numEdges1 = forwardEdges[5*(node1-minNode)];
						int numEdges2 = backEdges[5*(node3-minNode)];

						if(numEdges1 < 4 && numEdges2 < 4)
						{
							int pos1 = numEdges1+5*(node1-minNode)+1;
							int pos2 = numEdges2+5*(node3-minNode)+1;
							forwardEdges[pos1] = node3;
							backEdges[pos2] = node1;

							forwardEdges[5*(node1-minNode)]++;
							backEdges[5*(node3-minNode)]++;
						}	
					}
				}				
			}
		}
	}
}

//void removeContainments()
//Description: This function removes containment fragments from the graph
//Input: 
//Output:None
//Return:None
void removeContainments(Graph ** & graphs, int graphsLoaded, int * & contigs, int startPos, int count, int * & nodeMap, int minNode, int maxNode)
{
	for(int i = startPos; i < count && contigs[i] != -1; i++)
	{
		int currNode = nodeMap[contigs[i]];
		Graph * graph  = '\0';
		for(int j = 0; j < graphsLoaded; j++)
		{
			graph = graphs[j];
			if(graph->getOffset() <= currNode && currNode < graph->getOffset() + graph->getNumNodes())
				break;
		}

		int source = currNode-graph->getOffset();		
		for(int j = 0; j < graph->nodeDegree(source); j++)
		{
			int dest = graph->edgeDest(source, j);
			Graph * graph2  = '\0';
			for(int k = 0; k < graphsLoaded; k++)
			{
				graph2 = graphs[k];
				if(graph2->getOffset() <= dest && dest < graph->getOffset() + graph->getNumNodes())
					break;
			}

			if(dest >= minNode && dest <=  maxNode)
			{
				if(graph->sContained(source, j) && graph->dContained(source, j))
				{
					if(graph2->nodeDegree(dest-graph2->getOffset()) > graph->nodeDegree(source))
					{
						graph->maskNodeFromGraph(source);
						break;
					}else{

						if(graph2->nodeDegree(dest-graph2->getOffset()) == graph->nodeDegree(source) && currNode < dest)
						{
							graph->maskNodeFromGraph(source);
							break;

						}
					}

				}else{
					if(graph->sContained(source, j) && !graph->dContained(source, j))
					{
						graph->maskNodeFromGraph(source);
						break;
					} 
				}	
			}	
		}
	}
}

//void initHybSpec(string & workDir)
//Description: This function initiates the hybrid graph spectrum directory
//Input: workDir (string &): the name of the working directory
//Output:None
//Return:None
void initHybSpec(string & workDir)
{
	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Spectrum", workDir.c_str());
	int numFiles = countFiles(graphDir);

	for(int i = 0; i < numFiles; i++)
	{
		sprintf(graphDir, "%s/HybSpectrum/Graph%d", workDir.c_str(), i);
		makeDir(graphDir);

		sprintf(graphDir, "%s/Spectrum/Graph%d", workDir.c_str(), i);

		int numFiles = countFiles(graphDir) - 2;

		if(i > 0)
			numFiles--;

		long int numNodes = 0;

		Graph graph;
		char * fileName = new char [1000];
		sprintf(fileName, "%s/graph%d", graphDir, numFiles-1);
		FILE * pFile = fopen(fileName, "r");

		graph.read(pFile);
		fclose(pFile);

		numNodes = graph.getOffset() + graph.getNumNodes();

		int * map = new int [numNodes];
		for(int j = 0; j < numNodes; j++)
			map[j] = -1;


		sprintf(graphDir, "%s/HybSpectrum/Graph%d", workDir.c_str(), i);

		sprintf(fileName, "%s/tmpMap", graphDir);
		pFile = fopen(fileName, "w");
		fwrite(map, sizeof(int), numNodes, pFile);
		fclose(pFile);

		delete [] map;
		delete [] fileName;
	}

	delete [] graphDir;	
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

	char * sysCmd = new char [1000];
	sprintf(sysCmd, "rm -r %s", dirName);
	system(sysCmd);
	delete [] sysCmd;
}


//void recoverClusters()
//Description: This runs the the parent thread in a parallel program
//Input: workDir (string &): the working directory, minMerge (int):
//the minimum overlap length for merging, minDensity (double): The minimum inter-node
////graphLevel (int), level from which the clusters are being recovered
//density
//Output:None
//Return:None
int recoverHybClusters(string & workDir, int * & contigs, int graphIter) {

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
	contigs = new int [2 * numNodes];

	char * hybGraphDir = new char [1000];
	sprintf(hybGraphDir, "%s/HybSpectrum", workDir.c_str());

	int numGraphs = countFiles(graphDir);

	int count = 0;
	for(int i = numGraphs-1; i > 0; i--)
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

		sprintf(fileName, "%s/Graph%d/tmpMap", hybGraphDir, i-1);
		pFile = fopen(fileName, "r");

		fseek(pFile, 0, SEEK_END);
		int hybMapSize = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);

		int * hybMap = new int [hybMapSize];
		fread(hybMap, sizeof(int), hybMapSize, pFile);
		fclose(pFile);

		if(i == numGraphs-1)
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
			if(contigs[j] != -1)
			{
				tmp[count++] = buff[contigs[j]*2];

				if(i > graphIter && buff[contigs[j]*2+1] != -1)
					if(hybMap[buff[contigs[j]*2+1]] != hybMap[buff[contigs[j] * 2]])
						tmp[count++] = -1;

				if(buff[contigs[j]*2+1] != -1)
					tmp[count++] = buff[contigs[j]*2+1];
			}else{
				tmp[count++] = -1;
			}
		}

		delete [] contigs;

		contigs = tmp;

		delete [] buff;
		delete [] hybMap;
		delete [] fileName;
	}

	delete [] graphDir;
	return count;
}
//void recoverClusters()
//Description: This runs the the parent thread in a parallel program
//Input: workDir (string &): the working directory, minMerge (int):
//the minimum overlap length for merging, minDensity (double): The minimum inter-node
////graphLevel (int), level from which the clusters are being recovered
//density
//Output:None
//Return:None
int recoverClusters(string & workDir, int * & contigIndex, int * & contigs, int graphIter) {

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
	contigs = new int [2 * numNodes];

	int count = 0;
	for(int i = graphIter; i > 0; i--)
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

		if(i == graphIter)
		{
			for(int j = 0; j < buffSize; j+=2)
			{
				contigs[j] = contigIndex[j/2];
				contigs[j+1] = -1;
			}
			count = buffSize;
		}

		int * tmp = new int [2 * numNodes];

		int cHold = count; count = 0;
		for(int j = 0; j < cHold; j++)
		{
			if(contigs[j] != -1)
			{
				tmp[count++] = buff[contigs[j]*2];

				if(graphIter < i && buff[contigs[j]*2+1] != -1)
					tmp[count++] = -1;

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


//void recoverClusterIndex()
//Description: This recovers the indexes of the nodes at the current level
//to which the clusters belong 
//Input: workDir (string &): the working directory, minMerge (int):
//the minimum overlap length for merging, minDensity (double): The minimum inter-node
////graphLevel (int), level from which the clusters are being recovered
//density
//Output:None
//Return:None
int recoverClusterIndex(string & workDir, int * & contigIndex, int graphIter){

	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Spectrum/Graph%d", workDir.c_str(), graphIter);

	int numFiles = countFiles(graphDir) - 3;
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

	contigIndex = new int [numNodes];

	sprintf(graphDir, "%s/Spectrum", workDir.c_str());
	numFiles = countFiles(graphDir);

	int count = 0;
	if(graphIter == (numFiles-1))
	{
		for(int i = 0; i < numNodes; i++)
		{
			contigIndex[i] = i;
		}

		count = numNodes;
	}

	for(int i = numFiles-1; i > graphIter; i--)
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
			for(int j = 0; j < buffSize/2; j++)
			{
				contigIndex[j] = j;
			}
			count = buffSize/2;
		}

		int * tmp = new int [numNodes];

		int cHold = count; count = 0;
		for(int j = 0; j < cHold; j++)
		{
			tmp[count++] = buff[contigIndex[j]*2];

			if(buff[contigIndex[j]*2+1] != -1)
				tmp[count++] = buff[contigIndex[j]*2+1];
		}

		delete [] contigIndex;

		contigIndex = tmp;

		delete [] buff;
		delete [] fileName;
	}

	delete [] graphDir;
	return count;
}

