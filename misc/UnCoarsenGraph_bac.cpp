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
#include "Dictionary.h"
#include "Fragment_Index.h"
#include "Labels.h"
#include "Graph.h"
#include "BandedAlignment.h"
#include "Needleman_Wunsch.h"
#include "Consensus.h"
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
void removeContainments(Graph ** &, int, int * &, int, int, int * &, uint32_t * &, int, int);
void bestOverlaps(Graph ** &, int, int * &, int, int, int * &, uint32_t * &, uint32_t * &, int, int);
void transitiveReduce(Graph ** &, int, uint32_t * &, uint32_t * &, uint32_t * &, uint32_t * &, int * &, int, int, int * &, int, int);
int trimGraph(uint32_t * &, uint32_t * &, uint32_t * &, uint32_t * &, int, int, int, int * &);
int popBubbles(uint32_t * &, uint32_t * &, uint32_t * &, uint32_t * &, int, int, int, int * &);
int traverseGraph(uint32_t * &, uint32_t * &, uint32_t * &, uint32_t * &, uint32_t * &, int, int);
void makeEdges(string &, int, int, int * &);
void makeFinalEdges(string &, int, int, int * &, int * &, int * &);
void sortEdges(string &, int, int, int, int);
void mergeEdges(string &, int, int, int);
void makeGraph(string &, int, int, int, int);
void finalizeTasks(const int);
void cleanTmps(const char [], const char []);
void removeDir(const char []);
void moveDir(const char [], const char []);
void removeFile(const char []);
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
	string timeS = "--time"; //Report the time

	string workDir = "";
	int minPathLen = 4;
	int minBubbleLen = 4;
	bool sam = false;
	bool verbose = false;	
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

		//Display time
		if(argv[i] == timeS && reportTime != true)
		{
			reportTime = true;
			i--;
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
					char * iFileName  = new char [1000];
					char * oFileName = new char [1000];

					sprintf(iFileName, "%s/ContigsInfo/cons.%d", workDir.c_str(), status.MPI_SOURCE);
					sprintf(oFileName, "%s/ContigsInfo/Contigs.store", workDir.c_str());

					FILE * iFile = fopen(iFileName, "r");
					fseek(iFile, 0, SEEK_END);
					buffSize = ftell(iFile);
					fseek(iFile, 0, SEEK_SET);

					char * cBuff = new char[buffSize];
					fread(cBuff, sizeof(char), buffSize, iFile);
					fclose(iFile);

					removeFile(iFileName);

					FILE * oFile = fopen(oFileName, "a");

					fseek(oFile, 0, SEEK_END);
					int numPrevCons = ftell(oFile);
					fseek(oFile, 0, SEEK_SET);

					fwrite(cBuff, sizeof(char), buffSize, oFile);
					fclose(oFile);

					sprintf(iFileName, "%s/ContigsInfo/flags.%d", workDir.c_str(), status.MPI_SOURCE);
					sprintf(oFileName, "%s/ContigsInfo/Flags.store", workDir.c_str());

					iFile = fopen(iFileName, "r");
					fread(cBuff, sizeof(char), buffSize, iFile);
					fclose(iFile);

					removeFile(iFileName);

					oFile = fopen(oFileName, "a");
					fwrite(cBuff, sizeof(char), buffSize, oFile);
					fclose(oFile);

					delete [] cBuff;			

					sprintf(iFileName, "%s/ContigsInfo/contigBounds.%d", workDir.c_str(), status.MPI_SOURCE);
					sprintf(oFileName, "%s/ContigsInfo/contigBounds.store", workDir.c_str());

					iFile = fopen(iFileName, "r");
					fseek(iFile, 0, SEEK_END);
					buffSize = ftell(iFile)/sizeof(int);
					fseek(iFile, 0, SEEK_SET);

					buff = new int[buffSize];
					fread(buff, sizeof(int), buffSize, iFile);
					fclose(iFile);

					removeFile(iFileName);

					for(int j = 0; j < buffSize; j++)
					{
						buff[j]+=numPrevCons;
					}

					oFile = fopen(oFileName, "a");
					fwrite(buff, sizeof(int), buffSize, oFile);
					fclose(oFile);

					delete [] buff;

					sprintf(iFileName, "%s/ContigsInfo/readBounds.%d", workDir.c_str(), status.MPI_SOURCE);
					sprintf(oFileName, "%s/ContigsInfo/readBounds.store", workDir.c_str());

					iFile = fopen(iFileName, "r");
					fseek(iFile, 0, SEEK_END);
					buffSize = ftell(iFile)/sizeof(int);
					fseek(iFile, 0, SEEK_SET);

					buff = new int[buffSize];
					fread(buff, sizeof(int), buffSize, iFile);
					fclose(iFile);

					removeFile(iFileName);

					oFile = fopen(oFileName, "a");
					fwrite(buff, sizeof(int), buffSize, oFile);
					fclose(oFile);

					delete [] buff;

					delete [] iFileName; 
					delete [] oFileName;

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

			sprintf(fileName, "%s/Spectrum/Graph%d/nDen", workDir.c_str(), i);
			pFile = fopen(fileName, "r");			

			fseek(pFile, 0, SEEK_END);
			int size = ftell(pFile)/sizeof(int);

			int * nDen1 = new int [size]; 
			fseek(pFile, 0, SEEK_SET);
			fread(nDen1 ,sizeof(int), size, pFile);
			fclose(pFile);

			int * nDen2 = new int [relable];
			for(int j = 0; j < relable; j++)
				nDen2[j] = 0;	

			for(int j = 0; j < size; j++)
			{
				nDen2[map1[j]]+=nDen1[j];		
			}
			delete [] nDen1;

			sprintf(fileName, "%s/HybSpectrum/Graph%d/nDen", workDir.c_str(), i);
			pFile = fopen(fileName, "w");
			fwrite(nDen2, sizeof(int), relable, pFile);
			fclose(pFile);

			delete [] nDen2;

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
			sortEdges(workDir, 2, i, 1, nTasks);

			makeGraph(workDir, i, 1, relable, nTasks); 

			delete [] map1;
			delete [] fileName;
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

		sprintf(fileName, "%s/OvlGraph/finalMap", workDir.c_str());
		pFile = fopen(fileName, "r");
		fseek(pFile, 0, SEEK_END);
		int ovlSize = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);

		int * nMap = new int [ovlSize];
		fread(nMap, sizeof(int), ovlSize, pFile);  
		fclose(pFile);

		int * reLabelMap = new int [mapSize];
		for(int i = 0; i < mapSize; i++)
		{
			reLabelMap[i] = -1;
		}

		for(int i = 0; i < mapSize; i++)
		{
			if(map[i] == -1)
			{
				reLabelMap[nMap[i]] = i;
			}
		}

		char * singleTons = new char [9 * megabyte];
		char * flags = new char [9 * megabyte];

		int * rBounds = new int [9 * megabyte];
		int * cBounds = new int [9 * megabyte];

		int seqPos = 0;
		int rBoundsPos = 0;
		int cBoundsPos = 0;

		sprintf(fileName, "%s/ContigsInfo/Contigs.store", workDir.c_str());

		pFile = fopen(fileName, "r");
		fseek(pFile, 0, SEEK_END);
		int consOffset = ftell(pFile);
		fseek(pFile, 0, SEEK_SET);
		fclose(pFile);

		int currIndex = 0;			
		sprintf(fileName, "%s/Fragments/sequences.%d", workDir.c_str(), currIndex++);
		pFile = fopen(fileName, "r");

		unsigned long long int fragmentOffset;
		fread(&fragmentOffset, sizeof(unsigned long long int), 1, pFile);

		unsigned long long int indexOffset;
		fread(&indexOffset, sizeof(unsigned long long int), 1, pFile);

		unsigned long long int length;
		fread(&length, sizeof(unsigned long long int), 1, pFile);

		unsigned long long int numFragments;
		fread(&numFragments, sizeof(unsigned long long int), 1, pFile);

		fseek(pFile, 0, SEEK_SET);

		Fragment_Index index;
		index.reserve(numFragments, length, indexOffset, fragmentOffset);

		index.read(pFile);
		fclose(pFile);

		index.set();

		int relabel = clusterCount;
		for(int i = 0; i < mapSize; i++)
		{
			if(reLabelMap[i] != -1)
			{
				while(i >= index.getIndexOffset() + index.numFragments())
				{
					index.clear();

					sprintf(fileName, "%s/Fragments/sequences.%d", workDir.c_str(), currIndex++);

					pFile = fopen(fileName, "r");
					fread(&fragmentOffset, sizeof(unsigned long long int), 1, pFile);
					fread(&indexOffset, sizeof(unsigned long long int), 1, pFile);
					fread(&length, sizeof(unsigned long long int), 1, pFile);
					fread(&numFragments, sizeof(unsigned long long int), 1, pFile);

					fseek(pFile, 0, SEEK_SET);

					index.reserve(numFragments, length, indexOffset, fragmentOffset);
					index.read(pFile);

					index.set();

					fclose(pFile);

				}

				map[reLabelMap[i]] = relabel++;

				pair<long long int, long long int> bounds = index.indexBounds(i - index.getIndexOffset());


				int addLen = bounds.second - bounds.first;

				if(seqPos + addLen >= 9 * megabyte)
				{
					sprintf(fileName, "%s/ContigsInfo/Contigs.store", workDir.c_str());

					pFile = fopen(fileName, "a");
					fwrite(singleTons, sizeof(char), seqPos, pFile);
					fclose(pFile);

					sprintf(fileName, "%s/ContigsInfo/Flags.store", workDir.c_str());

					pFile = fopen(fileName, "a");
					fwrite(flags, sizeof(char), seqPos, pFile);
					fclose(pFile);

					seqPos = 0;
				}  

				if(rBoundsPos + 2 >= 9 * megabyte)
				{
					sprintf(fileName, "%s/ContigsInfo/readBounds.store", workDir.c_str());

					pFile = fopen(fileName, "a");
					fwrite(rBounds, sizeof(int), rBoundsPos, pFile);
					fclose(pFile);

					rBoundsPos = 0;
				}

				if(cBoundsPos + 1 >= 9 * megabyte)
				{
					sprintf(fileName, "%s/ContigsInfo/contigBounds.store", workDir.c_str());

					pFile = fopen(fileName, "a");
					fwrite(cBounds, sizeof(int), cBoundsPos, pFile);
					fclose(pFile);

					cBoundsPos = 0;
				}

				for(int j = 0; j < bounds.second - bounds.first; j++)
				{
					singleTons[seqPos] = index.at(bounds.first+j);
					flags[seqPos++] = '$'; 
				}

				rBounds[rBoundsPos++] = i;
				rBounds[rBoundsPos++] = 0;

				cBounds[cBoundsPos++] = consOffset;

				consOffset+=addLen; 
			}
		}

		if(seqPos > 0)
		{

			sprintf(fileName, "%s/ContigsInfo/Contigs.store", workDir.c_str());

			pFile = fopen(fileName, "a");
			fwrite(singleTons, sizeof(char), seqPos, pFile);
			fclose(pFile);

			sprintf(fileName, "%s/ContigsInfo/Flags.store", workDir.c_str());

			pFile = fopen(fileName, "a");
			fwrite(flags, sizeof(char), seqPos, pFile);
			fclose(pFile);

		}

		if(rBoundsPos > 0)
		{
			sprintf(fileName, "%s/ContigsInfo/readBounds.store", workDir.c_str());

			pFile = fopen(fileName, "a");
			fwrite(rBounds, sizeof(int), rBoundsPos, pFile);
			fclose(pFile);
		}

		if(cBounds > 0)
		{
			sprintf(fileName, "%s/ContigsInfo/contigBounds.store", workDir.c_str());

			pFile = fopen(fileName, "a");
			fwrite(cBounds, sizeof(int), cBoundsPos, pFile);
			fclose(pFile);
		}


		delete [] nMap; 
		delete [] reLabelMap;		
		delete [] singleTons;

		delete [] rBounds;
		delete [] cBounds;

		sprintf(fileName, "%s/HybSpectrum/Graph0/tmpMap", workDir.c_str());

		pFile = fopen(fileName, "w");
		fwrite(map, sizeof(int), mapSize, pFile);
		fclose(pFile);

		sprintf(fileName, "%s/Spectrum/Graph0/nDen", workDir.c_str());
		pFile = fopen(fileName, "r");			

		fseek(pFile, 0, SEEK_END);
		int size = ftell(pFile)/sizeof(int);
		int * nDen1 = new int [size]; 
		fseek(pFile, 0, SEEK_SET);

		fread(nDen1 ,sizeof(int), size, pFile);
		fclose(pFile);

		int * nDen2 = new int [relabel];
		for(int j = 0; j < relabel; j++)
			nDen2[j] = 0;	

		for(int j = 0; j < size; j++)
			nDen2[map[j]]+=nDen1[j];		

		delete [] nDen1;

		sprintf(fileName, "%s/HybSpectrum/Graph0/nDen", workDir.c_str());
		pFile = fopen(fileName, "w");
		fwrite(nDen2, sizeof(int), relabel, pFile);
		fclose(pFile);

		delete [] nDen2;	

		sprintf(fileName, "%s/OvlGraph/finalMap", workDir.c_str());
		pFile = fopen(fileName, "r");
		fseek(pFile, 0, SEEK_END);
		ovlSize = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);

		nMap = new int [ovlSize];
		fread(nMap, sizeof(int), ovlSize, pFile);  
		fclose(pFile);

		sprintf(fileName, "%s/HybSpectrum/deletedNodes", workDir.c_str());
		pFile = fopen(fileName, "r");

		fseek(pFile, 0, SEEK_END);
		ovlSize = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);

		int * dNodes = new int [ovlSize];
		fread(dNodes, sizeof(int), ovlSize, pFile);  
		fclose(pFile);

		makeFinalEdges(workDir, 0, 2, map, nMap, dNodes);

		delete [] nMap;
		delete [] dNodes;

		sortEdges(workDir, 0, 0, 2, nTasks);
		mergeEdges(workDir, 0, 2, nTasks);
		sortEdges(workDir, 2, 0, 2, nTasks);

		makeGraph(workDir, 0, 2, relabel, nTasks);

		rBounds = new int [size];
		for(int i = 0; i < size; i++)
		{
			rBounds[i] = -1;
		}

		sprintf(fileName, "%s/ContigsInfo/readBounds.store", workDir.c_str());

		FILE * rbFile = fopen(fileName, "r");
		fseek(rbFile, 0, SEEK_END);
		int rbSize = ftell(rbFile)/sizeof(int);
		fseek(rbFile, 0, SEEK_SET); 

		int * rbBuff = new int [rbSize];
		fread(rbBuff, sizeof(int), rbSize, rbFile);

		fclose(rbFile);

		for(int i = 0; i < rbSize; i+=2)
		{
			rBounds[rbBuff[i]] = rbBuff[i+1];
		}

		rbFile = fopen(fileName, "w");
		fwrite(rBounds, sizeof(int), size, rbFile);
		fclose(rbFile);

		delete [] rbBuff;
		delete [] rBounds;

		sprintf(fileName, "%s/ContigsInfo/contigBounds.store", workDir.c_str());

		FILE * cbFile = fopen(fileName, "r");
		fseek(cbFile, 0, SEEK_END);
		int cbSize = ftell(cbFile)/sizeof(int);
		fseek(cbFile, 0, SEEK_SET);

		int * cbBuff = new int [cbSize];
		fread(cbBuff, sizeof(int), cbSize, cbFile);

		fclose(cbFile);

		int currCbPos = 0;		

		char * cBuff = new char [100 * megabyte];
		char * fBuff = new char [100 * megabyte];

		sprintf(fileName, "%s/ContigsInfo/Contigs.store", workDir.c_str());
		FILE * cFile = fopen(fileName, "r");

		sprintf(fileName, "%s/ContigsInfo/Flags.store", workDir.c_str());
		FILE * fFile = fopen(fileName, "r");

		int numRead = 0; int sectionSize = 500; int runningTotalLen = 0; int runningTotalSections = 0; int numContigIndexes = 0;
		do{
			numRead = fread(cBuff, sizeof(char), 100*megabyte, cFile);
			fread(fBuff, sizeof(char), 100*megabyte, fFile);

			int counts = numRead;

			while(currCbPos < cbSize && cbBuff[currCbPos] < runningTotalLen+counts) currCbPos++;
			currCbPos--;

			if(counts == 100*megabyte && currCbPos > 0 && currCbPos != cbSize-1 && (cbBuff[currCbPos+1] - cbBuff[currCbPos]) < 100 * megabyte)
			{
				fseek(cFile, cbBuff[currCbPos]-(runningTotalLen+counts), SEEK_CUR);
				fseek(fFile, cbBuff[currCbPos]-(runningTotalLen+counts), SEEK_CUR);
				counts+=(cbBuff[currCbPos]- (runningTotalLen+counts));
			}	

			int numSections = counts/sectionSize + ((sectionSize-1) + counts%sectionSize)/sectionSize;
			Fragment_Index contigIndex;

			contigIndex.reserve(numSections, counts, runningTotalSections, runningTotalLen);
			runningTotalSections+=numSections; runningTotalLen+=counts;	

			for(int i = 0; i < counts; i+=sectionSize)
			{
				int addSize = min(sectionSize, (counts-i));
				contigIndex.append(cBuff+i, fBuff+i, addSize);
			}			

			contigIndex.set();	

			sprintf(fileName, "%s/ContigsInfo/Contigs.store.%d", workDir.c_str(), numContigIndexes++);

			FILE * oFile = fopen(fileName, "w");
			contigIndex.write(oFile);
			fclose(oFile);	

		}while(numRead == (100*megabyte));	

		delete [] cBuff;
		delete [] fBuff;
		delete [] cbBuff;

		sprintf(fileName, "%s/ContigsInfo/Contigs.store", workDir.c_str());
		removeFile(fileName);

		sprintf(fileName, "%s/ContigsInfo/Flags.store", workDir.c_str());
		removeFile(fileName);

		delete [] fileName;

		for(int i = numGraphs-1; i >= 0; i--)
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

				int * dNodes = new int [mSize];
				for(int i = 0; i < mSize; i++)
					dNodes[i] = 0;

				delete [] fileName;
				int currNode = 0;

				int gSize = 2; int graphsLoaded = 0; int currGraphSet = 0;
				Graph ** graphs = new Graph *[gSize];

				int iSize = 2; int indexesLoaded = 0; int currIndexSet = 0;
				Fragment_Index ** indexes = new Fragment_Index * [iSize];

				uint32_t * nodes = new uint32_t [1000];
				uint32_t * maxEdges = new uint32_t [2000];

				uint32_t * forwardEdges = new uint32_t [1000 * 1001]; //More Capacity
				uint32_t * backEdges = new uint32_t [1000 * 1001]; //More Capacity

				uint32_t * orderings = new uint32_t [2000];
				uint32_t * contig = new uint32_t [1000]; 
				uint32_t * sandBox = new uint32_t [1100];

				int consPos = 0; int cSize = 9 * megabyte;
				char * consensus = new char [cSize];

				char * flags = new char [cSize];

				int cbSize = 3 * megabyte; int cbPos = 0;
				int * cBounds = new int [cbSize];

				int rbSize = 3 * megabyte; int rbPos = 0;
				int * rBounds = new int [25 * megabyte];

				int pcSize = 6000; 
				int * parCons1 = new int [pcSize];
				int * parCons2 = new int [pcSize];

				int matxSize = 1000;
				int ** matrix = new int * [matxSize];
				for(int i = 0; i < matxSize; i++)
					matrix[i] = new int [1000];

				int prevContigLen = 0;

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
								delete [] nodes; nodes = new uint32_t [numNodes];
								delete [] maxEdges; maxEdges = new uint32_t[2*numNodes];
								delete [] forwardEdges; forwardEdges = new uint32_t[1001*numNodes];
								delete [] backEdges; backEdges = new uint32_t[1001*numNodes];

								delete [] orderings; orderings = new uint32_t[2*numNodes];
								delete [] contig; contig = new uint32_t[numNodes];
								delete [] sandBox, sandBox = new uint32_t[100+numNodes];
							}

							loadGraphs(workDir, graphs, gSize, graphsLoaded, currGraphSet, contigs, startPos, count, nodeMap);
							loadFragments(workDir, indexes, iSize, indexesLoaded, currIndexSet, contigs, startPos, count, nodeMap);

							for(int j = 0; j < numNodes; j++)
							{
								forwardEdges[1001*j] = 0;
								backEdges[1001*j] = 0;
								nodes[j] = 1;
								sandBox[100+j] = 0;
							}

							removeContainments(graphs, graphsLoaded, contigs, startPos, count, nodeMap, nodes, minNode, maxNode);
							transitiveReduce(graphs, graphsLoaded, nodes, forwardEdges, backEdges, sandBox, contigs, startPos, count, nodeMap, minNode, maxNode);
							int numTrimmed = trimGraph(nodes, forwardEdges, backEdges, sandBox, minNode, maxNode, minPathLen, dNodes);


							priority_queue<pair<int, int> > myQueue; map<int, int > prevVisit; //For Now
							int numPopped = popBubbles(nodes, forwardEdges, backEdges, sandBox, minNode, maxNode,  minPathLen, dNodes);
							float numContigs = traverseGraph(nodes, forwardEdges, backEdges, orderings, sandBox, minNode, maxNode);

							//ERASE THIS "<<endl;
							for(int j = 0; j < numNodes; j++)
							{
								for (int k = 0; k < forwardEdges[1001*j]; k++)
									cout<<j+minNode<<" "<<forwardEdges[1001*j+1+k]<<endl;
							}
							cout<<endl;

							for(int j = 0; j < numNodes; j++)
							{
								for (int k = 0; k < forwardEdges[1001*j]; k++)
									if(nodes[forwardEdges[1001*j+1+k]-minNode] && nodes[j])
									cout<<j+minNode<<" "<<forwardEdges[1001*j+1+k]<<endl;

							}

							cout<<"*"<<endl;
							//DONE ERASE 

							float maxContig = 0; float totalContigs = 0; int currContig = 0;
							int contigPos = 0;

							while(currContig < numContigs)
							{
								int num1 = orderings[contigPos];

								contigPos = contigPos+num1+1;		
								int num2 = orderings[contigPos];

								int contigLen = num1+num2;
								if(contigLen > maxContig)
									maxContig = contigLen;

								totalContigs+=contigLen;

								contigPos = contigPos+num2+1;		

								currContig++;
							}


							totalContigs+=(numTrimmed+numPopped);

							if(numContigs != 1 || maxContig/totalContigs < 0)
							{
								for(int j = startPos; j < count && contigs[j] != -1; j++)
								{
									dNodes[nodeMap[contigs[j]]] = 0;
								}
							}	

							cout<<numContigs<<"NUMBER "<<endl;
							if(numContigs == 1 && maxContig/totalContigs >= 0)
							{
								int num1 = orderings[0];
								int num2 = orderings[num1 + 1];

								int m = 1; contig[0] = 0;
								for(int j = num1 + num2 + 1; j > num1 + 1; j--)
								{
									contig[m++] = orderings[j];
									contig[0]++;
								} 

								for(int j = 1; j < num1 + 1; j++)
								{
									contig[m++] = orderings[j];
									contig[0]++;
								}												
								buff[buffSize++] = nIndex;

								if(cbPos + 1 >= cbSize)
								{
									char * fileName = new char [1000];
									sprintf(fileName, "%s/ContigsInfo/contigBounds.%d", workDir.c_str(), rank); 										
									FILE * cbFile = fopen(fileName, "a");
									fwrite(cBounds, sizeof(int), cbPos, cbFile);
									fclose(cbFile); 

									delete [] fileName;

									cbPos = 0;
								} 

								cBounds[cbPos++] = prevContigLen + consPos;

								Consensus myConsensus; int oPos = 0; int prevVal = 0;  
								while(myConsensus.getConsensus(contig, oPos, minNode, maxNode,  matrix,  matxSize , graphs, gSize, indexes, iSize, consensus, consPos, cSize, flags, rBounds, prevVal, rbPos, rbSize, parCons1, parCons2, pcSize) == false)
								{
									if(rbPos+2 >= rbSize)
									{
										char * fileName = new char [1000];
										sprintf(fileName, "%s/ContigsInfo/readBounds.%d", workDir.c_str(), rank); 										
										FILE * rbFile = fopen(fileName, "a");
										fwrite(rBounds, sizeof(int), rbPos, rbFile);
										fclose(rbFile); 

										delete [] fileName;

										rbPos = 0;

									}else{
										char * fileName = new char [1000];
										sprintf(fileName, "%s/ContigsInfo/cons.%d", workDir.c_str(), rank); 										
										FILE * consFile = fopen(fileName, "a");
										fwrite(consensus, sizeof(char), consPos, consFile);
										fclose(consFile); 

										sprintf(fileName, "%s/ContigsInfo/flags.%d", workDir.c_str(), rank); 										
										FILE * flagsFile = fopen(fileName, "a");
										fwrite(flags, sizeof(char), consPos, flagsFile);
										fclose(flagsFile); 

										delete [] fileName;

										prevContigLen+=consPos;

										consPos = 0;
									}	
								}
							}

							if(numNodes > 1000)
							{
								delete [] nodes; nodes = new uint32_t [1000];
								delete [] maxEdges; maxEdges = new uint32_t[2000];
								delete [] forwardEdges; forwardEdges = new uint32_t[1000 * 1001];
								delete [] backEdges; backEdges = new uint32_t[1000 * 1001];
								delete [] orderings; orderings = new uint32_t[2000];
								delete [] contig; contig = new uint32_t[1000];

								delete [] sandBox; sandBox = new uint32_t[numNodes+100];
							}

						}

						currNode++;
					}
				}

				if(cbPos > 0)
				{
					char * fileName = new char [1000];
					sprintf(fileName, "%s/ContigsInfo/contigBounds.%d", workDir.c_str(), rank); 									

					FILE * cbFile = fopen(fileName, "a");
					fwrite(cBounds, sizeof(int), cbPos, cbFile);
					fclose(cbFile); 

					delete [] fileName;
				}

				if(rbPos > 0)
				{
					char * fileName = new char [1000];
					sprintf(fileName, "%s/ContigsInfo/readBounds.%d", workDir.c_str(), rank); 										
					FILE * rbFile = fopen(fileName, "a");
					fwrite(rBounds, sizeof(int), rbPos, rbFile);
					fclose(rbFile); 

					delete [] fileName;
				}

				if(consPos > 0)
				{
					char * fileName = new char [1000];
					sprintf(fileName, "%s/ContigsInfo/cons.%d", workDir.c_str(), rank); 										
					FILE * consFile = fopen(fileName, "a");
					fwrite(consensus, sizeof(char), consPos, consFile);
					fclose(consFile); 

					sprintf(fileName, "%s/ContigsInfo/flags.%d", workDir.c_str(), rank); 										
					FILE * flagsFile = fopen(fileName, "a");
					fwrite(flags, sizeof(char), consPos, flagsFile);
					fclose(flagsFile); 

					delete [] fileName;
				}	



				delete [] consensus;
				delete [] flags;
				delete [] parCons1;
				delete [] parCons2;

				for(int i = 0; i < matxSize; i++)
					delete [] matrix [i];

				delete [] matrix;

				delete [] nodes;
				delete [] maxEdges;

				delete [] contigs;
				delete [] nodeMap;
				delete [] cIndex;

				delete [] cBounds;
				delete [] rBounds;

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

					char * fileName = new char [1000];
					sprintf(fileName, "%s/HybSpectrum/deletedNodes", workDir.c_str());

					pFile = fopen(fileName, "r");
					int * dNodesM = new int [mSize];
					fread(dNodesM, sizeof(int), mSize, pFile);
					fclose(pFile);

					for(int i = 0; i < mSize; i++)
					{
						if(dNodes[i] == 1)
							dNodesM[i] = 1;
					}


					fileName = new char [1000];
					sprintf(fileName, "%s/HybSpectrum/deletedNodes", workDir.c_str()); 										
					pFile = fopen(fileName, "w");
					fwrite(dNodesM, sizeof(int), mSize, pFile);
					fclose(pFile); 

					delete [] dNodesM;
					delete [] fileName;

					MPI_Send(&handshake, 1, MPI_INT, 0, workTag, MPI_COMM_WORLD);
				}

				delete [] buff;
				delete [] dNodes;

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

						cout<<"ADDING "<<arr[j+1]<<" "<<arr[j+2]<<endl;

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

	if(reportTime == true && rank == 0)
	{
		time_t rawtime;
		struct tm * timeinfo;
		time(&rawtime);
		timeinfo = localtime(&rawtime);
		cout<<"end time: "<<asctime(timeinfo)<<endl;
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

	sprintf(fileName, "%s/nDen", graphDir);
	if(fileExists(fileName))
	{
		sprintf(dest, "%s/nDen", tmpDir);
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
				numEdges++;
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

//void makeFinalEdges(string &, int, int * &)
//Description: This creates final edges for the collapsed graph
//Input: workDir (string &): the working directory, graphIter (int): The current
//graph iteration, nodeMap (int *): the mapping from the previous graph
//to the new one
//Output:None
//Return:None
void makeFinalEdges(string & workDir, int graphIter, int auxFiles, int * & nodeMapHyb, int * & nodeMapFinal, int * & dNodes) {

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
			if(!dNodes[nodeMapFinal[graph.edgeSource(j) + graph.getOffset()]] && !dNodes[nodeMapFinal[graph.edgeDest(j)]])
			{
				uint32_t source = nodeMapHyb[graph.edgeSource(j) + graph.getOffset()];
				uint32_t dest = nodeMapHyb[graph.edgeDest(j)];
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
					numEdges++;
					buff[count++] = source;
					buff[count++] = dest;
					buff[count++] = weight | ((sContained << 14) | (dContained << 13) | (rDovetail << 12));
				}

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
			for(int j = 0; j < forwardEdges[1001*i]; j++)
			{
				int n = forwardEdges[1001*i+j+1];
				if(nodes[n-minNode])
				{
					fIndex = j; fTally++;
				}
			}

			for(int j = 0; j < backEdges[1001*i]; j++)
			{
				int n = backEdges[1001*i+j+1];
				if(nodes[n-minNode])
				{
					bIndex = j; bTally++;
				}
			}

			if(fTally == 1)
			{		
				int nextNode = forwardEdges[1001*i+1+fIndex];
				fTally = 0; bTally = 0;
				for(int j = 0; j < forwardEdges[1001*(nextNode-minNode)]; j++)
				{
					int n = forwardEdges[1001*(nextNode-minNode)+j+1];
					if(nodes[n-minNode])
					{
						fIndex = j; fTally++;
					}
				}

				for(int j = 0; j < backEdges[1001*(nextNode-minNode)]; j++)
				{
					int n = backEdges[1001*(nextNode-minNode)+j+1];
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
					nextNode = forwardEdges[1001*(nextNode-minNode)+1+fIndex];

					fTally = 0; bTally = 0;
					for(int j = 0; j < forwardEdges[1001*(nextNode-minNode)]; j++)
					{
						int n = forwardEdges[1001*(nextNode-minNode)+j+1];
						if(nodes[n-minNode])
						{
							fIndex = j; fTally++;
						}
					}

					for(int j = 0; j < backEdges[1001*(nextNode-minNode)]; j++)
					{
						int n = backEdges[1001*(nextNode-minNode)+j+1];
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
			for(int j = 0; j < forwardEdges[1001*i]; j++)
			{
				int n = forwardEdges[1001*i+j+1];
				if(nodes[n-minNode])
				{
					fIndex = j; fTally++;
				}
			}

			for(int j = 0; j < backEdges[1001*i]; j++)
			{
				int n = backEdges[1001*i+j+1];
				if(nodes[n-minNode])
				{
					bIndex = j; bTally++;
				}
			}

			if(bTally == 1)
			{
				int nextNode = backEdges[1001*i+1+bIndex];
				fTally = 0; bTally = 0;
				for(int j = 0; j < forwardEdges[1001*(nextNode-minNode)]; j++)
				{
					int n = forwardEdges[1001*(nextNode-minNode)+j+1];
					if(nodes[n-minNode])
					{
						fIndex = j; fTally++;
					}
				}

				for(int j = 0; j < backEdges[1001*(nextNode-minNode)]; j++)
				{
					int n = backEdges[1001*(nextNode-minNode)+j+1];
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
					nextNode = backEdges[1001*(nextNode-minNode)+1+bIndex];

					fTally = 0; bTally = 0;
					for(int j = 0; j < forwardEdges[1001*(nextNode-minNode)]; j++)
					{
						int n = forwardEdges[1001*(nextNode-minNode)+j+1];
						if(nodes[n-minNode])
						{
							fIndex = j; fTally++;
						}
					}

					for(int j = 0; j < backEdges[1001*(nextNode-minNode)]; j++)
					{
						int n = backEdges[1001*(nextNode-minNode)+j+1];
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
int trimGraph(uint32_t * & nodes, uint32_t * & forwardEdges, uint32_t * & backEdges, uint32_t * & sandBox, int minNode, int maxNode, int pathLen, int * & dNodes)
{
	int numTrimmed = 0; int numNodesInGraph = 0;
	for(int i = 0; i < maxNode-minNode+1; i++)
		if(nodes[i])
			numNodesInGraph++;

	for(int i = 0; i < maxNode-minNode+1; i++)
	{
		if(nodes[i])
		{
			cout<<"TRIMMING "<<i+minNode<<endl;
			int numExplored = 0;
			int fIndex = 0; int bIndex = 0;
			int fTally = 0; int bTally = 0;
			for(int j = 0; j < forwardEdges[1001*i]; j++)
			{
				int n = forwardEdges[1001*i+j+1];
				if(nodes[n-minNode])
				{
					fIndex = j; fTally++;
				}
			}

			for(int j = 0; j < backEdges[1001*i]; j++)
			{
				int n = backEdges[1001*i+j+1];
				if(nodes[n-minNode])
				{
					bIndex = j; bTally++;
				}
			}

			cout<<bTally<<" "<<fTally<<endl;

			if(fTally == 0 && bTally == 1)
			{
				sandBox[numExplored++] = i + minNode;
				int nextNode = backEdges[1001*i+bIndex+1];
				fTally = 0; bTally = 0;

				for(int j = 0; j < forwardEdges[1001*(nextNode-minNode)]; j++)
				{
					int n = forwardEdges[1001*(nextNode-minNode)+j+1];
					if(nodes[n-minNode])
					{
						fIndex = j; fTally++;
					}
				}

				for(int j = 0; j < backEdges[1001*(nextNode-minNode)]; j++)
				{
					int n = backEdges[1001*(nextNode-minNode)+j+1];
					if(nodes[n-minNode])
					{
						bIndex = j; bTally++;
					}
				}


				while(fTally == 1 && bTally == 1)
				{
					sandBox[numExplored++] = nextNode;
					nextNode = backEdges[1001*(nextNode-minNode)+ bIndex + 1];

					fTally = 0; bTally = 0;
					for(int j = 0; j < forwardEdges[1001*(nextNode-minNode)]; j++)
					{
						int n = forwardEdges[1001*(nextNode-minNode)+j+1];
						if(nodes[n-minNode])
						{
							fIndex = j; fTally++;
						}
					}

					for(int j = 0; j < backEdges[1001*(nextNode-minNode)]; j++)
					{
						int n = backEdges[1001*(nextNode-minNode)+j+1];
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
				int nextNode = forwardEdges[1001*i+1];
				fTally = 0; bTally = 0;

				for(int j = 0; j < forwardEdges[1001*(nextNode-minNode)]; j++)
				{
					int n = forwardEdges[1001*(nextNode-minNode)+j+1];
					if(nodes[n-minNode])
					{
						fIndex = j; fTally++;
					}
				}

				for(int j = 0; j < backEdges[1001*(nextNode-minNode)]; j++)
				{
					int n = backEdges[1001*(nextNode-minNode)+j+1];
					if(nodes[n-minNode])
					{
						bIndex = j; bTally++;
					}
				}

				while(fTally == 1 && bTally  == 1)
				{
					sandBox[numExplored++] = nextNode;
					nextNode = forwardEdges[1001*(nextNode-minNode)+fIndex + 1];

					fTally = 0; bTally = 0;
					for(int j = 0; j < forwardEdges[1001*(nextNode-minNode)]; j++)
					{
						int n = forwardEdges[1001*(nextNode-minNode)+j+1];
						if(nodes[n-minNode])
						{
							fIndex = j; fTally++;
						}
					}

					for(int j = 0; j < backEdges[1001*(nextNode-minNode)]; j++)
					{
						int n = backEdges[1001*(nextNode-minNode)+j+1];
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
			if(numExplored <= pathLen && p1/g1 < .5)
			{
				for(int j = 0; j < numExplored; j++)
				{
					cout<<"DELETING "<<sandBox[j]<<endl;
					nodes[sandBox[j]-minNode] = 0;			
					dNodes[sandBox[j]] = 1;
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
int popBubbles(uint32_t * & nodes, uint32_t * & forwardEdges, uint32_t * & backEdges, uint32_t * & sandBox, int minNode, int maxNode, int pathLen, int * & dNodes)
{
	int numPopped = 0;
	for(int i = 0; i < maxNode - minNode + 1; i++)
	{
		if(nodes[i])
		{
			for(int x = 0; x < forwardEdges[1001*i]; x++)
			{
				if(nodes[forwardEdges[1001*i+x+1]- minNode] == 1)
				{ 
					for(int y = x+1; y < forwardEdges[1001*i]; y++)
					{
						if(nodes[forwardEdges[1001*i+y+1] - minNode] == 1)
						{
							int node1 = forwardEdges[i*1001+1+x];
							int node2 = forwardEdges[i*1001+1+y];

							if(nodes[node1-minNode] && nodes[node2-minNode])
							{
								int path1 = 0; int path2 = pathLen; 
								sandBox[path1++] = node1;
								sandBox[path2++] = node2;
								sandBox[100+(node1-minNode)] = 1;
								sandBox[100+(node2-minNode)] = 1;

								int fIndex1 = 0; int bIndex1 = 0;
								int fTally1 = 0; int bTally1 = 0;
								for(int j = 0; j < forwardEdges[1001*(node1-minNode)]; j++)
								{
									int n = forwardEdges[1001*(node1-minNode)+j+1];
									if(nodes[n-minNode])
									{
										fIndex1 = j; fTally1++;
									}
								}

								for(int j = 0; j < backEdges[1001*(node1-minNode)]; j++)
								{
									int n = backEdges[1001*(node1-minNode)+j+1];
									if(nodes[n-minNode])
									{
										bIndex1 = j; bTally1++;
									}
								}

								int fIndex2 = 0; int bIndex2 = 0;
								int fTally2 = 0; int bTally2 = 0;
								for(int j = 0; j < forwardEdges[1001*(node2-minNode)]; j++)
								{
									int n = forwardEdges[1001*(node2-minNode)+j+1];
									if(nodes[n-minNode])
									{
										fIndex2 = j; fTally2++;
									}
								}

								for(int j = 0; j < backEdges[1001*(node2-minNode)]; j++)
								{
									int n = backEdges[1001*(node2-minNode)+j+1];
									if(nodes[n-minNode])
									{
										bIndex2 = j; bTally2++;
									}
								}

								while((fTally1 == 1 && bTally1 == 1) || (fTally2 == 1 && bTally2 == 1))
								{
									if(fTally1 == 1 && bTally1 == 1 && path1 < pathLen)
									{
										node1 = forwardEdges[1001*(node1-minNode)+1+fIndex1];

										fTally1 = 0; bTally1 = 0;
										for(int j = 0; j < forwardEdges[1001*(node1-minNode)]; j++)
										{
											int n = forwardEdges[1001*(node1-minNode)+j+1];
											if(nodes[n-minNode])
											{
												fIndex1 = j; fTally1++;
											}
										}

										for(int j = 0; j < backEdges[1001*(node1-minNode)]; j++)
										{
											int n = backEdges[1001*(node1-minNode)+j+1];
											if(nodes[n-minNode])
											{
												bIndex1 = j; bTally1++;
											}
										}

										if(sandBox[100+(node1-minNode)] == 1) // changed from node2 to node1
										{
											int p1 = 0; 
											for(p1 = pathLen; p1 < path2; p1++)
											{
												if(sandBox[p1] == node1)
													break;
											}

											for(int j = p1-1; j >= pathLen; j--)
											{
												nodes[sandBox[j]-minNode] = 0;
												dNodes[sandBox[j]] = 1;
												numPopped++;
											}

											if(p1-1 == pathLen-1)
											{
												int tmpHold[forwardEdges[1001*i]];

												int numAdded = 0;
												for(int j = 0; j < forwardEdges[1001*i]; j++)
												{
													if(forwardEdges[1001*i + 1 + j] != node1)
													{
														tmpHold[numAdded++] = forwardEdges[1001*i + 1 + j];
													} 
												}

												for(int k = 0; k < numAdded; k++)
												{
													forwardEdges[1001*i + 1 + k] = tmpHold[k];
												}
												forwardEdges[1001*i] = numAdded;

												int tmpHold2[backEdges[1001*(node1 - minNode)]];

												numAdded = 0;
												for(int j = 0; j < backEdges[1001*(node1 - minNode)]; j++)
												{
													if(backEdges[1001*(node1 - minNode) + j + 1] != i + minNode)
													{
														tmpHold2[numAdded++] = backEdges[1001*(node1 - minNode) + j + 1];
													}
												}

												for(int k = 0; k < numAdded; k++)
												{
													backEdges[1001*(node1 - minNode) + k + 1] = tmpHold2[k]; 
												}

												backEdges[1001*(node1 - minNode)] = numAdded;
											}
											break;
										}

										sandBox[path1++] = node1;
										sandBox[100+(node1-minNode)] = 1;
									}

									if(fTally2 == 1 && bTally2 == 1 && path2 < (2*pathLen))
									{
										node2 = forwardEdges[1001*(node2-minNode)+1+fIndex2];

										fTally2 = 0; bTally2 = 0;
										for(int j = 0; j < forwardEdges[1001*(node2-minNode)]; j++)
										{
											int n = forwardEdges[1001*(node2-minNode)+j+1];
											if(nodes[n-minNode])
											{
												fIndex2 = j; fTally2++;
											}
										}

										for(int j = 0; j < backEdges[1001*(node2-minNode)]; j++)
										{
											int n = backEdges[1001*(node2-minNode)+j+1];
											if(nodes[n-minNode])
											{
												bIndex2 = j; bTally2++;
											}
										}

										if(sandBox[100+(node2-minNode)] == 1) //changed node1 to node2
										{
											int p1 = 0;
											for(p1 = 0; p1 < path1; p1++)
											{
												if(sandBox[p1] == node2)
												{
													break;
												}
											}

											for(int j = p1-1; j >= 0; j--)
											{
												nodes[sandBox[j]-minNode] = 0;
												dNodes[sandBox[j]] = 1;
												numPopped++;
											}

											if(p1-1 == -1)
											{
												int tmpHold[forwardEdges[1001*i]];

												int numAdded = 0;
												for(int j = 0; j < forwardEdges[1001*i]; j++)
												{
													if(forwardEdges[1001*i + 1 + j] != node2)
													{
														tmpHold[numAdded++] = forwardEdges[1001*i + 1 + j];
													} 
												}

												for(int k = 0; k < numAdded; k++)
												{
													forwardEdges[1001*i + 1 + k] = tmpHold[k];
												}
												forwardEdges[1001*i] = numAdded;

												int tmpHold2[backEdges[1001*(node2 - minNode)]];

												numAdded = 0;
												for(int j = 0; j < backEdges[1001*(node2 - minNode)]; j++)
												{
													if(backEdges[1001*(node2 - minNode) + j + 1] != i + minNode)
													{
														tmpHold2[numAdded++] = backEdges[1001*(node2 - minNode) + j + 1];
													}
												}

												for(int k = 0; k < numAdded; k++)
												{
													backEdges[1001*(node2 - minNode) + k + 1] = tmpHold2[k]; 
												}

												backEdges[1001*(node2 - minNode)] = numAdded;
											}

											break;
										}
										sandBox[path2++] = node2;
										sandBox[100+(node2-minNode)] = 1;									}

										if(path1 >= pathLen && path2 >= (2*pathLen))
											break;

										if(((path1 >= pathLen) && (fTally2 != 1 || bTally2 != 1)) || ((path2 >= 2*pathLen) && (fTally1 != 1 || bTally1 != 1)))
											break;
								}

								for(int j = 0; j < path1; j++)
								{
									sandBox[(sandBox[j]-minNode)+100] = 0;
								}

								for(int j = pathLen; j < path2; j++)
								{
									sandBox[(sandBox[j]-minNode)+100] = 0;
								}
							}
						}
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
void transitiveReduce(Graph ** & graphs, int graphsLoaded, uint32_t * & nodes, uint32_t * & forwardEdges, uint32_t * & backEdges, uint32_t * & sandBox, int * & contigs, int startPos, int count, int * & nodeMap, int minNode, int maxNode)
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

			for(int j = 0; j < graph1->nodeDegree(node1-graph1->getOffset()); j++)
			{
				uint32_t dest = graph1->edgeDest(node1-graph1->getOffset(), j);
				if(!graph1->rDovetail(node1-graph1->getOffset(), j) && dest >= minNode && dest <= maxNode && nodes[dest-minNode])
				{
					sandBox[100+(dest-minNode)] = 1;
				}

			}

			for(int j = 0; j < graph1->nodeDegree(node1-graph1->getOffset()); j++)
			{
				uint32_t node2 = graph1->edgeDest(node1-graph1->getOffset(), j);
				if(!graph1->rDovetail(node1-graph1->getOffset(), j) && node2 >= minNode && node2 <= maxNode && nodes[node2-minNode])
				{
					cout<<node1<<" "<<node2<<endl;
					Graph * graph2;
					for(int l = 0; l < graphsLoaded; l++)
					{
						graph2 = graphs[l];
						if(graph2->getOffset() <= node2 && node2 < graph2->getOffset() + graph2->getNumNodes())
							break;
					}

					for(int h = 0; h < graph2->nodeDegree(node2-graph2->getOffset()); h++)
					{
						uint32_t dest = graph2->edgeDest(node2-graph2->getOffset(), h);
						if(!graph2->rDovetail(node2-graph2->getOffset(), h) && dest >= minNode && dest <= maxNode && nodes[dest-minNode])
						{
							if(node1 == 54324) cout<<node2<<" ELIM "<<dest<<endl;
							sandBox[100+(dest-minNode)] = 0;
						}
					}
				}
			}
		}				

		for(int i = startPos; i < count && contigs[i] != -1; i++)
		{	
			int node2 = nodeMap[contigs[i]];

			int numEdges1 = forwardEdges[1001*(node1-minNode)];
			int numEdges2 = backEdges[1001*(node2-minNode)];
			
			if(node1 == 54324) cout<<numEdges1<<" "<<numEdges2<<" THE EDGES "<<node2<<" "<<sandBox[100+(node2-minNode)]<<endl;
			if(node1 != node2 && numEdges1 < 1000 && numEdges2 < 1000 && sandBox[100+(node2-minNode)])
			{
				int pos1 = numEdges1+1001*(node1-minNode)+1;
				int pos2 = numEdges2+1001*(node2-minNode)+1;
				forwardEdges[pos1] = node2;
				backEdges[pos2] = node1;

				forwardEdges[1001*(node1-minNode)]++;
				backEdges[1001*(node2-minNode)]++;

			}	

			sandBox[100+(node2-minNode)] = 0;
		}
	}
	cout<<"DONE "<<endl;
}


void bestOverlaps(Graph ** & graphs, int graphsLoaded, int * & contigs, int startPos, int count, int * & nodeMap, uint32_t * & nodes, uint32_t * & maxEdges, int minNode, int maxNode)
{
	for(int i = startPos; i < count && contigs[i] != -1; i++)
	{
		int currNode = nodeMap[contigs[i]];
		bool foundReverse = false; bool foundForward = false;

		if(nodes[currNode-minNode])
		{
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
				if(dest >= minNode && dest <= maxNode && nodes[dest-minNode])
				{
					if(!graph->sContained(source, j) && !graph->dContained(source, j))
					{
						if(!graph->rDovetail(source, j) && !foundForward)
						{ 
							int dest = graph->edgeDest(source, j);
							maxEdges[2*(currNode-minNode)] = dest;

							foundForward = true;
						}

						if(graph->rDovetail(source, j) && !foundReverse)
						{ 
							int dest = graph->edgeDest(source, j);
							maxEdges[2*(currNode-minNode)+1] = dest;

							foundReverse = true;
						}
					}
				}

				if(foundForward && foundReverse)
					break;
			}

		}

		if(!foundForward)
			maxEdges[2*(currNode-minNode)] = -1;

		if(!foundReverse)
			maxEdges[2*(currNode-minNode)+1] = -1;
	}
}

//void removeContainments()
//Description: This function removes containment fragments from the graph
//Input: 
//Output:None
//Return:None
void removeContainments(Graph ** & graphs, int graphsLoaded, int * & contigs, int startPos, int count, int * & nodeMap, uint32_t * & nodes, int minNode, int maxNode)
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
		//int dContained = 0; int doveTails = 0;
		for(int j = 0; j < graph->nodeDegree(source) && nodes[currNode-minNode] != 0; j++) //ADDED nodes[currNode-minNode] != 0
		{
			int dest = graph->edgeDest(source, j);
			if(dest >= minNode && dest <=  maxNode && nodes[dest-minNode] != 0) //ADDED nodes[currNode-minNode] != 0
			{
				Graph * graph2  = '\0';
				for(int k = 0; k < graphsLoaded; k++)
				{
					graph2 = graphs[k];
					if(graph2->getOffset() <= dest && dest < graph2->getOffset() + graph2->getNumNodes())
						break;
				}

				if(graph->sContained(source, j) && graph->dContained(source, j))
				{
					int ncTally1 = 0; int ncTally2 = 0;

					for(int h = 0; h < graph->nodeDegree(source); h++)
						if(!graph->dContained(source, h) && !graph->sContained(source, h))
							ncTally1++;

					for(int h = 0; h < graph2->nodeDegree(dest-graph2->getOffset()); h++)
						if(!graph2->dContained(dest-graph2->getOffset(), h) && !graph2->sContained(dest-graph2->getOffset(), h))
							ncTally2++;

					if(ncTally2 > ncTally1)
					{
						nodes[currNode-minNode] = 0;
						break;
					}else{
						if(ncTally1 > ncTally2)
						{
							nodes[dest-minNode] = 0;

						}else{
							if(ncTally2 == ncTally1 && currNode < dest)
							{
								nodes[currNode-minNode] = 0;
								break;
							}
						}
					}

					//	dContained = 0; doveTails = 0;
				}

				if(graph->sContained(source, j) && !graph->dContained(source, j))
				{
					nodes[currNode-minNode] = 0;
					break;
				}

				if(!graph->sContained(source, j) && graph->dContained(source, j))
				{
					nodes[dest-minNode] = 0;
				} 

				/*
				   if(!graph->sContained(source, j) && graph->dContained(source, j))
				   dContained++;

				   if(!graph->sContained(source, j) && !graph->dContained(source, j))	
				   doveTails++;

				   if(dContained > 3 || doveTails > 3)
				   break;
				   */
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

	char * fileName = new char [1000];
	sprintf(fileName, "%s/OvlGraph/finalMap", workDir.c_str());
	FILE * pFile = fopen(fileName, "r");
	fseek(pFile, 0, SEEK_END);
	int mSize = ftell(pFile)/sizeof(int);
	fclose(pFile);

	int * dNodes = new int [mSize];
	for(int i = 0; i < mSize; i++)
		dNodes[i] = 0;

	sprintf(fileName, "%s/HybSpectrum/deletedNodes", workDir.c_str());
	pFile = fopen(fileName, "w");
	fwrite(dNodes, sizeof(int), mSize, pFile);
	fclose(pFile);

	delete [] dNodes;
	delete [] fileName;

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

//void removeFile(const char [])
//Description: This function calls the linux rm command
//Input: fileName (const char [])
//Output:None
//Return:None
void removeFile(const char fileName []) {

	char * sysCmd = new char [1000];
	sprintf(sysCmd, "rm %s", fileName);
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

	char * contigDir = new char [1000];
	sprintf(contigDir, "%s/Contigs", workDir.c_str());

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

	int * cMap = new int [2 * numNodes];

	char * hybGraphDir = new char [1000];
	sprintf(hybGraphDir, "%s/HybSpectrum", workDir.c_str());

	int numGraphs = countFiles(graphDir);

	int count = 0; int sCount = 0; 
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

				if(graphIter == numGraphs-1)
					cMap[sCount++] = j/2;
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

		if(i == (graphIter + 1))
		{
			for(int j = 0; j < cHold; j++)
			{
				if(contigs[j] == -1)
				{
					cMap[sCount++] = hybMap[buff[contigs[(j-1)] * 2]];

					if(buff[contigs[(j-1)]*2+1] != -1 && (hybMap[buff[contigs[(j-1)]*2+1]] != hybMap[buff[contigs[(j-1)] * 2]]))
						cMap[sCount++] = hybMap[buff[contigs[(j-1)]*2+1]];
				}
			}
		}

		delete [] contigs;

		contigs = tmp;

		delete [] buff;
		delete [] hybMap;
		delete [] fileName;
	}

	int * cMap_rev = new int [sCount];
	for(int i = 0; i < sCount; i++)
	{
		cMap_rev[cMap[i]] = i; 
	}

	delete [] cMap;

	sprintf(contigDir, "%s/cMap%d", contigDir, graphIter);
	FILE * pFile = fopen(contigDir, "w");
	fwrite(cMap_rev, sizeof(int), sCount, pFile);

	fclose(pFile); 

	delete [] cMap_rev;
	delete [] contigDir;
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

