#include <iostream>
#include <sys/stat.h>
#include <map>
#include <fstream>
#include <string>
#include <algorithm>
#include <cstdio>
#include <stdint.h>
#include <vector>
#include <ctime>
#include <list>
#include <set>
#include <climits>
#include <cmath>
#include <queue>
#include <dirent.h>
#include "mpi.h"
using namespace std;
#include "Dictionary.h"
#include "Fragment_Index.h"
#include "Graph.h"

//Size unit
#define megabyte 1024000

char arr []  = {'A', 'C', 'T', 'G', '-'};

#define pathTag 1
#define endTag 0
#define help 911

void help_message(){

	cout<<"\n\n";
	cout<<"\t\t\t Path Finder Algorithm (Distributed)"<<endl;
	cout<<"\t\t\tWritten by Julia Warnke-Sommer\n"<<endl;
	cout<<"Command line arguments\n"<<endl;
	cout<<"--workDir :working directory :required"<<endl;
	cout<<"--distributedDir : the distributed directory "<<endl;
	cout<<"--graphIter : The current graph iteration"<<endl;
	cout<<"--minDensity : The minimum inter-cluster d"<<endl;
	cout<<"--percentMerged : The minimum percentage of nodes that have been merged "<<endl;
	cout<<"for the algorithm to continue merging"<<endl;
	cout<<"--verbose : verbose mode"<<endl;

	cout<<"\n\n"<<endl;

	cout<<"Exiting program"<<endl;
	exit(help);
}

int countFiles(const char []);
void transitiveReduce(Graph ** &, int, set<int> &, uint32_t * &, int);

int main(int argc, char * argv [])
{	
	int rank;
	int nTasks;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

	//User Input Flags
	string workDirS = "--workDir";               //Current Working Directory
	string minNodeWeightS = "--minNodeWeight";         //Minimum path length
	string minEdgeWeightS = "--minEdgeWeight";
	string verboseS = "--verbose";		     //verbose mode
	string distribDirS = "--DistributedGraph";   //Distributed Graph
	string timeS = "--time"; //Report the time

	string workDir = "";
	string distribDir = "";
	int minNodeWeight = 5;
	int minEdgeWeight = 2000;
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

		//The minimum node weight
		if(argv[i] == minNodeWeightS){
			minNodeWeight = atoi(argv[i+1]);
		}

		//The minimum edge weight
		if(argv[i] == minEdgeWeightS){
			minEdgeWeight = atoi(argv[i+1]);
		}

		//Verbose or not
		if(argv[i] == verboseS && verbose != true)
		{
			verbose = true;
			i--;
		}

		//Getting the distributed directory
		if(argv[i] == distribDirS)
		{
			distribDir = argv[i+1];		
		}

		//Display time
		if(argv[i] == timeS && reportTime != true)
		{
			reportTime = true;
			i--;
		}
	}

	int numPartitions = 0;

	char * distHybGraphDir = new char [1000];
	sprintf(distHybGraphDir, "%s/%s/Graph0", workDir.c_str(), distribDir.c_str());
	numPartitions = countFiles(distHybGraphDir);
	delete [] distHybGraphDir;

	if(rank == 0)
	{
		if(distribDir == "" || workDir == "")
		{
			help_message();
		}

		if(reportTime == true)
		{
			time_t rawtime;
			struct tm * timeinfo;
			time(&rawtime);
			timeinfo = localtime(&rawtime);
			cerr<<"start time: "<<asctime(timeinfo)<<endl;
		}

		if(numPartitions != (nTasks-1))
			help_message(); 

		char * pathDir = new char [1000];
		sprintf(pathDir, "%s/Paths", workDir.c_str());

		mkdir(pathDir, S_IRWXU);	

		delete [] pathDir;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if(rank != 0)
	{
		char * fileName = new char [1000];

		int part = log2(numPartitions); int section = rank-1;
		sprintf(fileName, "%s/%s/Graph0/Partition%d_%d", workDir.c_str(), distribDir.c_str(), part, section);

		FILE * pFile = fopen(fileName, "r");
		fseek(pFile, 0, SEEK_END);
		int size = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);

		int * partNodes = new int [size];
		fread(partNodes, sizeof(int), size, pFile);
		fclose(pFile);

		set<int> pNodes;
		for(int i = 0; i < size; i++)
		{
			pNodes.insert(partNodes[i]);
		}
		delete [] partNodes;

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

		set<int> visited; 
		int * outBuff = new int[10*megabyte];
		int currPos = 0;
		int totalPos = 0;

		int * outBuff2 = new int[10*megabyte]; 
		int currPos2 = 0;

		int * outBuff3 = new int[10* megabyte];
		int currPos3 = 0;

		vector<int> forwardPath;
		vector<int> backPath; 

		for(int i = 0; i < graphsLoaded; i++)
		{
			Graph * graph = graphs[i];
			for(int j = 0; j < graph->getNumNodes(); j++)
			{
				if(pNodes.find(j+graph->getOffset()) != pNodes.end() && visited.find(j+graph->getOffset()) == visited.end())
				{
					Graph * cGraph = graph;
					int currNode = j+cGraph->getOffset();
					visited.insert(currNode);
					int inDegree = 1; int outDegree = 0; int inEdge = 0; int outEdge = 0;
					forwardPath.push_back(currNode);

					for(int k = 0; k < cGraph->nodeDegree(currNode - cGraph->getOffset()); k++)
					{
						if(cGraph->edgeInGraph(currNode - cGraph->getOffset(), k) && !cGraph->rDovetail(currNode - cGraph->getOffset(), k))
						{
							outDegree++; outEdge = k;
						}
					}

					while(inDegree == 1 && outDegree == 1)
					{
						inDegree = 0; outDegree = 0;
						int prevNode = currNode;

						currNode = cGraph->edgeDest(currNode - cGraph->getOffset(), outEdge);
						for(int j = 0; j < graphsLoaded; j++)
						{
							cGraph = graphs[j];
							if(cGraph->getOffset() <= currNode && currNode < cGraph->getOffset() + cGraph->getNumNodes())
								break;
						}	

						for(int k = 0; k < cGraph->nodeDegree(currNode - cGraph->getOffset()); k++)
						{
							if(cGraph->edgeInGraph(currNode - cGraph->getOffset(), k) && !cGraph->rDovetail(currNode - cGraph->getOffset(), k))
							{
								outDegree++; outEdge = k;
							}

							if(cGraph->edgeInGraph(currNode - cGraph->getOffset(), k) && cGraph->rDovetail(currNode - cGraph->getOffset(), k))
							{
								inDegree++; inEdge = k;
							}
						}

						if(inDegree == 1)
						{
							if(visited.find(currNode) == visited.end())
							{
								visited.insert(currNode);
							}else{
								break;
							}

							if(pNodes.find(currNode) != pNodes.end())
							{
								forwardPath.push_back(currNode);
							}else{
								if(currPos3+2 >= 10 * megabyte)
								{
									char * fileName = new char [1000]; 
									sprintf(fileName, "%s/Paths/connections_%d", workDir.c_str(), rank-1);

									FILE * pFile = fopen(fileName, "a");
									fwrite(outBuff3, sizeof(int), currPos3, pFile);
									fclose(pFile);

									currPos3 = 0;

									delete [] fileName;
								}

								outBuff3[currPos3++] = prevNode; outBuff3[currPos3++] = currNode;
								break;
							}
						}
					}		

					inDegree = 0; outDegree = 1;	

					cGraph = graph;
					currNode = j+cGraph->getOffset();

					for(int k = 0; k < cGraph->nodeDegree(currNode - cGraph->getOffset()); k++)
					{
						if(cGraph->edgeInGraph(currNode - cGraph->getOffset(), k) && cGraph->rDovetail(currNode - cGraph->getOffset(), k))
						{
							inDegree++; inEdge = k;
						}
					}

					while(inDegree == 1 && outDegree == 1)
					{
						inDegree = 0; outDegree = 0;
						int prevNode = currNode;

						currNode = cGraph->edgeDest(currNode - cGraph->getOffset(), inEdge);
						for(int j = 0; j < graphsLoaded; j++)
						{
							cGraph = graphs[j];
							if(cGraph->getOffset() <= currNode && currNode < cGraph->getOffset() + cGraph->getNumNodes())
								break;
						}


						for(int k = 0; k < cGraph->nodeDegree(currNode - cGraph->getOffset()); k++)
						{
							if(cGraph->edgeInGraph(currNode - cGraph->getOffset(), k) && !cGraph->rDovetail(currNode - cGraph->getOffset(), k))
							{
								outDegree++; outEdge = k;
							}

							if(cGraph->edgeInGraph(currNode - cGraph->getOffset(), k) && cGraph->rDovetail(currNode - cGraph->getOffset(), k))
							{
								inDegree++; inEdge = k;
							}
						}


						if(outDegree == 1)
						{
							if(visited.find(currNode) == visited.end())
							{
								visited.insert(currNode);
							}else{
								break;
							}

							if(pNodes.find(currNode) != pNodes.end())
							{
								backPath.push_back(currNode);
							}else{

								if(currPos3+2 >= 10 * megabyte)
								{
									char * fileName = new char [1000]; 
									sprintf(fileName, "%s/Paths/connections_%d", workDir.c_str(), rank-1);

									FILE * pFile = fopen(fileName, "a");
									fwrite(outBuff3, sizeof(int), currPos3, pFile);
									fclose(pFile);

									currPos3 = 0;

									delete [] fileName;
								}

								outBuff3[currPos3++] = currNode; outBuff3[currPos3++] = prevNode;
								break;
							}
						}
					}

					if(currPos2 + 4 >= 10 * megabyte)
					{
						char * fileName = new char [1000];
						sprintf(fileName, "%s/Paths/index_%d", workDir.c_str(), rank-1);

						FILE * pFile = fopen(fileName, "a");
						fwrite(outBuff2, sizeof(int), currPos2, pFile);
						fclose(pFile);

						currPos2 = 0;

						delete [] fileName;
					}

					outBuff2[currPos2++] = forwardPath.back(); outBuff2[currPos2++] = totalPos;

					if(backPath.size() > 0)
					{
						outBuff2[currPos2++] = backPath.back(); outBuff2[currPos2++] = totalPos; 
					}

					if(currPos+forwardPath.size() + backPath.size() + 1 >= 10 * megabyte)
					{
						char * fileName = new char [1000];
						sprintf(fileName, "%s/Paths/path_%d", workDir.c_str(), rank-1);

						FILE * pFile = fopen(fileName, "a");
						fwrite(outBuff, sizeof(int), currPos, pFile);
						fclose(pFile);

						currPos = 0;

						delete [] fileName;
					}

					outBuff[currPos++] = forwardPath.size() + backPath.size(); totalPos++;

					for(vector<int>::reverse_iterator it = backPath.rbegin(); it != backPath.rend(); it++)
					{
						outBuff[currPos++] = (*it); totalPos++;
					}  	

					for(vector<int>::iterator it = forwardPath.begin(); it != forwardPath.end(); it++)
					{
						outBuff[currPos++] = (*it); totalPos++;
					}

					forwardPath.clear(); backPath.clear();
				}
			} 
		}

		if(currPos > 0)
		{
			char * fileName = new char [1000];
			sprintf(fileName, "%s/Paths/path_%d", workDir.c_str(), rank-1);

			FILE * pFile = fopen(fileName, "a");
			fwrite(outBuff, sizeof(int), currPos, pFile);
			fclose(pFile);
		}


		if(currPos2 > 0)
		{
			char * fileName = new char [1000];
			sprintf(fileName, "%s/Paths/index_%d", workDir.c_str(), rank-1);

			FILE * pFile = fopen(fileName, "a");
			fwrite(outBuff2, sizeof(int), currPos2, pFile);
			fclose(pFile);

		}


		if(currPos3 > 0)
		{
			char * fileName = new char [1000];
			sprintf(fileName, "%s/Paths/connections_%d", workDir.c_str(), rank-1);

			FILE * pFile = fopen(fileName, "a");
			fwrite(outBuff3, sizeof(int), currPos3, pFile);
			fclose(pFile);
		}

	

		int handshake = 0;
		MPI_Send(&handshake, 1, MPI_INT, 0, pathTag, MPI_COMM_WORLD);
	}

	if(rank == 0)
	{
		for(int i = 0; i < nTasks-1; i++)
		{
			MPI_Status status; int handshake = 0;
			MPI_Recv(&handshake, 1, MPI_INT, MPI_ANY_SOURCE, pathTag, MPI_COMM_WORLD, &status);
			cout<<"REcieved "<<endl;
		}


		int iSize = 0; int pSize = 0;
		for(int i = 0; i < nTasks-1; i++)
		{
			char * fileName = new char [1000];
			sprintf(fileName, "%s/Paths/path_%d", workDir.c_str(), i);

			FILE * pFile = fopen(fileName, "r");
			fseek(pFile, 0, SEEK_END);
			pSize+=ftell(pFile)/sizeof(int);				
			fclose(pFile);

			sprintf(fileName, "%s/Paths/index_%d", workDir.c_str(), i);
			pFile = fopen(fileName, "r");
			fseek(pFile, 0, SEEK_END);
			iSize+=ftell(pFile)/sizeof(int);
			fclose(pFile);

			delete [] fileName;		
		}

		int readPos1 = 0; int readPos2 = 0; 
		int * paths = new int[pSize];
		int * index = new int [iSize]; 

		map<int, int> links;
		map<int, int> links_rev;

		for(int i = 0; i < nTasks-1; i++)
		{
			char * fileName = new char [1000];
			sprintf(fileName, "%s/Paths/path_%d", workDir.c_str(), i);

			FILE * pFile = fopen(fileName, "r");
			fseek(pFile, 0, SEEK_END);
			int readSize1 = ftell(pFile)/sizeof(int);

			fseek(pFile, 0, SEEK_SET);

			fread(paths+readPos1, sizeof(int), readSize1, pFile);
			fclose(pFile);

			sprintf(fileName, "%s/Paths/index_%d", workDir.c_str(), i);
			pFile = fopen(fileName, "r");
			fseek(pFile, 0, SEEK_END);
			int readSize2 = ftell(pFile)/sizeof(int);

			fseek(pFile, 0, SEEK_SET);

			fread(index+readPos2, sizeof(int), readSize2, pFile);
			fclose(pFile);

			for(int j = 0; j < readSize2; j+=2)
			{
				index[readPos2+j+1]+=readPos1;
			}

			readPos1+=readSize1;
			readPos2+=readSize2;

			sprintf(fileName, "%s/Paths/connections_%d", workDir.c_str(), i);

			pFile = fopen(fileName, "r");
			fseek(pFile, 0, SEEK_END);
			int readSize3 = ftell(pFile)/sizeof(int);
			fseek(pFile, 0, SEEK_SET);

			int * tmpBuff = new int [readSize3];
			fread(tmpBuff, sizeof(int), readSize3, pFile);
			fclose(pFile);

			for(int j = 0; j < readSize3; j+=2)
			{
				links.insert(pair<int, int>(tmpBuff[j], tmpBuff[j+1]));
				links_rev.insert(pair<int, int>(tmpBuff[j+1], tmpBuff[j]));
			}
			delete [] tmpBuff;
			delete [] fileName;
		}

		map<int, int> locs;
		for(int i = 0; i < iSize; i+=2)
		{
			locs.insert(pair<int, int>(index[i], index[i+1]));
		} 
		delete [] index;

		int * outBuff = new int [pSize];
		int outPos = 0;
		int currPos = 0; 

		set<int> visited;

		while(currPos < pSize)
		{
			int cPathLen = paths[currPos];
			vector<int> backPath;

			int currNode = paths[currPos+1]; int nLoc = currPos;
			if(visited.find(nLoc) == visited.end())
			{
				while(links_rev.find(currNode) != links_rev.end() && visited.find(nLoc) == visited.end())
				{
					visited.insert(nLoc);
					backPath.push_back(nLoc);

					currNode = (*links_rev.find(currNode)).second;
					nLoc = (*locs.find(currNode)).second;
					currNode = paths[nLoc+1];
				}

				if(visited.find(nLoc) == visited.end())
				{
					visited.insert(nLoc);
					backPath.push_back(nLoc);				
				}

				int mergePathLen = 0;
				for(vector<int>::reverse_iterator it = backPath.rbegin(); it != backPath.rend(); it++)
				{
					int nLoc_t = (*it);
					for(int i = 0; i < paths[nLoc_t]; i++)
					{
						outBuff[outPos+1+mergePathLen] = paths[nLoc_t+i+1];
						mergePathLen++;
					}
				}

				vector<int> forwardPath;
				currNode = paths[currPos+paths[currPos]];
				if(links.find(currNode) != links.end())
				{
					currNode = (*links.find(currNode)).second; nLoc = (*locs.find(currNode)).second;
					currNode = paths[nLoc + paths[nLoc]];
					while(links.find(currNode) != links.end() && visited.find(nLoc) == visited.end())
					{
						visited.insert(nLoc);
						forwardPath.push_back(nLoc);

						currNode = (*links.find(currNode)).second;
						nLoc = (*locs.find(currNode)).second;
						currNode = paths[nLoc + paths[nLoc]];
					}				

					if(visited.find(nLoc) == visited.end())
					{
						visited.insert(nLoc);
						forwardPath.push_back(nLoc);				
					}
				}

				for(vector<int>::iterator it = forwardPath.begin(); it != forwardPath.end(); it++)
				{
					int nLoc = (*it);
					for(int i = 0; i < paths[nLoc]; i++)
					{
						outBuff[outPos+1+mergePathLen] = paths[nLoc+i+1];
						mergePathLen++;
					}
				}

				outBuff[outPos] = mergePathLen;
				outPos+=(mergePathLen+1);

			}

			currPos+=(paths[currPos] + 1);
		}

		for(int i = 0; i < nTasks-1; i++)
		{
			char * fileName = new char[1000];
			sprintf(fileName, "%s/Paths/path_%d", workDir.c_str(), i);

			remove(fileName);

			sprintf(fileName, "%s/Paths/index_%d", workDir.c_str(), i);	
			remove(fileName);

			sprintf(fileName, "%s/Paths/connections_%d", workDir.c_str(), i);
			remove(fileName);

			delete [] fileName;
		}

		char * fileName = new char[1000];
		sprintf(fileName, "%s/Paths/paths", workDir.c_str());

		FILE * pFile = fopen(fileName, "w");
		fwrite(outBuff, sizeof(int), outPos, pFile);
		fclose(pFile); 

		delete [] outBuff;
		
		if(reportTime == true)
		{
			time_t rawtime;
			struct tm * timeinfo;
			time(&rawtime);
			timeinfo = localtime(&rawtime);
			cerr<<"end time: "<<asctime(timeinfo)<<endl;
		}
	}	


	MPI_Finalize();

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

	return numFiles;
}
