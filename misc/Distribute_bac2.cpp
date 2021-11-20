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
#include "mpi.h"
using namespace std;
#include "Dictionary.h"
#include "Graph.h"

//MPI Tags
#define greedyTag 1
#define quitTag 2
#define multiKTag 3
#define help 911
#define bufOverflow 792
#define PartOutOfRange 2615

//The help message
void help_message(){

	cout<<"\n\n";
	cout<<"\t\t\t Merge and Traverse Graph Partitioning Algorithm"<<endl;
	cout<<"\t\t\tWritten by Julia Warnke-Sommer\n"<<endl;
	cout<<"Command line arguments\n"<<endl;
	cout<<"--workDir :working directory :required"<<endl;
	cout<<"--verbose : verbose mode"<<endl;
	cout<<"\n\n";

	cout<<"\n\n"<<endl;

	cout<<"Exiting program"<<endl;
	MPI_Abort(MPI_COMM_WORLD, help);
}

//Functions
int countFiles(const char []);
void initiatePartition(string);		
void greedyGraphGrow(string, int, int, int, int);
void projectPart(string, int, int, int, int);
void multiKernighanLin(string, int, int);
void kernighanLin(set<int> &, set<int> &, int &, int &, int * &, int, Graph ** &, int);
int getEdgeCut(string, Graph ** &, int,  set<int> &, set<int> &);


int main(int argc, char * argv [])
{
	//Set up the MPI
	int rank;
	int nTasks;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nTasks);

	string workDir = "NO_WORK_DIR";
	bool reportTime = false;
	int numPartitions = 2;
	int maxIter = 10;
	bool verbose = false;

	string workDirS = "--workDir"; //Current Working Directory
	string numPartitionsS = "--numPartitions"; //number of partitions to split graph into
	string maxIterS = "--maxIter";
	string verboseS = "--verbose";
	string timeS = "--time"; //Report the time

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
		if(argv[i] == numPartitionsS){
			numPartitions = atoi(argv[i+1]);


			if(numPartitions > 1024)
			{
				cerr<<"Please enter number less than 1024"<<endl;
				exit(PartOutOfRange);
			}

			int pExp = log2(numPartitions);
			numPartitions = pow(2, pExp);

			cout<<"Setting number of partitions to "<<numPartitions<<endl;
		}

		if(argv[i] == maxIterS){
			maxIter = atoi(argv[i]);
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


	//User didn't provided the working directory
	if(workDir == "NO_WORK_DIR")
	{
		help_message();
		exit(help);
	}


	if(rank == 0)
	{
		cerr<<"Graph Partitioning Algorithm"<<endl;
		cerr<<"Written by Julia Warnke-Sommer"<<endl;
		cerr<<endl;
		cerr<<"Graph Parameters Parameters"<<endl;
		cerr<<endl;
		cerr<<"Working directory : "<<workDir<<endl;
		cerr<<"Verbose mode : ";

		if(verbose) 
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

	if(rank == 0)
	{
		//Creating necessary directories
		int numGraphLevels = 0;
		{
			//Create the Distributed directory
			char * arr = new char [1000];
			sprintf(arr, "%s/DistributedGraph", workDir.c_str());
			mkdir(arr, S_IRWXU);

			sprintf(arr, "%s/DistributedGraph", workDir.c_str());
			mkdir(arr, S_IRWXU);

			char * graphDir = new char [1000];
			sprintf(graphDir, "%s/HybSpectrum", workDir.c_str());
			numGraphLevels = countFiles(graphDir);

			//Create the Distributed Graphs directories
			for(int i = 0; i < numGraphLevels; i++)
			{	
				sprintf(graphDir, "%s/DistributedGraph/Graph%d", workDir.c_str(), i);
				mkdir(graphDir, S_IRWXU);
			}

			delete [] graphDir;
			delete [] arr;
		}



		//Initial Greedy Graph Growing Partitioning
		{
			//This holds the jobs available 
			queue <pair<int, int> > roundAndSection;

			initiatePartition(workDir);		
			roundAndSection.push(pair<int, int>(0, 0));
			queue<int> availProcs;

			for(int i = 1; i < nTasks; i++)
			{
				availProcs.push(i);
			}

			int message [] = {0, 0, 0}; MPI_Status status; set<int> procsRecall; 
			while(!roundAndSection.empty())
			{
				while(!roundAndSection.empty() && !availProcs.empty())
				{
					pair<int, int> nextJob = roundAndSection.front(); roundAndSection.pop();
					int proc = availProcs.back(); availProcs.pop(); 
					procsRecall.insert(proc);

					message[0] = nextJob.first; message[1] = nextJob.second; message[2] = numGraphLevels-1; 
					MPI_Send(message, 3, MPI_INT, proc, greedyTag, MPI_COMM_WORLD);
				}	

				MPI_Recv(message, 3, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				cout<<"recieved "<<message[0]<<" "<<message[1]<<" "<<message[2]<<" FROM "<<status.MPI_SOURCE<<endl;
				availProcs.push(status.MPI_SOURCE); 

				procsRecall.erase(status.MPI_SOURCE);

				if(message[0]+1 <= log2(numPartitions))
				{
					cout<<"PUSH "<<message[0]+1<<" "<<message[1]*2<<endl;
					cout<<"PUSH "<<message[0]+1<<" "<<message[1]*2+1<<endl;
					roundAndSection.push(pair<int, int> (message[0]+1, message[1]*2));
					roundAndSection.push(pair<int, int> (message[0]+1, message[1]*2+1));
				}
			}

			for(set<int>::iterator it = procsRecall.begin(); it != procsRecall.end(); it++)
			{
				MPI_Recv(message, 3, MPI_INT, (*it), MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			}

			procsRecall.clear(); 
			MPI_Abort(MPI_COMM_WORLD, help);


			int cGraphLevel = 0;
			for(; cGraphLevel < numGraphLevels && !availProcs.empty(); cGraphLevel++)
			{		
				message[0] = cGraphLevel; int proc = availProcs.back(); availProcs.pop();
				MPI_Send(message, 3, MPI_INT, proc, multiKTag, MPI_COMM_WORLD);							

				procsRecall.insert(proc);
			}

			while(cGraphLevel < numGraphLevels)
			{
				MPI_Recv(message, 3, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				MPI_Send(message, 3, MPI_INT, status.MPI_SOURCE, multiKTag, MPI_COMM_WORLD);	

				cGraphLevel++;
			}

			for(set<int>::iterator it = procsRecall.begin(); it != procsRecall.end(); it++)
			{
				MPI_Recv(message, 3, MPI_INT, (*it), MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			}

			for(int i = 1; i < nTasks; i++)
			{
				int message [] = {0, 0, 0};
				MPI_Send(message, 3, MPI_INT, i, quitTag, MPI_COMM_WORLD);
			} 

			procsRecall.clear(); 
		}

	}else{

		while(true)
		{
			int message [] = {0, 0, 0}; MPI_Status status;
			MPI_Recv(message, 3, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);		

			if(status.MPI_TAG == greedyTag)
			{
				greedyGraphGrow(workDir, message[0], message[1], message[2], maxIter);	
				//projectPart(workDir, message[0], message[1], message[2], maxIter); 
				cout<<"RECEIVED "<<message[0]<<" "<<message[1]<<" "<<message[2]<<endl;

				MPI_Send(message, 3, MPI_INT, 0, greedyTag, MPI_COMM_WORLD);
			}

			if(status.MPI_TAG == multiKTag)
			{
				multiKernighanLin(workDir, message[0], maxIter);				
				MPI_Send(message, 3, MPI_INT, 0, multiKTag, MPI_COMM_WORLD);	
			}

			if(status.MPI_TAG == quitTag)
			{
				break;
			}
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

	//Log file stuff
	if(rank == 0)
	{
		cerr<<"Job is complete "<<endl;
		cerr<<"Thank you for running the Focus algorithm"<<endl;
	}

	MPI_Finalize();

	return 0;
}

int getEdgeCut(string workDir, Graph ** & graphs, int graphsLoaded,  set<int> & cPart1, set<int> & cPart2)
{
	int edgeCut = 0;
	for(set<int>::iterator it = cPart1.begin(); it != cPart1.end(); it++)
	{
		int node1 = (*it); 
		Graph * graph1;
		for(int l = 0; l < graphsLoaded; l++)
		{
			graph1 = graphs[l];
			if(graph1->getOffset() <= node1 && node1 < graph1->getOffset() + graph1->getNumNodes())
				break;
		}		

		for(int i = 0; i < graph1->nodeDegree(node1); i++)
		{
			int node2 = graph1->edgeDest(node1, i);
			if(cPart2.find(node2) != cPart2.end())
			{
				edgeCut+=graph1->edgeOvlLen(node1, i);	
			}
		}
	}

	return edgeCut;
}

//void greedyGraphGrow(int, int)
//Description: This function uses the greedy graph growing algorithm to
//create initial partitions in the graph
//Input: partition (int) : the current level of partitioning, section (int):
//the current section of the partition
//Output: None
//Return: None
void greedyGraphGrow(string workDir, int partition, int section, int graphLevel, int maxIter)  	
{	
	cout<<"HERE 1"<<endl;
	cout<<graphLevel<<" "<<maxIter<<endl;
	char * fileName = new char [1000];
	sprintf(fileName, "%s/DistributedGraph/Graph%d/Partition%d_%d", workDir.c_str(), graphLevel, partition, section);
	cout<<fileName<<endl;

	FILE * pFile = fopen(fileName, "r");

	fseek(pFile, 0, SEEK_END);
	int cpSize = ftell(pFile)/sizeof(int);
	fseek(pFile, 0, SEEK_SET);	

	int * currPart = new int [cpSize];
	fread(currPart, sizeof(int), cpSize, pFile);
	fclose(pFile);

	cout<<"HERE 1.5"<<endl;

	set< pair<int, int> > horizon;
	map<int, int> horizon_rev; 
	set<int> cPart1; set<int> cPart2;

	sprintf(fileName, "%s/HybSpectrum/Graph%d", workDir.c_str(), graphLevel);

	int numFiles = countFiles(fileName);
	int graphsLoaded = numFiles-2;

	Graph ** graphs = new Graph * [graphsLoaded];
	for(int i = 0; i < graphsLoaded; i++)
		graphs[i] = new Graph();

	cout<<"HERE 2"<<endl;

	for(int i = 0; i < graphsLoaded; i++)
	{
		FILE * pFile;

		char * fileName = new char [200];
		sprintf(fileName, "%s/HybSpectrum/Graph%d/graph%d", workDir.c_str(), graphLevel, i);

		pFile = fopen(fileName, "r");

		graphs[i]->read(pFile);

		fclose(pFile);
		delete [] fileName;
	}

	sprintf(fileName, "%s/HybSpectrum/Graph%d/nDen", workDir.c_str(), graphLevel);
	pFile = fopen(fileName, "r");

	fseek(pFile, 0, SEEK_END);
	int nDenSize = ftell(pFile)/sizeof(int);
	fseek(pFile, 0, SEEK_SET);

	int * nDen = new int [nDenSize];
	fread(nDen, sizeof(int), nDenSize, pFile);

	fclose(pFile);

	set<int> cPart;

	int totalWeight = 0;
	for(int i = 0; i < cpSize; i++)
	{
		int weight = nDen[currPart[i]];
		if(weight == 0) weight = 1;
		totalWeight+=weight;

		cPart.insert(currPart[i]);
	}

	int maxPartWeight = totalWeight/2; 
	int prevEdgeCut = INT_MAX; 

	int * ePart1 = '\0'; int * ePart2 = '\0'; int eSize1 = 0; int eSize2 = 0;

	cout<<"HERE 3"<<endl;

	for(int i = 0; i < 4; i++)
	{	
		int randStart = rand() % cpSize;
		int randEnd = randStart;

		int randNode = currPart[randStart++]; 
		horizon.insert(pair<int, int>(0, randNode));		
		horizon_rev.insert(pair<int, int>(randNode, 0));

		int part1Weight = 0; 
		cout<<"ITERATION "<<i<<endl;
		while(!horizon.empty() && part1Weight < maxPartWeight)
		{
			int node1 = (*horizon.rbegin()).second;

			set<pair<int, int> >::iterator erase = horizon.rbegin().base();
			set<pair<int, int> >::iterator test = erase;

			--test;

			horizon.erase(--erase); horizon_rev.erase(node1);

			Graph * graph1;
			for(int l = 0; l < graphsLoaded; l++)
			{
				graph1 = graphs[l];
				if(graph1->getOffset() <= node1 && node1 < graph1->getOffset() + graph1->getNumNodes())
					break;
			}

			int weight = nDen[node1]; if(weight == 0) weight = 1;
			if(weight + part1Weight <= maxPartWeight)
			{
				part1Weight+= weight; cPart1.insert(node1); 		
				for(int j = 0; j < graph1->nodeDegree(node1-graph1->getOffset()); j++)
				{
					int node2 = graph1->edgeDest(node1-graph1->getOffset(), j);	
					if(cPart.find(node2) != cPart.end() && cPart1.find(node2) == cPart1.end())
					{
						if(horizon_rev.find(node2) != horizon_rev.end())
						{
							int node2Weight = (*horizon_rev.find(node2)).second;
							set< pair<int, int> >::iterator foundIt = horizon.find(pair<int, int>(node2Weight, node2));

							int updateNode2Weight = (*foundIt).first + graph1->edgeOvlLen(node1-graph1->getOffset(), j); 
							pair<int, int> update = pair<int, int> (updateNode2Weight, node2);

							horizon.erase(foundIt); horizon.insert(update);

							(*horizon_rev.find(node2)).second = update.first;	
						}else{

							horizon_rev.insert(pair<int, int>(node2, graph1->edgeOvlLen(node1-graph1->getOffset(), j)));
							horizon.insert(pair<int, int>(graph1->edgeOvlLen(node1-graph1->getOffset(), j), node2));	
						}
					}
				}		
			}


			if(horizon.empty())
			{
				for(; randStart < cpSize && randStart != randEnd; randStart++)
				{
					int randNode = currPart[randStart];
					if(cPart1.find(randNode) == cPart1.end())
					{
						horizon.insert(pair<int, int>(0, randNode));		
						horizon_rev.insert(pair<int, int>(randNode, 0));
						break;
					}
				}

				if(randStart >= cpSize)
				{
					randStart = 0;
					for( ;randStart != randEnd; randStart++)
					{
						int randNode = currPart[randStart];
						if(cPart1.find(randNode) == cPart1.end())
						{
							horizon.insert(pair<int, int>(0, randNode));		
							horizon_rev.insert(pair<int, int>(randNode, 0));
							break;
						}
					}
				}
			}
		}

		int part2Weight = 0;
		for(int j = 0; j < cpSize; j++)
		{
			if(cPart1.find(currPart[j]) == cPart1.end())
			{
				cPart2.insert(currPart[j]);
				int weight = nDen[currPart[j]];

				if(weight == 0) weight = 1;
				part2Weight+=weight;
			}
		} 		

		cout<<cPart1.size()<<" "<<cPart2.size()<<endl;	
		cout<<part1Weight<<" "<<part2Weight<<endl;

		kernighanLin(cPart1, cPart2, part1Weight, part2Weight, nDen, maxIter, graphs, graphsLoaded);
		cout<<"RIGHT HERE "<<endl;
		int edgeCut = getEdgeCut(workDir, graphs, graphsLoaded, cPart1, cPart2); 
		cout<<"EDGE CUT "<<edgeCut<<endl;

		if(edgeCut < prevEdgeCut)
		{
			delete [] ePart1; delete [] ePart2;
			ePart1 = new int [cPart1.size()]; eSize1 = cPart1.size();
			ePart2 = new int [cPart2.size()]; eSize2 = cPart2.size();

			int j = 0; 
			for(set<int>::iterator it = cPart1.begin(); it != cPart1.end(); it++)
			{
				ePart1[j++] = (*it); 	
			}

			j = 0;
			for(set<int>::iterator it = cPart2.begin(); it != cPart2.end(); it++)
			{
				ePart2[j++] = (*it);
			}
		}

		cout<<"RIGHT HERE 2"<<endl;

		cPart1.clear(); cPart2.clear(); horizon.clear(); horizon_rev.clear();
	}

	sprintf(fileName, "%s/DistributedGraph/Graph%d/Partition%d_%d", workDir.c_str(), graphLevel, (partition+1), 2*section);
	pFile = fopen(fileName, "w");
	fwrite(ePart1, sizeof(int), cPart1.size(), pFile);
	fclose(pFile);

	sprintf(fileName, "%s/DistributedGraph/Graph%d/Partition%d_%d", workDir.c_str(), graphLevel, (partition+1), 2*section+1);
	pFile = fopen(fileName, "w");
	fwrite(ePart2, sizeof(int), cPart2.size(), pFile);
	fclose(pFile);

	sprintf(fileName, "%s/DistributedGraph/Graph%d/Partition%d_%d", workDir.c_str(), graphLevel, partition, section);

	delete [] ePart1; delete [] ePart2;

	remove(fileName);
}

//void projectPart(string & workDir, int, int, int, int)
////Description: This function projects the partitions on the 
////more granular graph levels 
////Input: partition (int) : the current level of partitioning, section (int):
////the current section of the partition
////Output: None
////Return: None
void projectPart(string workDir, int partition, int section, int graphLevel, int maxIter)
{
	char * fileName = new char[1000];
	sprintf(fileName, "%s/HybSpectrum/Graph%d/tmpMap", workDir.c_str(), graphLevel);		
	FILE * pFile = fopen(fileName, "r");
	fseek(pFile, 0, SEEK_END);

	int hybSize1 = ftell(pFile)/sizeof(int);

	fseek(pFile, 0, SEEK_SET);

	int * hybMap1 = new int [hybSize1];
	fread(hybMap1, sizeof(int), hybSize1, pFile);
	fclose(pFile);

	sprintf(fileName, "%s/HybSpectrum/Graph%d/tmpMap", workDir.c_str(), graphLevel-1);	
	pFile = fopen(fileName, "r"); 
	fseek(pFile, 0, SEEK_END);

	int hybSize2 = ftell(pFile)/sizeof(int);
	fseek(pFile, 0, SEEK_SET); 

	int * hybMap2 = new int [hybSize2];	
	fread(hybMap2,  sizeof(int), hybSize2, pFile);
	fclose(pFile);

	sprintf(fileName, "%s/Spectrum/Graph%d/map", workDir.c_str(), graphLevel);
	pFile = fopen(fileName, "r"); 
	fseek(pFile, 0, SEEK_END);

	int size = ftell(pFile)/sizeof(int);
	fseek(pFile, 0, SEEK_SET); 

	int * childMap = new int [size];

	fread(childMap, sizeof(int), size, pFile);	
	fclose(pFile);

	sprintf(fileName, "%s/DistributedGraph/Graph%d/Partition%d_%d", workDir.c_str(), graphLevel-1, partition, section);		
	pFile = fopen(fileName, "r");
	fseek(pFile, 0, SEEK_END);

	size = ftell(pFile)/sizeof(int);
	fseek(pFile, 0, SEEK_SET);

	int * prevPart_tmp = new int [size];
	fread(prevPart_tmp, sizeof(int), size, pFile);
	fclose(pFile);	

	set<int> prevPart;
	for(int i = 0; i < size; i++)
		prevPart.insert(prevPart_tmp[i++]);

	delete [] prevPart_tmp;

	sprintf(fileName, "%s/DistributedGraph/Graph%d/Partition%d_%d", workDir.c_str(), graphLevel, (partition+1), 2*section);		

	pFile = fopen(fileName, "r");       
	fseek(pFile, 0, SEEK_END);

	size = ftell(pFile)/sizeof(int);
	fseek(pFile, 0, SEEK_SET);

	int * cPart1_tmp = new int [size];
	fread(cPart1_tmp, sizeof(int), size, pFile);
	fclose(pFile);

	set<int> cPart;
	for(int i = 0; i < size; i++)
	{
		cPart.insert(cPart1_tmp[i]);
	}

	delete [] cPart1_tmp;

	set<int> pPart1;

	for(int i = 0; i < hybSize1; i++)
	{
		int currNode = hybMap1[i];
		if(cPart.find(currNode) != cPart.end())
		{
			if(prevPart.find(hybMap2[childMap[2*i]]) != prevPart.end())
				pPart1.insert(hybMap2[childMap[2*i]]);

			if(prevPart.find(hybMap2[childMap[2*i+1]]) != prevPart.end())
				if(hybMap2[childMap[2*i+1]] != hybMap2[childMap[2*i]])
					pPart1.insert(hybMap2[childMap[2*i+1]]); 	
		}			
	}

	cPart.clear();

	sprintf(fileName, "%s/DistributedGraph/Graph%d/Partition%d_%d", workDir.c_str(), graphLevel, (partition+1), 2*section);

	pFile = fopen(fileName, "r");
	fseek(pFile, 0, SEEK_END);

	size = ftell(pFile)/sizeof(int);
	fseek(pFile, 0, SEEK_SET);

	int * cPart2_tmp = new int [size];
	fread(cPart2_tmp, sizeof(int), size, pFile);
	fclose(pFile);

	for(int i = 0; i < size; i++)
	{
		cPart.insert(cPart2_tmp[i]);
	}

	delete [] cPart2_tmp;
	set<int> pPart2;

	for(int i = 0; i < hybSize1; i++)
	{
		int currNode = hybMap1[i];
		if(cPart.find(currNode) != cPart.end())
		{
			if(prevPart.find(hybMap2[childMap[2*i]]) != prevPart.end())
				pPart2.insert(hybMap2[childMap[2*i]]);

			if(prevPart.find(hybMap2[childMap[2*i+1]]) != prevPart.end())
				if(hybMap2[childMap[2*i+1]] != hybMap2[childMap[2*i]])
					pPart2.insert(hybMap2[childMap[2*i+1]]);               
		}
	}
	cPart.clear();		

	int numExtra = 0;	
	for(set<int>::iterator it = prevPart.begin(); it != prevPart.end(); it++)
	{
		if(pPart1.find((*it)) == pPart1.end() && pPart2.find((*it)) == pPart2.end())
		{
			if(numExtra%2)
			{
				pPart1.insert((*it));
			}else{
				pPart2.insert((*it)); 
			}
		}
	}	

	int * pPart1_tmp = new int [pPart1.size()];
	int currPos = 0;

	for(set<int>::iterator it = pPart1.begin(); it != pPart1.end(); it++)
	{
		pPart1_tmp[currPos++] = (*it);				
	}

	sprintf(fileName, "%s/DistributedGraph/Graph%d/Partition%d_%d", workDir.c_str(), graphLevel-1, (partition+1), 2*section);	

	pFile = fopen(fileName, "w");	
	fwrite(pPart1_tmp, sizeof(int), pPart1.size(), pFile);
	fclose(pFile);

	delete [] pPart1_tmp;

	int * pPart2_tmp = new int [pPart2.size()];
	currPos = 0;

	for(set<int>::iterator it = pPart2.begin(); it != pPart2.end(); it++)
	{
		pPart2_tmp[currPos++] = (*it);
	}

	sprintf(fileName, "%s/DistributedGraph/Graph%d/Partition%d_%d", workDir.c_str(), graphLevel-1, (partition+1), 2*section+1); 

	pFile = fopen(fileName, "w");
	fwrite(pPart2_tmp, sizeof(int), pPart2.size(), pFile);
	fclose(pFile);

	sprintf(fileName, "%s/DistributedGraph/Graph%d/Partition%d_%d", workDir.c_str(), graphLevel-1, partition, section);

	remove(fileName);

	delete [] fileName;
}


void multiKernighanLin(string workDir, int graphLevel, int maxIter)
{
	char * fileName = new char [1000];

	sprintf(fileName, "%s/HybSpectrum/Graph%d", workDir.c_str(), graphLevel);

	int numFiles = countFiles(fileName);
	int graphsLoaded = numFiles-2;

	Graph ** graphs = new Graph * [graphsLoaded];
	for(int i = 0; i < graphsLoaded; i++)
		graphs[i] = new Graph();

	for(int i = 0; i < graphsLoaded; i++)
	{
		FILE * pFile;

		char * fileName = new char [200];
		sprintf(fileName, "%s/HybSpectrum/Graph%d/graph%d", workDir.c_str(), graphLevel, i);

		pFile = fopen(fileName, "r");

		graphs[i]->read(pFile);

		fclose(pFile);
		delete [] fileName;
	}

	fileName = new char [1000];

	sprintf(fileName, "%s/DistributedGraph/Graph%d/", workDir.c_str(), graphLevel);
	int numSections = countFiles(fileName);

	int partition = log2(numSections);

	map<int, int> nodeToSection;

	for(int i = 0; i < numSections; i++)
	{
		sprintf(fileName, "%s/DistributedGraph/Graph%d/Partition%d_%d", workDir.c_str(), graphLevel, partition, i);
		FILE * pFile = fopen(fileName, "r");

		pFile = fopen(fileName, "r"); 
		fseek(pFile, 0, SEEK_END);

		int size = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET); 

		int * buff = new int [size];	
		fread(buff,  sizeof(int), size, pFile);
		fclose(pFile);

		for(int j = 0; j < size; j++)
		{
			nodeToSection.insert(pair<int, int>(buff[j], i));
		}
	}

	set<pair<int, pair<int, int> > > priority; map<int, pair<int, int> > revs_priority;
	for(map<int, int>::iterator it = nodeToSection.begin(); it != nodeToSection.end(); it++)
	{
		int scores[1024] = {0};
		int node1 = (*it).first;

		Graph * graph1;
		for(int l = 0; l < graphsLoaded; l++)
		{
			graph1 = graphs[l];
			if(graph1->getOffset() <= node1 && node1 < graph1->getOffset() + graph1->getNumNodes())
				break;
		}

		int max = 0; int maxIndex = (*it).second;

		for(int i = 0; i < graph1->nodeDegree(node1); i++)
		{
			int node2 = graph1->edgeDest(node1, i);
			int nSection = (*nodeToSection.find(node2)).second;					

			int weight = graph1->edgeOvlLen(node1, i);

			scores[nSection]+=weight;

			if(scores[nSection] >= max)
			{
				max = scores[nSection]; maxIndex = nSection;
			}
		}

		if(scores[maxIndex] - scores[(*it).second] > 0)
		{
			pair<int, int> locMax (node1, maxIndex);
			priority.insert(pair<int, pair<int, int> >(scores[maxIndex] - scores[(*it).second], locMax));		
			revs_priority.insert( pair<int, pair<int, int> > (node1, pair<int, int> (scores[maxIndex] - scores[(*it).second], maxIndex)));
		}
	}


	int numIter = 0; 
	while(numIter < maxIter)
	{	
		int negInARow = 0;  
		while(!priority.empty())
		{
			pair<int, pair<int, int> > update = (*priority.rbegin());
			set<pair<int, pair<int, int> > >::iterator erase = priority.rbegin().base();

			priority.erase(--erase);

			int node1 = update.second.first; 
			nodeToSection.erase(node1);
			nodeToSection.insert(pair<int, int>(node1, update.second.second));

			Graph * graph1;
			for(int l = 0; l < graphsLoaded; l++)
			{
				graph1 = graphs[l];
				if(graph1->getOffset() <= node1 && node1 < graph1->getOffset() + graph1->getNumNodes())
					break;
			}

			for(int j = 0; j < graph1->nodeDegree(j); j++)
			{
				int node2 = graph1->edgeDest(node1, j);
				int scores[1024] = {0};

				Graph * graph2;
				for(int l = 0; l < graphsLoaded; l++)
				{
					graph2 = graphs[l];
					if(graph2->getOffset() <= node2 && node2 < graph2->getOffset() + graph2->getNumNodes())
						break;
				}

				int max = 0; int maxIndex = (*nodeToSection.find(node2)).second;
				for(int k = 0; k < graph2->nodeDegree(node2); k++)
				{
					int node3 = graph2->edgeDest(node2, k);
					int nSection = (*nodeToSection.find(node3)).second;					

					int weight = graph2->edgeOvlLen(node2, k);

					scores[nSection]+=weight;

					if(scores[nSection] >= max)
					{
						max = scores[nSection]; maxIndex = nSection;
					}
				}

				if(revs_priority.find(node2) != revs_priority.end())
				{
					pair<int, pair<int, int> > foundIt = (*revs_priority.find(node2));
					pair<int, pair<int, int> > eraseIt (foundIt.second.first, pair<int, int>(foundIt.first, foundIt.second.second)); 
					priority.erase(eraseIt); revs_priority.erase(node2);
				}

				if(scores[maxIndex] - scores[(*nodeToSection.find(node2)).second] > 0)
				{
					pair<int, int> n (node2, maxIndex);
					pair<int, pair<int, int> > nUpdate (scores[maxIndex], n);

					priority.insert(nUpdate);

					pair<int, pair<int, int> > rnUpdate(node2, pair<int, int>(nUpdate.first, n.second));
					revs_priority.insert(rnUpdate); 
				}

			}
			numIter++;
		}
	}

	delete [] fileName;
}

//kernighanLin(set<int> &, set<int> &, Graph ** &, int)
//Description: This is the refinement algorithm
//Input: cPart1(set<int>): 1 Partition, cPart2 (set<int>): 2 Partition, graphs (Graph**): 
//The input graph, graphsLoaded (int): number of sections of graph loaded 
//Return:None
//Output:None 
void kernighanLin(set<int> & cPart1, set<int> & cPart2, int & cPart1Weight, int & cPart2Weight, int * & nDen, int maxIter, Graph ** & graphs, int graphsLoaded)
{
	set<pair<int, int> > c1Priority; set< pair<int, int> > c2Priority; 
	map<int, int> priority_rev; 

	cout<<" kernighanLin "<<endl;

	for(set<int>::iterator it = cPart1.begin(); it != cPart1.end(); it++)
	{
		int node1 = (*it);

		Graph * graph1;
		for(int l = 0; l < graphsLoaded; l++)
		{
			graph1 = graphs[l];
			if(graph1->getOffset() <= node1 && node1 < graph1->getOffset() + graph1->getNumNodes())
				break;
		}

		int c1Weight = 0; int c2Weight = 0;
		for(int k = 0; k < graph1->nodeDegree(node1 - graph1->getOffset()); k++)
		{
			int node2 = graph1->edgeDest(node1 - graph1->getOffset(), k);
			if(cPart1.find(node2) != cPart1.end())
			{
				c1Weight+=graph1->edgeOvlLen(node1 - graph1->getOffset(),k);
			}

			if(cPart2.find(node2) != cPart2.end())
			{
				c2Weight+=graph1->edgeOvlLen(node1 - graph1->getOffset(),k);
			}
		}	


		int weight = c2Weight - c1Weight;
		c1Priority.insert(pair<int, int>(weight, node1));
		priority_rev.insert(pair<int, int>(node1, weight));
	}

	for(set<int>::iterator it = cPart2.begin(); it != cPart2.end(); it++)
	{       
		int node1 = (*it);

		Graph * graph1;
		for(int l = 0; l < graphsLoaded; l++)
		{       
			graph1 = graphs[l];
			if(graph1->getOffset() <= node1 && node1 < graph1->getOffset() + graph1->getNumNodes())
				break;
		}


		int c1Weight = 0; int c2Weight = 0;
		for(int k = 0; k < graph1->nodeDegree(node1 - graph1->getOffset()); k++)
		{       
			int node2 = graph1->edgeDest(node1 - graph1->getOffset(), k);
			if(cPart1.find(node2) != cPart1.end())
			{       
				c1Weight+=graph1->edgeOvlLen(node1 - graph1->getOffset(),k);
			}

			if(cPart2.find(node2) != cPart2.end())
			{       
				c2Weight+=graph1->edgeOvlLen(node1 - graph1->getOffset(),k);
			}
		}

		int weight = c1Weight - c2Weight;
		c2Priority.insert(pair<int, int>(weight, node1));
		priority_rev.insert(pair<int, int>(node1, weight));
	}

	cout<<"LOADED Partitions"<<c1Priority.size()<<" "<<c2Priority.size()<<endl;

	int numIter = 0; bool maximumOpt = false;
	while(numIter < maxIter && maximumOpt != true)
	{	
		int negInARow = 0; vector< pair<int, int> > selectSwitch;  
		int maxd1 = 0; int maxd2 = 0; int maxCut = 0; int currCut = 0;

		while(!c1Priority.empty() && !c2Priority.empty() && negInARow < 50)
		{
			set<pair<int, int> >::iterator c1Priority_iterHold =  c1Priority.end(); 
			pair<int, int> coords; int maxGWeight = INT_MIN; int maxDWeight = INT_MIN+1; 

			while(maxDWeight > maxGWeight)
			{
				set<pair<int, int> >::iterator c1Priority_iter = c1Priority_iterHold;					
				--c1Priority_iter;

				maxDWeight = INT_MIN;
				for(set< pair<int, int> >::reverse_iterator c2Priority_iter = c2Priority.rbegin(); c2Priority_iter != c2Priority.rend() && c1Priority_iter != c1Priority.end(); c2Priority_iter++)
				{
					cout<<"IN HERE NOW "<<endl;
					pair<int, int> c1Next = (*c1Priority_iter);
					pair<int, int> c2Next = (*c2Priority_iter);

					cout<<c1Next.second<<" "<<c1Next.first<<endl;
					cout<<c2Next.first<<" "<<c2Next.second<<endl;

					int node1 = c1Next.second; int node2 = c2Next.second;

					int DWeight  = c1Next.first + c2Next.first;
					if(DWeight > maxDWeight)
						maxDWeight = DWeight;

					Graph * graph1;
					for(int l = 0; l < graphsLoaded; l++)
					{
						graph1 = graphs[l];
						if(graph1->getOffset() <= node1 && node1 < graph1->getOffset() + graph1->getNumNodes())
							break;
					}

	
					int cWeight = 0;
					for(int k = 0; k < graph1->nodeDegree(node1 - graph1->getOffset()); k++)
                			{
						if(graph1->edgeDest(node1 - graph1->getOffset(), k) == node2)
						{
							cWeight = graph1->edgeOvlLen(node1 - graph1->getOffset(), k);			
						}
					}

					int GWeight = c1Next.first + c2Next.first - 2 * cWeight;

					if(GWeight > maxGWeight)
					{
						maxGWeight = GWeight; 
 						coords.first = node1; coords.second = node2;
					}

					c1Priority_iter++;		
				}

				cout<<maxDWeight<<" "<<maxGWeight<<" MAX WEIGHTS "<<endl;
				c1Priority_iterHold--;
			}

			int tmpCPart1Weight = cPart1Weight - nDen[coords.first] + nDen[coords.second]; 
			int tmpCPart2Weight = cPart2Weight - nDen[coords.second] + nDen[coords.first];		

			int changeWeight = abs(tmpCPart1Weight - tmpCPart2Weight);
			int currWeight = abs(cPart1Weight - cPart2Weight);	

			int totalWeight = cPart1Weight + cPart2Weight;

			int padding = .03 * totalWeight;

			currCut+=maxGWeight; 
			if((currCut > maxCut && changeWeight <= padding) || (currCut == maxCut && changeWeight < currWeight))
			{
				cout<<tmpCPart1Weight<<" "<<tmpCPart2Weight<<" UPDATED WEIGHTS "<<endl;
				maxCut = currCut;
				negInARow = 0;
			}else{
				negInARow++; 
				cout<<"NEGATIVE!!"<<endl;
			}

			cPart1Weight = tmpCPart1Weight;  
			cPart2Weight = tmpCPart2Weight;

			selectSwitch.push_back(coords); 			
			int switchNode1 = coords.first;

			Graph * graphS;
			for(int l = 0; l < graphsLoaded; l++)
			{       
				graphS = graphs[l];
				if(graphS->getOffset() <= switchNode1 && switchNode1 < graphS->getOffset() + graphS->getNumNodes())
					break;
			}

			for(int j = 0; j < graphS->nodeDegree(switchNode1 - graphS->getOffset()); j++)
			{
				int updateNode = graphS->edgeDest(switchNode1 - graphS->getOffset(), j);			
				pair<int, int> update = (*priority_rev.find(updateNode));

				if(c1Priority.find(pair<int, int>(update.second, update.first)) != c1Priority.end() )
				{
					c1Priority.erase(pair<int, int>(update.second, update.first));
					update.second = update.second + graphS->edgeOvlLen(switchNode1 - graphS->getOffset(), j);
					c1Priority.insert(pair<int, int>(update.second, update.first));
				}	

				if(c2Priority.find(pair<int, int>(update.second, update.first)) != c2Priority.end())
				{       
					c2Priority.erase(pair<int, int>(update.second, update.first));
					update.second = update.second - graphS->edgeOvlLen(switchNode1 - graphS->getOffset(), j);
					c2Priority.insert(pair<int, int>(update.second, update.first));
				}
			}

			int switchNode2 = coords.second;
			for(int l = 0; l < graphsLoaded; l++)
			{       
				graphS = graphs[l];
				if(graphS->getOffset() <= switchNode2 && switchNode2 < graphS->getOffset() + graphS->getNumNodes())
					break;
			}

			for(int j = 0; j < graphS->nodeDegree(switchNode2 - graphS->getOffset()); j++)
			{
				int updateNode = graphS->edgeDest(switchNode2 - graphS->getOffset(), j);			
				pair<int, int> update = (*priority_rev.find(updateNode));

				if(c1Priority.find(pair<int, int>(update.second, update.first)) != c1Priority.end() )
				{
					c1Priority.erase(pair<int, int>(update.second, update.first));
					update.second = update.second - graphS->edgeOvlLen(switchNode2 - graphS->getOffset(), j);
					c1Priority.insert(pair<int, int>(update.second, update.first));
				}	

				if(c2Priority.find(pair<int, int>(update.second, update.first)) != c2Priority.end())
				{       
					c2Priority.erase(pair<int, int>(update.second, update.first));
					update.second = update.second + graphS->edgeOvlLen(switchNode2 - graphS->getOffset(), j);
					c2Priority.insert(pair<int, int>(update.second, update.first));
				}
			}
		}


		int currPos = 0;
		for(vector< pair<int, int> >::reverse_iterator it = selectSwitch.rbegin(); it !=  selectSwitch.rend(); it++)
		{
			if(currPos > negInARow) 
			{
				cPart1.erase((*it).first);  cPart1.insert((*it).second);
				cPart2.erase((*it).second); cPart2.insert((*it).first);	
			}else{
				cPart1Weight = cPart1Weight - nDen[(*it).second] + nDen[(*it).first];
				cPart2Weight = cPart2Weight - nDen[(*it).first] + nDen[(*it).second];
			} 
			currPos++;
		}

		if(maxCut == 0) maximumOpt = true;
		numIter++;
	}
}

//void initiatePartition(const char [])
////Description: This function initiatiates the first partition
////Input: workDir (string) : The working directory
////Output:None
////Return:None
void initiatePartition(string workDir)
{	
	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/HybSpectrum", workDir.c_str());

	int numGraphs = countFiles(graphDir);

	for(int i = 0; i < numGraphs; i++)
	{
		sprintf(graphDir, "%s/HybSpectrum/Graph%d", workDir.c_str(), i);

		char * fileName = new char [1000];

		sprintf(fileName, "%s/nDen", graphDir);

		FILE * pFile = fopen(fileName, "r");
		fseek(pFile, 0, SEEK_END);

		int graphSize = ftell(pFile)/sizeof(int);
		fclose(pFile);

		int * nPart = new int[graphSize];

		for(int j = 0; j < graphSize; j++)
		{
			nPart[j] = j;
		}

		sprintf(fileName, "%s/DistributedGraph/Graph%d/Partition0_0", workDir.c_str(), i);

		pFile = fopen(fileName, "w");

		fwrite(nPart, sizeof(int), graphSize, pFile);

		fclose(pFile);

		delete [] nPart;
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

