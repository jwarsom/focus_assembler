#include <iostream>
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
#include "Labels.h"
#include "Graph.h"

//Size unit
#define megabyte 1024000

char arr []  = {'A', 'C', 'T', 'G', '-'};

#define reduceTag 1
#define trimTag 2
#define popTag 3
#define endTag 0
#define help 911

void help_message(){

	cout<<"\n\n";
	cout<<"\t\t\t Trimming Algorithm"<<endl;
	cout<<"\t\t\tWritten by Julia Warnke\n"<<endl;
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
int trimGraph(Graph ** &, int, set<int> &, int * &, int);
int popBubbles(Graph ** &, int, set<int> &, int * &, int);

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
	string pathLenS = "--pathLen"; //The maximum path length
	string bubbleLenS = "--bubbleLen"; //The maximum bubble length
	string timeS = "--time"; //Report the time

	string workDir = "";
	string distribDir = "";
	int minNodeWeight = 5;
	int minEdgeWeight = 2000;
	int pathLen = 5;
	int bubbleLen = 5;
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

		//Getting the distributed directory
		if(argv[i] == pathLenS)
		{       
			pathLen = atoi(argv[i+1]);
		}

		//Getting the distributed directory
		if(argv[i] == bubbleLenS)
		{       
			bubbleLen = atoi(argv[i+1]);
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

		char * hybGraphDir = new char [1000];
		sprintf(hybGraphDir, "%s/HybSpectrum/Graph0", workDir.c_str());

		int numGraphs = countFiles(hybGraphDir) - 2;

		for(int j = 0; j < numGraphs; j++)
		{
			Graph graph; 
			char * fileName = new char [1000];
			sprintf(fileName, "%s/HybSpectrum/Graph0/graph%d", workDir.c_str(), j);
			FILE * pFile = fopen(fileName, "r");
			graph.read(pFile);
			fclose(pFile);

			for(int i = 0; i < graph.getNumNodes(); i++)
			{
				graph.setNodeInGraph(i);
			}

			for(int i = 0; i < graph.getNumEdges(); i++)
			{
				graph.setEdgeInGraph(i);
				if(graph.sContained(i))
				{
					graph.maskNodeFromGraph(graph.edgeSource(i));
				}
			}

			pFile = fopen(fileName, "w");
			graph.write(pFile);
			fclose(pFile);

			delete [] fileName;
		}

		delete [] hybGraphDir;
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

		char * contigInfoDir = new char [1000];
		sprintf(contigInfoDir, "%s/ContigsInfo", workDir.c_str());

		int numContigs = countFiles(contigInfoDir) - 2;

		char * graphDir = new char [1000];
		sprintf(graphDir, "%s/HybSpectrum/Graph0", workDir.c_str());

		int numGraphs = countFiles(graphDir)-2;

		sprintf(fileName, "%s/HybSpectrum/Graph0/nDen", workDir.c_str());

		pFile = fopen(fileName, "r");
		fseek(pFile, 0, SEEK_END);
		int nDenSize = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);

		int * nDen = new int [nDenSize];

		fread(nDen, sizeof(int), nDenSize, pFile);
		fclose(pFile);

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

		int numNodes = 0; 
		for(int i = 0; i < graphsLoaded; i++)
		{
			numNodes+=graphs[i]->getNumNodes();
		}

		uint32_t * sandBox = new uint32_t [numNodes];
		cout<<"TRANSITIVE "<<rank<<endl;
		transitiveReduce(graphs, graphsLoaded, pNodes, sandBox, numNodes);
		delete [] sandBox;
		cout<<"DONE TRANSITIVE "<<rank<<endl;

		int handshake = 0; MPI_Status status;

		MPI_Send(&handshake, 1, MPI_INT, 0, reduceTag, MPI_COMM_WORLD);
		MPI_Recv(&handshake, 1, MPI_INT, 0, reduceTag, MPI_COMM_WORLD, &status);


		int erasedEdges = 0;
		for(int i = 0; i < graphsLoaded; i++)
		{
			sprintf(fileName, "%s/HybSpectrum/Graph0/graph%d", workDir.c_str(), i);

			FILE * pFile = fopen(fileName, "r");
			Graph uGraph;
			uGraph.read(pFile);
			fclose(pFile);

			for(int j = 0; j < uGraph.getNumEdges(); j++)
			{
				if(!uGraph.edgeInGraph(j))
				{
					graphs[i]->maskEdgeFromGraph(j);
					erasedEdges++;
				}
			}

			pFile = fopen(fileName, "w");
			graphs[i]->write(pFile);

			fclose(pFile);
		}


		MPI_Send(&handshake, 1, MPI_INT, 0, reduceTag, MPI_COMM_WORLD);

		for(int i = 0; i < graphsLoaded; i++)
			delete graphs[i];

		MPI_Recv(&handshake, 1, MPI_INT, 0, trimTag, MPI_COMM_WORLD, &status);

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

		cout<<"START TRIM "<<endl;
		int numTrimmed =  trimGraph(graphs, graphsLoaded, pNodes, nDen, pathLen);		
		cout<<"NUM TRIMMED "<<numTrimmed<<endl;
		cout<<"DONE TRIM "<<endl;

		MPI_Send(&handshake, 1, MPI_INT, 0, trimTag, MPI_COMM_WORLD);
		MPI_Recv(&handshake, 1, MPI_INT, 0, trimTag, MPI_COMM_WORLD, &status);

		for(int i = 0; i < graphsLoaded; i++)
		{
			sprintf(fileName, "%s/HybSpectrum/Graph0/graph%d", workDir.c_str(), i);

			FILE * pFile = fopen(fileName, "r");
			Graph uGraph;
			uGraph.read(pFile);
			fclose(pFile);

			for(int j = 0; j < uGraph.getNumEdges(); j++)
			{
				if(!uGraph.edgeInGraph(j))
				{
					graphs[i]->maskEdgeFromGraph(j);
					erasedEdges++;
				}

				if(!uGraph.nodeInGraph(uGraph.edgeSource(j)))
				{
					graphs[i]->maskNodeFromGraph(uGraph.edgeSource(j));
				}
			}

			pFile = fopen(fileName, "w");
			graphs[i]->write(pFile);

			fclose(pFile);
		}

		MPI_Send(&handshake, 1, MPI_INT, 0, trimTag, MPI_COMM_WORLD);

		for(int i = 0; i < graphsLoaded; i++)
			delete graphs[i];

		MPI_Recv(&handshake, 1, MPI_INT, 0, popTag, MPI_COMM_WORLD, &status);	

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
		
		cout<<"POPPING BUBBLES "<<endl;
		int numPopped = popBubbles(graphs, graphsLoaded, pNodes, nDen, bubbleLen);	
		cout<<"DONE POP "<<endl;

		MPI_Send(&handshake, 1, MPI_INT, 0, popTag, MPI_COMM_WORLD);
		MPI_Recv(&handshake, 1, MPI_INT, 0, popTag, MPI_COMM_WORLD, &status);

		for(int i = 0; i < graphsLoaded; i++)
		{
			sprintf(fileName, "%s/HybSpectrum/Graph0/graph%d", workDir.c_str(), i);

			FILE * pFile = fopen(fileName, "r");
			Graph uGraph;
			uGraph.read(pFile);
			fclose(pFile);

			for(int j = 0; j < uGraph.getNumEdges(); j++)
			{
				if(!uGraph.edgeInGraph(j))
				{
					graphs[i]->maskEdgeFromGraph(j);
					erasedEdges++;
				}

				if(!uGraph.nodeInGraph(uGraph.edgeSource(j)))
				{
					graphs[i]->maskNodeFromGraph(uGraph.edgeSource(j));
				}
			}

			pFile = fopen(fileName, "w");
			graphs[i]->write(pFile);

			fclose(pFile);
		}

		MPI_Send(&handshake, 1, MPI_INT, 0, popTag, MPI_COMM_WORLD);

		for(int i = 0; i < graphsLoaded; i++)
			delete graphs[i];

		delete [] graphs;	

		delete [] nDen;
		delete [] graphDir; 
		delete [] contigInfoDir;	
	}

	if(rank == 0)
	{
		MPI_Status status; int handshake = 0;
		for(int i = 0; i < nTasks - 1; i++)
		{
			MPI_Recv(&handshake, 1, MPI_INT, MPI_ANY_SOURCE, reduceTag, MPI_COMM_WORLD, &status);
			MPI_Send(&handshake, 1, MPI_INT, status.MPI_SOURCE, reduceTag, MPI_COMM_WORLD);
			MPI_Recv(&handshake, 1, MPI_INT, status.MPI_SOURCE, reduceTag, MPI_COMM_WORLD, &status);
		}

		for(int i = 1; i < nTasks; i++)
		{
			MPI_Send(&handshake, 1, MPI_INT, i, trimTag, MPI_COMM_WORLD);
		}
		
		for(int i = 0; i < nTasks - 1; i++)
		{
			MPI_Recv(&handshake, 1, MPI_INT, MPI_ANY_SOURCE, trimTag, MPI_COMM_WORLD, &status);
			MPI_Send(&handshake, 1, MPI_INT, status.MPI_SOURCE, trimTag, MPI_COMM_WORLD);
			MPI_Recv(&handshake, 1, MPI_INT, status.MPI_SOURCE, trimTag, MPI_COMM_WORLD, &status);
		}

		for(int i = 1; i < nTasks; i++)
		{
			MPI_Send(&handshake, 1, MPI_INT, i, popTag, MPI_COMM_WORLD);
		}

		for(int i = 0; i < nTasks - 1; i++)
		{
			MPI_Recv(&handshake, 1, MPI_INT, MPI_ANY_SOURCE, popTag, MPI_COMM_WORLD, &status);
			MPI_Send(&handshake, 1, MPI_INT, status.MPI_SOURCE, popTag, MPI_COMM_WORLD);
			MPI_Recv(&handshake, 1, MPI_INT, status.MPI_SOURCE, popTag, MPI_COMM_WORLD, &status);
		}

		char * fileName = new char [1000];
		sprintf(fileName, "%s/OvlGraph", workDir.c_str());
		cout<<fileName<<endl;
		int numFiles = countFiles(fileName);
		int graphsLoaded = numFiles-1; //-2

		Graph ** graphs = new Graph * [graphsLoaded];
		for(int i = 0; i < graphsLoaded; i++)
			graphs[i] = new Graph();

		cout<<numFiles<<endl;

		for(int i = 0; i < graphsLoaded; i++)
		{
			FILE * pFile;

			char * fileName = new char [200];
			sprintf(fileName, "%s/OvlGraph/graph%d", workDir.c_str(), i);
			cout<<fileName<<endl;

			pFile = fopen(fileName, "r");

			graphs[i]->read(pFile);

			fclose(pFile);
			delete [] fileName;
		}

		for(int i = 0; i < graphsLoaded; i++)
		{
			Graph * graph = graphs[i];
			for(int j = 0; j < graph->getNumNodes(); j++)
			{
				if(graph->nodeInGraph(j))
				{
					for(int k = 0; k < graph->nodeDegree(j); k++)
					{
						if(graph->edgeInGraph(j, k))
						{
							int nNode = graph->edgeDest(j, k);
							Graph * graph2;
							for(int l = 0; l < graphsLoaded; l++)
							{
								graph2 = graphs[l];
								if(nNode >= graph2->getOffset() && nNode < graph2->getOffset() + graph2->getNumNodes())
									break; 
							}

							if(graph2->nodeInGraph(nNode-graph2->getOffset()))
							{
								cout<<j+graph->getOffset()<<" "<<nNode<<" "<<graph->edgeOvlLen(j, k)<<endl;
							}
						}
					}
				}				
			}
		}

		delete [] fileName;

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

//void popBubbles()
//Description: This function makes the graph
//Input: cMap (map<uint32_t, uint32_t>): the map from containments, graphArr (Graph **), pointers to 
//the graphs, numGraphs (int): the number of graphs that are loaded, graph (Graph &) the graph we are 
//creating.  
//Output:None
//Return:None
int popBubbles(Graph ** & graphs, int graphsLoaded, set<int> & pNodes, int * & nDen, int bubbleLen)
{
	int numPopped = 0;
	priority_queue< pair<int, int> >  myQueue;
	map<int, int> backTrack;

	for(set<int>::iterator it = pNodes.begin(); it != pNodes.end(); it++)
	{
		int node1 = (*it);
		Graph * graph1;
		for(int l = 0; l < graphsLoaded; l++)
		{
			graph1 = graphs[l];
			if(graph1->getOffset() <= node1 && node1 < graph1->getOffset() + graph1->getNumNodes())
				break;
		}

		if(graph1->nodeInGraph(node1-graph1->getOffset()) )
		{

			int tally = 0;	
			for(int j = 0; j < graph1->nodeDegree(node1-graph1->getOffset()); j++)
			{
				int nNode = graph1->edgeDest(node1-graph1->getOffset(), j);
				Graph * graph2;
				for(int l = 0; l < graphsLoaded; l++)
				{       
					graph2 = graphs[l];
					if(graph2->getOffset() <= nNode && nNode < graph2->getOffset() + graph2->getNumNodes())
						break;
				}

				if(graph1->edgeInGraph(node1-graph1->getOffset(), j) && graph2->nodeInGraph(nNode-graph2->getOffset()))
					tally++;	

			}

			if(tally > 1)
			{
				backTrack.insert(pair<int, int>(node1, -1));
				for(int j = 0; j < graph1->nodeDegree(node1-graph1->getOffset()); j++)
				{
					int nNode = graph1->edgeDest(node1-graph1->getOffset());
					Graph * graph2;
					for(int l = 0; l < graphsLoaded; l++)
					{
						graph2 = graphs[l];
						if(graph2->getOffset() <= nNode && nNode < graph2->getOffset() + graph2->getNumNodes())
							break;
					}

					if(graph1->edgeInGraph(node1-graph1->getOffset(), j) && graph2->nodeInGraph(nNode-graph2->getOffset()))
					{
						myQueue.push(pair<int, int> (nDen[nNode], nNode));
						backTrack.insert(pair<int, int>(nNode, node1));
					}
				}

				while(!myQueue.empty())
				{
					pair<int, int> next = myQueue.top(); myQueue.pop();
					int currNode = next.second;
					int prevNode = (*backTrack.find(currNode)).second;

					if(next.first - nDen[currNode] > (3*bubbleLen))
						break;

					Graph * graph1;
					for(int l = 0; l < graphsLoaded; l++)
					{
						graph1 = graphs[l];
						if(graph1->getOffset() <= currNode && currNode < graph1->getOffset() + graph1->getNumNodes())
							break;
					}

					for(int j = 0; j < graph1->nodeDegree(currNode-graph1->getOffset()); j++)
					{
						int nextNode = graph1->edgeDest(currNode-graph1->getOffset(), j);

						Graph * graph2;
						for(int l = 0; l < graphsLoaded; l++)
						{
							graph2 = graphs[l];
							if(graph2->getOffset() <= nextNode && nextNode < graph2->getOffset() + graph2->getNumNodes())
								break;
						}	

						if(graph1->edgeInGraph(currNode-graph1->getOffset(), j) &&  graph2->nodeInGraph(nextNode-graph2->getOffset()) && nextNode != prevNode)
						{
							if(backTrack.find(nextNode) == backTrack.end())
							{
								int pathWeight = next.first+nDen[nextNode];						

								backTrack.insert(pair<int, int>(nextNode, currNode));
								next.first = pathWeight; next.second = nextNode;		
								myQueue.push(next);

							}else{
								next = (*backTrack.find(nextNode));

								int cNode = next.first; int pNode = next.second;
								cNode = pNode; pNode = (*backTrack.find(cNode)).second;					

								vector<int> toErase; 

								while(pNode != -1)
								{
									Graph * graph3;
									for(int l = 0; l < graphsLoaded; l++)
									{
										graph3 = graphs[l];
										if(graph3->getOffset() <= cNode && cNode < graph3->getOffset() + graph3->getNumNodes())
											break;
									}

									int tally2 = 0;
									for(int l = 0; l < graph3->nodeDegree(cNode-graph3->getOffset()); l++)
									{
										int nNode = graph3->edgeDest(cNode-graph3->getOffset(), j);
										Graph * graph4;
										for(int l = 0; l < graphsLoaded; l++)
										{
											graph4 = graphs[l];
											if(graph4->getOffset() <= nNode && nNode < graph4->getOffset() + graph4->getNumNodes())
												break;
										}

										if(graph4->nodeInGraph(nNode-graph4->getOffset()) && graph3->edgeInGraph(cNode-graph3->getOffset(), j) && nNode != pNode)
											tally2++;
									}

									if(tally2 > 1)
										break;

									toErase.push_back(cNode);
									backTrack.erase(cNode);

									cNode = pNode;
									pNode = (*backTrack.find(cNode)).second;
								}


								for(vector<int>::iterator it = toErase.begin(); it != toErase.end(); it++)
								{
									int eNode = (*it);

									Graph * graph3;
									for(int l = 0; l < graphsLoaded; l++)
									{       
										graph3 = graphs[l];
										if(graph3->getOffset() <= eNode && eNode < graph3->getOffset() + graph3->getNumNodes())
											break;
									}									


									graph3->maskNodeFromGraph(eNode-graph3->getOffset());
									numPopped++;
								}

								if(cNode == next.second)
								{
									Graph * graph3;
									for(int l = 0; l < graphsLoaded; l++)
									{
										graph3 = graphs[l];
										if(graph3->getOffset() <= cNode && cNode < graph3->getOffset() + graph3->getNumNodes())
											break;
									}

									for(int l = 0; l < graph3->nodeDegree(cNode-graph3->getOffset()); l++)
									{
										int nNode = graph3->edgeDest(cNode-graph3->getOffset(), j);

										if(nNode == nextNode)
										{
											graph3->maskEdgeFromGraph(cNode-graph3->getOffset(), j);
											break;
										} 
									}

									for(int l = 0; l < graphsLoaded; l++)
									{
										graph3 = graphs[l];
										if(graph3->getOffset() <= nextNode && nextNode < graph3->getOffset() + graph3->getNumNodes())
											break;
									}

									for(int l = 0; l < graph3->nodeDegree(nextNode-graph3->getOffset()); l++)
									{
										int nNode = graph3->edgeDest(nextNode-graph3->getOffset(), j);
										if(nNode == cNode)
										{
											graph3->maskEdgeFromGraph(nextNode-graph3->getOffset(), j);
											break;
										}
									}
								}

								backTrack.erase(next.first);
								next.second = currNode; 
								backTrack.insert(next);
							}
						}
					}				
				}
			}

		}
		
		backTrack.clear(); 
	
		while(!myQueue.empty())
			myQueue.pop();
	}	

	return numPopped;
}

//void trimGraph()
//Description: This function trims the overlap graph
//Input: graph (graph &): the overlap graph, pathLen (int): The length of the ends
//that can be trimmed.
//Output:None
//Return:None
int trimGraph( Graph ** & graphs, int graphsLoaded, set<int> & pNodes, int * & nDen, int pathLen)
{
	int numTrimmed = 0;
	for(set<int>::iterator it = pNodes.begin(); it != pNodes.end(); it++)
	{
		int node1 = (*it);
		Graph * graph1;
		for(int l = 0; l < graphsLoaded; l++)
		{
			graph1 = graphs[l];
			if(graph1->getOffset() <= node1 && node1 < graph1->getOffset() + graph1->getNumNodes())
				break;
		}

		if(graph1->nodeInGraph(node1-graph1->getOffset()))
		{
			int eIndex; int tally = 0;
			for(int j = 0; j < graph1->nodeDegree(node1-graph1->getOffset()); j++)
			{
				int node2 = graph1->edgeDest(node1-graph1->getOffset(), j);
				Graph * graph2;
				for(int l = 0; l < graphsLoaded; l++)
				{
					graph2 = graphs[l];
					if(graph2->getOffset() <= node2 && node2 < graph2->getOffset() + graph2->getNumNodes())
						break;
				}

				if(graph1->edgeInGraph(node1-graph1->getOffset(), j) &&  graph2->nodeInGraph(node2-graph2->getOffset()))
				{
					eIndex = j; tally++;
				}
			}

			int nextNode = node1; int pathWeight = 0;
			vector<int> toErase;

			if(tally == 0)
			{
				toErase.push_back(nextNode);
				pathWeight+=nDen[nextNode];
			}

			while(tally == 1)
			{
				toErase.push_back(nextNode);
				pathWeight+=nDen[nextNode];

				Graph * graph2;
				for(int l = 0; l < graphsLoaded; l++)
				{
					graph2 = graphs[l];
					if(graph2->getOffset() <= nextNode && nextNode < graph2->getOffset() + graph2->getNumNodes())
						break;
				}

				int prevNode = nextNode;
				nextNode = graph2->edgeDest(nextNode-graph2->getOffset(), eIndex);

				tally = 0; eIndex = 0;

				for(int j = 0; j < graph2->nodeDegree(nextNode-graph2->getOffset()); j++)
				{
					int n = graph2->edgeDest(nextNode-graph2->getOffset(), j);
					if(n != prevNode)
					{
						Graph * graph3;
						for(int l = 0; l < graphsLoaded; l++)
						{
							graph3 = graphs[l];
							if(graph3->getOffset() <= n && n < graph3->getOffset() + graph3->getNumNodes())
								break;
						}

						if(graph2->edgeInGraph(nextNode-graph2->getOffset(), j) && graph3->nodeInGraph(n-graph3->getOffset()))
						{
							eIndex = j; tally++;
						}
					}
				}

				if(pathWeight > pathLen)
					break;
			}

			if(tally == 0 && nextNode != node1)
			{
				toErase.push_back(nextNode);
				pathWeight+=nDen[nextNode];
			}

			if(pathWeight <= pathLen)
			{
				for(int i = 0; i < toErase.size(); i++)
				{
					int currNode = toErase.at(i);

					Graph * graph1;
					for(int l = 0; l < graphsLoaded; l++)
					{
						graph1 = graphs[l];
						if(graph1->getOffset() <= currNode && currNode < graph1->getOffset() + graph1->getNumNodes())
							break;
					}

					graph1->maskNodeFromGraph(currNode-graph1->getOffset());

					numTrimmed++;
				}
		
			}		
		}
	}

	return numTrimmed;
}


//void transitiveReduce()
//Description: This function transitively reduces a graph
//Input: 
//Output:None
//Return:None
void transitiveReduce(Graph ** & graphs, int graphsLoaded, set<int> & pNodes, uint32_t * & sandBox, int numNodes)
{
	for(set<int>::iterator it = pNodes.begin(); it != pNodes.end(); it++)
	{
		uint32_t node1 = (*it);
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
			if(!graph1->rDovetail(node1-graph1->getOffset(), j))
				sandBox[dest] = 1;
		}

		for(int j = 0; j < graph1->nodeDegree(node1-graph1->getOffset()); j++)
		{
			uint32_t node2 = graph1->edgeDest(node1-graph1->getOffset(), j);
			if(graph1->nodeInGraph(node2) && !graph1->rDovetail(node1-graph1->getOffset(), j))
			{
				Graph * graph2;
				for(int l = 0; l < graphsLoaded; l++)
				{
					graph2 = graphs[l];
					if(graph2->getOffset() <= node2 && node2 < graph2->getOffset() + graph2->getNumNodes())
						break;
				}

				bool foundLongest = false; int pos1 = 0;
				for(int h = 0; h < graph2->nodeDegree(node2-graph2->getOffset()); h++)
				{
					uint32_t dest = graph2->edgeDest(node2-graph2->getOffset(), h);
					if(!graph2->rDovetail(node2-graph2->getOffset(), h))
						sandBox[dest] = 0;
				}

			}
		}

		for(int j = 0; j < graph1->nodeDegree(node1-graph1->getOffset()); j++)
		{
			uint32_t node2 = graph1->edgeDest(node1-graph1->getOffset(), j);
			if(!sandBox[node2] && !graph1->rDovetail(node1-graph1->getOffset(), j))
			{
				graph1->maskEdgeFromGraph(node1-graph1->getOffset(), j);	
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
					if(dest == node1)
					{
						graph2->maskEdgeFromGraph(node2-graph2->getOffset(), h);
						break;
					}
				}	
			}
		}

	}
}
