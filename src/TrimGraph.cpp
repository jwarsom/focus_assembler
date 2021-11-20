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

	cout<<"\n\n"<<endl;

	cout<<"Exiting program"<<endl;
	exit(help);
}

int countFiles(const char []);
void transitiveReduce(Graph ** &, int, uint32_t * &, int);

int main(int argc, char * argv [])
{	
	//User Input Flags
	string workDirS = "--workDir";               //Current Working Directory
	string minNodeWeightS = "--minNodeWeight";         //Minimum path length
	string minEdgeWeightS = "--minEdgeWeight";
	string verboseS = "--verbose";		     //verbose mode
	string timeS = "--time"; //Report the time

	string workDir = "";
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

		//Display time
		if(argv[i] == timeS && reportTime != true)
		{
			reportTime = true;
			i--;
		}
	}

	if(reportTime == true)
	{
		time_t rawtime;
		struct tm * timeinfo;
		time(&rawtime);
		timeinfo = localtime(&rawtime);
		cerr<<"start time: "<<asctime(timeinfo)<<endl;
	}

	char * contigInfoDir = new char [1000];
	sprintf(contigInfoDir, "%s/ContigsInfo", workDir.c_str());

	int numContigs = countFiles(contigInfoDir) - 2;

	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/HybSpectrum/Graph0", workDir.c_str());

	int numGraphs = countFiles(graphDir)-2;

	for(int j = 0; j < numGraphs; j++)
	{
		Graph graph; 
		char * fileName = new char [1000];
		sprintf(fileName, "%s/HybSpectrum/Graph0/graph%d", workDir.c_str(), j);
		FILE * pFile = fopen(fileName, "r");
		graph.read(pFile);
		fclose(pFile);

		for(int i = 0; i < graph.getNumEdges(); i++)
		{
			graph.setEdgeInGraph(i);
		}

		for(int i = 0; i < graph.getNumNodes(); i++)
		{
			graph.setNodeInGraph(i);
		}

		pFile = fopen(fileName, "w");
		graph.write(pFile);
		fclose(pFile);

		delete [] fileName;
	}

	char * fileName = new char [1000];

	sprintf(fileName, "%s/HybSpectrum/Graph0/nDen", workDir.c_str());

	FILE * pFile = fopen(fileName, "r");
	fseek(pFile, 0, SEEK_END);
	int nDenSize = ftell(pFile)/sizeof(int);
	fseek(pFile, 0, SEEK_SET);

	int * nDen = new int [nDenSize];

	fread(nDen, sizeof(int), nDenSize, pFile);
	fclose(pFile);

	cerr<<"NUMBER OF NODES "<<nDenSize<<endl;

	sprintf(fileName, "%s/ContigsInfo/contigBounds.store", workDir.c_str());

	pFile = fopen(fileName, "r");
	fseek(pFile, 0, SEEK_END);
	int cBoundsSize = ftell(pFile)/sizeof(long long int);
	fseek(pFile, 0, SEEK_SET);	

	long long int * cBounds = new long long int [cBoundsSize];

	fread(cBounds, sizeof(long long int), cBoundsSize, pFile);
	fclose(pFile);

	sprintf(fileName, "%s/ContigsInfo/readBounds.store", workDir.c_str());

	pFile = fopen(fileName, "r");
	fseek(pFile, 0, SEEK_END);
	int rBoundsSize = ftell(pFile)/sizeof(long long int);
	fseek(pFile, 0, SEEK_SET);	

	long long int * rBounds = new long long int [rBoundsSize];

	fread(rBounds, sizeof(long long int), rBoundsSize, pFile);
	fclose(pFile);

	int ** matrix = new int * [101];

	for(int i = 0; i < 101; i++)
		matrix[i] = new int[500];

	cerr<<"SIZE OF "<<nDenSize<<endl;
	int * eNodes = new int [nDenSize];
	for(int i = 0; i < nDenSize; i++)
	{
		eNodes[i] = 0;
	}


	int prevContig = 0; int prevGraph = 0; int nodesErased = 0;
	for(int i = 0; i < numContigs; i++)
	{
		sprintf(fileName, "%s/ContigsInfo/Contigs.store.%d", workDir.c_str(), i);

		pFile = fopen(fileName, "r");
		Fragment_Index index1;

		unsigned long long int fragmentOffset;
		fread(&fragmentOffset, sizeof(unsigned long long int), 1, pFile);

		unsigned long long int indexOffset;
		fread(&indexOffset, sizeof(unsigned long long int), 1, pFile);

		unsigned long long int length;
		fread(&length, sizeof(unsigned long long int), 1, pFile);

		unsigned long long int numFragments;
		fread(&numFragments, sizeof(unsigned long long int), 1, pFile);

		fseek(pFile, 0, SEEK_SET);

		index1.reserve(numFragments, length, indexOffset, fragmentOffset);

		index1.read(pFile);
		fclose(pFile);

		index1.set();

		int currContig = prevContig; int currGraph = prevGraph;
		for(int j = 0; j < numContigs /*i+1*/; j++)
		{
			sprintf(fileName, "%s/ContigsInfo/Contigs.store.%d", workDir.c_str(), i);

			pFile = fopen(fileName, "r");
			Fragment_Index index2;

			fread(&fragmentOffset, sizeof(unsigned long long int), 1, pFile);
			fread(&indexOffset, sizeof(unsigned long long int), 1, pFile);
			fread(&length, sizeof(unsigned long long int), 1, pFile);
			fread(&numFragments, sizeof(unsigned long long int), 1, pFile);

			fseek(pFile, 0, SEEK_SET);

			index2.reserve(numFragments, length, indexOffset, fragmentOffset);
			index2.read(pFile);
			fclose(pFile);

			index2.set();

			currContig = prevContig; currGraph = prevGraph;
			
			Graph graph; 
			sprintf(fileName, "%s/HybSpectrum/Graph0/graph%d", workDir.c_str(), currGraph);

			pFile = fopen(fileName, "r");
			graph.read(pFile);
			fclose(pFile);

			while(currContig < cBoundsSize && cBounds[currContig] >= index1.getFragmentOffset() && cBounds[currContig] < index1.getFragmentOffset() + index1.length())
			{
				if(i == j && nDen[currContig] != 0)
				{
					char name[200];
					sprintf(name, ">contig%d %d", currContig, nDen[currContig]);
					int end = 0;
					if(currContig == cBoundsSize-1){ end = index1.length();}
					else{
						end = cBounds[currContig+1]-index1.getFragmentOffset();
						if(end > index1.length()) end = index1.length();
					}
					cout<<name<<endl;
					for(int m = cBounds[currContig]-index1.getFragmentOffset(); m < end; m++)
						cout<<index1.at(m);
					cout<<endl;
				}

				if(currContig >= graph.getNumNodes() + graph.getOffset())
				{
					sprintf(fileName, "%s/HybSpectrum/Graph0/graph%d", workDir.c_str(), currGraph);
					pFile = fopen(fileName, "w");
					graph.write(pFile);
					fclose(pFile);

					currGraph++;			
					sprintf(fileName, "%s/HybSpectrum/Graph0/graph%d", workDir.c_str(), currGraph);
					pFile = fopen(fileName, "r");
					graph.read(pFile);
					fclose(pFile);
				}
					
				if(nDen[currContig] != 0 && nDen[currContig] <= minNodeWeight && graph.nodeInGraph(currContig-graph.getOffset()))
				{
					for(int k = 0; k < graph.nodeDegree(currContig-graph.getOffset()) && !eNodes[currContig]; k++)
					{
						int destContig = graph.edgeDest(currContig-graph.getOffset(), k);
						if(nDen[destContig] >= nDen[currContig] &&  cBounds[destContig] >= index2.getFragmentOffset() && cBounds[destContig] < index2.getFragmentOffset() + index2.length())	
						{
							int beginX = cBounds[currContig] - index1.getFragmentOffset(); int endX;
							if(currContig == cBoundsSize-1)
							{
								endX = index1.length();
							}else{
								endX = cBounds[currContig+1] - index1.getFragmentOffset();	

								if(endX > index1.length()) endX = index1.length();
							} 

							int beginY = cBounds[destContig] - index2.getFragmentOffset(); int endY;				
							if(destContig == cBoundsSize-1)
							{
								endY = index2.length();
							}else{

								endY = cBounds[destContig+1] - index2.getFragmentOffset();	

								if(endY > index2.length()) endY = index2.length();
							} 

							for(int m = beginY; m < endY; m+=50)
							{
								int x_len = 50; int y_len = min(100, endY-m);
								Needleman_Wunsch<Fragment_Index, string> needle(x_len, y_len, beginX, m);											
								needle.Align(matrix, index1, index2, 0, 0);

								float misAligned = needle.misAligned() + needle.xAlignmentStart() + (x_len - needle.xAlignmentEnd());
								float alignLen = needle.alignLen() + needle.xAlignmentStart() + (x_len - needle.xAlignmentEnd());

								bool notContained = false;
								if(misAligned/alignLen <= .10)
								{
									int alignLoc = m + needle.yAlignmentEnd()+1;	
									for(int n = beginX+50; n < endX && (endX-n) > 10; n+=50)
									{
										int xAlignmentStart = n;
										int yAlignmentStart = max(beginY, alignLoc - 10);

										int xAlignmentEnd = min(n+50, endX);
										int yAlignEnd = min(alignLoc+(xAlignmentEnd - xAlignmentStart)+25, endY);	

										x_len = xAlignmentEnd - xAlignmentStart;
										y_len = yAlignEnd - yAlignmentStart;

										Needleman_Wunsch<Fragment_Index, string> needle2(x_len, y_len, xAlignmentStart, yAlignmentStart);
										needle2.Align(matrix, index1, index2, 0, 0);

										float a1 = needle2.misAligned() + needle2.xAlignmentStart() + (x_len - needle2.xAlignmentEnd()-1);
										float a2 = needle2.alignLen() + needle2.xAlignmentStart() + (x_len - needle2.xAlignmentEnd()-1);

										a1+=abs(((yAlignmentStart+needle2.yAlignmentStart()) - alignLoc));
										a2+=abs(((yAlignmentStart+needle2.yAlignmentStart()) - alignLoc));

										if(a1/a2 > .30)
										{
											misAligned+=(a1 + (endX-xAlignmentEnd));
											alignLen+=(a2+ (endX-xAlignmentEnd));
											break;
										}

										misAligned+=a1;
										alignLen+=a2;	

										alignLoc = yAlignmentStart + needle2.yAlignmentEnd()+1;
										if(alignLoc >= yAlignEnd)
										{
											notContained = true; break;
										}
									}							

									if(misAligned/alignLen <= .10 && !notContained)
									{
										eNodes[currContig] = 1;
										graph.maskNodeFromGraph(currContig-graph.getOffset());
										nodesErased++;	
										break;
									}
								} 
							}				
						}
					} 
				}	

				currContig++;				
			}

			sprintf(fileName, "%s/HybSpectrum/Graph0/graph%d", workDir.c_str(), currGraph);
			pFile = fopen(fileName, "w");
			graph.write(pFile);
			fclose(pFile);
		}	
		prevContig = currContig; prevGraph = currGraph; 
	}	


	cerr<<nodesErased<<" ERASED NODES "<<endl;
	delete [] cBounds;
	delete [] rBounds;

	for(int i = 0; i < 101; i++)
		delete [] matrix [i];


	delete [] matrix;

	for(int i = 0; i < numGraphs; i++)
	{
		sprintf(fileName, "%s/HybSpectrum/Graph0/graph%d", workDir.c_str(), i);
		pFile = fopen(fileName, "r");
		cerr<<fileName<<endl;

		Graph graph;
		graph.read(pFile);
		fclose(pFile);

		for(int j = 0; j < graph.getNumEdges(); j++)
		{
			cerr<<graph.edgeSource(j)<<" "<<graph.edgeDest(j)<<endl;

			if(eNodes[graph.edgeSource(j)] == 1 || eNodes[graph.edgeDest(j)] == 1)
			{
				graph.maskEdgeFromGraph(j);
			}
		}


		pFile = fopen(fileName, "w");
		graph.write(pFile);
		fclose(pFile);
	}


	int total = 0;
	for(int i = 0; i < nDenSize; i++)
		if(eNodes[i] == 1)
			total++;

	cerr<<"ERASED "<<total<<endl;	

	delete [] eNodes;

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

	sprintf(fileName, "%s/Contigs/cMap0", workDir.c_str());

	pFile = fopen(fileName, "r");
	int * cMap = new int [numNodes];
	fread(cMap, sizeof(int), numNodes, pFile);
	fclose(pFile);

	int tMasked = 0; 

	cerr<<"Transitive "<<endl;
	uint32_t * sandBox = new uint32_t [numNodes];
	transitiveReduce(graphs, graphsLoaded, sandBox, numNodes);
	delete [] sandBox;

	cerr<<"Transitive DONE "<<endl;

	for(int i = 0; i < graphsLoaded; i++)
	{
		sprintf(fileName, "%s/HybSpectrum/Graph0/graph%d", workDir.c_str(), i);

		pFile = fopen(fileName, "w");

		graphs[i]->write(pFile);

		fclose(pFile);
	}

	ofstream of1;  ofstream of2; ofstream of3;

	sprintf(fileName, "%s_graph.txt", workDir.c_str());

	of1.open(fileName);
	sprintf(fileName, "%s_eWeights.txt", workDir.c_str());
	of2.open(fileName);
	sprintf(fileName, "%s_nWeights.txt", workDir.c_str());
	of3.open(fileName);

	for(int i = 0; i < graphsLoaded; i++)
	{
		for(int j = 0; j < graphs[i]->getNumEdges(); j++)
		{
			if(graphs[i]->edgeInGraph(j))
			{
				of1<<graphs[i]->edgeSource(j)+graphs[i]->getOffset()<<" "<<graphs[i]->edgeDest(j)<<endl;
			}
		}
	}

	int totalEdges = 0;
	for(int i = 0; i < graphsLoaded; i++)
	{
		for(int j = 0; j < graphs[i]->getNumEdges(); j++)
		{
			if(graphs[i]->edgeInGraph(j))
			{
				of2<<graphs[i]->edgeOvlLen(j)<<endl;
			}

			totalEdges++;
		}
	}

	cerr<<"Total Edges "<<totalEdges<<endl;

	for(int i = 0; i < numNodes; i++)
	{
		of3<<nDen[i]<<endl;
	}

	delete [] fileName; 
	delete [] graphDir; 
	delete [] contigInfoDir;	
	delete [] cMap;
	delete [] nDen;

	if(reportTime == true)
	{
		time_t rawtime;
		struct tm * timeinfo;
		time(&rawtime);
		timeinfo = localtime(&rawtime);
		cerr<<"end time: "<<asctime(timeinfo)<<endl;
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

	return numFiles;
}

//void transitiveReduce()
//Description: This function transitively reduces a graph
//Input: 
//Output:None
//Return:None
void transitiveReduce(Graph ** & graphs, int graphsLoaded, uint32_t * & sandBox, int numNodes)
{
	for(int i = 0; i < numNodes; i++)
	{
		uint32_t node1 = i;
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
