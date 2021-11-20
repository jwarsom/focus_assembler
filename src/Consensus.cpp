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
#include "MappingValues.h"
#include "Dictionary.h"
#include "Graph.h"
#include "Fragment_Index.h"
#include "BandedAlignment.h"
#include "Needleman_Wunsch.h"
#include "Labels.h"
#include "Consensus.h"


//Consensus()
//Description:Default Constructor
//Input: None
//Ouput:None
//Return:None
Consensus::Consensus(void)
{

}

//~Consensus()
//Description:Default Destructor
//Input: None
//Ouput:None
//Return:None
Consensus:: ~Consensus(void)
{

}


//addSequences(int [], int, int * [], int &)
//Description:This adds a sequence to the consensus
//Input: cSeq (an int array): Represents the current consensus sequences, cSeqCutOff(int):
//How much of the sequence we should add
//Ouput:None
//Return:None
void Consensus::addSequence(int cSeq [], int cSeqCutOff, char consensus [], int & currConsPos, char flags [])
{
	for(int i = 0; i < 6*cSeqCutOff; i+=6)
	{
		float max = 0; int maxPos = 0; float total = 0; int nextMax = 0;
		for(int j = 0; j < 6; j++)
		{
			if(cSeq[i+j] >= max)
			{
				nextMax = max; max = cSeq[i+j]; maxPos = j; 
			}else{
				if(cSeq[i+j] > nextMax)
					nextMax = cSeq[i+j];
			}				

			total+=cSeq[i+j];
		}

		if(maxPos > 1)
		{
			if(nextMax/total >= .20 && total > 5)
			{
				flags[currConsPos] = '@';
			}else{
				flags[currConsPos] = '$';
			}
			consensus[currConsPos++] = inverseMapVals[maxPos-2];
		}else{
			if(currConsPos > 0)
			{
				if(nextMax/total >= .20 && total > 5)
					flags[currConsPos-1] = '@';
			}
		}
	}

}


//getGraph(Graph * [], int, Graph * &, int) 
//Description:This retrieves the current graph that the node belongs to 
//Input: graphs (Graph * []): The array of graphs, gSize (int): The size of the graphs, graph (Graph * &), currNode (int): The current node we are searching for
//Ouput:None
//Return:None
void Consensus::getGraph(Graph * graphs [], int gSize, Graph * & graph, int currNode)
{
	for(int i = 0; i < gSize; i++)
	{
		graph = graphs[i];

		if(graph->getOffset() <= currNode && currNode < graph->getOffset() + graph->getNumNodes())
			break;
	}
}

//getFragment_Index(Fragment_Index * [], int, Fragment_Index * &, int) 
//Description:This retrieves the current fragment index that the node belongs to 
//Input: indexes (Fragment_Index * []): The array of fragment indexes, iSize (int): The size of the indexes, index1 (Fragment_Index * &), currNode (int): The current node we are searching for
//Ouput:None
//Return:None
void Consensus::getFragment_Index(Fragment_Index * indexes [], int iSize, Fragment_Index * & index, int currNode) {

	for(int i = 0; i < iSize; i++)
	{
		index = indexes[i];

		if(index->getIndexOffset() <= currNode && currNode < index->getIndexOffset() + index->numFragments())
			break;
	}
}

//addSequences(int [], int, int * [], int &)
//Description:This adds a sequence to the consensus
//Input: cSeq (an int array): Represents the current consensus sequences, cSeqCutOff(int):
//How much of the sequence we should add
//Ouput:None
//Return:None
bool Consensus::getConsensus(uint32_t contig [], int & contigPos, int minNode, int maxNode,  int * matrix [], int mSize, Graph * graphs [], int gSize, Fragment_Index * indexes [], int iSize, char  consensus [], int & consensusPos, int consensusSize, char flags [], long long int rBounds [], int & prevVal, long long int & bPos, long long int bSize, int tmpHold1 [], int tmpHold2 [], int tSize)
{
	string alignStr; int contigSize = contig[0];
	for(int i = contigPos; i < contigSize; i++)
	{
		int currNode = contig[i+1];
		Graph * graph; Fragment_Index * index1; Fragment_Index * index2;

		getGraph(graphs, gSize, graph, currNode);
		getFragment_Index(indexes, iSize, index1, currNode);

		int nodeIndex = currNode - graph->getOffset();
		pair<unsigned long long int, unsigned long long int> bounds = index1->indexBounds(currNode - index1->getIndexOffset());

		if(i == contigPos)
		{
			for(int j = bounds.first; j < bounds.second; j++)
			{
				alignStr.push_back(index1->at(j));
			} 

			for(int j = 0; j < tSize; j++)
			{
				tmpHold1[j] = 0; tmpHold2[j] = 0;
			}

			for(int j = 0; j < alignStr.length(); j++)
			{
				tmpHold1[6*j + mapVals2[alignStr.at(j)]] = 1;
			}
		}

		for(int j = 0; j < graph->nodeDegree(nodeIndex); j++)
		{
			if(graph->dContained(nodeIndex, j) && graph->edgeDest(nodeIndex, j) >= minNode && graph->edgeDest(nodeIndex, j) <= maxNode)
			{

				int cNode = graph->edgeDest(nodeIndex, j);
				getFragment_Index(indexes, iSize, index2, cNode);

				pair<unsigned long long int, unsigned long long int> bounds2 = index2->indexBounds(cNode - index2->getIndexOffset());
				
				int strLen = alignStr.length();
				Needleman_Wunsch<string, Fragment_Index>  needle(alignStr.length(), bounds2.second-bounds2.first, 0, bounds2.first);
			
				if(alignStr.length() < mSize && bounds2.second-bounds2.first  < mSize)
				{
					alignStr = needle.Align(matrix, alignStr, *index2, 0, 0);
				}else{
					int ** largeMatrix = new int * [alignStr.length()];
					for(int k = 0; k < alignStr.length(); k++)
					{
						largeMatrix[k] = new int [bounds2.second-bounds2.first];
					}	

					alignStr = needle.Align(largeMatrix, alignStr, *index2, 0, 0);

					for(int k = 0; k < alignStr.length(); k++)
					{
						delete [] largeMatrix[k];
					}	

					delete [] largeMatrix;
				}

				int addedGap = 0; int start = needle.yAlignmentStart(); int end = (alignStr.length()/2 - (bounds2.second-bounds2.first - needle.yAlignmentEnd()-1)) /*- (strLen - needle.xAlignmentEnd()-1)*/; //FIGURE THIS OUT 

				for(int k = start; k < end; k++)
				{
					if(alignStr.at(k) == '-' && tmpHold1[6*(k-start) + 1  - 6 * addedGap] == 0)		
					{
						int numCover = 0;
						for(int m = 0; m < 6; m++)
						{
							tmpHold2[6*(k-start)+m] = 0;
							numCover+=tmpHold1[6*(k-start) + m - 6 * addedGap];
						}
						tmpHold2[6*(k-start)+1] = 1;

						tmpHold2[6*(k-start)] = numCover;
						tmpHold2[6*(k-start) + mapVals2[alignStr.at(alignStr.length()/2 + k)]] += 1;
						addedGap++; 

					}else{		

						if((k-start) >= needle.xAlignmentStart())
						{
							tmpHold1[6*(k-start) +  mapVals2[alignStr.at(alignStr.length()/2+k)] - 6 * addedGap] += 1; 	
						}


						for(int m = 0; m < 6; m++)
						{
							tmpHold2[6*(k-start) + m] = tmpHold1[6*(k-start) + m - 6 * addedGap]; 
						}
					}
				} 	

				alignStr = alignStr.substr(start, end-start);

				int * tmp = tmpHold1;
				tmpHold1 = tmpHold2; 

				tmpHold2 = tmp;
			}		

		}

		if(i < contigSize - 1)
		{
			int cNode = contig[i+2];
			getFragment_Index(indexes, iSize, index2, cNode);

			pair<unsigned long long int, unsigned long long int> bounds2 = index2->indexBounds(cNode - index2->getIndexOffset());

			int strLen = alignStr.length();
			Needleman_Wunsch<string, Fragment_Index>  needle(alignStr.length(), bounds2.second-bounds2.first, 0, bounds2.first);

			if(alignStr.length() < mSize && bounds2.second-bounds2.first < mSize)
			{
				alignStr = needle.Align(matrix, alignStr, *index2, 0, 0);
			}else{
				int ** largeMatrix = new int * [alignStr.length()];
				for(int k = 0; k < alignStr.length(); k++)
				{
					largeMatrix[k] = new int [bounds2.second-bounds2.first];
				}	

				alignStr = needle.Align(largeMatrix, alignStr, *index2, 0, 0);

				for(int k = 0; k < alignStr.length(); k++)
				{
					delete [] largeMatrix[k];
				}	

				delete [] largeMatrix;
			}

			int addedGap = 0; int start = needle.yAlignmentStart(); int end = (alignStr.length()/2 - (bounds2.second-bounds2.first - needle.yAlignmentEnd()-1)) - (strLen - needle.xAlignmentEnd() -1) ;

			for(int k = start; k < end; k++)
			{
				if(alignStr.at(k) == '-' && tmpHold1[6*(k-start) + 1  - 6 * addedGap] == 0)		
				{
					int numCover = 0;
					for(int m = 0; m < 6; m++)
					{
						tmpHold2[6*(k-start)+ m] = 0;
						numCover+=tmpHold1[6*(k-start) + m - 6 * addedGap];
					}
					tmpHold2[6*(k-start)+1] = 1;


					tmpHold2[6*(k-start)] = numCover;
					tmpHold2[6*(k-start) + mapVals2[alignStr.at(alignStr.length()/2 + k)]] += 1;
					addedGap++; 

				}else{		

					if((k-start) >= needle.xAlignmentStart())
					{
						tmpHold1[6*(k-start) +  mapVals2[alignStr.at(alignStr.length()/2+k)] - 6 * addedGap] += 1; 	
					}

					for(int m = 0; m < 6; m++)
					{
						tmpHold2[6*(k-start) + m] = tmpHold1[6*(k-start) + m - 6 * addedGap]; 
					}
				}
			} 	

			int * tmp = tmpHold1;
			tmpHold1 = tmpHold2; 
			tmpHold2 = tmp;

			if(consensusPos + needle.xAlignmentEnd() < consensusSize &&  bPos + 2 < bSize) 	
			{
				addSequence(tmpHold1, needle.xAlignmentStart(), consensus, consensusPos, flags);

				rBounds[bPos++] = cNode; rBounds[bPos++] =  prevVal;
				prevVal = prevVal + needle.xAlignmentEnd();

			}else{
				contigPos = i;
				return false;
			}

			int addInt = 0;
			for(int k = needle.xAlignmentStart() + needle.yAlignmentStart(); k < end; k++)
			{
				tmpHold2[6*(addInt)] = tmpHold1[6*k];
				tmpHold2[6*(addInt)+1] = 0;

				for(int m = 2; m < 6; m++)
				{
					tmpHold2[6*(addInt) + m] = tmpHold1[6*k+m]; 
				}

				if(alignStr.at(alignStr.length()/2 + k) == '-')
					tmpHold2[6*(addInt)+1] = 1;

				addInt++;
			}

			for(int k = needle.yAlignmentEnd() + 1; k < bounds2.second - bounds2.first; k++)
			{
				for(int m = 0; m < 6; m++)
					tmpHold2[6*(addInt) + m] = 0; 

				tmpHold2[6*(addInt) + mapVals2[index2->at(bounds2.first+k)]] =  1; 
	
				addInt++;
			}


			tmp = tmpHold1;
			tmpHold1 = tmpHold2; 
			tmpHold2 = tmp;
			
			alignStr = alignStr.substr(alignStr.length()/2 + needle.xAlignmentStart() + needle.yAlignmentStart());

		}else{

			if(consensusPos + alignStr.size() < consensusSize &&  bPos + 2 < bSize) 	
			{
				addSequence(tmpHold1, alignStr.size(), consensus, consensusPos, flags);

				rBounds[bPos++] = currNode; rBounds[bPos++] =  prevVal;
				prevVal = prevVal + alignStr.size();

			}else{
				contigPos = i;
				return false;
			}

		}
	}		
	return true;
}


