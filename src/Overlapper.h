/*
 * Overlapper.h
 *
 *      Author: Julia
 */

#ifndef OVERLAPPER_H_
#define OVERLAPPER_H_
#include <queue>
#include<vector>
#include <math.h>
#include <map>
#include<set>
#include<algorithm>
using namespace std;
#include "Chain.h"
#include "SuffixArray.h"
#include "Fragment_Index.h"
#include "FM_Index.h"
#include "Seed.h"

/* This is a templated class that actually controls that overlapping process*/
/* This class uses a compressed suffix array that has been passed to it along with
 * a FM Index to find overlaps in the fragments. It writes the overlaps to an array
 * which is then written to a file
 */

template <class T>
class Overlapper {
	public:
		Overlapper(SuffixArray<T> &, FM_Index<T> &, T &, T &, int * &, int * &, const int, const int, const int, const int, const bool, const bool);
		virtual ~Overlapper();
		void findOverlaps(const char [], bool);
	private:
		void getAlignments(unsigned long long int, int * [], int [], unsigned long long int [], int [], uint32_t [], FILE * &, int &, int);
		bool passedQgram(Chain &, pair<unsigned long long int, unsigned long long int> &, pair<unsigned long long int, unsigned long long int> &, int [], const unsigned int);
		pair<long int, long int> findBounds(const char, pair<long int, long int> &);
		void collectHits(pair<unsigned long long int, unsigned long long int> &, set<Seed> &, int [], unsigned long long int [], const unsigned int);
		int encode(T * &, const unsigned long long int, const int, const unsigned int);
		void alignSeqs(Chain &, pair<unsigned long long int, unsigned long long int> &, pair<unsigned long long int, unsigned long long int> &, pair<float, float> &, pair<int, int> &, pair<int, int> &, int * [], const unsigned int);

		//The data structure pointers
		SuffixArray<T> * suffix;
		FM_Index<T> * fmIndex;
		T * qIndex;
		T * rIndex;
		int * cOvlq;
		int * cOvlr;

		//The params
		int kmerSize;
		int stepSize;
		int identity;
		int length;
		bool containment;
		bool alignRepsOnly;
};

//Overlapper(FM_Index &, hierarchical_CompressedArray &, T &, T &, const int, const int, const int, const int)
//Description: This is the constructor
//Input: FM_Index & f, the fm index, hierarchical_CompressedArray & c, the compressed array, T & qIndex, T & rIndex, the string indexes, const int k, kmerSize, const int s, stepSize,
//const int i, the minimum identity, const int l, the minimum length of the overlap.
//Output:None
//Return:None
	template <class T>
Overlapper<T>::Overlapper(SuffixArray<T> & c, FM_Index<T> & f, T & q, T & r, int * & coq, int * & cor, const int k, const int s, const int i, const int l,  const bool cntnmnt, const bool repsOnly)
	:suffix(&c), fmIndex(&f), qIndex(&q), rIndex(&r), cOvlq(coq), cOvlr(cor), kmerSize(k), stepSize(s), identity(i), length(l) , containment(cntnmnt), alignRepsOnly(repsOnly){;;;}

	//~Overlapper()
	//Description: This is the destructor
	//Input:None
	//Output:None
	//Return:None
	template <class T>
	Overlapper<T>::~Overlapper() {}

	//void findOverlaps(const char [], const char [], const char [])
	//Description: This function finds the overlaps between the forward and reverse sequences
	//Input: char array, overlapFile; the file that the overlap information is stored into, const char readLenFile [], the file that the read length information goes into
	//const char readQualFile []; the file that will contain the quality values
	//Output: None
	//Return:None
	template <class T>
void Overlapper<T>::findOverlaps(const char o [],  bool verbose)
{
	//Start up the qGram
	int * qGram = new int[2048];
	int ** matrix = new int * [1001];
	for(int i = 0; i < 1001; i++)
		matrix[i] = new int[1001];

	uint32_t * overlaps = new uint32_t[25 * megabyte];
	int oTrack = 0;

	int * hitStats = new int [1000];
	unsigned long long int * hitIndexes = new unsigned long long int [1000];

	int currFile = 0;
	char * overlapFile = new char[500];
	sprintf(overlapFile, "%s_%d", o, currFile);

	//TRACK HOW MANY MB
	FILE * oFile1 = fopen(overlapFile, "a");
	for(unsigned int i = 0; i < qIndex->numFragments(); i++)
	{
		if(verbose && i%1000 == 0)
			cout<<"Aligned "<<i<<" of "<<qIndex->numFragments()<<" sequences"<<endl;
	
		if((i+qIndex->getIndexOffset() == cOvlq[3*((i+qIndex->getIndexOffset()) - cOvlq[0]) + 1]) || (!containment && !alignRepsOnly && cOvlq[3*((i+qIndex->getIndexOffset()) - cOvlq[0]) + 3] >= 90))
			getAlignments(i, matrix, hitStats, hitIndexes, qGram, overlaps, oFile1, oTrack, 0);

		//getAlignments(i, matrix, qGram, overlaps, oFile1, oTrack, 1); //No more reverse
	}

	if(verbose)
		cout<<"Aligned "<<qIndex->numFragments()<<" of "<<qIndex->numFragments()<<" sequences"<<endl;

	if(oTrack > 0)
	{
		fwrite(overlaps, sizeof(uint32_t), oTrack, oFile1);
		fclose(oFile1);
	}

	delete [] qGram;
	for(int i = 0; i < 1001; i++)
	{
		delete [] matrix[i];
	}

	delete [] overlapFile;
	delete [] hitStats;
	delete [] hitIndexes;
	delete [] matrix;
	delete [] overlaps;
}

//void getAlignments(int, int [], int [], int [], int [], int, int, int, int)
//Description: This function returns the alignments of various fragments
//Input: qFragmentIndex: the index of the query sequence, qGrams: The array for storing the qGrams, overlaps: the array that will store the overlap info
//readLengths: the array that stores the readLength values. qualVals: the array that stores the quality, the next three integers are tracking values
//int reverse: Determines if the query is the reverse sequence or not
//Output:None
//Return: None
	template <class T>
void Overlapper<T>::getAlignments(unsigned long long  qFragmentIndex, int * matrix[], int hitStats [], unsigned long long int hitIndexes [], int qGram [], uint32_t overlaps [], FILE * & oFile1, int & oTrack, int reverse)
{
	set<Seed> hits;
	pair<unsigned long long int, unsigned long long int> qBounds = qIndex->indexBounds(qFragmentIndex, reverse);
	fill(qGram, qGram + 2048, -1);

	//Collect the hits that the query has in the reference
	collectHits(qBounds, hits, hitStats, hitIndexes, reverse);

	set<Seed>::iterator it = hits.begin();

	//Iterate through the hits
	while(it != hits.end())
	{
		unsigned long long int rFragmentIndex = rIndex->positionToIndex((*it).coord1Y);
		pair<unsigned long long int, unsigned long long int> rBounds = rIndex->indexBounds(rFragmentIndex);
		Chain chain;

		while(it != hits.end() && (*it).coord1Y < rBounds.second)
		{
			if((*it).coord2Y < rBounds.second)
			{
				chain.addSeed((*it));
			}
			it++;
		}

		if(cOvlq[3*((qFragmentIndex+qIndex->getIndexOffset()) - cOvlq[0]) + 1] != cOvlr[3*((rFragmentIndex+rIndex->getIndexOffset()) - cOvlr[0]) + 1] && (cOvlr[3*((rFragmentIndex+rIndex->getIndexOffset()) - cOvlr[0]) + 1] == rFragmentIndex+rIndex->getIndexOffset() || (cOvlr[3*((rFragmentIndex+rIndex->getIndexOffset()) - cOvlr[0]) + 3] >= 90 && !alignRepsOnly) || containment))
		{
			chain.chainSeeds();
			if(passedQgram(chain, qBounds, rBounds, qGram, reverse))
			{
				pair<float, float> alignScore; pair<int, int> qAlignment; pair<int, int> rAlignment; 
				alignSeqs(chain, qBounds, rBounds, alignScore, qAlignment, rAlignment, matrix, reverse);

				if(((alignScore.first - alignScore.second)/alignScore.first)*100 >= identity && alignScore.first >= length)
				{
					if(containment)
					{
						float rLen = rBounds.second-rBounds.first; float qLen = qBounds.second-qBounds.first;
						float cInterval = .03;

						uint32_t rContained = (alignScore.first/rLen + cInterval);			

						if((qFragmentIndex+qIndex->getIndexOffset() < rFragmentIndex+rIndex->getIndexOffset()) && rContained)
						{
							int cRatio = 100 * (rLen/qLen);
							if(((alignScore.first - alignScore.second)/alignScore.first)*100 >= cOvlr[3*((rFragmentIndex+rIndex->getIndexOffset()) - cOvlr[0]) + 2])
							{
								cOvlr[3*((rFragmentIndex+rIndex->getIndexOffset()) - cOvlr[0]) + 1] = qFragmentIndex+qIndex->getIndexOffset();
								cOvlr[3*((rFragmentIndex+rIndex->getIndexOffset()) - cOvlr[0]) + 2] = ((alignScore.first - alignScore.second)/alignScore.first)*100; 
								cOvlr[3*((rFragmentIndex+rIndex->getIndexOffset()) - cOvlr[0]) + 3] = cRatio; 
							}
						}
					}else{
						if((oTrack+6) > 25*megabyte)
						{
							fwrite(overlaps, sizeof(uint32_t), oTrack, oFile1); oTrack = 0;
						}

						overlaps[oTrack++] = cOvlq[3*((qFragmentIndex+qIndex->getIndexOffset()) - cOvlq[0]) + 1];
						overlaps[oTrack++] = cOvlr[3*((rFragmentIndex+rIndex->getIndexOffset()) - cOvlr[0]) + 1];
						
						float ovlInterval = 10; float idenInterval = .5; float idenLowerBound = 70; float cInterval = .03;
						float qLen = qBounds.second-qBounds.first; float rLen = rBounds.second-rBounds.first;

						uint32_t ovlLen = (alignScore.first)/ovlInterval;
						uint32_t ovlIden = (((alignScore.first - alignScore.second)/alignScore.first) * 100 - idenLowerBound)/idenInterval;

						//qContained = 1, rContained = 1 if the alignment length is at 97% of qFragment's, rFragment's length respectively

						uint32_t qContained = (alignScore.first/qLen + cInterval);
						uint32_t rContained = (alignScore.first/rLen + cInterval); 

						uint32_t rDovetail = 0;
						if((qAlignment.first < rAlignment.first) && !(qContained || rContained))
							rDovetail = 1;

						overlaps[oTrack] = 0;
						overlaps[oTrack++] = ((ovlLen << 22) | (ovlIden << 16) | (reverse << 15) | (qContained << 14) | (rContained << 13) | (rDovetail << 12));

						//Adding the back edge as well
						overlaps[oTrack++] = cOvlr[3*((rFragmentIndex+rIndex->getIndexOffset()) - cOvlr[0]) + 1];
						overlaps[oTrack++] = cOvlq[3*((qFragmentIndex+qIndex->getIndexOffset()) - cOvlq[0]) + 1];

						qContained = (alignScore.first/rLen + cInterval);
						rContained = (alignScore.first/qLen + cInterval);

						rDovetail = 0;
						if((rAlignment.first < qAlignment.first) && !(qContained || rContained))
							rDovetail = 1;

						overlaps[oTrack] = 0;
						overlaps[oTrack++] = ((ovlLen << 22) | (ovlIden << 16) | (reverse << 15) | (qContained << 14) | (rContained << 13) | (rDovetail << 12));
					}

				}

			}
		}
	}
}

//void collectHits(pair<unsigned long long int, unsigned long long int> &, set <Seed> &, const unsigned int	)
//Description: This function collects the hits that a fragment has in an index set
//Input: pair<unsigned long long int, unsigned long long int> & qBounds, the bounds of the query sequence in the fragment index, set<Seed> the set of seeds
//const unsigned int reverse, if reverse = 1, then the query segment is the reverse complement
//Output:None
//Return:None
	template <class T>
void Overlapper<T>::collectHits(pair<unsigned long long int, unsigned long long int> & qBounds, set<Seed> & hits, int hitStats [], unsigned long long int hitIndexes [], const unsigned int reverse)
{
	int repArr [64];int numKmers = 0; 
	memset(repArr, 0, 64 * sizeof(int));
	pair<unsigned long long int, unsigned long long > suffixBounds;

	if(qBounds.second - qBounds.first > 1000)
	{
		delete [] hitStats; delete [] hitIndexes;
		hitStats = new int [qBounds.second - qBounds.first];
		hitIndexes = new unsigned long long int [qBounds.second - qBounds.first];
	}

	for(unsigned long long int i = qBounds.first; i < qBounds.second-kmerSize; i+=stepSize)
	{
		bool goodKmer = true; float qValsTotal = 0;
		for(int j = 0; j < kmerSize; j++)
		{
			int tripleVal = ((4*4) * mapVals[qIndex->at(i+j)]) + (4 * mapVals[qIndex->at(i+j+1)]) + mapVals[qIndex->at(i+j+2)];
			repArr[tripleVal]++;

			if(repArr[tripleVal] >= 5)
			{
				goodKmer = false; break;
			}
			qValsTotal+=qIndex->qualAt(i+j);
		}


		if(qValsTotal/kmerSize < 25)
			goodKmer = false;

		memset(repArr, 0, 64 * sizeof(int));

		if(goodKmer)
		{
			hitIndexes[numKmers] = i;

			suffixBounds = fmIndex->kmerBounds(*qIndex, i, kmerSize, reverse);
			hitStats[numKmers++] = suffixBounds.second - suffixBounds.first; 
		}
	}

	sort(hitStats, hitStats+numKmers);
	unsigned int  outlierBounds = 0;

	if(numKmers >= 4)
	{
		int upperQuart = .75 * numKmers;
		int lowerQuart = .25 * numKmers;
		int innerRange = hitStats[upperQuart-1] - hitStats[lowerQuart-1];
		outlierBounds = upperQuart + innerRange;
	}

	if(outlierBounds == 0)
		outlierBounds = 200;

	for(int i = 0; i < numKmers; i++)
	{
		suffixBounds = fmIndex->kmerBounds(*qIndex, hitIndexes[i], kmerSize, reverse);
		if(suffixBounds.second - suffixBounds.first <= outlierBounds)
		{
			for(unsigned long long int j = suffixBounds.first; j < suffixBounds.second; j++)
			{
				unsigned long long int indexPosition = suffix->SA_to_SI(j);
				unsigned long long int lowIndex = rIndex->positionToIndex(indexPosition); //change
				unsigned long long int highIndex = rIndex->positionToIndex(indexPosition+kmerSize-1); //change

				if(lowIndex == highIndex)
				{
					Seed seed(hitIndexes[i], hitIndexes[i]+kmerSize-1-stepSize, indexPosition, indexPosition+kmerSize-1-stepSize);
					set<Seed>::iterator upperBound = hits.upper_bound(seed);
					bool merged = false;

					//Merging the seeds
					while(upperBound != hits.begin())
					{
						upperBound--;
						if((*upperBound).coord2X == hitIndexes[i]+kmerSize-stepSize-1 && (*upperBound).coord2Y == indexPosition+kmerSize-stepSize-1)
						{
							unsigned long long int lowIndex = rIndex->positionToIndex((*upperBound).coord1Y); //change
							if(lowIndex == highIndex)
							{
								Seed mergeSeed((*upperBound).coord1X, hitIndexes[i]+kmerSize-1, (*upperBound).coord1Y,  indexPosition+kmerSize-1);
								hits.erase(upperBound);
								hits.insert(mergeSeed);
								merged = true;
								break;
							}
						}

						if((*upperBound).coord2Y < indexPosition+kmerSize-stepSize-1)
							break;
					}

					if(!merged)
					{
						seed.coord2X = hitIndexes[i]+kmerSize-1; seed.coord2Y = indexPosition+kmerSize-1;
						seed.weight = seed.coord2X - seed.coord1X+1;
						hits.insert(seed);
					}
				}
			}

		}
	}

	if(qBounds.second - qBounds.first > 1000)
	{
		delete [] hitStats; delete [] hitIndexes;
		hitStats = new int [1000];
		hitIndexes = new unsigned long long int [1000];
	}
}

//bool passedQgram(Chain & chain, pair<unsigned long long int, unsigned long long int> &, pair<unsigned long long int, unsigned long long int> &, int [], const unsigned int )
//Description: This function determines if there are sufficient qGram matches to consider a chain of seeds a candidate for alignment
//Input: Chain & chain, the chain of seeds, pair<unsigned long long int, unsigned long long int> qBounds, the query sequences' bounds, pair<unsigned long long int, unsigned long long int> rBounds, the reference sequence's bounds
//int qGram [], the matrix that holds the qGram counts, const unsigned int reverse, if reverse = 1, then the query sequence is the reverse complement
//Output:None
//Return: bool, if the chain of seeds passed the qGram filter
	template <class T>
bool Overlapper<T>::passedQgram(Chain & chain, pair<unsigned long long int, unsigned long long int> & qBounds, pair<unsigned long long int, unsigned long long int> & rBounds, int qGram [], const unsigned int reverse)
{
	int maxChain = chain.longestChain();
	int qGramSize = 4;

	int lastHit = maxChain;
	int firstHit = lastHit;
	int nextHit = firstHit;
	while(nextHit != -1)
	{
		firstHit = nextHit;
		nextHit = chain.predecessor(firstHit);
	}

	unsigned long long int qBegin  = chain.at(firstHit).coord1X;
	unsigned long long int qEnd = chain.at(lastHit).coord2X;
	unsigned long long int rBegin  = chain.at(firstHit).coord1Y;
	unsigned long long int rEnd = chain.at(lastHit).coord2Y;

	int qGramTotal = 0;
	int qDiag = 0; int rDiag = 0;

	//Checking for qgrams in a window of 50 bps
	for(unsigned int i = 0; i < qEnd - qBegin; i+=50)
	{
		for(unsigned int j = 0; j < 50 && (j+i < qEnd-qGramSize-qBegin || (j+i < 50 && j+i < qBounds.second-qBegin-qGramSize)); j++)
		{
			qGram[encode(qIndex, i+j+qBegin, qGramSize, reverse)] = i+j;
			qDiag++;
		}

		for(unsigned int j = 0; j < 50 && (i+j < rEnd-qGramSize-rBegin || (j+i < 50 && i+j < rBounds.second-rBegin-qGramSize)); j++)
		{
			int pos = i+j;
			int match = qGram[encode(rIndex, i+j+rBegin, qGramSize, 0)];

			if(match != -1 && pow(match-pos, 2) < 1600)
				qGramTotal+=1;
			rDiag++;
		}
	}

	float norm = 100;
	float allowedError = 1-(identity - 5)/norm;
	float alignmentDiagonal = min(qDiag, rDiag);
	int estQgrams = alignmentDiagonal - qGramSize + 1 - (allowedError * alignmentDiagonal) * qGramSize;

	if(qGramTotal >= estQgrams)
	{

		return true;
	}

	return false;
}

//int encode(T * &, const unsigned long long int, const int, const unsigned int)
//Description: This function encodes a kmer in a sequence
//Input: T *, the fragment index, const long long int pos, the position in the index, const int len, the length of the kmer, const unsigned int reverse
//Output:None
//Return: int, the value of the kmer
template <class T>
int Overlapper<T>::encode(T * & index, const unsigned long long int pos, const int len, const unsigned int reverse){
	int rtrnVal = 0; int factor = len-1;
	for(int i = 0; i < len; i++)
	{
		rtrnVal+=mapVals[index->at(pos+i, reverse)]*pow(4, factor--);
	}
	return rtrnVal;
}

//pair<int, pair<int, int> alignSeqs(Chain & chain, pair<unsigned long long int, unsigned long long int> &, pair<unsigned long long int, unsigned long long int> &, int * [], const unsigned int )
//Description: This function actually aligns two sequences
//Input: Chain & chain, the seed chain, pair<unsigned long long int, unsigned long long int> & qBounds, the bounds of the query sequence, pair<unsigned long long int, unsigned long long int> & rBounds, the reference sequences bounds, int * matrix, the score matrix
//const unsigned int reverse, if reverse = 1, then the query sequence is the reverse complement
//Output:None
//Return: None
template <class T>
void Overlapper<T>::alignSeqs(Chain & chain, pair<unsigned long long int, unsigned long long int> & qBounds, pair<unsigned long long int, unsigned long long int> & rBounds, pair<float, float> & alignmentScore, pair<int, int> & qAlignment, pair<int, int> & rAlignment,  int * matrix [], const unsigned int reverse){

	int alignmentWeight = 0; int alignmentMisMatch = 0;
	int currChain = chain.longestChain();

	unsigned long long int xBegin = chain.at(currChain).coord2X+1; int xLen = qBounds.second - xBegin;
	unsigned long long int yBegin = chain.at(currChain).coord2Y+1; int yLen = rBounds.second - yBegin;

	float divisor = 100;
	float range = .10;
	float multFactor = max(range, 1-identity/divisor);

	int bandWidth = 2* min(multFactor*xLen, multFactor*yLen);
	bandWidth = max(bandWidth, 10);

	{
		BandedAlignment<Fragment_Index> align(xBegin, yBegin, xLen, yLen, bandWidth, true);
		if(xLen <= 1000 && yLen <= 1000)
		{
			alignmentScore = align.align(matrix, *qIndex, *rIndex, reverse, true);

		}else{

			int size = max(xLen, yLen);
			int ** largeMatrix = new int * [size+1];
			for(int i = 0; i < size+1; i++)
				largeMatrix[i] = new int[size+1];

			alignmentScore = align.align(largeMatrix, *qIndex, *rIndex, reverse, true);

			for(int i = 0; i < size+1; i++)
				delete [] largeMatrix[i];

			delete [] largeMatrix;
		}

		alignmentWeight+=(alignmentScore.first+chain.at(currChain).weight);
		alignmentMisMatch+=alignmentScore.second;
		qAlignment.second = align.getXEnd() + (xBegin - qBounds.first);
		rAlignment.second = align.getYEnd() + (yBegin - rBounds.first);
	}

	int nextChain = chain.predecessor(currChain);
	while(nextChain != -1)
	{
		xBegin = chain.at(nextChain).coord2X+1; xLen = chain.at(currChain).coord1X - xBegin;
		yBegin = chain.at(nextChain).coord2Y+1; yLen = chain.at(currChain).coord1Y - yBegin;

		bandWidth = 2* max(multFactor*xLen, multFactor*yLen);
		bandWidth = max(bandWidth, 10);

		BandedAlignment<Fragment_Index> align(xBegin, yBegin, xLen, yLen, bandWidth, false);

		pair<int, int> alignmentScore;
		if(xLen <= 1000 && yLen <= 1000)
		{
			alignmentScore = align.align(matrix, *qIndex, *rIndex, reverse);

		}else{
			int size = max(xLen, yLen);
			int ** largeMatrix = new int * [size];
			for(int i = 0; i < size; i++)
				largeMatrix[i] = new int[size];

			alignmentScore = align.align(largeMatrix, *qIndex, *rIndex, reverse);

			for(int i = 0; i < size; i++)
				delete [] largeMatrix[i];

			delete [] largeMatrix;
		}

		alignmentWeight+=(alignmentScore.first+chain.at(nextChain).weight);
		alignmentMisMatch+=alignmentScore.second;
		currChain = nextChain;
		nextChain = chain.predecessor(currChain);
	}

	xBegin = qBounds.first; xLen = chain.at(currChain).coord1X - xBegin;
	yBegin = rBounds.first; yLen = chain.at(currChain).coord1Y - yBegin;

	bandWidth = 2* min(multFactor*xLen, multFactor*yLen);
	bandWidth = max(bandWidth, 10);

	{
		BandedAlignment<Fragment_Index> align(xBegin, yBegin, xLen, yLen, bandWidth, true);

		if(xLen <= 1000 && yLen <= 1000)
		{
			alignmentScore = align.alignReverse(matrix, *qIndex, *rIndex, reverse, true);

		}else{
			int size = max(xLen, yLen);
			int ** largeMatrix = new int * [size];
			for(int i = 0; i < size; i++)
				largeMatrix[i] = new int[size];

			alignmentScore = align.alignReverse(largeMatrix, *qIndex, *rIndex, reverse, true);

			for(int i = 0; i < size; i++)
				delete [] largeMatrix[i];

			delete [] largeMatrix;
		}

		alignmentWeight+=alignmentScore.first;
		alignmentMisMatch+=alignmentScore.second;
		qAlignment.first = align.getXStart();
		rAlignment.first = align.getYStart();

	}

	alignmentScore.first = alignmentWeight; alignmentScore.second = alignmentMisMatch;
}

#endif /* OVERLAPPER_H_ */
