/*
 * FM_Index.h
 *
 *      Author: Julia
 */

#ifndef FM_INDEX_H_
#define FM_INDEX_H_
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include "math.h"
using namespace std;

template <class T>
class FM_Index {
	public:
		FM_Index();
		virtual ~FM_Index();
		void calcBWT(const SuffixArray<T> &, const T &, const unsigned long long int);
		char atBWT(const int) const;
		long int numOcc(const char, const int) const;
		long int numLexSmaller(const int) const;
		pair<unsigned long long int, unsigned long long int> kmerBounds(T &, const unsigned long long int, const int, const unsigned int reverse = 0);
	private:
		void encode(const char, const int);
		void decode(char &, const int) const;
		unsigned char  * bwt;
		long int ** occuranceATCG;
		long int * tallyATCG;
		long int * lexSmallerCounts;
		long int length;
		long int endLoc;
};

//FM_Index()
//Description: This is the class constructor
//Input:None
//Output:None
//Return:None
	template <class T>
FM_Index<T>::FM_Index()
	:bwt(nullptr), occuranceATCG(nullptr), tallyATCG(nullptr), lexSmallerCounts(nullptr), length(0), endLoc(0){;;;}

//~FM_Index
//Description:This is the class destructor
//Input:None
//Output:None
//Return:None
template <class T>
FM_Index<T>::~FM_Index() {
	delete [] bwt;
	for(int i = 0; i < length/128+1; i++)
	{
		delete [] occuranceATCG[i];
	}

	delete [] occuranceATCG;
	delete [] tallyATCG;
	delete [] lexSmallerCounts;
}

//calcBWT
//Description:This function calculates the BWT string from the suffix array
//and genome or fragment index
//Input:suffix (SuffixArray &), the suffix array, index ( the fragment or genome index)
//Output:None
//Return:None
template <class T>
void FM_Index<T>::calcBWT(const SuffixArray<T> & suffix, const T & index, const unsigned long long int numChars){

		long int size = numChars+1;
		length = size;
		long int cSize = (size+3)/4;

		bwt = new unsigned char[cSize];
		memset(bwt, 0, cSize);
		occuranceATCG = new long int * [size/128+1];

		for(int i = 0; i < size/128+1; i++)
		{
			occuranceATCG[i] = new long int[4];
			memset(occuranceATCG[i], 0, sizeof(long int) * 4);
		}

		tallyATCG =  new long int[4];
		memset(tallyATCG, 0, sizeof(long int) * 4);
		lexSmallerCounts = new long int [5];
		memset(lexSmallerCounts, 0, sizeof(long int) * 5);
		lexSmallerCounts[0] = 1;

		endLoc = suffix.ISA_to_SA(0);
		encode('A', endLoc);

		for(int i = 0; i < size; i++)
		{
			if(i != endLoc)
			{
				tallyATCG[mapVals[index.at(suffix.SA_to_SI(i)-1)]]++;
				encode(index.at(suffix.SA_to_SI(i)-1), i);
			}

			if((i+1)%128 == 0)
			{
				occuranceATCG[(i+1)/128][0] = tallyATCG[0];
				occuranceATCG[(i+1)/128][1] = tallyATCG[1];
				occuranceATCG[(i+1)/128][2] = tallyATCG[2];
				occuranceATCG[(i+1)/128][3] = tallyATCG[3];
			}
		}

		for(int i = 1; i < 5; i++)
		{
			lexSmallerCounts[i] = lexSmallerCounts[i-1]+tallyATCG[i-1];
		}
}

//kmerBounds(T &, unsigned long long int, int)
//Description: This function finds kmer bounds in the hierarchical array
//Input: T &, the index that contains the kmer, unsigned long long int start, the beginning of the kmer, int length, the length of the kmer
//Output: None
//Return: pair<unsigned long long int, unsigned long long int>, the bounds of the kmer in the suffix array
	template <class T>
pair<unsigned long long int, unsigned long long int> FM_Index<T>::kmerBounds(T & index, const unsigned long long int start, const int kmerLen, const unsigned int reverse)
{
	unsigned long long int currPos = kmerLen-1; char c = index.at(start+currPos, reverse);
	unsigned long long int first = numLexSmaller(mapVals[c]); unsigned long long int last = numLexSmaller(mapVals[c] + 1)-1;

	while(first <= last && currPos >= 1)
	{
		c = index.at(start+currPos-1, reverse);
		first = numLexSmaller(mapVals[c]) + numOcc(c, first-1);
		last = numLexSmaller(mapVals[c]) + numOcc(c, last)-1;
		currPos = currPos-1;
	}

	return pair<unsigned long long int, unsigned long long int> (first, last+1);
}

//atBWT(const int)
//Description:This function returns a nucleotide in the BWT string
//Input:i (an int), the position of the nucleotide
//Output:None
//Return:char, the nucleotide at position i in the BWT string
template <class T>
	char FM_Index<T>::atBWT(const int i) const {
		if(i == endLoc)
			return '$';

		char nucleotide;
		decode(nucleotide, i);
		return nucleotide;
}

//numOcc(const char, const int) const
//Description: This function determines the number of the occurances
//of a character c in the BWT string in the range BWT[0, i]. Note that
//this range includes the boundary values 0 and i
//Input:c (a char) the character we are counting, i (an int), the upper bound
//of the range that we are searching
//Output:None
//Return: long int, the number of occurances of char c in the range[0, i]
template <class T>
long int FM_Index<T>::numOcc(const char c, const int i) const {

	int charType = mapVals[c];
	long int count = occuranceATCG[(i+1)/128][charType];
	int leftOver = (i+1)%128;
	int startPos = (i+1)-leftOver;

	for(int j = startPos; j < i+1; j++)
	{
		if(atBWT(j) == c) //this used to have if(j != endPos) .. I don't see why I needed that
		{
			count++;
		}
	}
	return count;
}

//numLexSmaller(const char)
//Description: This function returns the number of characters that are lexographically
//smaller than the character c in the BWT string
//Input: c, the char that we are using to tabulate the lexographically smaller characters
//Output:None
//Return: long int, the number of occurrences of characters that are lexographically smaller
//than c in the BWT string
template <class T>
long int FM_Index<T>::numLexSmaller(const int c) const {

	return lexSmallerCounts[c];
}

//encode(const char, const int)
//Description: This function encodes a nucleotide into a 2-bit representation
//Input:nucleotide (a char), i (const int), the position of the nucleotide
//Output:None
//Return:None
template <class T>
void FM_Index<T>::encode(const char nucleotide, const int i) {
	unsigned char c = mapVals[nucleotide];
	c = c << (6 - 2 * (i%4));
	bwt[i/4] = bwt[i/4] | c;
}

//decode(char &, const int) const
//Description:This function decodes a 2-bit representation of a nucleotide
//Input: nucleotide (a char) the value will be copied into this nucleotide, i (const int)
//the position of the nucleotide
//Output:None
//Return:None
template <class T>
void FM_Index<T>::decode(char & nucleotide, const int i) const{
	nucleotide = inverseMapVals[((bwt[i/4] >> (6 - 2 * (i%4))) & 3)];
}

#endif /* FM_INDEX_H_ */
