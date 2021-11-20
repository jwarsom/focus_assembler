/*
 * Fragment_Index.cpp
 *
 *      Author: Julia
 */

#include<string>
#include<iostream>
#include<math.h>
#include<utility>
#include<cstdio>
#include<stdint.h>
#include<cstring>
using namespace std;
#include "MappingValues.h"
#include "Dictionary.h"
#include "Fragment_Index.h"

#define resizeVal 10485760

/* Class Description */
// This Fragment_Index class indexes and stores Fragments
// It also stores the Fragment Names and 0 - 64 qual vals

//Fragment_Index(const int)
//Description: Default constructor, initializes everything to null except for the tagLen
//Input:t (an int) The length of the fragment labels
//Output:None
//Return:None
Fragment_Index::Fragment_Index()
	:boundsData(nullptr), fragments(nullptr), indexOffset(0), fragmentOffset(0), totalFragments(0), indexLength(0){;;;}


//~Fragment_Index(void)
//Description: Class destructor: deletes the fragments and labels, bounds is deleted by the dictionary
//Input:None
//Output:None
//Return:None
Fragment_Index::~Fragment_Index(void) {
	delete [] fragments;
}

//( const char &, const char &, const int )
//Description: This function will encode a nucleotide and qual into a
//single char
//Input: nucleotide ( a char) the nucleotide, qual (a char) the quality value, i (an int) the position of the nucleotide
//Output:none
//Return: None
void Fragment_Index::encode(const char & nucleotide, const char & qual, const int i) {

	int subscript = nucleotide;
	unsigned char appVal = mapVals[subscript];
	fragments[i] = (appVal | ((qual-33) << 2));
}

//decode(char &, char &, const int) const
//Description: This function decodes a combined nucleotide and quality value
//Input: nucleotide (a char), qual (a char), i (const int) the position of the nucleotide
//Output:None
//Return:None
void Fragment_Index::decode(char & nucleotide, char & qual, const int i) const {

	nucleotide = inverseMapVals[fragments[i] & 3];
	qual = fragments[i] >> 2;
}

//Description: This function returns a quality value at position i of the index
//Input:  i (const int) the position of the nucleotide
//Output:None
//Return: Teh quality value at i
int Fragment_Index::qualAt(const unsigned long long int i, const unsigned int reverse) const {
	return ((fragments[i] >> 2) - reverse * (fragments[i] >> 2)) | (reverse * fragments[indexLength-i-1] >> 2);
}

//at(const int) const
//Description: This function returns the character at an index
//position in the dataset
//Input: i (an int), the position of the character in the index
//Output:none
//Return:None
char Fragment_Index::at(const unsigned long long int i, const unsigned int reverse) const {

	return ((inverseMapVals[(fragments[i] & 3)] - reverse * inverseMapVals[(fragments[i] & 3)]) | (reverse * inverseComplementMapVals[(fragments[indexLength-i-1] & 3)]));
}

//uint64_t ithWords(const unsigned long long int, const unsigned int)
//Description: This function implicitly divides the fragment index into words of size l and returns the nth word
//Input:n, (const unsigned long long int) the rank of the word we are trying to find, l, (const unsigned int) the word size
//Output:None
//Return: a uint64_t the encoded value of the nth word
uint64_t Fragment_Index::ithWord(const unsigned long long int n, const unsigned int l, const unsigned int reverse){

		uint64_t rtrnVal = 0; int factor = l-1;
		for(unsigned int i = 0; i < l; i++)
		{
			int subscript = at(n*l+i, reverse);	
			rtrnVal+=mapVals[subscript]*pow(4, factor--);
		}

		return rtrnVal;
}

//length(void)
//Description: This function returns the length of the index
//ie. how many characters it contains
//Input:none
//Output:None
//Return: int , the size of the index
unsigned long long int Fragment_Index::length(void) const {
	return indexLength;
}

//numFragments(void)
//Description: This function returns the number of fragments contained in
//the index
//Input:None
//Output:None
//Return: int, the number of fragments in the index
unsigned long long int Fragment_Index::numFragments(void) const {
	return totalFragments;
}

//append(const char * &, const char * &, const char * &)
//Description: This function appends fragment and its quality values to the index
//Input: dna ( char *), quals (char *), tag (char *) the label of the fragment
//Output: None
//Return: none
void Fragment_Index::append(const char * & dna, const char * & quals) {

	int size = strlen(dna);
	for(int i = 0; i < size; i++)
		encode(dna[i], quals[i], indexLength++);

	boundsData[(indexLength+fragmentOffset%64)/64+1] = boundsData[(indexLength+fragmentOffset%64)/64+1] | (appValL >> (indexLength+fragmentOffset)%64);

	totalFragments++;
}

//append(const char * , const char * , const char * )
//Description: This function appends fragment and its quality values to the index
//Input: dna ( char *), quals (char *), tag (char *) the label of the fragment
//Output: None
//Return: none
void Fragment_Index::append(const char * dna, const char * quals) {

	int size = strlen(dna);

	for(int i = 0; i < size; i++)
	{
		encode(dna[i], quals[i], indexLength++);
	}

	boundsData[(indexLength+fragmentOffset%64)/64+1] = boundsData[(indexLength+fragmentOffset%64)/64+1] | (appValL >> (indexLength+fragmentOffset)%64);

	totalFragments++;

}


//append(const char * , const char * , unsigned int )
//Description: This function appends fragment and its quality values to the index
//Input: dna ( char *), quals (char *), size (int) the size of the fragment
//Output: None
//Return: none
void Fragment_Index::append(const char * dna, const char * quals, unsigned int size) {

	for(unsigned int i = 0; i < size; i++)
	{
		encode(dna[i], quals[i], indexLength++);
	}

	boundsData[(indexLength+fragmentOffset%64)/64+1] = boundsData[(indexLength+fragmentOffset%64)/64+1] | (appValL >> (indexLength+fragmentOffset)%64);

	totalFragments++;

}

//append(const char * , const char * &, unsigned int )
//Description: This function appends fragment and its quality values to the index
//Input: dna ( char *), quals (char *), size, (int) the size of the fragment
//Output: None
//Return: none
void Fragment_Index::append(const char * & dna, const char * & quals, unsigned int size) {

	for(unsigned int i = 0; i < size; i++)
		encode(dna[i], quals[i], indexLength++);

	boundsData[(indexLength+fragmentOffset%64)/64+1] = boundsData[(indexLength+fragmentOffset%64)/64+1] | (appValL >> (indexLength+fragmentOffset)%64);

	totalFragments++;

}

//positionToIndex(const int) const
//Description: This function converts a position in the index to
//its corresponding fragment's index
//Input: pos (an int), the position in the index
//Output:None
//Return: int, the index of the corresponding fragment
unsigned long long int Fragment_Index::positionToIndex(const unsigned long long int pos, const unsigned int reverse) {

	if(!reverse)
	{
		return bounds.rank1(pos+1+fragmentOffset%64) - 1;
	}
	else
	{
		return bounds.rank1(indexLength - pos+1+fragmentOffset%64 - 1) - 1;
	}
}

//pair<int, int> indexBounds(const int i) const
//Description: This function returns the bounds of a
//fragment in the index
//Input: i (an int), the index of the fragment
//Output:None
//Return: a pair of long ints, the bounds of the fragment
pair<unsigned long long int, unsigned long long int> Fragment_Index::indexBounds(const unsigned long long int i, const unsigned int reverse) {

	if(!reverse)
	{
		return pair<unsigned long long int, unsigned long long int>(bounds.select1(i)-fragmentOffset%64, bounds.select1(i+1)-fragmentOffset%64);

	}else{
		return pair<unsigned long long int, unsigned long long int>(indexLength - (bounds.select1(i+1)-fragmentOffset%64), indexLength-(bounds.select1(i)-fragmentOffset%64));
	}
}

//reserve(const int, const int, const int)
//Description: this function reserves space for the fragments, bounds, and labels arrays
//Input: n (an int) the number of fragments to reserve, f (an int) the length of the fragment
//Ouput: None
//Return: None
void Fragment_Index::reserve(const unsigned long long int n, const unsigned long long int f, const unsigned long long int i_offset, const unsigned long long int f_offset) {
	
	clear();

	boundsData = new uint64_t[((f+1)+63+f_offset%64)/64+1];
	fragments = new unsigned char[f];
	memset(fragments, 0, f);
	memset(boundsData, 0, (((f+1)+63+f_offset%64)/64+1)*sizeof(uint64_t));

	boundsData[1] = appValL >> f_offset%64;

	indexOffset = i_offset;
	fragmentOffset = f_offset;
}

//long int getIndexOffset(void)
//Description: this function gets the indexOffset
//Input: None
//Output: None
//Return: a long int, the indexOffset
unsigned long long int Fragment_Index::getIndexOffset(void) const {
	return indexOffset;
}

//unsigned long long int getFragmentOffset(void) const
//Description: this function get the fragmentOffset
//Input: None
//Output: None
//Return: a long int, the fragmentOffset
unsigned long long int Fragment_Index::getFragmentOffset(void) const {
	return fragmentOffset;
}

//void set(void)
//Description:This function sets the dictionary
//This function must be called before the fragment array can be used
//Input:None
//Output:Mome
//Return: None
void Fragment_Index::set(void) {
	bounds.setBits(boundsData, totalFragments+1, indexLength+1+fragmentOffset%64);
	bounds.setRanks();  bounds.setSelect1();
}

//void clean(void)
//Description: this function clears the fragment index
//Input: None
//Output: None
//Return: None
void Fragment_Index::clear(void) {
	bounds.clear();
	delete [] fragments;
	boundsData = nullptr;
	fragments = nullptr;

	indexOffset = 0; 
        fragmentOffset = 0; 
	totalFragments = 0;
	indexLength = 0;
}

//void write(FILE * & pFile)
//Description: This function writes the fragment index to a file
//Input: pFile, a FILE pointer
//Output: None
//Return: None
void Fragment_Index::write(FILE * & pFile){

	fwrite(&fragmentOffset, sizeof(unsigned long long int), 1, pFile);
	fwrite(&indexOffset, sizeof(unsigned long long int), 1, pFile);

	fwrite(&indexLength, sizeof(unsigned long long int), 1, pFile);
	fwrite(&totalFragments, sizeof(unsigned long long int), 1, pFile);

	fwrite(fragments, sizeof(char), indexLength, pFile);

	unsigned long long int boundsLength = ((indexLength+1)+63+fragmentOffset%64)/64;
	fwrite(&boundsLength, sizeof(unsigned long long int), 1, pFile);
	fwrite(boundsData+1, sizeof(uint64_t), boundsLength, pFile);
}

//void read(FILE * & pFile)
//Description: This function reads a fragment index from file
//Input: pFile, a FILE pointer
//Output: None
//Return: None
void Fragment_Index::read(FILE * & pFile){

	unsigned long long int intervalFragmentOffset;
	unsigned long long int intervalIndexOffset;
	fread(&intervalFragmentOffset, sizeof(unsigned long long int), 1, pFile);
	fread(&intervalIndexOffset, sizeof(unsigned long long int), 1, pFile);

	unsigned long long int intervalFragmentLength;
	fread(&intervalFragmentLength, sizeof(unsigned long long int), 1, pFile);

	unsigned long long int intervalTotalFragments;
	fread(&intervalTotalFragments, sizeof(unsigned long long int), 1, pFile);

	unsigned long long int startPos = intervalFragmentOffset - fragmentOffset;
	fread(fragments+startPos, sizeof(unsigned char), intervalFragmentLength, pFile);

	unsigned long long int intervalBoundsLength;
	fread(&intervalBoundsLength, sizeof(unsigned long long int), 1, pFile);

	uint64_t mergePoint;
	fread(&mergePoint, sizeof(uint64_t), 1, pFile);

	startPos = (intervalFragmentOffset - fragmentOffset + fragmentOffset%64)/64+1;
	boundsData[startPos] = boundsData[startPos] | mergePoint;
	fread(boundsData+startPos+1, sizeof(uint64_t), intervalBoundsLength-1, pFile);
    
	totalFragments+=intervalTotalFragments;
	indexLength+=intervalFragmentLength;
}
