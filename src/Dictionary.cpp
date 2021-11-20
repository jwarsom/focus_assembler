/*
 * Dictionary.cpp
 *
 *      Author: Julia
 */
#include <iostream>
#include <stdint.h>
#include <cstdio>
#include <cstring>
using namespace std;
#include "math.h"
#include "stdint.h"
#include "math.h"
#include "Dictionary.h"


/* This is a succinct dictionary data structure that supports rank and select in O(1) time
 * while still using space that is very close to the information theoretic limit
 * The concept of the succinct data structure was introduced by Jacobson to encode bit vectors
 */
/* This provides a two level indexing system for the dictionary
 * This data structure is modeled after the one introduced by
 * Sebastiano Vigna for 64 bit systems: Ref: Broadword Implementation of Rank/Select Queries
 * Not theoretically optimal by much more practical than other algorithms
 * Runs much faster
 */

#define offset 2
#define data_offset 1
#define binCount 8192

static int cArr[] = {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};


//Dictionary()
//Description: This is the constructor for the dictionary
//Input:None
//Ouput:None
//Return:None
Dictionary::Dictionary()
	:data(nullptr), bounds(nullptr), capacity(0), numMembers(0), universeSize(0), ones_inventory(nullptr), ones_subInventory(nullptr), ones_spill(nullptr), ones_bin(0), ones_numBins(0), ones_wordsPerSubBin(0), ones_numSubBins(0), ones_numSpill(0),
	zeros_inventory(nullptr), zeros_subInventory(nullptr), zeros_spill(nullptr), zeros_bin(0), zeros_numBins(0), zeros_wordsPerSubBin(0), zeros_numSubBins(0), zeros_numSpill(0){;;;}

	//~Dictionary()
	//Description: This is the destructor for the dictionary
	//Input:None
	//Output:None
	//Return:None
	Dictionary::~Dictionary() {
		delete [] data;
		delete [] bounds;
		delete [] ones_inventory;
		delete [] ones_subInventory;
		delete [] ones_spill;
		delete [] zeros_inventory;
		delete [] zeros_subInventory;
		delete [] zeros_spill;
	}


//setBits(uint64_t * &, const long long int, const long long int)
//Description:This function sets the bit members in the dictionary
//Input: d (uint64_t * &): a pointer to the bit array, m (const long long int): the
//number of members in the dictionary, n(const long long int): the size of the universe
//Output:None
//Return:None
void Dictionary::setBits(uint64_t * & d, const long long int m, const long long int n) {


	delete [] data;
	numMembers = m; universeSize = n; ones_numBins = 0; zeros_numBins = 0;
	capacity = (n+63)/64+1;

	ones_bin = max(floor(binCount * (long double) m /(long double) n), /*64.0*/(long double) 72.0); //CHANGE MINS 
	zeros_bin = max(floor(binCount * (long double) (n-m) /(long double) n), /*64.0*/(long double) 72.0); //CHANGE MIS 
	data = d;

	double max_words_per_sub_bin = 3;
	ones_wordsPerSubBin = min(max_words_per_sub_bin, floor(log(m/ones_bin+2)/log(2)));
	zeros_wordsPerSubBin = min(max_words_per_sub_bin, floor(log((n-m)/zeros_bin+2)/log(2))); //change from ones_bin to zeros_bin

	if(ones_wordsPerSubBin == 0) ones_wordsPerSubBin = 1;
	if(zeros_wordsPerSubBin == 0) zeros_wordsPerSubBin = 1; 

	ones_bin-=ones_bin%(4*ones_wordsPerSubBin);   
	zeros_bin-=zeros_bin%(4*zeros_wordsPerSubBin);    

	ones_numSubBins = 0; zeros_numSubBins = 0;
}


/* All of the accessor functions */

//getCapacity(void) const
//Description: Returns the dictionary size (in this case, the number of uint64_t)
//Input:None
//Output:None
//Return: The size of the bits arry ( in uint64_t)
long long int Dictionary::getCapacity(void) const {
	return capacity;
}

//long long int getNumMembers(void)
//Description: return the number of members in the
//dictionary
//Input:None
//Output:None
//Return:None
long long int Dictionary::getNumMembers(void) const {
	return numMembers;
}

//long long int getUniverseSize(void)
//Description: return the size of the dictionary universe
//Input:None
//Output:None
//Return:None
long long int Dictionary::getUniverseSize(void) const {
	return universeSize;
}

//uint64_t * getData(void)
//Description: returns the bit data
//Input: None
//Output:None
//Return: uint64_t * (data)
uint64_t * Dictionary::getData(void) const {
	return data;
}

//clear(void)
//Description: This function clears the dictionary data
//Input:None
//Output:None
//Return:None
void Dictionary::clear(void){

	delete [] data;
	delete [] bounds;
	delete [] ones_inventory;
	delete [] ones_subInventory;
	delete [] ones_spill;
	delete [] zeros_inventory;
	delete [] zeros_subInventory;
	delete [] zeros_spill;
	data = nullptr; bounds = nullptr; ones_inventory = nullptr;
	ones_subInventory = nullptr; ones_spill = nullptr; zeros_inventory = nullptr;
	zeros_subInventory = nullptr; zeros_spill = nullptr;
	ones_bin = 0; ones_numBins = 0; ones_wordsPerSubBin = 0; ones_numSubBins = 0;
	ones_numSpill = 0;
	zeros_bin = 0; zeros_numBins = 0; zeros_wordsPerSubBin = 0; zeros_numSubBins = 0;
	zeros_numSpill = 0; capacity = 0; numMembers = 0; universeSize = 0;
}

//setRanks(void)
//Description:This function populates the rank data structure
//Input:None
//Output:None
//Return:None
void Dictionary::setRanks(void){

	delete [] bounds;
	long int bounds_cap = (capacity+7-1)/8;
	bounds_cap = bounds_cap * 2;
	bounds = new uint64_t[bounds_cap+offset];
	memset(bounds, 0, sizeof(uint64_t) * offset);

	for(long long int size = 0; size < capacity-1; size++)
	{
		//Get the locations in the bounds array
		int64_t  div1 = size/8;
		int64_t  mod1 = size%8 - 1;
		int64_t  mod2 = (size-1)%8 - 1;

		uint64_t x = bounds[(2* div1 + offset) - 2 * (mod1 >> 60 & 8)/8];
		x += bounds[(2*div1 + offset) - 2 * (mod1 >> 60 & 8)/8 + 1] >> ( mod2  + ( mod2 >> 60 & 8)) * 9 & 0x1FF;

		//sideways addition
		uint64_t s = sidewaysAddition(data[size]);
		x+=s;

		//Update the bounds
		bounds[2*div1+offset] = bounds[2*div1+offset] & (bounds[2*div1+offset] * (1-((mod1 >> 60 & 8)/8)) );
		bounds[2*div1+offset] = bounds[2*div1+offset] | (x * (mod1 >> 60 & 8)/8);
		bounds[2*div1+offset+1] = bounds[2*div1+offset+1] & (bounds[2*div1+offset+1] * (1-((mod1 >> 60 & 8)/8)));

		x = x - bounds[2*div1 + offset];
		bounds[2*div1+offset+1] = bounds[2*div1+offset+1] | x << (9 * mod1);
	}

}

//bool isSet(const unsigned long long int)
//Description: This returns true if a bit is set in the dictionary
//Input: n, (an unsigned long long int) the position of the bit
//Output:None
//Return: a bool, true if the bit is set
bool Dictionary::isSet(const unsigned long long int n) const{
	uint64_t x = data[n/64+1];
	return x & (0x8000000000000000 >> n%64);
}


//setSelect1(void)
//Description:This function populates that select structure that allows
//for the selection of a ones bit with a given ones rank
//Input:None
//Output:None
//Return:None
void Dictionary::setSelect1(void) {

	delete [] ones_inventory;
	delete [] ones_subInventory;
	delete [] ones_spill;

	ones_inventory = new uint64_t[numMembers/ones_bin + 2];
	ones_inventory[0] = 0; ones_numBins++; ones_inventory[ones_numBins] = 0;


	for(long long int size = 0; size < capacity-1; size++)
	{
		//Get the locations in the bounds array
		int64_t  div1 = size/8;
		int64_t  mod1 = size%8 - 1;
		int64_t  mod2 = (size-1)%8 - 1;

		uint64_t x = bounds[(2* div1 + offset) - 2 * (mod1 >> 60 & 8)/8];
		x += bounds[(2*div1 + offset) - 2 * (mod1 >> 60 & 8)/8 + 1] >> ( mod2  + ( mod2 >> 60 & 8)) * 9 & 0x1FF;

		//sideways addition
		uint64_t s = sidewaysAddition(data[size]); uint64_t l = sidewaysAddition(data[size+1]);
		x+=s;

		//Track the position of the nth ranked bit for the select search
		ones_inventory[0] = (ones_inventory[0] | (data[size+1] >> (63-computeBitPos(data[size+1], 1))) * ((size * 64 + computeBitPos(data[size+1], 1))) * (1-(x >> (63 - computeBitPos(x, 1))))); //added Paren around 63 - computeBit
		ones_inventory[ones_numBins] = (1-(x+l)/(ones_bin*ones_numBins+1))* ones_inventory[ones_numBins] | (x+l)/(ones_bin*ones_numBins+1) * ((64 * size) + computeBitPos(data[size+1], (ones_bin*ones_numBins)-x+1));

		//Figure this out note
		if(((ones_inventory[ones_numBins] - ones_inventory[ones_numBins-1])/(1ULL << 16))*((x+l)/(ones_bin*ones_numBins+1)))
			ones_numSpill+=1;

		ones_numBins+=(x+l)/(ones_bin*ones_numBins+1);
		ones_inventory[ones_numBins] = ((size * 64)+computeBitPos(data[size+1], sidewaysAddition(data[size+1])));
	}

	ones_inventory[ones_numBins]+=1;

	//Figure this out
	if((ones_inventory[ones_numBins] - ones_inventory[ones_numBins-1])/(1ULL << 16))
		ones_numSpill+=1;

	ones_numSubBins = (ones_numBins+1) * ones_wordsPerSubBin;
	ones_subInventory = new uint64_t [ones_numSubBins];
	ones_spill = new uint64_t [max((ones_numSpill*ones_bin), 1LL)];
	memset(ones_subInventory, 0, ones_numSubBins * sizeof(uint64_t));
	memset(ones_spill, 0, max((ones_numSpill*ones_bin), 1LL) * sizeof(uint64_t));

	uint64_t spillSize = 0;

	for(long long int i = 1; i < ones_numBins+1; i++)
	{
		if((ones_inventory[i] - ones_inventory[i-1])/(1ULL << 16))
		{
			ones_subInventory[(i-1) * ones_wordsPerSubBin] = spillSize;
			for(unsigned long long int j = ones_inventory[i-1]+1; j < ones_inventory[i]-i/ones_numBins; j++)
			{
				ones_spill[spillSize] = (j - ones_inventory[i-1]) * (rank1(j+1)-rank1(j));
				spillSize+=(rank1(j+1)-rank1(j));
			}
			//check first for errors 
			ones_spill[spillSize-1+i/ones_numBins] = ones_spill[spillSize-1+i/ones_numBins] | ((ones_inventory[i]-1 - ones_inventory[i-1]) * (i/ones_numBins));

		}else{

			uint64_t subBinSize = (ones_bin)/(4*ones_wordsPerSubBin);
			int subInventorySize = 1; unsigned int bitCount = 0; uint16_t * p16;
			p16 = (uint16_t *)&ones_subInventory[(i-1) * ones_wordsPerSubBin];
			for(unsigned int j = 1; j < ones_inventory[i] - ones_inventory[i-1]-i/ones_numBins; j++) //SEE WHAT HAPPENS
			{
				bitCount+=(rank1(ones_inventory[i-1]+j+1)-rank1(ones_inventory[i-1]+j));
				p16[subInventorySize-1] = (bitCount/(subInventorySize*subBinSize))*j;
				subInventorySize+=bitCount/(subInventorySize*subBinSize);
			}

			if(i/ones_numBins)
			{
				bitCount++;
				p16[subInventorySize-1] = (bitCount/(subInventorySize*subBinSize))*((ones_inventory[i]-1) - ones_inventory[i-1]);
				subInventorySize+=bitCount/(subInventorySize*subBinSize);
			}
		}
	}
}

//setSelect0(void)
//Description:This function populates that select structure that allows
//for the selection of a zero bit with a given zero rank
//Input:None
//Output:None
//Return:None
void Dictionary::setSelect0(void) {

	delete [] zeros_inventory;
	delete [] zeros_subInventory;
	delete [] zeros_spill;

	zeros_inventory = new uint64_t[(universeSize-numMembers)/zeros_bin + 2];
	zeros_inventory[0] = 0; zeros_numBins++; zeros_inventory[zeros_numBins] = 0;
	data[0] = 0xffffffffffffffff;

	for(long long int size = 0; size < capacity-1; size++)
	{
		//Get the locations in the bounds array
		int64_t  div1 = size/8;
		int64_t  mod1 = size%8 - 1;
		int64_t  mod2 = (size-1)%8 - 1;

		uint64_t x = bounds[(2* div1 + offset) - 2 * (mod1 >> 60 & 8)/8];
		x += bounds[(2*div1 + offset) - 2 * (mod1 >> 60 & 8)/8 + 1] >> ( mod2  + ( mod2 >> 60 & 8)) * 9 & 0x1FF;

		if(size != 0)
			x = 64 * (size-1) - x;

		//sideways addition
		uint64_t s = sidewaysAddition(~data[size]);
		uint64_t l = sidewaysAddition(~data[size+1]);
		x+=s;

		//Track the position of the nth ranked bit for the select search
		zeros_inventory[0] = (zeros_inventory[0] | (~data[size+1] >> (63-computeBitPos(~data[size+1], 1))) * ((size * 64 + computeBitPos(~data[size+1], 1))) * (1-(x >> (63 - computeBitPos(x, 1))))); //add paren around 63 - computeBit
		zeros_inventory[zeros_numBins] = (1-(x+l)/(zeros_bin*zeros_numBins+1))* zeros_inventory[zeros_numBins] | (x+l)/(zeros_bin*zeros_numBins+1) * ((64 * size) + computeBitPos(~data[size+1], (zeros_bin*zeros_numBins)-x+1));
		zeros_numSpill+=((zeros_inventory[zeros_numBins] - zeros_inventory[zeros_numBins-1])/(1ULL << 16))*((x+l)/(zeros_bin*zeros_numBins+1));
		zeros_numBins+=(x+l)/(zeros_bin*zeros_numBins+1);
		zeros_inventory[zeros_numBins] = ((size * 64)+computeBitPos(~data[size+1], sidewaysAddition(~data[size+1])));
	}

	zeros_inventory[zeros_numBins]+=1;
	zeros_numSpill+=(zeros_inventory[zeros_numBins] - zeros_inventory[zeros_numBins-1])/(1ULL << 16);
	zeros_numSubBins = (zeros_numBins+1) * zeros_wordsPerSubBin;
	zeros_subInventory = new uint64_t [zeros_numSubBins];
	zeros_spill = new uint64_t [max((zeros_numSpill*zeros_bin), 1LL)];
	memset(zeros_subInventory, 0, zeros_numSubBins * sizeof(uint64_t));
	memset(zeros_spill, 0, zeros_numSpill * sizeof(uint64_t));
	uint64_t spillSize = 0;

	for(long long int i = 1; i < zeros_numBins+1; i++)
	{
		if((zeros_inventory[i] - zeros_inventory[i-1])/(1ULL << 16))
		{
			zeros_subInventory[(i-1) * zeros_wordsPerSubBin] = spillSize;
			for(unsigned long long int j = zeros_inventory[i-1]+1; j < zeros_inventory[i]-i/zeros_numBins; j++)
			{
				zeros_spill[spillSize] = (j - zeros_inventory[i-1]) * (rank0(j+1)-rank0(j));
				spillSize+=(rank0(j+1)-rank0(j));
			}
			//check here for errors 
			zeros_spill[spillSize-1+i/zeros_numBins] = zeros_spill[spillSize-1+i/zeros_numBins] | ((zeros_inventory[i]-zeros_inventory[i-1]-1) * (i/zeros_numBins));

		}else{
			uint64_t subBinSize = (zeros_bin)/(4*zeros_wordsPerSubBin);
			int subInventorySize = 1; unsigned int bitCount = 0; uint16_t * p16;
			p16 = (uint16_t *)&zeros_subInventory[(i-1) * zeros_wordsPerSubBin];
			for(unsigned int j = 1; j < zeros_inventory[i] - zeros_inventory[i-1]-i/zeros_numBins; j++)
			{
				bitCount+=(rank0(zeros_inventory[i-1]+j+1)-rank0(zeros_inventory[i-1]+j));
				p16[subInventorySize-1] = (bitCount/(subInventorySize*subBinSize))*j;
				subInventorySize+=bitCount/(subInventorySize*subBinSize);
			}

			if(i/zeros_numBins)
			{
				bitCount++;
				p16[subInventorySize-1] = (bitCount/(subInventorySize*subBinSize))*((zeros_inventory[i]-1) - zeros_inventory[i-1]);
				subInventorySize+=bitCount/(subInventorySize*subBinSize);
			}
		}
	}

	data[0] = 0;
}

//computeBitPos(const uint64_t v, unsigned int r)
//Description: This function computes the position of a bit with rank r r:[1-64]
//Input: (a const uint64_t) this is the 64 bit containing the bit that we are looking for
//r: an int, the position of the bit
//Output:None
//Return: an int, the position of the bit
int Dictionary::computeBitPos( const uint64_t v,  unsigned int r ) const {

	uint64_t s, a, b, c, d, t;

	a =  v - ((v >> 1) & ~0UL/3);
	b = (a & ~0UL/5) + ((a >> 2) & ~0UL/5);
	c = (b + (b >> 4)) & ~0UL/0x11;
	d = (c + (c >> 8)) & ~0UL/0x101;
	t = (d >> 32) + (d >> 48);

	s  = 64;
	s -= ((t - r) & 256) >> 3; r -= (t & ((t - r) >> 8));
	t  = (d >> (s - 16)) & 0xff;
	s -= ((t - r) & 256) >> 4; r -= (t & ((t - r) >> 8));
	t  = (c >> (s - 8)) & 0xf;
	s -= ((t - r) & 256) >> 5; r -= (t & ((t - r) >> 8));
	t  = (b >> (s - 4)) & 0x7;
	s -= ((t - r) & 256) >> 6; r -= (t & ((t - r) >> 8));
	t  = (a >> (s - 2)) & 0x3;
	s -= ((t - r) & 256) >> 7; r -= (t & ((t - r) >> 8));
	t  = (v >> (s - 1)) & 0x1;
	s -= ((t - r) & 256) >> 8;
	s = 64 - s;
	return s;
}

//rank1(const unsigned long long int n)
//Description: This function returns the ones rank of a bit inside of the dictionary
//rank(p) = sum (b < p) : b = 1;
//Input:n (an unsigned long long int)
//Return:A unsigned long long int, the rank of the nth bit
unsigned long long int Dictionary::rank1(const unsigned long long int n) const {


	uint64_t x = bounds[2*(n/512) + offset];
	uint64_t s = bounds[2*(n/512) + offset+1] >> ((n/64)%8-1 + (((n/64)%8-1) >> 60 & 8)) * 9 & 0x1FF;

	uint64_t l = 0xffffffffffffffff;
	l = l << (63 - n % 64); //add Paren Here arond 63 to 64
	l = l << 1;

	uint64_t p = data[n/64+data_offset] & l;
	p = p - ((p >> 1) & 0x5555555555555555);
	p = (p & 0x3333333333333333) + ((p >> 2) & 0x3333333333333333);
	p = (p + (p >> 4)) & 0x0f0f0f0f0f0f0f0f;
	p = (p * 0x0101010101010101) >> 56;

	return(x+s+p);
}

//rank1(const unsigned long long int n)
//Description: This function returns the zero rank of a bit inside of the dictionary
//rank(p) = sum (b < p) : b = 1;
//Input:n (an unsigned long long int)
//Return:A unsigned long long int, the zero rank of the nth bit
unsigned long long int Dictionary::rank0(const unsigned long long int n) const {
	return n - rank1(n);
}

//Select1(const unsigned long long)
//Description: This function returns the first bit with rank n
//This function is a guided binary search, there are better ways to implement this
//I can look into this later
//Input: n, the rank of the bit I am trying to find
//Return: unsigned long long int, the position of the bit ranked n
unsigned long long int Dictionary::select1(const unsigned long long int n) const {

	//Upper Bounds
	unsigned long long int low = ones_inventory[n/ones_bin];
	unsigned long long int high = ones_inventory[n/ones_bin + 1];

	if((high-low)/(1ULL << 16))
	{

		if(n%ones_bin == 0)
			return low;

		uint64_t * p64;
		p64 = (uint64_t *)&ones_subInventory[(n/ones_bin) * ones_wordsPerSubBin];
		low+=ones_spill[p64[0]+n%ones_bin-1];

	}
	else{
		uint16_t * p16;
		p16 = (uint16_t *)&ones_subInventory[(n/ones_bin) * ones_wordsPerSubBin];
		uint64_t ones_subBinSize = (ones_bin)/(4*ones_wordsPerSubBin);
		low+=(p16[(n%ones_bin)/ones_subBinSize-cArr[(n%ones_bin)/ones_subBinSize]] * cArr[(n%ones_bin)/ones_subBinSize]);
		long long int size = low/64;
		long long unsigned int currSum = rank1(size*64);

		for(; currSum < n+1; size++)
		{
			currSum+=sidewaysAddition(data[size+data_offset]);
		}


		size--;
		low = 64 * size + computeBitPos(data[size+data_offset], n - rank1(size*64)+1);

	}
	return low;
}


//Select0(const unsigned long long)
//Description: This function returns the first bit with rank n (zeros)
//Simple search
//Input: n, the rank of the bit I am trying to find
//Return: unsigned long long int, the position of the bit ranked n
unsigned long long int Dictionary::select0(const unsigned long long int n) const{

	//Upper Bounds
	unsigned long long int low = zeros_inventory[n/zeros_bin];
	unsigned long long int high = zeros_inventory[n/zeros_bin + 1];

	if((high-low)/(1ULL << 16))
	{
		uint64_t * p64;
		p64 = (uint64_t *)&zeros_subInventory[(n/zeros_bin) * zeros_wordsPerSubBin];
		low+=(zeros_spill[p64[0]+n%zeros_bin-(1-((n%zeros_bin-1) >> 60 & 8)/8)] * (1-((n%zeros_bin-1) >> 60 & 8)/8));
	}
	else{

		uint16_t * p16;
		p16 = (uint16_t *)&zeros_subInventory[(n/zeros_bin) * zeros_wordsPerSubBin];
		uint64_t zeros_subBinSize = (zeros_bin)/(4*zeros_wordsPerSubBin);
		low+=(p16[(n%zeros_bin)/zeros_subBinSize-cArr[n%zeros_bin/zeros_subBinSize]] * cArr[n%zeros_bin/zeros_subBinSize]);
		long long int size = low/64;
		long long unsigned int currSum = rank0(size*64);
		for(; currSum < n+1; size++)
		{
			currSum+=sidewaysAddition(~data[size+data_offset]);
		}
		size--;
		low = 64 * size + computeBitPos(~data[size+data_offset], n - rank0(size*64)+1);
	}

	return low;
}
