/*
 * Dictionary.h
 *
 *      Author: Julia
 */

#ifndef DICTIONARY_H_
#define DICTIONARY_H_

class Dictionary {
public:
	Dictionary();
	virtual ~Dictionary();
	Dictionary & operator = (const Dictionary &);
	void clear(void);
	void setBits(uint64_t * &, const long long int, const long long int);
	void setRanks(void);
	void setSelect1(void);
	void setSelect0(void);
	unsigned long long int rank1(const unsigned long long int) const;
	unsigned long long int rank0(const unsigned long long int) const;
	unsigned long long int select1(const unsigned long long int) const;
	unsigned long long int select0(const unsigned long long int) const;
	int computeBitPos(const uint64_t, unsigned int) const;
	bool isSet(const unsigned long long int) const;
	//All the accession functions
	long long int getCapacity(void) const;
	long long int getNumMembers(void) const;
	long long int getUniverseSize(void) const;
	uint64_t * getData(void) const;


	inline uint64_t sidewaysAddition(const uint64_t x) const
	{
		uint64_t s = x - ((x >> 1) & 0x5555555555555555);
		s = (s & 0x3333333333333333) + ((s >> 2) & 0x3333333333333333);
		s = (s + (s >> 4)) & 0x0f0f0f0f0f0f0f0f;
		s = (s * 0x0101010101010101) >> 56;
		return s;
	}

private:
	inline uint64_t cmp1(uint64_t x, uint64_t y){
		return ((((y | 0x8080808080808080) - (x & ~0x8080808080808080)) | (x ^ y) ) ^ (x & ~y)) & 0x8080808080808080;
	}

	inline uint64_t cmp2(uint64_t x, uint64_t y){
		return (((y | 0x8080808080808080) - (x & ~0x8080808080808080)) ^ x ^ y ) & 0x8080808080808080;
	}

	inline uint64_t cmp3(uint64_t x, uint64_t y){
		return ((((x | 0x8080808080808080) - 0x01010101010101010) | x)) & 0x8080808080808080;
	}

	uint64_t * data;
	uint64_t * bounds;
	long long int capacity;
	long long int numMembers;
	long long int universeSize;

	//Ones inventory
	uint64_t * ones_inventory;
	uint64_t * ones_subInventory;
	uint64_t * ones_spill;
	long long int ones_bin;
	long long int ones_numBins;
	long long int ones_wordsPerSubBin;
	long long int ones_numSubBins;
	long long int ones_numSpill;

	//Zeroes inventory
	uint64_t * zeros_inventory;
	uint64_t * zeros_subInventory;
	uint64_t * zeros_spill;
	long long int zeros_bin;
	long long int zeros_numBins;
	long long int zeros_wordsPerSubBin;
	long long int zeros_numSubBins;
	long long int zeros_numSpill;
};

#endif /* DICTIONARY_H_ */
