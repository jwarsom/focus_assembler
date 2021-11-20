/*
 * Fragment_Index.h
 *
 *      Author: Julia
 */

#ifndef FRAGMENT_INDEX_H_
#define FRAGMENT_INDEX_H_

class Fragment_Index {
public:

	//Dummy initializer/destructor functions
	Fragment_Index();
	virtual ~Fragment_Index(void);

	//accession functions
	char at(const unsigned long long int, const unsigned int reverse = 0) const;
	int qualAt(unsigned long long int, const unsigned int reverse = 0) const;
	uint64_t ithWord(const unsigned long long int, unsigned int const, const unsigned int reverse = 0);
	pair<unsigned long long int, unsigned long long int> indexBounds(const unsigned long long int, const unsigned int reverse = 0);
	unsigned long long int positionToIndex(const unsigned long long int, const unsigned int reverse = 0);
	unsigned long long int length(void) const;
	unsigned long long int numFragments(void) const;
	unsigned long long int getIndexOffset(void) const;
	unsigned long long int getFragmentOffset(void) const;

	//appending functions
	void append(const char * &, const char * &);
	void append(const char *, const char *);
	void append(const char * &, const char * &, unsigned int);
	void append(const char *, const char *, unsigned int);

	//Internal settings
	void reserve(const unsigned long long int, const unsigned long long int, const unsigned long long int, const unsigned long long int);
	void clear(void);
	void set(void);

	//File reading and writing
	void read(FILE * &);
	void write(FILE * &);

private:

	//encode stuff
	void encode(const char &, const char &, const int);
	void decode(char &, char &, const int) const;

	//The major data structures

	//The bounds
	Dictionary bounds;
	uint64_t * boundsData;

	//Major data structure
	unsigned char * fragments;

	//length tracking values
	unsigned long long int indexOffset; //This tracks the number of fragments in previous datasets
 	unsigned long long int fragmentOffset; //This tracks the lengths of the fragments in previous datasets
	unsigned long int totalFragments;
	unsigned long long int indexLength;
};

#endif /* FRAGMENT_INDEX_H_ */
