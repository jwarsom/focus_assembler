/*
 * MergeFragments.h
 *
 *      Author: Julia
 */

#ifndef MERGEFRAGMENTS_H_
#define MERGEFRAGMENTS_H_

/* This class merges two files */


/* This class merges two sets of sorted files
 */
class MergeFragments {

	public:
		MergeFragments();
		~MergeFragments();
		void merge(int,  int,  char [],  char [],  char [] ,  char []);
	private:
		void writeFinalFragments(Fragment_Index &, Labels &, int, char * & , int, int, int, char []);
		void makeIndexes(char [], char [], char [], int , int);
		void writeBuff(char * &,  int,  char [],  int);
		void appendBuff(char * &,  int,  char [],  int);
		void appendFragment(char * &, int &, Fragment_Index &, pair<unsigned long long int, unsigned long long int> &, string &);
		void readData( char [],  char [], Fragment_Index &, Labels &, int);
		bool fileExists(const char []);

		bool fileExists( char []);
};

#endif /* MERGEFRAGMENTS_H_ */
