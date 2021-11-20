/*
 * Merge.h
 *
 *      Author: Julia
 */

#ifndef MERGEEDGES_H_
#define MERGEEDGES_H_

/* This class merges two files */

//This compartor compares the sources and
//destinations of two edges
struct fileCmp
{
	public:
		fileCmp(uint32_t * & f)
			:F(f){;;;}

		bool operator() (const int & lhs, const int & rhs) const
		{
			if(F[30 + lhs] > F[30 + rhs])
				return true;
			if(F[30 + lhs] == F[30 + rhs])
				return F[40 + lhs] > F[40 + rhs];

			return false;
		}

	private:
		uint32_t * F;
};


/* This class merges two sets of sorted files
*/
class MergeEdges {

public:
	MergeEdges(const int);
	~MergeEdges();
	void merge(const int, const int, const char [], const char []);
private:
	void addEdge(uint32_t [], int &, uint32_t, uint32_t, uint32_t);
	void getNextEdge(priority_queue<int, vector<int>, fileCmp> &, const char [], uint32_t [], const int, const int, const int);
	FILE * initiateQueue(priority_queue<int, vector<int>, fileCmp> &, uint32_t [], const char [], const int, const int , const int );
	void readFile(uint32_t [], const char [], const int, const int, const int);
	void writeFile(uint32_t [], const int, const char [], const int);
	void appendFile(uint32_t [], const int, const char [], const int);
	bool fileExists(const char [], const int);

	FILE ** pFile; //The file pointers used to access files
	uint32_t * fileInfo; // This contains other information related to the files

	//0 - 9: The current file sets
	//10 - 19: The current position in the buffer
	//20 - 29: The number of elements read from a file
	//30 - 39: The source nodes of the current edge read from a file
	//40- 49: The dests nodes of the current edge read from a file
	//50 - 59: The weights of the current edge read from a file

	int mergeType; //Merging the weights = 1; Removing duplicates = 0

};

#endif /* MERGEEDGES_H_ */
