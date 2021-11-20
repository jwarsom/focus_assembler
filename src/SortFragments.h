/*
 * SortFragments.h
 *
 *      Author: Julia
 */

#ifndef SORTFRAGMENTS_H_
#define SORTFRAGMENTS_H_

/*
 * This class sorts a fragment index according to 
 * an ordering obtained from the overlap graph.
 * 
 */


//This compartor compares the rank of two fragments
//in the reordered overlap graph

//This compartor compares rank of the fragments in the
//new overlap graph
struct fragCmp
{
	public:
		fragCmp(int * & n)
			:N(n){;;;}

		bool operator() (const unsigned long long & lhs, const unsigned long long & rhs) const
		{
			if(N[lhs] < N[rhs])
				return true;
			return false;
		}

	private:
		int * N; //the node map
};

class SortFragments{

	public:
		SortFragments();
		void sort(const char [], const char [], const char [], const char [],  int * &);
	private:
		bool fileExists(const char []);
};

#endif /* SORTFRAGMENTS_H_ */
