/*
 * Sort.h
 *
 *      Author: Julia
 */

#ifndef SORTEDGES_H_
#define SORTEDGES_H_

/*
 * This class sorts a directory with files of edges.
 * Edges are sorted by node labels or by weights
 */


class SortEdges{

	public:
		SortEdges(const int);
		void sort_edges(const char [], const char []);
	private:
		bool fileExists(const char []);
		int sortType; //0 for sorting by Node Labels //1 for Sorting by Edge Weight
};

#endif /* SORT_H_ */
