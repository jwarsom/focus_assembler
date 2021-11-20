/*
 * Graph.h
 *
 *      Author: Julia
 */

#ifndef GRAPH_H_
#define GRAPH_H_

/* This is the graph data structure
 * It is an edge list.
 * This graph is static, but it can mask and unmask edges
 * Nodes are indexed by a dictionary structure: see the dictionary class for more info
 * O(1) Node Access by Edge Index,
 * Log(node_degree) Edge Access by Node Labels
 * O(1) Edge Access by Index
 *
 * This graph is strict on input requirements
 * The input list must be sorted by Edge Labels
 */

class Graph {

public:

	Graph(void);
	~Graph(void);
	void clear(void);

	//Reserving and Adding to the Graph
	void reserveNumEdges(const long long int);
	void reserveNumNodes(const long long int);
	void addSilentNodes(const int);
	void newEdgeSet(void);
	void addEdge(const uint32_t, const uint32_t);
	void set(void);
	void setHeavyGraph(void);
	void setOffset(long long int);

	//Edge Accession by Edge index
	int edgeOvlLen(const long long int) const;
	float edgeOvlIden(const long long int) const;
	bool edgeInGraph(const long long int) const;
	bool isReverse(const long long int) const;
	bool sContained(const long long int) const;
	bool dContained(const long long int) const;
	bool rDovetail(const long long int) const;
	bool matePair(const long long int) const;
	uint32_t edgeSource(const long long int) const;
	uint32_t edgeDest(const long long int) const;

	//Edge Accession by Node index and rank of edge in list
	int edgeOvlLen(const uint32_t, const int) const;
	float edgeOvlIden(const uint32_t, const int) const;
	bool edgeInGraph(const uint32_t, const int) const;
	bool isReverse(const uint32_t, const int) const;
	bool sContained(const uint32_t, const int) const;
	bool dContained(const uint32_t, const int) const;
	bool rDovetail(const uint32_t, const int) const;
	bool matePair(const uint32_t, const int) const;
	uint32_t edgeSource(const uint32_t, const int) const;
	uint32_t edgeDest(const uint32_t, const int) const;

	//Set Function for Edges
	void setEdgeInGraph(const long long int); //By Edge index
	void setEdgeInGraph(const uint32_t, const int); //By Node index
	void maskEdgeFromGraph(const long long int); //By Edge index
	void maskEdgeFromGraph(const uint32_t, const int); //By Node index

	//Accession and Set functions for Nodes
	int nodeDegree(const uint32_t) const;
	bool nodeInGraph(const uint32_t) const;
	void setNodeInGraph(const uint32_t);
	void maskNodeFromGraph(const uint32_t);

	//Graph info
	long long int getNumNodes(void) const;
	long long int getNumEdges(void) const;

	long long int getOffset(void) const;

	//Reading and Writing
	void write(FILE * &);
	void read(FILE * &);

private:

	//This is a dictionary of all the nodes
	//Allows for Edge information to be retrieved from the edge list
	//O(1) Node Search O(1) Node Degree
	//Log(Node_Degree) Edge Retrieval
	Dictionary EdgeIndex;
	uint64_t * Edges;

	Dictionary NodeIndex; //This indexes all nodes including those with no edges
			      //Nodes with edges = 1 bit, Nodes with no edges  = 0 (singleton)

	uint64_t * Nodes;	 //The Node Information

	//This contains the destination nodes of an edge
	uint32_t * DestNodes;

	//This contains edge weights
	uint16_t * EdgeWeights;

	//Contain flags such as containment info, masking info, etc.
	uint64_t * EdgeFlags;
	uint64_t * NodeFlags;

	//Number of nodes and edges and offset if graph is a subset
	long long int numEdges;
	long long int numNodes;
	long long int numEdgeSets;
	long long int offset;

	//Number of masked nodes and edges
	long long int maskedEdges;
	long long int maskedNodes;

	bool isSet;
	int heavyEdge; // this determines whether or not the full 16 bit Edge Weight is used for the overlap weight
};

#endif /* GRAPH_H_ */
