#include <iostream>
#include <stdint.h>
#include <cstring>
#include <cstdio>
using namespace std;
#include "MappingValues.h"
#include "Dictionary.h"
#include "Graph.h"

//(reverse << 15) | (qContained << 14) | (rContained << 13) | (matePairs << 12) | (edgeInGraph << 11)

//Graph()
//Description:Default Constructor
//Input: None
//Ouput:None
//Return:None
Graph::Graph(void)
:Edges ('\0'), Nodes('\0'), DestNodes('\0'), EdgeWeights('\0'), EdgeFlags('\0'), NodeFlags('\0'), numEdges(0), numNodes(0), numEdgeSets(0), offset(0), maskedEdges(0), maskedNodes(0), isSet(false), heavyEdge(0){;;;}

//~Graph()
//Description:Default Destructor
//Input: None
//Ouput:None
//Return:None
Graph::~Graph(void)
{
	delete [] DestNodes;
	delete [] EdgeWeights;
	delete [] EdgeFlags;
	delete [] NodeFlags;

	if(!isSet)
	{
		delete [] Nodes;
		delete [] Edges;
	}
}

//setOffset(long long int)
//Description:This sets the offset
//Input: o, the offset
//Ouput:None
//Return:None
void Graph::setOffset(long long int o) {
	offset = o;
}

//setOffset(void)
//Description:This sets the offset
//Input: None
//Ouput:None
//Return:None
long long int Graph::getOffset(void) const {
	return offset;
}

//clear()
//Description:Default clears the graph
//Input: None
//Ouput:None
//Return:None
void Graph::clear(void)
{
	delete [] DestNodes;
	delete [] EdgeWeights;
	delete [] EdgeFlags;
	delete [] NodeFlags;

	if(!isSet)
	{
		delete [] Nodes;
		delete [] Edges;
	}

	NodeIndex.clear();
	EdgeIndex.clear();

	DestNodes = '\0'; EdgeWeights = '\0';
	EdgeFlags = '\0'; NodeFlags = '\0';
	Nodes = '\0'; Edges = '\0';

	numEdges = 0; numNodes = 0; numEdgeSets = 0;
	maskedNodes = 0; maskedEdges = 0;
	isSet = false; heavyEdge = 0;
}

//Reserving and Adding to the Graph

//reserveNumEdges(const long long int n)
//Description: This function reserves edges in the graph
//Input: n, a const long long int (the number of edges)
//Ouput:None
//Return:None
void Graph::reserveNumEdges(const long long int n)
{
	Edges = new uint64_t [(n+64)/64+1];
	DestNodes = new uint32_t [n];
	EdgeWeights = new uint16_t[n];
	EdgeFlags = new uint64_t [(n+63)/64 * 6]; //fix this

	memset(Edges, 0, ((n+64)/64+1) * sizeof(uint64_t));
	memset(DestNodes, 0, n * sizeof(uint32_t));
	memset(EdgeWeights, 0, n * sizeof(uint16_t));
	memset(EdgeFlags, 0, ((n+63)/64 * 6) * sizeof(uint64_t));
}

//reserveNumEdges(const long long int n)
//Description: This function reserves nodes in the graph
//Input: n, a const long long int (the number of edges)
//Ouput:None
//Return:None
void Graph::reserveNumNodes(const long long int n){


	Nodes = new uint64_t[(n+63)/64+1];
	NodeFlags = new uint64_t [(n+63)/64];

	memset(Nodes, 0, ((n+63)/64+1) * sizeof(uint64_t));
	memset(NodeFlags, 0, ((n+63)/64) * sizeof(uint64_t));
}

//void set(void)
//Description: This function sets the graph after all of the edges have been added
//Input:None
//Output:None
//Return:None
void Graph::set(void) {

	long long int indexPos = numEdges/64 + 1;

	if(Edges != '\0')
	{
		Edges[indexPos] = (Edges[indexPos] | appValL >> numEdges%64);
		EdgeIndex.setBits(Edges, numEdgeSets+1, numEdges+1);
		EdgeIndex.setRanks(); EdgeIndex.setSelect1();
	}

	NodeIndex.setBits(Nodes, numEdgeSets, numNodes);
	NodeIndex.setRanks(); NodeIndex.setSelect1();

	isSet = true;
}

//void setHeavyGraph(void) const
//Description: This allows the graph to store one edge weight
//in 16 bits (I use this for the multilevel merging)
//Uses a flag switch
//Input:None
//Output:None
//Return:None
void Graph::setHeavyGraph(void) {
	heavyEdge = 1;
}

//addSilentNodes(const int n)
//Description: This function adds a node to the graph
//Input:n, an int, the number of "silent" 0 degree nodes
//Output:None
//Return:None
void Graph::addSilentNodes(const int n){
	numNodes+=n;
}

//newEdgeSet(void)
//Description: This function adds a new edge set to the graph
//Input:none
//Output:None
//Return:None
void Graph::newEdgeSet(void) {

	long long int indexPos = numEdges/64 + 1;
	Edges[indexPos] = (Edges[indexPos] | appValL >> numEdges%64);
	numEdgeSets++;

	indexPos = numNodes/64+1;
	Nodes[indexPos] = (Nodes[indexPos] | appValL >> numNodes%64);
	numNodes++;
}

//addEdge(const uint32_t, const uint32_t)
//Description: This function add an edge to the graph
//Input:destNode, a uint32_t (the destination node of the edge)
//Input:edgeInfo, a uint32 (contains the edgeOvl
//Output:None
//Return:None
void Graph::addEdge(const uint32_t destNode, const uint32_t edgeInfo) {

	DestNodes[numEdges] = destNode;
	EdgeWeights[numEdges] = edgeInfo >> 16;

	uint64_t eFlags = edgeInfo << 16;
	eFlags = eFlags << 32;	

	long long int indexPos = (numEdges * 6)/64;
	EdgeFlags[indexPos] = EdgeFlags[indexPos] | eFlags >> ((numEdges * 6)%64);

	eFlags = eFlags >> 59;
	indexPos = (numEdges * 6 + 4)/64;
	EdgeFlags[indexPos] = EdgeFlags[indexPos] | eFlags << (63 - ((numEdges * 6 + 4)%64));
	numEdges++;
}

//Edge Accession by Edge index

//edgeOvlLen(const long long int) const;
//Description: This function accesses an edge by its index
//Input: n, a const long long int, the edge index
//Output:None
//Return: int, the overlap length of the fragments
//Return: represented by the end nodes of the edge
int Graph::edgeOvlLen(const long long int n) const {

	int ovlInterval = 10;
	uint16_t edgeInfo = EdgeWeights[n];
	return ovlInterval * (edgeInfo >> (6-heavyEdge*6));
}

//edgeOvlIden(const long long int) const;
//Description: This function accesses an edge by its index
//Input: n, a const long long int, the edge index
//Output:None
//Return: int, the overlap length of the fragments
//Return: represented by the end nodes of the edge
float Graph::edgeOvlIden(const long long int n) const {

	float idenInterval = .5; float idenLowerBound = 70;
	uint16_t edgeInfo = EdgeWeights[n];
	edgeInfo = edgeInfo << 10;
	
	float score = edgeInfo >> 10;	
	return idenInterval * score + idenLowerBound;
}

//edgeInGraph(const long long int) const;
//Description: This function accesses an edge by its index
//Input: n, a const long long int, the edge index
//Output:None
//Return: int, the overlap length of the fragments
//Return: represented by the end nodes of the edge
bool Graph::edgeInGraph(const long long int n) const {

	long long int indexPos = (n*6+5)/64;
	uint64_t eFlags = EdgeFlags[indexPos];

	eFlags = eFlags << (n*6+5)%64;
	return !(eFlags >> 63);
}

//isReversed(const long long int) const;
//Description: This function accesses an edge by its index
//Input: n, a const long long int, the edge index
//Output:None
//Return: bool, if the edge is a reverse edge
bool Graph::isReverse(const long long int n) const {

	long long int indexPos = (n*6)/64;
	uint64_t eFlags = EdgeFlags[indexPos];

	eFlags = eFlags << (n*6)%64;
	return eFlags >> 63;
}

//qContained(const long long int) const;
//Description: This function accesses an edge by its index
//Input: n, a const long long int, the edge index
//Output:None
//Return: bool, if the query is contained in the reference, query = source, reference = dest
bool Graph::sContained(const long long int n) const {

	long long int indexPos = (n*6 + 1)/64;
	uint64_t eFlags = EdgeFlags[indexPos];

	eFlags = eFlags << (n*6 + 1)%64;
	return eFlags >> 63;
}

//rConstrained(const long long int) const;
//Description: This function accesses an edge by its index
//Input: n, a const long long int, the edge index
//Output:None
//Return: bool, if the reference is contained in the query, query = source, reference = dest
bool Graph::dContained(const long long int n) const {
	long long int indexPos = (n*6 + 2)/64;
	uint64_t eFlags = EdgeFlags[indexPos];

	eFlags = eFlags << (n*6 + 2)%64;
	return eFlags >> 63;
}

//rDovetail(const long long int) const;
//Description: This function accesses an edge by its index
//Input: n, a const long long int, the edge index
//Output:None
//Return: bool, if the source fragment is to the right of the dest fragment
bool Graph::rDovetail(const long long int n) const {

	long long int indexPos = (n*6 + 3)/64;
	uint64_t eFlags = EdgeFlags[indexPos];

	eFlags = eFlags << (n*6 + 3)%64;
	return eFlags >> 63;
}

//rConstrained(const long long int) const;
//Description: This function accesses an edge by its index
//Input: n, a const long long int, the edge index
//Output:None
//Return: bool, if the edge represents a mate pair
bool Graph::matePair(const long long int n) const {

	long long int indexPos = (n*6 + 4)/64;
	uint64_t eFlags = EdgeFlags[indexPos];

	eFlags = eFlags << (n*6 + 4)%64;
	return eFlags >> 63;
}

//uint32_t edgeSource(const long long int) const;
//Description: This function accesses an edge by its index
//Input: n, a const long long int, the edge index
//Output:None
//Return: uint32_t, the label of the source node of the edge
uint32_t Graph::edgeSource(const long long int n) const {
	return NodeIndex.select1(EdgeIndex.rank1(n+1) - 1);
}

//uint32_t edgeDest(const long long int) const;
//Description: This function accesses an edge by its index
//Input: n, a const long long int, the edge index
//Output:None
//Return: uint32_t, the label of the destinatin node of the edge
uint32_t Graph::edgeDest(const long long int n) const {
	return DestNodes[n];
}

//Edge Accession by Node indexes

//edgeOvlLen(const long long int) const;
//Description: Uses Binary search to access an edge by its node labels
//Input: source, uint32_t: the source node, dest, uint32_t: the dest node
//Output:None
//Return: int, the overlap length of the fragments
//Return: represented by the end nodes of the edge
int Graph::edgeOvlLen(const uint32_t source, const int edgeRank) const {

	long long int nodeLoc = EdgeIndex.select1(NodeIndex.rank1(source));
	return edgeOvlLen(nodeLoc+edgeRank);
}

//edgeIdenLen(const long long int) const;
//Description: Uses Binary search to access an edge by its node labels
//Input: source, uint32_t: the source node, dest, uint32_t: the dest node
//Output:None
//Return: float, the overlap identity of the fragments
//Return: represented by the end nodes of the edge
float Graph::edgeOvlIden(const uint32_t source, const int edgeRank) const {

	long long int nodeLoc = EdgeIndex.select1(NodeIndex.rank1(source));
	return edgeOvlIden(nodeLoc+edgeRank);
}

//edgeInGraph(const long long int) const;
//Description: Uses Binary search to access an edge by its node labels
//Input: source, uint32_t: the source node, dest, uint32_t: the dest node
//Output:None
//Return: bool, if the edge is in the graph
bool Graph::edgeInGraph(const uint32_t source, const int edgeRank) const {

	long long int nodeLoc = EdgeIndex.select1(NodeIndex.rank1(source));
	return edgeInGraph(nodeLoc+edgeRank);
}

//isReverse(const long long int) const;
//Description: Uses Binary search to access an edge by its node labels
//Input: source, uint32_t: the source node, dest, uint32_t: the dest node
//Output:None
//Return: bool, if the edge is a reverse edge
bool Graph::isReverse(const uint32_t source, const int edgeRank) const {

	long long int nodeLoc = EdgeIndex.select1(NodeIndex.rank1(source));
	return isReverse(nodeLoc+edgeRank);
}

//edgeInGraph(const long long int) const;
//Description: Uses Binary search to access an edge by its node labels
//Input: source, uint32_t: the source node, dest, uint32_t: the dest node
//Output:None
//Return: bool, if the query sequence is contained
bool Graph::sContained(const uint32_t source, const int edgeRank) const {

	long long int nodeLoc = EdgeIndex.select1(NodeIndex.rank1(source));
	return sContained(nodeLoc+edgeRank);
}

//rContained(const long long int) const;
//Description: Uses Binary search to access an edge by its node labels
//Input: source, uint32_t: the source node, dest, uint32_t: the dest node
//Output:None
//Return: bool, if the query sequence is contained
bool Graph::dContained(const uint32_t source, const int edgeRank) const {

	long long int nodeLoc = EdgeIndex.select1(NodeIndex.rank1(source));
	return dContained(nodeLoc+edgeRank);
}

//rDovetail(const long long int) const;
//Description: This function accesses an edge by its node labels
//Input: n, a const long long int, the edge index
//Output:None
//Return: bool, if the source fragment is to the right of the dest fragment
bool Graph::rDovetail(const uint32_t source, const int edgeRank) const {

	long long int nodeLoc = EdgeIndex.select1(NodeIndex.rank1(source));
	return rDovetail(nodeLoc+edgeRank);
}

//matePair(const long long int) const;
//Description: Uses Binary search to access an edge by its node labels
//Input: source, uint32_t: the source node, dest, uint32_t: the dest node
//Output:None
//Return: bool, if the edge is a matepair
bool Graph::matePair(const uint32_t source, const int edgeRank) const {

	long long int nodeLoc = EdgeIndex.select1(NodeIndex.rank1(source));
	return matePair(nodeLoc+edgeRank);
}

//edgeSource(const long long int) const;
//Description: Uses Binary search to access an edge by its node labels
//Input: source, uint32_t: the source node, dest, uint32_t: the dest node
//Output:None
//Return: uint32_t the source label
uint32_t Graph::edgeSource(const uint32_t source, const int edgeRank) const {

	long long int nodeLoc = EdgeIndex.select1(NodeIndex.rank1(source));
	return edgeSource(nodeLoc + edgeRank);
}

//edgeDest(const long long int) const;
//Description: Uses Binary search to access an edge by its node labels
//Input: source, uint32_t: the source node, dest, uint32_t: the dest node
//Output:None
//Return: uint32_t the dest label
uint32_t Graph::edgeDest(const uint32_t source, const int edgeRank) const {

	long long int nodeLoc = EdgeIndex.select1(NodeIndex.rank1(source));
	return edgeDest(nodeLoc + edgeRank);
}

//Set Function for Edges

//setEdgeInGraph(const long long int);
//Description: This sets an edge into the graph
//Input: n, a const long long int, the edge Index
//Output:None
//Return: None
void Graph::setEdgeInGraph(const long long int n) {

	long long int indexPos = (n*6 + 5)/64;
	EdgeFlags[indexPos] = EdgeFlags[indexPos] & ~(appValL >> (n* 6 + 5)%64);
}


//setEdgeInGraph(const uint32_t, const uint32_t);
//Description: This sets an edge into the graph
//Input: source, a uint32_t, dest, a uint32_t
//Output:None
//Return: None
void Graph::setEdgeInGraph(const uint32_t source, const int edgeRank) {

		long long int nodeLoc = EdgeIndex.select1(NodeIndex.rank1(source));
		setEdgeInGraph(nodeLoc + edgeRank);
}

//maskEdgeFromGraph(long long int);
//Description: This masks an edge in the graph
//Input: n, a long long int, the index of the edge
//Output:None
//Return:None
void Graph::maskEdgeFromGraph(long long int n) {

	long long int indexPos = (n*6 + 5)/64;
	EdgeFlags[indexPos] = EdgeFlags[indexPos] | (appValL >> (n*6 + 5)%64);
}

//maskEdgeFromGraph(long long int)
//Description: This masks an edge in the graph
//Input:source, a uint32_t (the source node), dest, a uint32_t (the dest node)
//Output:None
//Return:None
void Graph::maskEdgeFromGraph(uint32_t source, const int edgeRank) {

	long long int nodeLoc = EdgeIndex.select1(NodeIndex.rank1(source));
	maskEdgeFromGraph(nodeLoc + edgeRank);
}

//int nodeDegree (uint32_t)
//Description: This returns the degree of a node
//Input:n, a uint32_t, the node index
//Output:None
//Return: int, the node's degree
int Graph::nodeDegree(const uint32_t n) const {
   
	if(NodeIndex.isSet(n))
		return (EdgeIndex.select1(NodeIndex.rank1(n)+1) - EdgeIndex.select1(NodeIndex.rank1(n)));
	else
		return 0;
}

//bool nodeInGraph(uint32_t n)
//Description: This returns true if the node is in the graph
//Input:n, a uint32_t, the node index
//Output:None
//Return:true if the node is in the graph
bool Graph::nodeInGraph(const uint32_t n) const {
	long long int indexPos = n/64;
	return !(( appValL >> n%64 ) & NodeFlags[indexPos]);
}

//void setNodeInGraph(uint32_t n)
//Description: This returns true if the node is in the graph
//Input:n, a uint32_t, the node index
//Output:None
//Return:None
void Graph::setNodeInGraph(const uint32_t n) {
	long long int indexPos = n/64;
	NodeFlags[indexPos] = NodeFlags[indexPos] & ~(appValL >> n%64);
}

//void setNodeInGraph(uint32_t n)
//Description: This returns true if the node is in the graph
//Input:n, a uint32_t, the node index
//Output:None
//Return:None
void Graph::maskNodeFromGraph(const uint32_t n) {
	long long int indexPos = n/64;
	NodeFlags[indexPos] = NodeFlags[indexPos] | (appValL >> n%64);
}

//Graph information

//long long int numNodes(void) const
//Description: Returns the number of nodes in a graph
//Input:None
//Output:None
//Return:A long long int, the number of nodes in the graph
long long int Graph::getNumNodes(void) const {
	return numNodes;
}

//long long int numEdges(void) const
//Description: Returns the number of edges in a graph
//Input:None
//Output:None
//Return:A long long int, the number of edges in the graph
long long int Graph::getNumEdges(void) const {
	return numEdges;
}

//void write(FILE * & pFile)
//Description: This function writes the graph to a file
//Input: pFile, a FILE pointer
//Output: None
//Return: None
void Graph::write(FILE * & pFile) {

	fwrite(&numEdges, sizeof(long long int), 1, pFile);
	fwrite(&numNodes, sizeof(long long int), 1, pFile);
	fwrite(&numEdgeSets, sizeof(long long int), 1, pFile);
	fwrite(&maskedEdges, sizeof(long long int), 1, pFile);
	fwrite(&maskedNodes, sizeof(long long int), 1, pFile);
	fwrite(&offset, sizeof(long long int), 1, pFile);
	fwrite(&heavyEdge, sizeof(int), 1, pFile);

	fwrite(Edges, sizeof(uint64_t), (numEdges+64)/64+1, pFile);
	fwrite(Nodes, sizeof(uint64_t), (numNodes+63)/64+1, pFile);
	fwrite(DestNodes, sizeof(uint32_t), numEdges, pFile);
	fwrite(EdgeWeights, sizeof(uint16_t), numEdges, pFile);
	fwrite(EdgeFlags, sizeof(uint64_t), (numEdges+63)/64 * 6, pFile);
	fwrite(NodeFlags, sizeof(uint64_t), (numNodes+63)/64, pFile);
}

//void write(FILE * & pFile)
//Description: This function writes the fragment index to a file
//Input: pFile, a FILE pointer
//Output: None
//Return: None
void Graph::read(FILE * & pFile) {

	clear();

	fread(&numEdges, sizeof(long long int), 1, pFile);
	fread(&numNodes, sizeof(long long int), 1, pFile);
	fread(&numEdgeSets, sizeof(long long int), 1, pFile);
	fread(&maskedEdges, sizeof(long long int), 1, pFile);
	fread(&maskedNodes, sizeof(long long int), 1, pFile);
	fread(&offset, sizeof(long long int), 1, pFile);
	fread(&heavyEdge, sizeof(int), 1, pFile);

	reserveNumEdges(numEdges);
	reserveNumNodes(numNodes);

	fread(Edges, sizeof(uint64_t), (numEdges+64)/64+1, pFile);
	fread(Nodes, sizeof(uint64_t), (numNodes+63)/64+1, pFile);
	fread(DestNodes, sizeof(uint32_t), numEdges, pFile);
	fread(EdgeWeights, sizeof(uint16_t), numEdges, pFile);
	fread(EdgeFlags, sizeof(uint64_t), (numEdges+63)/64 * 6, pFile);
	fread(NodeFlags, sizeof(uint64_t), (numNodes+63)/64, pFile);

	set();
}
