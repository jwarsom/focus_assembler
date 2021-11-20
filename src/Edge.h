/*
 * Edge.h
 *
 *      Author: Julia
 */


#ifndef EDGE_H_
#define EDGE_H_


class Edge {
public:

	//Constructors
	Edge();
	Edge (const Edge &);
	Edge (uint32_t, uint32_t, int, int);

	//Comparison operators
	Edge & operator = (const Edge &);
	bool operator < (const Edge &) const;
	bool operator > (const Edge &) const;
	bool operator == (const Edge &) const;

	//Accession operators
	uint32_t getVertex1(void) const;
	uint32_t getVertex2(void) const;
	int getLength(void) const;
	int getIdentity(void) const;

	//Set functions
	void setEdgeValues(uint32_t, uint32_t, int, int);

	//Destructors
	virtual ~Edge();
private:

	//Private members
	uint32_t vertex1;
	uint32_t vertex2;
	int length;
	int iden;
};

#endif
