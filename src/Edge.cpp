/*
 * Edge.h
 *
 *      Author: Julia
 */

/* This class hold edge information. It is a simple class;
 * it only consists of three variables -- the two edge labels and the edge weight
 * as a int*/

#include <iostream>
using namespace std;
#include "Edge.h"


//Constructor
//Description: This is the generic class constructor
//Input:None
//Output:None
//Return:None
Edge::Edge()
:vertex1(0), vertex2(0), length(0), iden(0){;;;}

//Constructor
//Description: This is a constructor that sets the edge values
//Input:v1 an uint64_t, v2 an uint64_t, weight a int
//Output:None
//Return:None
Edge::Edge (uint32_t v1, uint32_t v2, int l, int i)
:vertex1(v1), vertex2(v2), length(l), iden(i) {;;;}

//Copy constructor
//Description: This is the copy constructor
//Input:An edge
//Output:None
//Return:None
Edge::Edge (const Edge & e){
	vertex1 = e.getVertex1();
	vertex2 = e.getVertex2();
	length = e.getLength();
	iden = e.getIdentity();
}

//void setEdgeValues(uint32_t, uint32_t, int)
//Description: This is the copy constructor
//Input:An edge
//Output:None
//Return:None
void Edge::setEdgeValues(uint32_t v1, uint32_t v2, int l, int i)
{
	vertex1 = v1; vertex2 = v2; length = l; iden = i;
}

//Operator =
//Description: This function copies another edges value
//Input:e, an Edge
//Output:None
//Return:None
Edge & Edge::operator = (const Edge & e){
	vertex1 = e.getVertex1();
	vertex2 = e.getVertex2();
	length = e.getLength();
	iden = e.getIdentity();
	return * this;
}

//Operator <
//Description: Less than
//Input: e, an Edge
//Output:None
//Return:bool
bool Edge::operator < (const Edge & e) const {
	return vertex1 < e.getVertex1();
}

//Operator >
//Description: Greater than
//Input: e, an Edge
//Output:None
//Return:bool
bool Edge::operator > (const Edge & e) const {
	return e < *this;
}

//Operator ==
//Description: Greater than
//Input: e, an Edge
//Output:None
//Return:bool
bool Edge::operator == (const Edge & e) const {
	return !(e < *this) && !(*this < e);
}

//getVertex1(void)
//Description: returns the value of the first vertex
//Input: none
//Output:none
//Return: a uint64_t, the value of vertex1
uint32_t Edge::getVertex1(void) const {
	return vertex1;
}

//getVertex2(void)
//Description: returns the value of the first vertex
//Input: none
//Output:none
//Return: a uint64_t, the value of vertex2
uint32_t Edge::getVertex2(void) const {
	return vertex2;
}

//getLength(void)
//Description: returns the ovl length
//Input:none
//Output:none
//Return: a int, the value of the weight
int Edge::getLength(void) const {
	return length;
}

//getIdentity(void)
//Description: returns the identity
//Input:none
//Output:none
//Return: a int, the value of the weight
int Edge::getIdentity(void) const {
	return iden;
}


//Destructors
//Description: This is a dummy destructor
//Return: None
Edge::~Edge() {
	//NOTHING TO DO
}
