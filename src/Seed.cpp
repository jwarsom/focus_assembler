/*
 * Seed.cpp
 *
 *      Author: Julia
 */

#include "Seed.h"

//Seed(int, int, int, int, int)
//Description:This is the constructor for the class
//Input:x1 (an int) the position in the fragment, y1 (an int) the position in the genome or fragment, x2 (an int) y2 (an int)
//Output:None
//Return:None
Seed::Seed(unsigned long long int x1, unsigned long long int x2, unsigned long long int y1, unsigned long long int y2)
:weight(x2-x1+1), coord1X(x1), coord2X(x2), coord1Y(y1), coord2Y(y2){;;;}

//Seed(const Seed & s)
//Description: copy constructor
//Input:None
//Output:None
//Return:None
Seed::Seed(const Seed & s)
{
	weight = s.weight;
	coord1X = s.coord1X;
	coord1Y = s.coord1Y;
	coord2X = s.coord2X;
	coord2Y = s.coord2Y;
}

//Seed()
//Description: Default description
//Input:None
//Output:None
//Return:None
Seed::Seed()
: weight(0), coord1X(0), coord2X(0), coord1Y(0), coord2Y(0){;;;}

//~Seed()
//Class Destructor
//Input:None
//Output:None
//Return:None
Seed::~Seed()
{

}

//Seed & operator = (const Seed & )
//Description: This is the assignment operator
//Input: s (a Seed)
//Output:None
//Return:None
Seed & Seed::operator = (const Seed & s)
{
	weight = s.weight;
	coord1X = s.coord1X;
	coord1Y = s.coord1Y;
	coord2X = s.coord2X;
	coord2Y = s.coord2Y;
	return * this;
}

//bool & operator < (const Seed &)
//Description: This is a less than operator
//Input: Seed & s, the seed we are being compared to
//Output:None
//Return: bool, if the seed is smaller
bool Seed::operator < (const Seed & s) const
{
	if(coord2Y < s.coord2Y)
		return true;

	if(coord2Y == s.coord2Y)
	{
		if(coord2X < s.coord2X)
			return true;
	}

	return false;
}
