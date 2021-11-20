/*
 * Chain.cpp
 *
 *      Author: Julia
 */

#include<iostream>
#include<map>
#include<vector>
#include<climits>
using namespace std;
#include "Seed.h"
#include "Chain.h"


Chain::Chain()
:maxChain(-1), prevX(INT_MAX), prevY(INT_MAX){;;;}

Chain::~Chain() {

}

Chain::Chain (Chain & c)
{
	for(int i = 0; i < c.size(); i++)
	{
		Link link;
		link.seed = c.at(i);
		link.chainWeight = link.seed.weight;
		links.push_back(link);
	}
	maxChain = c.longestChain();
	prevX = c.getPrevX();	
	prevY = c.getPrevY();	
}

int Chain::getPrevX() const
{
	return prevX;
}

int Chain::getPrevY() const
{
	return prevY;
}

Chain & Chain::operator = (Chain & c)
{
	for(int i = 0;  i < c.size(); i++)
	{
		Link link;
		link.seed = c.at(i);
		link.chainWeight = link.seed.weight;
		links.push_back(link);
	}
	maxChain = c.longestChain();
	prevX = c.getPrevX();	
	prevY = c.getPrevY();	
	return * this;
}

Seed & Chain::at(const int index)
{
	return links.at(index).seed;
}

int Chain::predecessor(const int n) const
{
	return links.at(n).predecessor;
}

int Chain::weight(const int n ) const
{
	return links.at(n).chainWeight;
}

int Chain::size(void) const
{
	return links.size();
}

int Chain::longestChain(void) const
{
	return maxChain;
}

int Chain::chainWeight(void) const
{
	return links.at(maxChain).chainWeight;
}

void Chain::addSeed(const Seed & s1)
{
	 Link link;
     link.seed = s1;
     link.chainWeight = s1.weight;
     links.push_back(link);
}

void Chain::chainSeeds(void)
{
	int maxWeight = 0;
	multimap<int, int> X;
	int count = 0;
	for(vector<Link>::iterator it = links.begin(); it != links.end(); it++)
	{
		X.insert(pair<int, int>((*it).seed.coord1X, count));
		X.insert(pair<int, int>((*it).seed.coord2X, count++));
	}
	multimap<pair<int, int>, int > Y;
	for(multimap<int, int>::iterator it = X.begin(); it != X.end(); it++)
	{
		int index = (*it).second;
		unsigned int coord = (*it).first;
		if(coord == links.at(index).seed.coord1X)
		{
			if(!Y.empty())
			{
				multimap<pair<int, int>, int >::iterator found;
				pair<int, int> key(-1*links.at(index).seed.coord1Y, -1*links.at(index).chainWeight);
				found = Y.upper_bound(key);
				if(found != Y.end())
				{
					links.at(index).predecessor = (*found).second;
					links.at(index).chainWeight+=(-1*(*found).first.second);
				}
			}
		}else
		{
			if(!conflicts(index, Y))
			{
				pair<int, int> insert(-1*links.at(index).seed.coord2Y, -1*links.at(index).chainWeight);
				deleteConflicts(index, Y);
				Y.insert(pair<pair<int, int>, int >(insert, index));
				if(maxWeight <= links.at(index).chainWeight)
				{
					maxWeight = links.at(index).chainWeight;
					maxChain = index;
				}

			}
		}
	}
}

void Chain::deleteConflicts(const int index, multimap<pair<int, int>, int> & Y){

	multimap<pair<int, int>, int >::iterator found;
	pair<int, int> key(-1*links.at(index).seed.coord2Y, -1*links.at(index).chainWeight);
	found = Y.upper_bound(key);
	if(found != Y.begin())
	{
		found--;
		if(links.at((*found).second).chainWeight <= links.at(index).chainWeight)
			Y.erase(found);
	}
}

bool Chain::conflicts(const int index, multimap<pair<int, int>, int> & Y)
{
	int closest = 0;
	if(!Y.empty())
	{
		multimap<pair<int, int>, int >::iterator found;
		found =  Y.lower_bound(pair<int, int>(-1*links.at(index).seed.coord2Y, -1*links.at(index).chainWeight));
		if(found != Y.end())
		{
			closest = (*found).second;
			return links.at(closest).chainWeight >= links.at(index).chainWeight;
		}
	}
	return closest;
}
