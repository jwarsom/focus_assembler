/*
 * Chain.h
 *
 *      Author: Julia
 */

#ifndef CHAIN_H_
#define CHAIN_H_

class Chain {
public:
	Chain();
	~Chain();
	Chain (Chain &);
	Chain & operator = (Chain &);
	void addSeed(const Seed &);
	void chainSeeds(void);
	Seed & at(const int);
	int longestChain(void) const;
	int chainWeight(void) const;
	int getPrevX() const;
	int getPrevY() const;
	int predecessor(const int) const;
	int weight(const int) const;
	int size(void) const;
	int MUM(void) const;
private:
	bool conflicts(const int, multimap<pair<int, int>, int > & );
	void deleteConflicts(const int, multimap<pair<int, int>, int> &);

	class Link{
	public:
		Link()
		:chainWeight(0), predecessor(-1){;;;}
		Link(const Link & l)
		{
			seed = l.seed;
			chainWeight = l.chainWeight;
			predecessor = l.predecessor;
		}
		Link & operator = (const Link & l)
		{
			seed = l.seed;
			chainWeight = l.chainWeight;
			predecessor = l.predecessor;
			return * this;
		}
		Seed seed;
		int chainWeight;
		int predecessor;
	};
	vector<Link> links;
	int maxChain;
	int prevX;
	int prevY;
};

#endif /* CHAIN_H_ */
