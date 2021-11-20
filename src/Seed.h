/*
 * Seed.h
 *
 *  Created on: Jun 10, 2011
 *      Author: Julia
 */

#ifndef SEED_H_
#define SEED_H_

class Seed {
	public:
		Seed(unsigned long long int, unsigned long long int, unsigned long long int, unsigned long long int);
		Seed(const Seed &);
		Seed();
		virtual ~Seed();
		Seed & operator = (const Seed &);
		bool operator < (const Seed &) const;
		int weight;
		unsigned long long int coord1X;
		unsigned long long int coord2X;
		unsigned long long int coord1Y;
		unsigned long long int coord2Y;
};

#endif /* SEED_H_ */
