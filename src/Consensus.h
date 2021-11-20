/*
 * Consensus.h
 *
 *      Author: Julia
 */

#ifndef CONSENSUS_H_
#define CONSENSUS_H_

class Consensus {
public:
	Consensus();	
	virtual ~Consensus();
	bool getConsensus(uint32_t [], int &, int, int, int * [], int, Graph * [], int, Fragment_Index * [], int, char [], int &, int, char [], long long int [], int &, long long int &, long long int, int [], int [], int);
private:
	void addSequence(int [], int, char [], int &, char []);
	void getGraph(Graph * [], int, Graph * &, int);
	void getFragment_Index(Fragment_Index * [], int, Fragment_Index * &, int);

};

#endif /* DICTIONARY_H_ */
