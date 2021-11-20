/*
 * Labels.h
 *
 *      Author: Julia
 */

#ifndef LABELS_H_
#define LABELS_H_

class Labels {
public:
	Labels();
	virtual ~Labels();
	void reserve(const unsigned int, const unsigned long long int, const unsigned long long int);
	void clear(void);

	//Append Labels
	void append(const char * &);
	void append(const char * );
	void append(const char *, const int);

	//Accession
	string labelAt(const unsigned long long int);
	unsigned int getTagLen(void) const;
	unsigned long long int getNumLabels(void) const;
	unsigned long long int getOffset(void) const;
	
	//Read/Write to file
	void write(FILE * &);
	void read(FILE * &);

private:
	unsigned int tagLen;
	unsigned long long int numLabels;
	unsigned long long int offset;
	char * data;
};

#endif /* LABELS_H_ */
