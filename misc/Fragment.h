/*
 * Fragment.h
 *
 *  Created on: Jun 5, 2011
 *      Author: Julia
 */

#ifndef FRAGMENT_H_
#define FRAGMENT_H_

class Fragment {
public:
	Fragment();
	Fragment(const char *, const char *);
	Fragment(const char *);
	Fragment(const char * &, const char * &);
	Fragment(const char * &);
	Fragment(const Fragment &);
	Fragment(const unsigned char * , const int);
	Fragment & operator = (const Fragment &);
	bool operator < (const Fragment &) const;
	bool operator > (const Fragment &) const;
	bool operator == (const Fragment &) const;
	virtual ~Fragment();
	char at(const int) const;
	char qual(const int) const;
	const unsigned char * c_str(void) const;
	int size(void) const;
private:
	void encode(const char &, const char &, const int);
	void decode(char &, char &, const int) const;
	unsigned char * str;
	int length;
};

#endif /* FRAGMENT_H_ */
