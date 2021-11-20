/*
 * Fragment.cpp
 *
 *  Created on: Jun 5, 2011
 *      Author: Julia
 */
#include<cstdio>
#include<stdint.h>
#include<math.h>
#include<cstring>
using namespace std;
#include "Fragment.h"
#define dummyQual 53
const char nucleotides[] = { 'A', 'C', 'T', 'G'};

/* Class Description */
// This Fragment class stores fragment info
// from next generation sequencing technologies
// (454, illumina, etc).
// It stores nucleotide values (A, T, C, G).
// Non-nucleotide values such as 'N' are assigned a random
// nucleotide
// It is able to store Phred quality values as well
// A nucleotide and quality value is encoded into
// a single char to save space.


//Fragment()
//Description:Default Constructor
//Input: None
//Ouput:None
//Return:None
Fragment::Fragment(void)
	:str('\0'), length(0) {;;;}

//Fragment(const char *, const char *)
//Description: Constructor -- encodes a DNA fragment and quals
//Input: dna (a char str), quals, (a char str)
//Output:None
//Return:None
Fragment::Fragment(const char * dna, const char * quals) {

	length = strlen(dna);
	str = new unsigned char[length];
	memset(str, 0, length);
	for(int i = 0; i < length; i++)
		encode(dna[i], quals[i], i);
}

//Fragment(const char *)
//Description: Constructor -- encodes a DNA fragment
//Input: dna (a char str)
//Output:None
//Return:None
Fragment::Fragment(const char * dna) {

	length = strlen(dna);
	str = new unsigned char[length];
	memset(str, 0, length);
	for(int i = 0; i < length; i++)
		encode(dna[i], dummyQual, i);
}

//Fragment(const char * &, const char * &)
//Description: Constructor -- encodes a DNA fragment and quals
//Input: dna (a char str ref), quals, (a char str ref)
//Output:None
//Return:None
Fragment::Fragment(const char * & dna, const char * & quals) {

	length = strlen(dna);
	str = new unsigned char[length];
	memset(str, 0, length);
	for(int i = 0; i < length; i++)
		encode(dna[i], quals[i], i);
}

//Fragment(const char * &)
//Description: Constructor -- encodes a DNA fragment
//Input: dna (a char str ref)
//Output:None
//Return:None
Fragment::Fragment(const char * & dna) {

	length = strlen(dna);
	str = new unsigned char[length];
	memset(str, 0, length);
	for(int i = 0; i < length; i++)
		encode(dna[i], dummyQual, i);
}

//Fragment(const Fragment & f)
//Description: Copy constructor
//Input: A Fragment object
//Output: None
//Return: None
Fragment::Fragment(const Fragment & f) {
	 length = f.size();
	 str = new unsigned char[length];
	 memcpy(str, f.c_str(), length);
}

//Fragment(const char * &, const int)
//Description: Constructor: This constructs a fragment from
//an encoded char string.
//Input: source (const char * &) The original encoded string
//len (an int): The number of characters to be copied
//Output: None
//Return:None
Fragment::Fragment(const unsigned char * source, const int len) {
	str = new unsigned char[len];
	length = len;
	memcpy(str, source, length);
}

//~Fragment(void)
//Decription: Destructor -- deletes str
//Input: None
//Output: None
//Return:None
Fragment::~Fragment(void) {

	delete [] str;
}

//Fragment::operator =
//Description: Assignment operator
//Input: A Fragment object
//Output: None
//Return: A Fragment object
Fragment & Fragment::operator = (const Fragment & f) {
	length = f.size();
	str = new unsigned char[length];
	memcpy(str, f.c_str(), length);
	return * this;
}

//Fragment::operator <
//Description: Less than operator
//Input: A Fragment object
//Output: None
//Return: Bool
bool Fragment::operator < (const Fragment & f) const {
	int i = 0;
	while(f.at(i) == at(i++) && i < length && i < f.size());
	if(f.at(i-1) == at(i-1))
	{
		return length < f.size();
	}
	return at(i-1) < f.at(i-1);
}

//Fragment::operator >
//Description: Greater than operator
//Input: A Fragment object
//Output: None
//Return: Bool
bool Fragment::operator > (const Fragment & f) const {
	int i = 0;
	while(f.at(i) == at(i++) && i < length && i < f.size());
	if(f.at(i-1) == at(i-1))
	{
		return length > f.size();
	}
	return at(i-1) > f.at(i-1);
}

//Fragment::operator ==
//Description: Equal than operator
//Input: A Fragment object
//Output: None
//Return: Bool
bool Fragment::operator == (const Fragment & f ) const {

	if(length != f.size())
		return false;

	int i = 0;
	while(f.at(i) == at(i++) && i < length);

	return i == length;
}

//at(const int) const
//Description: Returns the nucleotide at position i
//of the fragment
//Input: i (an int) -- The position of the nucleotide
//Output:None
//Return: A char -- the nucleotide
char Fragment::at(const int i) const {

	char nucleotide; char qual; // The quality value and nucleotide
	decode(nucleotide, qual, i);
	return nucleotide;
}

//qual(const int) const
//Description: Returns the quality value at position i
//of the fragment
//Input: i (an int) -- The position of the quality value
//Output:None
//Return: A char -- the quality value
char Fragment::qual(const int i) const {

	char nucleotide; char qual; // The quality value and nucleotide
	decode(nucleotide, qual, i);
	return qual;
}

//c_str(void) const
//Description: Returns encodes c_str
//Input: None
//Output:None
//Return: The encoded c_str
const unsigned char * Fragment::c_str(void) const {
	return str;
}

//size(void) const
//Description: Returns the size of the fragment
//Input: None
//Output:None
//Return: The length of the fragment
int Fragment::size(void) const {
	return length;
}

//encode(const char &, const char &, int) const
//Description: This function encodes a nucleotide and qual value into a
//single char and places it into position i in str
//Input: nucleotide (a char), qual (a char)
//Output:None
//Return:None
//Note: A useful equation for remembering Phred encoding --
//       Phred qual = qualityChar (what we are given) - 33;
//		 Since we are storing this quality value in 6 bits
//       We will subtract 33 from the original qualityChar
//Note: I show an example of encoding 'C' and quality 32 down below (binary values)
void Fragment::encode(const char & nucleotide, const char & qual, const int i){

	str[i] = nucleotide & 6;  //str[i] = 'C' ^ 6  = 01000011 ^ 00000110 = 00000010
	str[i] = str[i] >> 1;    // 00000011 >> 1 = 00000001
	str[i] = str[i] << 6;    // 00000001 << 6 = 01000000
	str[i] = str[i] | ((qual-33) & 63); // 01000000 | ( 00100000 (32) & 00111111 (63) ) = 01100000 = 01 ('C') & 100000 (32)
}

//dencode(char &, char &, int) const
//Description: This function decodes a compact char into a nucleotide and qual
//Input: nucleotide (a char), qual (a  char), i (an int)
//Output:None
//Return:None
//Note: A const global array is used to decode nucleotides
void Fragment::decode(char & nucleotide, char & qual, const int i) const{

	nucleotide = nucleotides[str[i] >> 6];
	qual = str[i] & 63;
}
