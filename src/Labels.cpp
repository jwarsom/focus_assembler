/*
 * Labels.cpp
 *
 *      Author: Julia
 */

#include <iostream>
#include<string>
#include <cstring>
#include <cstdio>
using namespace std;
#include "Labels.h"


/* This class holds the labels and descriptions of the reads */

//Labels()
//Description: This is the constructor for the labels class
//Input:None
//Output:None
//Return:None
Labels::Labels()
	:tagLen(0), numLabels(0), offset(0), data(nullptr){;;;}

	//~Labels()
	//Description: This is the destructor for the class
	//Input:None
	//Output:None
	//Return: None
	Labels::~Labels() {
		delete [] data;
	}

//void reserve(const unsigned int, const unsigned long long int, const unsigned long long int)
//Description: This function reserve space for the labels
//Input: t, const unsigned int -- the sequence tag length,
//n, const unsigned long long int -- the number of label
//Output:None
//Return: None
void Labels::reserve(const unsigned int t, const unsigned long long int n, const unsigned long long int o) {
	
	clear();

	tagLen = t; offset = o; 
	data = new char[n*(tagLen+1)];
}


//void clear()
//Description: This function clears the labels
//Input: None
//Output:None
//Return: None
void Labels::clear(void){
	delete [] data; 
	tagLen = 0; numLabels = 0; offset = 0; data = nullptr;
}

//void append(const char * &)
//Description: This function appends a label to the data structure
//Input: l, const char * & -- the label to be appended
//Output:None
//Return: None
void Labels::append(const char * & l){

	const char * end = strchr(l,'\n'); 

	if(end == '\0')			//If there is no newline
		end = strchr(l, '\0');

	unsigned int cpyLen = end-l;
	if(cpyLen > tagLen)
		cpyLen = tagLen;

	memcpy(data+((tagLen+1)*numLabels), l, cpyLen);

	data[((tagLen+1)*numLabels++)+cpyLen] = '\0';
}

//void append(const char * &)
//Description: This function appends a label to the data structure
//Input: l, const char * & -- the label to be appended
//Output:None
//Return: None
void Labels::append(const char * l) {

	const char * end = strchr(l, '\n'); 

	if(end == '\0')               //If there is no newline   
		end = strchr(l, '\0');

	unsigned int cpyLen = end-l;
	if(cpyLen > tagLen)
		cpyLen = tagLen;

	memcpy(data+((tagLen+1)*numLabels), l, cpyLen);
	data[((tagLen+1)*numLabels++)+cpyLen] = '\0';
}

//void append(const char *, int)
//Description: This function appends a label to the data structure
//Input: l, const char * & -- the label to be appended; cpyLen, 
//const int -- the size of the label
//Output:None
//Return: None
void Labels::append(const char * l, const int cpyLen) {

	memcpy(data+((tagLen+1)*numLabels), l, cpyLen);
	data[((tagLen+1)*numLabels++)+cpyLen] = '\0';
}

//string labelAt(const unsigned long long int n)
//Description: This function returns a label at the ith position
//Input: n, const unsigned long long int
//Output: None
//Return: string -- the label
string Labels::labelAt(const unsigned long long int n){

	return string(data + (n * (tagLen+1)));
}

//unsigned int getTagLen(void) const
//Description: This function returns the tag length
//Input:None
//Output:None
//Return:unsigned int -- the tag length
unsigned int Labels::getTagLen(void) const {
	return tagLen;
}

//unsigned long long int getNumLabels(void) const
//Description: This function returns the number of 
//labels that the Labels container contains
//Input:None
//Output:None
//Return: unsigned long long int -- the number of labels
unsigned long long int Labels::getNumLabels(void) const {
	return numLabels;
}

//unsigned long long int getOffset(void) const
//Description: This function returns the offset of the label set
//Input:None
//Output:None
//Return: unsigned long long int -- the offset
unsigned long long int Labels::getOffset(void) const {
	return offset;
}

//void write(FILE * &)
//Description: This function writes the labels to a file
//Input: FILE * & -- the file pointer
//Output:None
//Return: None
void Labels::write(FILE * & pFile){

	fwrite(&tagLen, sizeof(unsigned int), 1, pFile);
	fwrite(&offset, sizeof(unsigned long long int), 1, pFile);
	fwrite(&numLabels, sizeof(unsigned long long int), 1, pFile);
	fwrite(data, sizeof(char), numLabels * (tagLen+1), pFile);
}

//void write(FILE * &)
//Description: This function reads the labels from a file
//Input: FILE * & -- the file pointer
//Output:None
//Return: None
void Labels::read(FILE * & pFile){

	unsigned long long int intervalOffset; unsigned long long int intervalNumLabels;

	fread(&tagLen, sizeof(unsigned int), 1, pFile);
	fread(&intervalOffset, sizeof(unsigned long long int), 1, pFile);
	fread(&intervalNumLabels, sizeof(unsigned long long int), 1, pFile);
	
	fread(data+(numLabels*(tagLen+1)), sizeof(char), intervalNumLabels * (tagLen+1), pFile);
	numLabels+=(intervalNumLabels);
}
