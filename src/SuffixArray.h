/*
 * SuffixArray.h
 *
 *      Author: Julia
 */

#ifndef SUFFIXARRAY_H_
#define SUFFIXARRAY_H_

#include<iostream>
#include<time.h>
#include<string>
#include <stdio.h>
#include <stdlib.h>
#include<assert.h>
#include<fstream>
#include<math.h>
#include<map>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;


//Suffix Array Class
//This class uses Sadakane's algorithm (n log n) algorithm to produce the suffix array
//The class is templated so that it is able to accept both genome indexes and
//fragment indexes

template <class T>
class SuffixArray{
	public:
		SuffixArray();
		~SuffixArray();
		void buildArray(T &, const unsigned long long int, const unsigned long long int, const int numIter = INT_MAX);
		void write(FILE * &);
		void read(FILE * & );
		unsigned long long int SA_to_SI(const unsigned long long int) const;
		unsigned long long int ISA_to_SA(const unsigned long long int) const;
	private:
		void radixSort_LS(T &, const unsigned long long int, const unsigned long long int);
		void ternary_split_quicksort_LS(const unsigned long long int, const unsigned long long int, const unsigned long long int);
		unsigned long long int SA_l;
		unsigned long long int * SA;
		unsigned long long int * ISA;
};

//SuffixArray()
//Description:This is the constructor for the suffix array
//Input:None
//Output:None
//Return:None
	template <class T>
SuffixArray<T>::SuffixArray()
	:SA_l(0), SA(nullptr), ISA(nullptr){;;;}

	//buildArray(T &);
	//Description: This function builds the suffix array using Sadakane's algorithm
	//Input: T (Either a Genome or Fragment index)
	//Output:None
	//Return:None
	template <class T>
void SuffixArray<T>::buildArray(T & index, const unsigned long long int begin, const unsigned long long int numChars, const int numIter)
{
	SA_l = numChars+1;
	SA = new unsigned long long int [numChars+1];
	ISA = new unsigned long long int [numChars+1];
	radixSort_LS(index, begin, numChars);

	bool iteration; unsigned long long int h = 1; int n = 0;
	do{
		iteration = false;
		for(unsigned long long int i = 0; i < numChars+1; i = ISA[SA[i]]+1)
		{
			if(i < ISA[SA[i]])
			{
				ternary_split_quicksort_LS(i, ISA[SA[i]], h);
				iteration = true;
			}
		}
		h = h*2;
	}while(iteration && n++ < numIter);
}

//write(FILE * &)
//Description:This function writes the suffix array to
//Input:FILE *, a pointer to the file
//Output:None
//Return:None
	template <class T>
void SuffixArray<T>::write(FILE * & out)
{
	fwrite(&SA_l, sizeof(unsigned long long int), 1, out);
	fwrite(SA, sizeof(unsigned long long int), SA_l, out);
	fwrite(ISA, sizeof(unsigned long long int), SA_l, out);
}

//read(FILE * &)
//Description:This function reads a suffix array data structure
//Input:FILE * &, a pointer to the file
//Output:None
//Return:none
	template <class T>
void SuffixArray<T>::read(FILE * & in)
{
	delete [] SA;
	delete [] ISA;

	fread(&SA_l, sizeof(unsigned long long int), 1, in);
	SA = new unsigned long long int[SA_l];
	fread(SA, sizeof(unsigned long long int), SA_l, in);
	ISA = new unsigned long long int[SA_l];
	fread(ISA, sizeof(unsigned long long int), SA_l, in);
}

//SA_to_SI(const int) const
//Description:This function converts a position in the suffix array to
//the position in the index
//Input:n ( an int) the position in the SA
//Output:None
//Return:None
template <class T>
unsigned long long int SuffixArray<T>::SA_to_SI(const unsigned long long int n) const
{
	return SA[n];
}

//ISA_to_SA(const long int)
//Description: This function returns the value of the ISA at position n
//Input:N, (a long int) the position in the ISA
//Output:None
//Return:None
template <class T>
unsigned long long int SuffixArray<T>::ISA_to_SA(const unsigned long long int n) const
{
	return ISA[n];
}

//SuffixArray()
//Description:This is a class destructor
//Input:None
//Output:None
//Return:None
	template <class T>
SuffixArray<T>::~SuffixArray()
{
	delete [] SA;
	delete [] ISA;
}

//radixSort(T & )
//Description:This function performs a radix sort
//Input:T (template class ) Genome or a Fragment index, const unsigned long long int begin, const unsigned long long int numChars, const int stepSize
//Output:None
//Return:None
	template <class T>
void SuffixArray<T>::radixSort_LS(T & index, const unsigned long long int begin, const unsigned long long int numChars)
{
	unsigned int m = 5;

	//Setting the values for the $ character
	int * count = new int[m];  memset(count, 0, sizeof(int)*m); ISA[numChars] = 0; SA[0] = numChars;  count[0] = 1;

	for(unsigned long long int i = 0; i < numChars; i++)
	{
		count[mapVals[index.at(begin+i)]+1]++;
	}

	for(unsigned long long int i = 1; i < m; i++)
	{
		count[i]+=count[i-1];
	}

	for(long long int i = numChars-1; i >= 0; i--)
	{
		ISA[i] = count[mapVals[index.at(begin+i)]+1]-1;
	}

	for(long long int i = numChars-1; i >= 0; i--)
	{
		unsigned long long int charType = mapVals[index.at(begin+i)]+1;
		SA[count[charType]-1] = i;
		count[charType]--;
	}

	delete [] count;
}


//ternary_split_quicksort_LS(int l, int r, int h)
//Description:This function performs a quicksort on a range of the SA
//Input: l (an int), r (an int), h, (an int), the right and left bounds of
//the interval, h is the iteration number
//Output:None
//Return:None
	template <class T>
void SuffixArray<T>::ternary_split_quicksort_LS(const unsigned long long int l, const unsigned long long int r, const unsigned long long int h)
{
	if(l > r) return;
	if(l == r){
		ISA[SA[l]] = l; return;
	}
	srand(time(NULL));
	int random = rand()%(r-l+1);
	unsigned long long int v = ISA[SA[l+ random]+h];
	unsigned long long int i = l; unsigned int long long mi = l; unsigned long long int j = r; unsigned long long int mj = r; unsigned long long int tmp;
	for(;;)
	{
		for(; (i <= j && ISA[SA[i]+h] <= v); i++)
			if(ISA[SA[i]+h] == v)
			{
				tmp = SA[i]; SA[i] = SA[mi]; SA[mi] = tmp; mi++;
			}
		for(; i<= j && v <= ISA[SA[j]+h]; j--)
			if(ISA[SA[j]+h] == v)
			{
				tmp = SA[j]; SA[j] = SA[mj]; SA[mj] = tmp; mj--;
			}
		if(i > j) break;
		tmp = SA[i]; SA[i] = SA[j]; SA[j] = tmp; i++; j--;
	}
	for(mi--, i--; l<= mi; mi--, i--){
		tmp = SA[i];  SA[i] = SA[mi]; SA[mi] = tmp;
	}
	for(mj++, j++; mj <= r; mj++, j++){
		tmp = SA[j]; SA[j] = SA[mj]; SA[mj] = tmp;
	}

	ternary_split_quicksort_LS(l, i, h);
	for(unsigned long long int k = i+1; k < j; k++) ISA[SA[k]] = j-1;
	ternary_split_quicksort_LS(j, r, h);
}

#endif /* SUFFIXARRAY_H_ */
