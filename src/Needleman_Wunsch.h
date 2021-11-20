/*
 * Needlename_Wunsch.cpp
 *
 *      Author: Julia
 */

#include<iostream>
#include<string>
#include<algorithm>
#include<cstring>
using namespace std;

/* This is an algorithm for performing alignment. */

#define gapPenalty -2
#define matchScore 1
#define misMatch -1

template <class T1, class T2>
class Needleman_Wunsch {
public:
	Needleman_Wunsch(const int, const int, const int, const int);
	virtual ~Needleman_Wunsch();
	string Align(int **, T1 &, T1 &, const int, const int);
	string Align(int **, T1 &, T2 &, const int, const int);
	int misAligned() const;
	int alignLen() const;
	int yAlignmentStart() const;
	int xAlignmentStart() const;
	int yAlignmentEnd() const;
	int xAlignmentEnd() const;
	int yGapsTotal() const;
	int xGapsTotal() const;
private:
	int x_bound;
	int y_bound;
	int x_offset;
	int y_offset;
	int misAlign;
	int x_gaps;
	int y_gaps;
	int length;
	int x_start;
	int y_start;
	int x_end;
	int y_end;
};


template <class T1, class T2>
Needleman_Wunsch<T1, T2>::Needleman_Wunsch(const int x_len, const int y_len, const int x_off, const int y_off)
	:x_bound(x_len), y_bound(y_len), x_offset(x_off), y_offset(y_off), misAlign(0), x_gaps(0), y_gaps(0), length(0), x_start(x_len), y_start(y_len), x_end(0), y_end(0){;;;}


template <class T1, class T2>
string Needleman_Wunsch<T1, T2>::Align(int * score [], T1 & x_str, T1 & y_str, const int forwardGaps, const int backwardGaps)
{
	if(y_bound <=  0 || x_bound <= 0)
	{
		if(backwardGaps != 0 && forwardGaps != 0)
		{
			misAlign+=y_bound;
			misAlign+=x_bound;
			y_gaps+=x_bound;
			x_gaps+=y_bound;
			length+=(y_bound+x_bound);
		}
		string rtrn = "";
		return rtrn;
	}

	for(int x = 0; x < x_bound+1; x++) score[x][0] =  x * gapPenalty * forwardGaps;
	for(int y = 0; y < y_bound+1; y++) score[0][y] =  y * gapPenalty * forwardGaps;
	int oneScore;

	for(int x = 1; x < x_bound+1; x++){
		for(int y = 1; y < y_bound+1; y++){
			if(x_str.at(x-1+x_offset) == y_str.at(y-1+y_offset)) oneScore = matchScore;
			else oneScore = misMatch;
			score[x][y] = max(score[x-1][y-1]+oneScore, score[x-1][y] + gapPenalty);
			score[x][y] = max(score[x][y-1] + gapPenalty, score[x][y]);
		}
	}

	char alignment[2][(y_bound+x_bound)+1];
	int i = 0; int x = x_bound; int y = y_bound;

	if(backwardGaps == 0)
	{
		int max = INT_MIN;
		for(int k = 1; k < x_bound+1; k++)
		{
			if(score[k][y_bound] > max)
			{
				max = score[k][y_bound];
				x = k;
				y = y_bound;
			}
		}

		for(int k = 1; k < y_bound+1; k++)
		{
			if(score[x_bound][k] > max)
			{
				max = score[x_bound][k];
				x = x_bound;
				y = k;
			}
		}	

		for(int k = x_bound; k > x; k--)
		{   
			alignment[0][(y_bound+x_bound)-i-1] = x_str.at(k-1+x_offset);
			alignment[1][(y_bound+x_bound)-i-1] = '-';
			i++;
		}

		for(int k = y_bound; k > y; k--)
		{   
			alignment[0][(y_bound+x_bound)-i-1] = '-';
			alignment[1][(y_bound+x_bound)-i-1] = y_str.at(y-1+y_offset);
			i++;
		}

	}

	x_end = x-1;
	y_end = y-1;

	for(; !(x == 0 || y == 0); i++)
	{
		if(0 < x && 0 < y && x_str.at(x-1+x_offset) == y_str.at(y-1+y_offset))
			oneScore = matchScore;
		else{  oneScore = misMatch; }
		if(0 < x && 0 < y && score[x][y] == (score[x-1][y-1] + oneScore)){
			if(oneScore == misMatch)
				misAlign++;


			alignment[0][(y_bound+x_bound)-i-1] = x_str.at(x-1+x_offset);
			alignment[1][(y_bound+x_bound)-i-1] = y_str.at(y-1+y_offset);
			y = y-1; x = x-1;
		}else
		{
			if(0 < x && 0 < y && score[x][y] == score[x-1][y]+gapPenalty)
			{
				alignment[0][(y_bound+x_bound)-i-1] = x_str.at(x-1+x_offset);
				alignment[1][(y_bound+x_bound)-i-1] = '-';
				x = x-1;
				misAlign++;
				y_gaps++;

			}else
			{
				if(0 < x && 0 < y && score[x][y] == score[x][y-1]+gapPenalty)
					alignment[0][(y_bound+x_bound)-i-1]= '-';
				alignment[1][(y_bound+x_bound)-i-1] = y_str.at(y-1+y_offset);
				y = y-1;
				misAlign++;
				x_gaps++;
			}

		}

		length++;
	}

	y_start = y;
	x_start = x;

	if(forwardGaps == 1)
	{
		while(x != 0)
		{
			alignment[0][(y_bound+x_bound)-i-1] = x_str.at(x-1+x_offset);
			alignment[1][(y_bound+x_bound)-i-1] = '-';
			x--;
			i++;
			misAlign++;
			y_gaps++;
			
			length++;

		}

		while(y != 0)
		{
			alignment[0][(y_bound+x_bound)-i-1] = '-';
			alignment[1][(y_bound+x_bound)-i-1] = y_str.at(y-1+y_offset);
			y--;
			i++;
			misAlign++;
			x_gaps++;

			length++;
		}

		y_start = y;
		x_start = x;
	}else{

		while(x != 0)
		{
			alignment[0][(y_bound+x_bound)-i-1] = x_str.at(x-1+x_offset);
			alignment[1][(y_bound+x_bound)-i-1] = '-';
			x--;
			i++;
		}

		while(y != 0)
		{
			alignment[0][(y_bound+x_bound)-i-1] = '-';
			alignment[1][(y_bound+x_bound)-i-1] = y_str.at(y-1+y_offset);
			y--;
			i++;
		}
	}

	char arr[2*i+1];
	memcpy(arr, alignment[0]+(y_bound+x_bound)-i, i);  //change by one
	memcpy(arr+i, alignment[1]+(y_bound+x_bound)-i, i); //change by one
	arr[2*i] = '\0';

	return string (arr);
}

template <class T1, class T2>
string Needleman_Wunsch<T1, T2>::Align(int * score [], T1 & x_str, T2 & y_str, const int forwardGaps, const int backwardGaps)
{
	if(y_bound <=  0 || x_bound <= 0)
	{
		if(backwardGaps != 0 && forwardGaps != 0)
		{
			misAlign+=y_bound;
			misAlign+=x_bound;
			y_gaps+=x_bound;
			x_gaps+=y_bound;

			length+=(x_bound+y_bound);
		}
		string rtrn = "";
		return rtrn;
	}

	for(int x = 0; x < x_bound+1; x++) score[x][0] =  x * gapPenalty * forwardGaps;
	for(int y = 0; y < y_bound+1; y++) score[0][y] =  y * gapPenalty * forwardGaps;
	int oneScore;
	for(int x = 1; x < x_bound+1; x++){
		for(int y = 1; y < y_bound+1; y++){
			if(x_str.at(x-1+x_offset) == y_str.at(y-1+y_offset)) oneScore = matchScore;
			else oneScore = misMatch;
			score[x][y] = max(score[x-1][y-1]+oneScore, score[x-1][y] + gapPenalty);
			score[x][y] = max(score[x][y-1] + gapPenalty, score[x][y]);
		}
	}


	char alignment[2][(y_bound+x_bound)+1];
	int i = 0; int x = x_bound; int y = y_bound;

	if(backwardGaps == 0)
	{
		int max = INT_MIN;
		for(int k = 1; k < x_bound+1; k++)
		{
			if(score[k][y_bound] > max)
			{
				max = score[k][y_bound];
				x = k;
				y = y_bound;
			}
		}

		for(int k = 1; k < y_bound+1; k++)
		{
			if(score[x_bound][k] > max)
			{
				max = score[x_bound][k];
				x = x_bound;
				y = k;
			}
		}	

		for(int k = x_bound; k > x; k--)
		{   
			alignment[0][(y_bound+x_bound)-i-1] = x_str.at(k-1+x_offset);
			alignment[1][(y_bound+x_bound)-i-1] = '-';
			i++;
		}

		for(int k = y_bound; k > y; k--)
		{   
			alignment[0][(y_bound+x_bound)-i-1] = '-';
			alignment[1][(y_bound+x_bound)-i-1] = y_str.at(k-1+y_offset);
			i++;
		}

		length++;
	}

	x_end = x-1;
	y_end = y-1;

	for(; !(x == 0 || y == 0); i++)
	{
		if(0 < x && 0 < y && x_str.at(x-1+x_offset) == y_str.at(y-1+y_offset))
			oneScore = matchScore;
		else{  oneScore = misMatch; }
		if(0 < x && 0 < y && score[x][y] == (score[x-1][y-1] + oneScore)){
			if(oneScore == misMatch)
				misAlign++;
			alignment[0][(y_bound+x_bound)-i-1] = x_str.at(x-1+x_offset);
			alignment[1][(y_bound+x_bound)-i-1] = y_str.at(y-1+y_offset);
			y = y-1; x = x-1;
		}else
		{
			if(0 < x && 0 < y && score[x][y] == score[x-1][y]+gapPenalty)
			{
				alignment[0][(y_bound+x_bound)-i-1] = x_str.at(x-1+x_offset);
				alignment[1][(y_bound+x_bound)-i-1] = '-';
				x = x-1;
				misAlign++;
				y_gaps++;
			}else
			{
				if(0 < x && 0 < y && score[x][y] == score[x][y-1]+gapPenalty)
					alignment[0][(y_bound+x_bound)-i-1]= '-';
				alignment[1][(y_bound+x_bound)-i-1] = y_str.at(y-1+y_offset);
				y = y-1;
				misAlign++;
				x_gaps++;
			}

		}
	}

	y_start = y;
	x_start = x;

	if(forwardGaps == 1)
	{
		while(x != 0)
		{
			alignment[0][(y_bound+x_bound)-i-1] = x_str.at(x-1+x_offset);
			alignment[1][(y_bound+x_bound)-i-1] = '-';
			x--;
			i++;

			misAlign++;
			y_gaps++;

			length++;
		}

		while(y != 0)
		{
			alignment[0][(y_bound+x_bound)-i-1] = '-';
			alignment[1][(y_bound+x_bound)-i-1] = y_str.at(y-1+y_offset);
			y--;
			i++;

			misAlign++;
			x_gaps++;
			
			length++;
		}

		y_start = y;
		x_start = x;
	}else{

		while(x != 0)
		{
			alignment[0][(y_bound+x_bound)-i-1] = x_str.at(x-1+x_offset);
			alignment[1][(y_bound+x_bound)-i-1] = '-';
			x--;
			i++;
		}

		while(y != 0)
		{
			alignment[0][(y_bound+x_bound)-i-1] = '-';
			alignment[1][(y_bound+x_bound)-i-1] = y_str.at(y-1+y_offset);
			y--;
			i++;
		}
	}


	char arr[2*i+1];
	memcpy(arr, alignment[0]+(y_bound+x_bound)-i, i);  //change by one
	memcpy(arr+i, alignment[1]+(y_bound+x_bound)-i, i); //change by one
	arr[2*i] = '\0';

	return string (arr);
}

	template <class T1, class T2>
int Needleman_Wunsch<T1, T2>::misAligned() const
{
	return misAlign;
}

	template <class T1, class T2>
int Needleman_Wunsch<T1, T2>::alignLen() const
{
	return length;
}

	template <class T1, class T2>
int Needleman_Wunsch<T1, T2>::yGapsTotal() const
{
	return y_gaps;
}

	template <class T1, class T2>
int Needleman_Wunsch<T1, T2>::xGapsTotal() const
{
	return x_gaps;
}

	template <class T1, class T2>
int Needleman_Wunsch<T1, T2>::yAlignmentStart() const
{
	return y_start;
}

	template <class T1, class T2>
int Needleman_Wunsch<T1, T2>::xAlignmentStart() const
{
	return x_start;
}

	template <class T1, class T2>
int Needleman_Wunsch<T1, T2>::yAlignmentEnd() const
{
	return y_end;

}

	template <class T1, class T2>
int Needleman_Wunsch<T1, T2>::xAlignmentEnd() const
{
	return x_end;
}

	template <class T1, class T2>
Needleman_Wunsch<T1, T2>::~Needleman_Wunsch() {

}

