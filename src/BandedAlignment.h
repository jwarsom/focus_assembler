/*
 * BandedAlignment.h
 *
 *      Author: Julia
 */


#ifndef BANDEDALIGNMENT_H_
#define BANDEDALIGNMENT_H_


/* This is an algorithm for performing alignment This algorithm
 * performs a banded alignment This algorithm is meant to quickly
 * align regions between seeds. This alignment algorithm will handle
 * sequences of all sizes provided with the appropriate matrix */

template <class T>
class BandedAlignment {
	public:
		BandedAlignment(unsigned long long int, unsigned long long int, int, int, int, bool);
		pair<int, int> align(int * [], T &, T &, const unsigned int reverse = 0, const bool gapPenalty = false);
		pair<int, int> alignReverse(int * [], T &, T &, const unsigned int reverse = 0, const bool gapPenalty = false);
		int getXStart(void) const;
		int getXEnd(void) const;
		int getYStart(void) const;
		int getYEnd(void) const;
		pair<int, int> traceBackTrim(int * [], T &, T &, const unsigned int, const unsigned int, const unsigned int reverse = 0);
		pair<int, int> traceBackReverseTrim(int * [], T &, T &, const unsigned int, const unsigned int, const unsigned int reverse = 0);
	private:
		pair<int, int> traceBack(int * [], T &, T &, const unsigned int reverse = 0);
		pair<int, int> traceBackReverse(int * [], T &, T &, const unsigned int reverse = 0);

		//Where we start the alignment
		unsigned long long int beginX;
		unsigned long long int beginY;

		//The lengths of the alignment
		int alignLenX;
		int alignLenY;

		//The band length
		int band;

		//Are gaps allowed in the alignment
		bool allowGaps;

		//Tracks gaps and misAlignment
		int misAlign;
		int xGaps;
		int yGaps;

		//Tracks alignment start/end
		int xStart;
		int xEnd;
		int yStart;
		int yEnd;
};


//BandedAlignment(unsigned long long int, unsigned long long int, int, int, int, bool)
//Description: This is the constructor for the class
//Input: unsigned long long int bX, the start int position in the first index, unsigned long long int bY, the
//starting position in the second index, int lX, the length of the string to be aligned, int lY the length of the
//second string to be aligned, int k, the band length, bool g, if gaps are allowed or not
//Output:None
//Retrun:None
	template <class T>
BandedAlignment<T>::BandedAlignment(unsigned long long int bX, unsigned long long int bY, int lX, int lY, int k, bool g)
	:beginX(bX), beginY(bY), alignLenX(lX), alignLenY(lY), band(k), allowGaps(g), misAlign(0), xGaps(0), yGaps(0), xStart(0), xEnd(0), yStart(0), yEnd(0) {;;;}


	//pair<int, int> align(int * [], T &, T &)
	//Description: This function aligns two sequences
	//Input: int * [], matrix, the alignment matrix, T & index1, T & index2
	//Output:None
	//Return: pair<int, int> the length of the alignment and the mismatches in the alignment
	template <class T>
	pair<int, int> BandedAlignment<T>::align(int * matrix [], T & index1, T & index2, const unsigned int reverse, const bool gapPenalty){

		if(alignLenY <=  0 || alignLenX <= 0)
		{
			if(allowGaps) return pair<int, int>(0, 0);
			return pair<int, int>(abs(alignLenY+alignLenX), abs(alignLenY+alignLenX));
		}
		
		int gp = 0;
	
		if(gapPenalty)
			gp = -1;

		//Initiate the alignment matrix with 0s at the boundaries
		for(int i = 0; i < band+1 && i < alignLenY+1; i++) matrix[0][i] = gp * i;
		for(int i = 0; i < band+1 && i < alignLenX+1; i++) matrix[i][0] = gp * i;

		//Initiate the bounds of the k-band with an extremely low value
		for(int i = 0; i < alignLenX+1 && i + band+1 < alignLenY+1; i++) matrix[i][i+band+1] = -10000;
		for(int i = 0; i < alignLenY+1 && i + band+1 < alignLenX+1; i++) matrix[i+band+1][i] = -10000;

		//The actual alignment
		int xCoordBound = min(alignLenX+1, alignLenY+band+1);

		for(int i = 1; i < xCoordBound; i++)
		{
			//Bounds of the bands in the diagonal
			int yCoordBound1 = 1 > i-band ? 1 : i-band;
			int yCoordBound2 = alignLenY+1 <  i+band+1 ? alignLenY+1 : i+band+1;

			for(int j = yCoordBound1; j < yCoordBound2; j++)
			{
				//The gap scores
				int verticalVal = matrix[i-1][j]-2;
				int horizontalVal = matrix[i][j-1]-2;

				//If there is a match
				int match = -1;
				if(index1.at(beginX+i-1, reverse) == index2.at(beginY+j-1)) match = 1;

				int diagonalVal = matrix[i-1][j-1] + match;

				if(diagonalVal >= horizontalVal && diagonalVal >= verticalVal)
				{
					matrix[i][j] = diagonalVal;
				}else
				{
					if(horizontalVal >= verticalVal)
					{
						matrix[i][j] = horizontalVal;
					}else
					{
						matrix[i][j] = verticalVal;
					}
				}
			}

		}
		return traceBack(matrix, index1, index2, reverse);
	}


//pair<int, int> alignReverse(int * [], T &, T &)
//Description: This function aligns two sequences in reverse
//Input: int * [], matrix, the alignment matrix, T & index1, T & index2
//Output:None
//Return: pair<int, int> the length of the alignment and the mismatches in the alignment
template <class T>
pair<int, int> BandedAlignment<T>::alignReverse(int * matrix [], T & index1, T & index2, const unsigned int reverse, bool gapPenalty){


	if(alignLenY <=  0 || alignLenX <= 0)
	{
		yEnd = alignLenY; xEnd = alignLenX; yStart = alignLenY; xStart = alignLenX; 
		if(allowGaps) return pair<int, int>(0, 0);
		return pair<int, int>(abs(alignLenY+alignLenX), abs(alignLenY+alignLenX));
	}
	
	int gp = 0;

	if(gapPenalty)
		gp = -1;
	
	//Initiate the alignment matrix with 0s at the boundaries
	for(int i = 0; i < band+1 && i < alignLenY+1; i++) matrix[0][i] = gp*i;
	for(int i = 0; i < band+1 && i < alignLenX+1; i++) matrix[i][0] = gp*i;

	//Intitiate the bounds of the k-band with an extremely low value
	for(int i = 0; i < alignLenX+1 && i + band+1 < alignLenY+1; i++) matrix[i][i+band+1] = -10000;
	for(int i = 0; i < alignLenY+1 && i + band+1 < alignLenX+1; i++) matrix[i+band+1][i] = -10000;

	//The actual alignment
	int xCoordBound = min(alignLenX+1, alignLenY+band+1);
	for(int i = 1; i < xCoordBound; i++)
	{
		//Bounds of the bands in the diagonal
		int yCoordBound1 = 1 > i-band ? 1 : i-band;
		int yCoordBound2 = alignLenY+1 <  i+band+1 ? alignLenY+1 : i+band+1;

		for(int j = yCoordBound1; j < yCoordBound2; j++)
		{
			//The gap scores
			int verticalVal = matrix[i-1][j]-2;
			int horizontalVal = matrix[i][j-1]-2;

			//If there is a match
			int match = -1;
			if(index1.at(beginX+alignLenX-i, reverse) == index2.at(beginY+alignLenY-j)) match = 1;

			int diagonalVal = matrix[i-1][j-1] + match;

			if(diagonalVal >= horizontalVal && diagonalVal >= verticalVal)
			{
				matrix[i][j] = diagonalVal;
			}else
			{
				if(horizontalVal >= verticalVal)
				{
					matrix[i][j] = horizontalVal;
				}else
				{
					matrix[i][j] = verticalVal;
				}
			}
		}
	}
	return traceBackReverse(matrix, index1, index2, reverse);
}


//pair<int, int> traceBack(int * [])
//Description: This function traces back an alignment
//Input: int * [], matrix, the alignment matrix
//Output:None
//Return: pair<int, int> the alignmnent length and mismatches
template <class T>
pair<int, int> BandedAlignment<T>::traceBack(int * matrix [], T & index1, T & index2, const unsigned int reverse){

	//Going to find the maximum
	int currentMax = INT_MIN;
	int xCoord = 0;  int yCoord = 0;

	int begin = 1 >  alignLenY-band ? 1 : alignLenY-band;
	int end = min(alignLenX+1, alignLenY+band+1);

	pair<int, int> alignScore(0, 0);

	for(int i = begin; i < end; i++)
	{
		if(matrix[i][alignLenY] >= currentMax)
		{
			xCoord = i; yCoord = alignLenY; currentMax = matrix[i][alignLenY];
		}
	}

	begin = (1 > alignLenX-band) ? 1 : alignLenX-band;
	end = min(alignLenY+1, alignLenX+band+1);

	for(int i = begin; i < end; i++)
	{
		if(matrix[alignLenX][i] >= currentMax)
		{
			xCoord = alignLenX; yCoord = i; currentMax = matrix[alignLenX][i];
		}
	}

	xEnd = xCoord-1;
	yEnd = yCoord-1;

	if(!allowGaps)
	{
		misAlign+=(alignLenX-xCoord + alignLenY-yCoord);
		xGaps+=(alignLenY-yCoord);
		yGaps+=(alignLenX-xCoord);
		alignScore.first+=(xGaps+yGaps);
		alignScore.second+=(xGaps+yGaps);
	}

	while(xCoord != 0 && yCoord != 0)
	{
		if(matrix[xCoord-1][yCoord-1] == matrix[xCoord][yCoord]-1 && index1.at(beginX+xCoord-1, reverse) == index2.at(beginY+yCoord-1))
		{
			xCoord = xCoord-1; yCoord = yCoord-1;
		}else
		{
			if(matrix[xCoord-1][yCoord-1] == matrix[xCoord][yCoord]+1)
			{
				xCoord = xCoord-1; yCoord = yCoord-1; misAlign++;
			}else
			{
				if(matrix[xCoord-1][yCoord] == matrix[xCoord][yCoord]+2)
				{
					xCoord = xCoord-1; misAlign++; yGaps++;
				}else
				{
					yCoord = yCoord-1; misAlign++; xGaps++;
				}
			}

			alignScore.second++;
		}

		alignScore.first++;
	}

	xStart = xCoord;
	yStart = yCoord;

	misAlign+=(xCoord + yCoord);
	xGaps+=(yCoord);
	yGaps+=(xCoord);

	alignScore.first+=(xCoord+yCoord);
	alignScore.second+=(xCoord+yCoord);

	return alignScore;
}



//pair<int, int> traceBackReverse(int * [])
//Description: This function traces back an alignment
//Input: int * [], matrix, the alignment matrix
//Output:None
//Return: pair<int, int> the alignmnent length and mismatches
template <class T>
pair<int, int> BandedAlignment<T>::traceBackReverse(int * matrix [], T & index1, T & index2, const unsigned int reverse){

	//Going to find the maximum
	int currentMax = INT_MIN;
	int xCoord = 0;  int yCoord = 0;

	int begin = 1 >  alignLenY-band ? 1 : alignLenY-band;
	int end = min(alignLenX+1, alignLenY+band+1);

	for(int i = begin; i < end; i++)
	{
		if(matrix[i][alignLenY] >= currentMax)
		{
			xCoord = i; yCoord = alignLenY; currentMax = matrix[i][alignLenY];
		}
	}

	begin = (1 > alignLenX-band) ? 1 : alignLenX-band;
	end = min(alignLenY+1, alignLenX+band+1);

	for(int i = begin; i < end; i++)
	{
		if(matrix[alignLenX][i] >= currentMax)
		{
			xCoord = alignLenX; yCoord = i; currentMax = matrix[alignLenX][i];
		}
	}
	

	xStart = alignLenX - xCoord;
	yStart = alignLenY - yCoord;

	if(!allowGaps)
	{
		misAlign+=(alignLenX-xCoord + alignLenY-yCoord);
		xGaps+=(alignLenY-yCoord);
		yGaps+=(alignLenX-xCoord);
	}

	while(xCoord != 0 && yCoord != 0)
	{
		if(matrix[xCoord-1][yCoord-1] == matrix[xCoord][yCoord]-1 && index1.at(beginX+alignLenX-xCoord, reverse) == index2.at(beginY+alignLenY-yCoord))
		{
			xCoord = xCoord-1; yCoord = yCoord-1;
		}else
		{
			if(matrix[xCoord-1][yCoord-1] == matrix[xCoord][yCoord]+1)
			{
				xCoord = xCoord-1; yCoord = yCoord-1; misAlign++;
			}else
			{
				if(matrix[xCoord-1][yCoord] == matrix[xCoord][yCoord]+2)
				{
					xCoord = xCoord-1; misAlign++; yGaps++;
				}else
				{
					yCoord = yCoord-1; misAlign++; xGaps++;
				}
			}
		}
	}

	xEnd = alignLenX - xCoord - 1;
	yEnd = alignLenY - yCoord - 1;

	misAlign+=(xCoord + yCoord);

	xGaps+=(yCoord);
	yGaps+=(xCoord);

	return pair<int, int>(min(alignLenX+xGaps,alignLenY+yGaps), misAlign);
}

//pair<int, int> traceBack(int * [], T &, T &, const unsigned int, const unsigned int)
//Description: traces back the alignment from a certain point on the x fragment
//Input:int * [], (matrix): the alignment matrix, T & (index1): the x sequence, T & (index2):
//the y sequence, const unsigned int (trimLen): the trimming length of the x sequence,
//const unsigned int (reverse): is the sequence a reverse complement. 
//Output:None
//Return:pair<int, int>: the alignment length and alignment mismatches
	template <class T>
pair<int, int> BandedAlignment<T>::traceBackTrim(int * matrix [], T & index1, T & index2, const unsigned int trimLenX, const unsigned int trimLenY, const unsigned int reverse)
{
	alignLenX = alignLenX - trimLenX;
	alignLenY = alignLenY - trimLenY;
	misAlign = 0; xGaps = 0; yGaps = 0;

	//Going to find the maximum
	int currentMax = INT_MIN;
	int xCoord = 0;  int yCoord = 0;

	int begin = 1 >  alignLenY-band ? 1 : alignLenY-band;
	int end = min(alignLenX+1, alignLenY+band+1);

	pair<int, int> alignScore(0, 0);

	if(alignLenX <= 0 || alignLenY <= 0)
	{
		yStart = yEnd = xStart = xEnd = 0;
		alignLenX = alignLenX + trimLenX;
		alignLenY = alignLenY + trimLenY;
		return alignScore;
	}

	for(int i = begin; i < end; i++)
	{
		if(matrix[i][alignLenY] >= currentMax)
		{
			xCoord = i; yCoord = alignLenY; currentMax = matrix[i][alignLenY];
		}
	}

	begin = (1 > alignLenX-band) ? 1 : alignLenX-band;
	end = min(alignLenY+1, alignLenX+band+1);

	for(int i = begin; i < end; i++)
	{
		if(matrix[alignLenX][i] >= currentMax)
		{
			xCoord = alignLenX; yCoord = i; currentMax = matrix[alignLenX][i];
		}
	}

	xEnd = xCoord-1;
	yEnd = yCoord-1;

	if(!allowGaps)
	{
		misAlign+=(alignLenX-xCoord + alignLenY-yCoord);
		xGaps+=(alignLenY-yCoord);
		yGaps+=(alignLenX-xCoord);
		alignScore.first+=(xGaps+yGaps);
		alignScore.second+=(xGaps+yGaps);
	}

	while(xCoord != 0 && yCoord != 0)
	{
		if(matrix[xCoord-1][yCoord-1] == matrix[xCoord][yCoord]-1 && index1.at(beginX+xCoord-1, reverse) == index2.at(beginY+yCoord-1))
		{
			xCoord = xCoord-1; yCoord = yCoord-1;
		}else
		{
			if(matrix[xCoord-1][yCoord-1] == matrix[xCoord][yCoord]+1)
			{
				xCoord = xCoord-1; yCoord = yCoord-1; misAlign++;
			}else
			{
				if(matrix[xCoord-1][yCoord] == matrix[xCoord][yCoord]+2)
				{
					xCoord = xCoord-1; misAlign++; yGaps++;
				}else
				{
					yCoord = yCoord-1; misAlign++; xGaps++;
				}
			}

			alignScore.second++;
		}

		alignScore.first++;
	}

	xStart = xCoord;
	yStart = yCoord;

	misAlign+=(xCoord + yCoord);
	xGaps+=(yCoord);
	yGaps+=(xCoord);

	alignScore.first+=(xCoord+yCoord);
	alignScore.second+=(xCoord+yCoord);

	alignLenX = alignLenX + trimLenX;
	alignLenY = alignLenY + trimLenY;

	return alignScore;
}

//pair<int, int> traceBackReverse(int * [])
//Description: This function traces back an alignment
//Input: int * [], matrix, the alignment matrix
//Output:None
//Return: pair<int, int> the alignmnent length and mismatches
template <class T>
pair<int, int> BandedAlignment<T>::traceBackReverseTrim(int * matrix [], T & index1, T & index2, const unsigned int trimLenX, const unsigned int trimLenY,  const unsigned int reverse){

	alignLenX = alignLenX - trimLenX;
	alignLenY = alignLenY - trimLenY;

	//Going to find the maximum
	int currentMax = INT_MIN;
	int xCoord = 0;  int yCoord = 0;

	misAlign = 0; xGaps = 0; yGaps = 0;

	if(alignLenX <= 0 || alignLenY <= 0)
	{
		yStart = yEnd = xStart = xEnd = 0;
		alignLenX = alignLenX + trimLenX;
		alignLenY = alignLenY + trimLenY;
		return pair<int, int> (0,0);
	}

	int begin = 1 >  alignLenY-band ? 1 : alignLenY-band;
	int end = min(alignLenX+1, alignLenY+band+1);

	for(int i = begin; i < end; i++)
	{
		if(matrix[i][alignLenY] >= currentMax)
		{
			xCoord = i; yCoord = alignLenY; currentMax = matrix[i][alignLenY];
		}
	}

	begin = (1 > alignLenX-band) ? 1 : alignLenX-band;
	end = min(alignLenY+1, alignLenX+band+1);

	for(int i = begin; i < end; i++)
	{
		if(matrix[alignLenX][i] >= currentMax)
		{
			xCoord = alignLenX; yCoord = i; currentMax = matrix[alignLenX][i];
		}
	}

	xStart = alignLenX - xCoord;
	yStart = alignLenY - yCoord;

	if(!allowGaps)
	{
		misAlign+=(alignLenX-xCoord + alignLenY-yCoord);
		xGaps+=(alignLenY-yCoord);
		yGaps+=(alignLenX-xCoord);
	}

	while(xCoord != 0 && yCoord != 0)
	{
		if(matrix[xCoord-1][yCoord-1] == matrix[xCoord][yCoord]-1 && index1.at((beginX+trimLenX)+alignLenX-xCoord, reverse) == index2.at((beginY+trimLenY)+alignLenY-yCoord))
		{
			xCoord = xCoord-1; yCoord = yCoord-1;
		}else
		{
			if(matrix[xCoord-1][yCoord-1] == matrix[xCoord][yCoord]+1)
			{
				xCoord = xCoord-1; yCoord = yCoord-1; misAlign++;
			}else
			{
				if(matrix[xCoord-1][yCoord] == matrix[xCoord][yCoord]+2)
				{
					xCoord = xCoord-1; misAlign++; yGaps++;
				}else
				{
					yCoord = yCoord-1; misAlign++; xGaps++;
				}
			}
		}
	}

	xEnd = alignLenX - xCoord - 1;
	yEnd = alignLenY - yCoord - 1;

	misAlign+=(xCoord + yCoord);

	xGaps+=(yCoord);
	yGaps+=(xCoord);

	pair<int, int> alignScore (min(alignLenX+xGaps,alignLenY+yGaps), misAlign);

	alignLenX = alignLenX + trimLenX;
	alignLenY = alignLenY + trimLenY;

	return alignScore;
}

//int getXStart(void) const
//Description: Returns the start of the alignment
//on the x fragment
//Input:None
//Output:None
//Return:(int):the start position
template <class T>
int BandedAlignment<T>::getXStart(void) const {
	return xStart;
}

//int getXEnd(void) const
//Description: Returns the end of the alignment
//on the x fragment
//Input:None
//Output:None
//Return:(int):the end position
template <class T>
int BandedAlignment<T>::getXEnd(void) const {
	return xEnd;
}

//int getYStart(void) const
//Description: Returns the start of the alignment
//on the y fragment
//Input:None
//Output:None
//Return:(int):the start position
template <class T>
int BandedAlignment<T>::getYStart(void) const {
	return yStart;
}

//int getYEnd(void) const
//Description: Returns the end of the alignment
//on the y fragment
//Input:None
//Output:None
//Return:(int):the end position
template <class T>
int BandedAlignment<T>::getYEnd(void) const {
	return yEnd;
}

#endif /* BANDEDALIGNMENT_H_ */




