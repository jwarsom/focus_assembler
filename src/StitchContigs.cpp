#include <iostream>
#include <sys/stat.h>
#include <map>
#include <fstream>
#include <string>
#include <algorithm>
#include <cstdio>
#include <stdint.h>
#include <vector>
#include <ctime>
#include <list>
#include <set>
#include <climits>
#include <cmath>
#include <queue>
#include <dirent.h>
using namespace std;
#include "Dictionary.h"
#include "Fragment_Index.h"
#include "Needleman_Wunsch.h"

//Size unit
#define megabyte 1024000

char arr []  = {'A', 'C', 'T', 'G', '-'};

#define pathTag 1
#define endTag 0
#define help 911

void help_message(){

	cout<<"\n\n";
	cout<<"\t\t\t Contig Stitcher"<<endl;
	cout<<"\t\t\tWritten by Julia Warnke-Sommer\n"<<endl;
	cout<<"Command line arguments\n"<<endl;
	cout<<"--workDir :working directory :required"<<endl;
	cout<<"for the algorithm to continue merging"<<endl;
	cout<<"--verbose : verbose mode"<<endl;

	cout<<"\n\n"<<endl;

	cout<<"Exiting program"<<endl;
	exit(help);
}

int countFiles(const char []);

int main(int argc, char * argv [])
{
	//User Input Flags
	string workDirS = "--workDir";               //Current Working Directory
	string timeS = "--time"; //Report the time

	string workDir = "";
	bool reportTime = false;

	for(int i = 1; i < argc; i+=2)
	{
		//The working directory
		if(argv[i] == workDirS)
		{
			workDir = argv[i+1];
			unsigned int pos = workDir.find_last_of("/");
			if(pos == workDir.size()-1)
				workDir.erase(pos);
		}

		//Display time
		if(argv[i] == timeS && reportTime != true)
		{
			reportTime = true;
			i--;
		}
	}

	if(reportTime == true)
	{
		time_t rawtime;
		struct tm * timeinfo;
		time(&rawtime);
		timeinfo = localtime(&rawtime);
		cerr<<"start time: "<<asctime(timeinfo)<<endl;
	}

	char * contigDir = new char [1000];
	sprintf(contigDir, "%s/ContigsInfo", workDir.c_str());

	char * fileName = new char [1000];

	sprintf(fileName, "%s/ContigsInfo/contigBounds.store", workDir.c_str());

	FILE * pFile = fopen(fileName, "r");
	fseek(pFile, 0, SEEK_END); 
	int cSize = ftell(pFile)/sizeof(long long int);
	fseek(pFile, 0, SEEK_SET);

	long long int * contigBounds = new long long int[cSize];
	fread(contigBounds, sizeof(long long int), cSize, pFile);
	fclose(pFile);
	
	sprintf(fileName, "%s/Paths/paths", workDir.c_str());
	pFile = fopen(fileName, "r");
	fseek(pFile, 0, SEEK_END);

	int pSize = ftell(pFile)/sizeof(int);

	fseek(pFile, 0, SEEK_SET);

	int * paths = new int [pSize];
	fread(paths, sizeof(int), pSize, pFile);

	fclose(pFile);

	int numContigsLoaded = countFiles(contigDir) - 2;
	Fragment_Index ** contigs = new Fragment_Index * [numContigsLoaded]; 

	for(int i = 0; i < numContigsLoaded; i++)
	{
		contigs[i] = new Fragment_Index(); 
		sprintf(fileName, "%s/ContigsInfo/Contigs.store.%d", workDir.c_str(), i);

		pFile = fopen(fileName, "r");

		unsigned long long int fragmentOffset;
		fread(&fragmentOffset, sizeof(unsigned long long int), 1, pFile);

		unsigned long long int indexOffset;
		fread(&indexOffset, sizeof(unsigned long long int), 1, pFile);

		unsigned long long int length;
		fread(&length, sizeof(unsigned long long int), 1, pFile);

		unsigned long long int numFragments;
		fread(&numFragments, sizeof(unsigned long long int), 1, pFile);

		fseek(pFile, 0, SEEK_SET);

		contigs[i]->reserve(numFragments, length, indexOffset, fragmentOffset);

		contigs[i]->read(pFile);
		fclose(pFile);

		contigs[i]->set();
	}

	int ** matrix = new int * [200];

	for(int j = 0; j < 200; j++)
		matrix[j] = new int [200];

	int currPos = 0; int currCount = 0;
	while(currPos < pSize)
	{
		int pathLen = paths[currPos];
		cout<<">Contig "<<currCount++<<endl;
		for(int i = 0; i < pathLen-1; i++)
		{
			int currContig = paths[currPos + 1 + i];
			long long int contigPos1 = contigBounds[currContig];

			Fragment_Index * contig1;
			for(int j = 0; j < numContigsLoaded; j++)
			{
				if(contigPos1 >= contigs[j]->getFragmentOffset() && contigPos1 < contigs[j]->getFragmentOffset() + contigs[j]->length())	
				{
					contig1 = contigs[j];
					break;
				}
			}

			int nextContig = paths[currPos + 2 + i];
			long long int contigPos2 = contigBounds[nextContig];

			Fragment_Index * contig2;
			for(int j = 0; j < numContigsLoaded; j++)
			{
				if(contigPos2 >= contigs[j]->getFragmentOffset() && contigPos2 < contigs[j]->getFragmentOffset() + contigs[j]->length())
				{
					contig2 = contigs[j];
					break;
				}
			}

			long long int beginX = contigBounds[currContig]; long long int endX = 0;

			if(currContig + 1  >= cSize)
			{
				endX = contig1->length() + contig1->getFragmentOffset();
			}else{	
				endX = contigBounds[currContig + 1];
			}

			long long int beginY = contigBounds[nextContig]; long long int endY = 0; 

			if(nextContig + 1  >= cSize)
			{
				endY = contig2->length() + contig2->getFragmentOffset();
			}else{
				endY = contigBounds[nextContig + 1];
			}

			long long int bSecY = beginY; long long int eSecY = beginY + 50;
			if(eSecY > endY) eSecY = endY;

			for(int j = 0; j < endX-beginX-50; j+=50)
			{
				long long int bSecX = beginX+j; 
				long long int eSecX = beginX+j+100;

				if(eSecX > endX) eSecX = endX;

				Needleman_Wunsch<Fragment_Index, string>  needle(eSecX-bSecX, eSecY-bSecY, bSecX-contig1->getFragmentOffset(), bSecY-contig2->getFragmentOffset());	
				needle.Align(matrix, *contig1, *contig2, 0, 0);

				float misAligned = needle.misAligned(); float length = needle.alignLen();
				if(misAligned/length <= .15)
				{
					int alignStart = needle.xAlignmentStart();		
					for(int k = 0; k < (bSecX-beginX)+alignStart; k++)
						cout<<contig1->at(beginX-contig1->getFragmentOffset()+k);

					break;
				}
			}
		}

		int currContig = paths[currPos + pathLen];
		long long int begin = contigBounds[currContig]; long long int end = 0;

		Fragment_Index * contig;
		for(int j = 0; j < numContigsLoaded; j++)
		{
			if(begin >= contigs[j]->getFragmentOffset() && begin < contigs[j]->getFragmentOffset() + contigs[j]->length())	
			{
				contig = contigs[j];
				break;
			}
		}


		if(currContig + 1  > cSize)
		{
			end = contig->length() + contig->getFragmentOffset();
		}else{
			end = contigBounds[currContig + 1];
		}	


		for(int i = 0; i < end-begin; i++)
			cout<<contig->at(begin-contig->getFragmentOffset()+i);

		cout<<endl;

		currPos+=(pathLen+1);
	}

	for(int i = 0; i <  numContigsLoaded; i++)
	{
		delete contigs [i];
	}

	delete [] contigs;

	for(int i = 0; i < 200; i++)
		delete [] matrix [i];

	delete [] matrix;

	if(reportTime == true)
	{
		time_t rawtime;
		struct tm * timeinfo;
		time(&rawtime);
		timeinfo = localtime(&rawtime);
		cerr<<"end time: "<<asctime(timeinfo)<<endl;
	}

	return 0;
}

//void countFiles(const char [])
//Description: This function counts files in a directory
//Input: dir (cosnt char []) : The directory where we are counting the files
//Output:None
//Return:None
int countFiles(const char dir [])
{
	DIR *pDir = opendir(dir);
	struct dirent * pDirent;

	//Count Files
	int numFiles = 0;
	while ((pDirent = readdir(pDir)) != NULL)
		if(strcmp(pDirent->d_name, ".") != 0 && strcmp(pDirent->d_name, "..") != 0)
			numFiles++;

	closedir(pDir);

	return numFiles;
}
