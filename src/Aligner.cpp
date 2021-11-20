/*
 * Aligner.cpp
 *
 *      Author: Julia
 */

#include<iostream>
#include<vector>
#include<stdint.h>
#include<climits>
#include<cstring>
#include<dirent.h>
using namespace std;
#include "Dictionary.h"
#include "Fragment_Index.h"
#include "MappingValues.h"
#include "SuffixArray.h"
#include "Seed.h"
#include "Chain.h"
#include "BandedAlignment.h"
#include "Overlapper.h"


//The help message
void help_message(){

	cerr<<"\n\n";
	cerr<<"\t\t\t Fragment aligner"<<endl;
	cerr<<"\t\t\tWritten by Julia Warnke\n"<<endl;
	cerr<<"Command line arguments\n"<<endl;
	cerr<<"--workDir :working directory :default dat"<<endl;
	cerr<<"--queryJob : the index of the query job that will ran by the aligner"<<endl;
	cerr<<"--refJob : the index of the reference job that will ran by the aligner"<<endl;
	cerr<<"--minIden :the minimum identity needed to consider an overlap valid"<<endl;
	cerr<<"--minOverlap :the minimum overlap needed to consider an overlap valid"<<endl;
	cerr<<"--kmerSize :the size of the kmers used for seeding overlaps"<<endl;
	cerr<<"--stepSize :the distance between the starting positions of the kmers used to seed alignments"<<endl;
	cerr<<"--numSets : the number of sets per job"<<endl;
	cerr<<"--maxMem :the maximum memory that is available to the aligner"<<endl;
	cerr<<"--verbose : verbose mode"<<endl;
	cerr<<"\n\n";

	cerr<<"\n\n"<<endl;

	cerr<<"Exiting program"<<endl;
	exit(help);
}


//bool fileExists( char fileName [] )
//Description: This function checks to see if a file exists
//Input: fileName (char []) the name of the file
//Output:None
//Return:None
bool fileExists(const char fileName []) {

	FILE * pFile = fopen(fileName, "r");
	if(pFile != '\0')
	{
		fclose(pFile);
		return true;
	}
	return false;
}

//readOffsets(FILE * &, unsigned long long int &, unsigned long long int &)
//Description: This function reads the offsets of a fragment data set
//Input:None
//Output:None
//Return:None
void readOffsets(FILE * & pFile, unsigned long long int & fragmentOffset, unsigned long long int & indexOffset)
{
	fread(&fragmentOffset, sizeof(unsigned long long int), 1, pFile);
	fread(&indexOffset, sizeof(unsigned long long int), 1, pFile);
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

//readLengths(FILE * &, unsigned long long int &, unsigned long long int &)
//Description: This function reads in the lengths of the fragment data set
//Input:None
//Output:None
//Return:None
void readLengths(FILE * & pFile, unsigned long long int & totalFragmentLength, unsigned long long int & numFragments)
{
	fseek(pFile, 2*sizeof(unsigned long long int), SEEK_SET);

	unsigned long long int intervalLength;
	fread(&intervalLength, sizeof(unsigned long long int), 1, pFile);

	totalFragmentLength+=intervalLength;

	unsigned long long int numIntervalFragments;
	fread(&numIntervalFragments, sizeof(unsigned long long int), 1, pFile);
	numFragments+=numIntervalFragments;

}

int main(int argc, char * argv [])
{
	//Param options and defaults
	string workDir = "NO_WORK_DIR";
	int queryJob = -1;
	int refJob = -1;
	int numSets = -1;
	int minIden = 90;
	int minOverlap = 50;
	int kmerSize = 22;
	int stepSize = 11;
	float maxMem = 1024;
	bool verbose = false;
	bool containment = false;
	bool alignRepsOnly = false;

	//user options
	string workDirS = "--workDir";
	string queryJobS = "--queryJob";
	string refJobS = "--refJob";
	string numSetsS = "--numSets";
	string minIdenS = "--minIden";
	string minOverlapS = "--minOverlap";
	string kmerSizeS =  "--kmerSize";
	string stepSizeS = "--stepSize";
	string maxMemS = "--maxMem";
	string verboseS = "--verbose";
	string containmentS = "--containment";
	string alignRepsOnlyS = "--alignRepsOnly";

	//Getting the commandline options
	for(int i = 1; i < argc; i+=2)
	{
		//The working directory
		if(argv[i] == workDirS)
		{
			workDir = argv[i+1];
			unsigned int pos = workDir.find_last_of("/");
			if(pos == workDir.size()-1)
				workDir.erase(pos);

			if(workDir.size() > 200)
			{
				cout<<"Please keep the working directory name less than 200 characters in length"<<endl;
				exit(bufOverflow);
			}
		}

		//The query job that is to be run
		if(argv[i] == queryJobS)
		{
			queryJob = atoi(argv[i+1]);
		}

		//The reference job that is to be run
		if(argv[i] == refJobS)
		{
			refJob = atoi(argv[i+1]);
		}

		//The sets of sequences allocated to each job
		if(argv[i] == numSetsS)
		{
			numSets = atoi(argv[i+1]);
		}

		//The minimum identity
		if(argv[i] == minIdenS)
		{
			minIden = atoi(argv[i+1]);
		}

		//The minimum overlap
		if(argv[i] == minOverlapS){

			minOverlap = atoi(argv[i+1]);
		}

		//The kmerSize
		if(argv[i] == kmerSizeS)
		{
			kmerSize = atoi(argv[i+1]);
		}

		//The kmer step size
		if(argv[i] == stepSizeS){

			stepSize = atoi(argv[i+1]);
		}

		//The maximum memory available to the system
		if(argv[i] == maxMemS){
			maxMem = atoi(argv[i+1]);
		}

		if(argv[i] == verboseS)
		{
			verbose = true;
			i--;
		}

		if(argv[i] == containmentS)
		{
			containment = true;
			i--;
		}

		if(argv[i] == alignRepsOnlyS)
		{
			alignRepsOnly = true;
			i--;
		}
		
	}

	if(queryJob == -1 || workDir == "NO_WORK_DIR" || refJob == -1)
	{
		help_message();
	}

	if(verbose)
	{
		//Log file stuff
		cerr<<"\n\n\t\tFragment Aligner Parameters"<<endl;
		cerr<<endl;
		cerr<<endl;
		cerr<<"Parameters: "<<endl;
		cerr<<"Working directory: "<<workDir<<endl;
		cerr<<"Query job: "<<queryJob<<endl;
		cerr<<"Reference job: "<<refJob<<endl;
		cerr<<"Number of sets per job: "<<numSets<<endl;
		cerr<<"Minimum identity: "<<minIden<<endl;
		cerr<<"Minimum overlap: "<<minOverlap<<endl;
		cerr<<"Kmer size: "<<kmerSize<<endl;
		cerr<<"Step size: "<<stepSize<<endl;
		cerr<<"Maximum memory availiable (mb): "<<maxMem<<endl;

		//Begin time of the assembly
		time_t rawtime;
		struct tm * timeinfo;
		time(&rawtime);
		timeinfo = localtime(&rawtime);
		cerr<<"Begin time: "<<asctime(timeinfo)<<endl;

	}

	//The number of bytes allocated to the various sets
	int DataSetAlloc = 10;

	//The allocation is proportional to the memory specified by the user
	DataSetAlloc = DataSetAlloc * maxMem/1024;


	int * cOvlq; int * cOvlr; int cOvlrSize = 0;
	{
		char * ovlInfoDir = new char [1000];
		char * fileName = new char [1000];

		sprintf(ovlInfoDir, "%s/OvlInfo/", workDir.c_str());
		sprintf(fileName, "%s/sequences.job%d.info", ovlInfoDir, queryJob);

		FILE * pFile = fopen(fileName, "r");
		fseek(pFile, 0, SEEK_END);
		int mSize = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);

		cOvlq = new int [mSize];
		fread(cOvlq, sizeof(int), mSize, pFile);
		fclose(pFile);

		if(queryJob != refJob)
		{
			sprintf(fileName, "%s/sequences.job%d.info", ovlInfoDir, refJob);

			pFile = fopen(fileName, "r");    
			fseek(pFile, 0, SEEK_END);
			mSize = ftell(pFile)/sizeof(int);
			fseek(pFile, 0, SEEK_SET);

			cOvlrSize = mSize;

			cOvlr = new int [mSize];
			fread(cOvlr, sizeof(int), mSize, pFile);
			fclose(pFile);

		}else{

			cOvlr = cOvlq;
			cOvlrSize = mSize;
		}
	
		delete [] fileName;
		delete [] ovlInfoDir;
	}

	cout<<"START "<<queryJob<<" "<<refJob<<endl;


	int rCurrSet = 0;
	while(rCurrSet < numSets)
	{
		//IN THIS SECTION THE PROGRAM IS COLLECTING THE REFERECE SEQUENCES
		char fileName [300];
		sprintf(fileName, "%s/Fragments/sequences.job%d.set%d", workDir.c_str(), refJob, rCurrSet);
		FILE * pFile  = fopen(fileName, "r");

		unsigned long long int numFragments = 0; unsigned long long int totalFragmentLength = 0; unsigned long long int fragmentOffset = 0; unsigned long long int indexOffset = 0;
		readOffsets(pFile, fragmentOffset, indexOffset);
		fclose(pFile);

		int tmpSet = rCurrSet;
		for(; tmpSet < numSets && totalFragmentLength/(megabyte*DataSetAlloc) < 1; tmpSet++)
		{
			sprintf(fileName, "%s/Fragments/sequences.job%d.set%d", workDir.c_str(), refJob, tmpSet);
			pFile  = fopen(fileName, "r");
			readLengths(pFile, totalFragmentLength, numFragments);
			fclose(pFile);
		}

		Fragment_Index rIndex;
		rIndex.reserve(numFragments, totalFragmentLength, indexOffset, fragmentOffset);

		for(; rCurrSet < tmpSet; rCurrSet++)
		{
			sprintf(fileName, "%s/Fragments/sequences.job%d.set%d", workDir.c_str(), refJob, rCurrSet);
			pFile = fopen(fileName, "r");
			if(verbose)
				cerr<<"READING REFERENCE FILE: "<<fileName<<endl;
			rIndex.read(pFile);
			fclose(pFile);
		}

		rIndex.set();

		//pair<int, int> bounds = rIndex.indexBounds(rIndex.numFragments()-1);
		//DONE COLLECTING THE REFERENCE SEQUENCES

		if(verbose)
			cout<<"\nBUILDING  SUFFIX ARRAY"<<endl;

		SuffixArray<Fragment_Index> suffix;
		if(queryJob == 0)
		{

			char * arr = new char [1000];
			sprintf(arr, "%s/Arrays/array%d_%d", workDir.c_str(), refJob, rCurrSet);

			if(!fileExists(arr))
			{
				suffix.buildArray(rIndex, 0, rIndex.length());
				FILE * pFile = fopen(arr, "w");
				suffix.write(pFile);
				fclose(pFile);
			}else{

				FILE * pFile = fopen(arr, "r");

				suffix.read(pFile);
				fclose(pFile);
				delete [] arr;
			}

		}else{
			char * arr = new char [1000];
			sprintf(arr, "%s/Arrays/array%d_%d", workDir.c_str(), refJob, rCurrSet);
			FILE * pFile = fopen(arr, "r");

			suffix.read(pFile);
			fclose(pFile);
			delete [] arr;
		}

		FM_Index<Fragment_Index> fmIndex;
		fmIndex.calcBWT(suffix, rIndex, rIndex.length());

		if(verbose)
			cout<<"DONE BUILDING SUFFIX ARRAY\n"<<endl;

		//IN THIS SECTION THE PROGRAM IS COLLECTING THE QUERY SEQUENCES
		int qCurrSet = 0;
		while(qCurrSet < numSets)
		{
			char fileName [300];
			sprintf(fileName, "%s/Fragments/sequences.job%d.set%d", workDir.c_str(), queryJob, qCurrSet);
			FILE * pFile  = fopen(fileName, "r");

			numFragments = 0; totalFragmentLength = 0; fragmentOffset = 0; indexOffset = 0;
			readOffsets(pFile, fragmentOffset, indexOffset);
			fclose(pFile);

			int tmpSet = qCurrSet;
			for(; tmpSet < numSets && totalFragmentLength/(megabyte*DataSetAlloc) < 1; tmpSet++)
			{
				sprintf(fileName, "%s/Fragments/sequences.job%d.set%d", workDir.c_str(), queryJob, tmpSet);
				pFile  = fopen(fileName, "r");
				readLengths(pFile, totalFragmentLength, numFragments);
				fclose(pFile);
			}

			Fragment_Index qIndex;
			qIndex.reserve(numFragments, totalFragmentLength, indexOffset, fragmentOffset);

			for(; qCurrSet < tmpSet; qCurrSet++)
			{
				sprintf(fileName, "%s/Fragments/sequences.job%d.set%d", workDir.c_str(), queryJob, qCurrSet);
				pFile = fopen(fileName, "r");
				if(verbose)
					cerr<<"READING QUERY FILE: "<<fileName<<endl;
				qIndex.read(pFile);
				fclose(pFile);
			}

			qIndex.set();

			char * overlapFile = new char [500];
			sprintf(overlapFile, "%s/Overlaps/ovl%d_%d", workDir.c_str(), queryJob, refJob /*qCurrSet, rCurrSet*/);
			Overlapper<Fragment_Index> overlap(suffix, fmIndex, qIndex, rIndex, cOvlq, cOvlr, kmerSize, stepSize, minIden, minOverlap, containment, alignRepsOnly);

			overlap.findOverlaps(overlapFile, verbose);
			delete [] overlapFile;
		}
	}

	if(containment)
	{
		char * ovlInfoDir = new char [1000];
		char * fileName = new char [1000];

		sprintf(ovlInfoDir, "%s/OvlInfo/", workDir.c_str());
		sprintf(fileName, "%s/sequences.job%d.info", ovlInfoDir, refJob);

		FILE * pFile = fopen(fileName, "w");

		fwrite(cOvlr, sizeof(int), cOvlrSize, pFile);
		fclose(pFile);

		delete [] ovlInfoDir;
		delete [] fileName;
	}

	delete [] cOvlq;

	if(queryJob != refJob)
		delete [] cOvlr;

	if(verbose)
	{
		//Get the time
		time_t rawtime;
		struct tm * timeinfo;
		time(&rawtime);
		timeinfo = localtime(&rawtime);
		cerr<<"\nEnd time: "<<asctime(timeinfo)<<endl;
		cerr<<endl;
		cerr<<endl;
	}

	cout<<" END"<<queryJob<<" "<<refJob<<endl;

	return 0;
}

