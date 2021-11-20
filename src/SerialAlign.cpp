/*
 * SerialAlign.cpp
 *
 *      Author: Julia
 */

#include<iostream>
#include<vector>
#include<stdint.h>
#include<climits>
#include<cstring>
using namespace std;
#include <dirent.h>
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
	cerr<<"\t\t\t Serial Fragment Aligner"<<endl;
	cerr<<"\t\t\tWritten by Julia Warnke\n"<<endl;
	cerr<<"Command line arguments\n"<<endl;
	cerr<<"--workDir :working directory :default dat"<<endl;
	cerr<<"--minIden :the minimum identity needed to consider an overlap valid"<<endl;
	cerr<<"--minOverlap :the minimum overlap needed to consider an overlap valid"<<endl;
	cerr<<"--kmerSize :the size of the kmers used for seeding overlaps"<<endl;
	cerr<<"--stepSize :the distance between the starting positions of the kmers used to seed alignments"<<endl;
	cerr<<"--maxMem :the maximum memory that is available to the aligner"<<endl;
	cerr<<"--verbose : verbose mode"<<endl;
	cerr<<"--trimLen: the length to be trimmed off of the ends of reads if necessary"<<endl;
	cerr<<"\n\n";

	cerr<<"\n\n"<<endl;

	cerr<<"Exiting program"<<endl;
	exit(help);
}

int main(int argc, char * argv [])
{
	//Param options and defaults
	string workDir = "NO_WORK_DIR";
	int minIden = 90;
	int minOverlap = 50;
	int kmerSize = 22;
	int stepSize = 11;
	float maxMem = 1024;
	bool verbose = false;
	int trimLen = 0;

	//user options
	string workDirS = "--workDir";
	string minIdenS = "--minIden";
	string minOverlapS = "--minOverlap";
	string kmerSizeS =  "--kmerSize";
	string stepSizeS = "--stepSize";
	string maxMemS = "--maxMem";
	string verboseS = "--verbose";
	string trimLenS = "--trimLen";

	//Getting the command line options
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

		//The length to trim on reads
		if(argv[i] == trimLenS){
			trimLen = atoi(argv[i+1]);
		}

		if(argv[i] == verboseS)
		{
			verbose = true;
			i--;
		}
	}

	struct dirent *pDirent;
	DIR *pDir;
	int numSets = 0;

	//count the subsets
	char * cmdStr = new char[1000];
	sprintf(cmdStr, "%s/Fragments", workDir.c_str());

	pDir = opendir(cmdStr);
	if(pDir != NULL)
	{
		while ((pDirent = readdir(pDir)) != NULL)
		{
			if(strcmp(pDirent->d_name, ".") !=0 && strcmp(pDirent->d_name, "..") != 0)
			{
				numSets++;
			}
		}

		closedir (pDir);

	}else{

		cerr<<"Could not open work directory"<<endl;
		help_message();
	}

	string appLoc = argv[0];
	int endLoc = appLoc.find_last_of('/');
	appLoc.erase(endLoc);

	//Make the overlaps directory
	sprintf(cmdStr, "%s/Overlaps", workDir.c_str());
	pDir = opendir(cmdStr);
	if(pDir != NULL)
	{
		closedir (pDir);
		sprintf(cmdStr, "rm -r %s/Overlaps", workDir.c_str());
		system(cmdStr);
	}
	sprintf(cmdStr, "mkdir %s/Overlaps", workDir.c_str());
	system(cmdStr);

	if(!verbose)
		sprintf(cmdStr, "%s/align --workDir %s --queryJob 0 --numSets %d --refJob 0 --minIden %d --minOverlap %d --kmerSize %d --stepSize %d --maxMem %f --trimLen %d", appLoc.c_str(), workDir.c_str(), numSets, minIden, minOverlap, kmerSize, stepSize, maxMem, trimLen);
	else
		sprintf(cmdStr, "%s/align --workDir %s --queryJob 0 --numSets %d --refJob 0 --minIden %d --minOverlap %d --kmerSize %d --stepSize %d --maxMem %f --trimLen %d --verbose", appLoc.c_str(), workDir.c_str(), numSets, minIden, minOverlap, kmerSize, stepSize, maxMem, trimLen);

	//Run the program
	system(cmdStr);

	delete [] cmdStr;
	return 0;
}
