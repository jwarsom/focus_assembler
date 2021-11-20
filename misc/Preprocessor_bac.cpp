/*
 * Preprocessor.cpp
 *
 *  Created on: Aug 3, 2011
 *      Author: Julia
 */

#include<iostream>
#include<vector>
#include<dirent.h>
#include<time.h>
#include<fstream>
#include<stdint.h>
#include<cassert>
#include<cstdlib>
#include<cstring>
#include<math.h>
#include<string>
using namespace std;
#include "MappingValues.h"
#include "Dictionary.h"
#include "Fragment_Index.h"
#include "Labels.h"
#include "SortFragments.h"
#include "MergeFragments.h"

//Exit signals
#define help 911
#define bufOverflow 792
#define inputProblem 716

void help_message(){

	cerr<<"\n\n";
	cerr<<"\t\t\t Fragment preprocessor"<<endl;
	cerr<<"\t\t\tWritten by Julia Warnke\n"<<endl;
	cerr<<"Command line arguments\n"<<endl;
	cerr<<"--workDir :working directory :default dat"<<endl;
	cerr<<"--numJobs :this divides the reads into sets of jobs that can be distributed among different processors"<<endl;
	cerr<<"--singleReads : input file for single reads (fastq or fasta)."<<endl;
	cerr<<"--seqTagLen : length of fragment name retained :default 10 characters"<<endl;
	cerr<<"--minReadLen : minimum read length : default 75 characters"<<endl;
	cerr<<"--Trim5 : Trim at the 5' end : default 0 characters"<<endl;
	cerr<<"--Trim3 : Trim at the 3' end : default 0 characters"<<endl; 
	cerr<<"\n\n";

	cerr<<"If quality values are provided they must be named as x.qual, where x is the full path name of the single read file"<<endl;
	cerr<<"\n\n"<<endl;

	cerr<<"Exiting program"<<endl;
	exit(help);
}

//Function declarations
void createDirectories(string &);
bool fileExists(const char []);
int  countFiles(const char []);
void makeDir(const char []);
void moveDir(const char [], const char []);
void removeDir(const char []);

int main(int argc, char * argv [])
{

	//Param options and defaults
	string workDir = "Merge_and_Traverse_Assembly";
	unsigned int  numJobs = 1;
	string singleReads = "NO FILE";
	int seqTagLen = 10;
	int windowLen = 20;
	int minQualVal = 25;
	int minReadLen = 75;
	int trim5 = 0;
	int trim3 = 0; 

	//user options
	string workDirS = "--workDir";
	string numJobsS = "--numJobs";
	string singleReadsS = "--singleReads";
	string seqTagLenS =  "--seqTagLen";
	string windowLenS = "--trimWindowLen";
	string minQualValS  = "--minQualVal";
	string minReadLenS = "--minReadLen";
	string trim5S = "--trim5";
	string trim3S = "--trim3";

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

		//The number of jobs
		if(argv[i] == numJobsS)
		{
			numJobs = atoi(argv[i+1]);
		}

		//The single reads
		if(argv[i] == singleReadsS)
		{
			singleReads = argv[i+1];
		}

		//The tag length
		if(argv[i] == seqTagLenS){

			seqTagLen = atoi(argv[i+1]);
		}

		//The window length used for trimming
		if(argv[i] == windowLenS){
			windowLen = atoi(argv[i+1]);
		}

		//The minimum quality value for trimming
		if(argv[i] == minQualValS){
			minQualVal = atoi(argv[i+1]);
		}

		//The minimum read length
		if(argv[i] == minReadLenS){
			minReadLen = atoi(argv[i+1]);
		}

		//The minimum read length
		if(argv[i] == trim3S){
			trim3 = atoi(argv[i+1]);
		}

		//The minimum read length
		if(argv[i] == trim5S){
			trim5 = atoi(argv[i+1]);
		}
	}

	//User forgot a file
	if(argc == 1 || singleReads == "NO FILE")
	{
		help_message();
	}

	//Create directories
	createDirectories(workDir);

	//Log file stuff
	{
		cerr<<"\n\n\t\tFragment Preprocessor Parameters"<<endl;
		cerr<<endl;
		cerr<<endl;
		cerr<<"Parameters: "<<endl;
		cerr<<"Working directory: "<<workDir<<endl;
		cerr<<"Number of jobs: "<<numJobs<<endl;
		cerr<<"Input read file: "<<singleReads<<endl;
		cerr<<"Maximum length of sequence names "<<seqTagLen<<endl;
		cerr<<"Minimum read length "<<minReadLen<<endl;
		cerr<<"Window size for quality trimming "<<windowLen<<endl;
		cerr<<"Minimum quality value "<<minQualVal<<endl;

		char arr[300];
		sprintf(arr, "%s.qual", singleReads.c_str());
		FILE * pFile = fopen(arr, "r");
		cerr<<"Quality file provided:";
		if(pFile)
		{
			cerr<<" yes"<<endl;
			fclose(pFile);
		}else
		{
			cerr<<" no"<<endl;
		}

		//Begin time of the assembly
		time_t rawtime;
		struct tm * timeinfo;
		time(&rawtime);
		timeinfo = localtime(&rawtime);
		cerr<<"Begin time: "<<asctime(timeinfo)<<endl;
	}

	//The input file name
	char inReadFileName[300];
	char inQualFileName[300];
	sprintf(inReadFileName, "%s", singleReads.c_str());
	sprintf(inQualFileName, "%s.qual", singleReads.c_str());
	FILE * r_inFile = fopen(inReadFileName, "r");
	FILE * q_inFile = fopen(inQualFileName, "r");

	char fmtReadFileName[300];
	char fmtQualFileName[300]; 

	string inReadFileNameS(inReadFileName);

	int cutPoint = inReadFileNameS.find_last_of(".");
	sprintf(fmtReadFileName, "%s_formatted.fasta", inReadFileNameS.substr(0, cutPoint).c_str());
	sprintf(fmtQualFileName, "%s_formatted.fasta.qual", inReadFileNameS.substr(0, cutPoint).c_str());

	//Trimming and generating reverse complements  
	{
		if(!q_inFile)
		{
			FILE * r_formatted = fopen(fmtReadFileName, "w");
			char * inReads = new char [11 * megabyte + 1];
			char * outReads = new char [11 * megabyte];
			unsigned long long int counts = 0; bool reachedEnd = false;
			unsigned long long int currPos = 0;

			do{
				counts = fread(inReads, sizeof(char), 11*megabyte, r_inFile);
				inReads[counts] = '\0';

				char * curr = inReads;
				char * next = strchr(curr+1, '>');

				while(next != '\0')
				{
					if(currPos + 2*(next - curr)+4 > 11 * megabyte)
					{
						fwrite(outReads, sizeof(char), currPos, r_formatted);
						currPos = 0;
					}

					char * newline = strchr(curr, '\n');

					int stopPoint = (next-newline-2);
					while(stopPoint > 0 && newline[stopPoint] == 'N') stopPoint--;  

					if(stopPoint - (trim3 + trim5) >= minReadLen)
					{
						memcpy(outReads+currPos, curr, newline-curr);
						currPos+=(newline-curr);

						outReads[currPos] = '\t'; //CHANGED FROM \n
						currPos++;

						memcpy(outReads+currPos, (newline+1)+trim5, (stopPoint)-(trim3+trim5));
						currPos+=(stopPoint)-(trim3+trim5);

						outReads[currPos] = '\n';
						currPos++;

						memcpy(outReads+currPos, curr, newline-curr);
						currPos+=(newline-curr);

						outReads[currPos] = '_'; currPos++; 	
						outReads[currPos] = 'r'; currPos++; 	
						outReads[currPos] = 'e'; currPos++; 	
						outReads[currPos] = 'v'; currPos++; 	

						outReads[currPos] = '\t'; //CHANGED FROM \n
						currPos++;

						int startPoint = trim5;

						for(int i = stopPoint-trim3; i > startPoint; i--)
						{
							outReads[currPos + (stopPoint-i)] = inverseComplementMapVals[mapVals[newline[i]]]; 
						}	
						currPos+=(stopPoint)-(trim3+trim5);

						outReads[currPos] = '\n';
						currPos++;

					}

					curr = next; next = strchr(curr+1, '>');
					if(counts != (11*megabyte) && !reachedEnd && next == '\0')
					{
						next = strchr(curr, '\0');
						reachedEnd = true;
					}
				}


				char * end = strchr(curr, '\0');
				fseek(r_inFile, curr-end, SEEK_CUR);

			}while(counts == (11*megabyte)); //Read through until the end of the file

			if(currPos  > 0)
			{
				fwrite(outReads, sizeof(char), currPos, r_formatted);
				currPos = 0;
			}

			delete [] inReads; delete [] outReads;
			fclose(r_formatted);

		}else{
			FILE * r_formatted = fopen(fmtReadFileName, "w");
			FILE * q_formatted = fopen(fmtQualFileName, "w");

			unsigned long long int currPos = 0; unsigned long long int counts = 0;
			bool reachedEnd = false;
			char * inReads = new char [11 * megabyte + 1];
			char * outReads = new char [11 * megabyte];
			char * inQuals = new char [11 * megabyte + 1];
			char * outQuals = new char [11 * megabyte];

			do{
				counts = fread(inReads, sizeof(char), 11*megabyte, r_inFile);
				assert(counts == fread(inQuals, sizeof(char), 11*megabyte, q_inFile));
				inReads[counts] = '\0'; inQuals[counts] = '\0';

				char * curr = inReads;
				char * currQuals = inQuals;

				char * next = strchr(curr+1, '>');
				char * nextQuals = strchr(currQuals, '\n');
				nextQuals = (strchr(nextQuals+1, '\n')+1);

				while(next != '\0')
				{
					if(currPos + 2*(next - curr) + 4 > 11 * megabyte)
					{
						fwrite(outReads, sizeof(char), currPos, r_formatted);
						fwrite(outQuals, sizeof(char), currPos, q_formatted);
						currPos = 0;
					}

					char * newline = strchr(curr, '\n');
					char * newLineQuals = strchr(currQuals, '\n');

					int stopPoint = next-newline-2;
					while(stopPoint > 0 && newline[stopPoint] == 'N') stopPoint--;

					stopPoint = max(0, stopPoint - trim3);
					for(int i = stopPoint; (i - windowLen) > 0; i--)
					{
						float totalQualVals = 0; stopPoint = i;
						for(int j = 0; j < windowLen; j++)
						{
							totalQualVals+=(*(newLineQuals+(i-j)) - 33);
						}

						float avgQualVal = totalQualVals/windowLen;

						if(avgQualVal >= minQualVal)
						{
							break;
						}
					}				

					int startPoint = trim5;
					for(; startPoint + windowLen < stopPoint; startPoint++)
					{
						float totalQualVals = 0; 
						for(int j = 0; j < windowLen; j++)
						{
							totalQualVals+=(*(newLineQuals+(startPoint+j)) - 33);
						}

						float avgQualVal = totalQualVals/windowLen;

						if(avgQualVal >= minQualVal)
						{
							break;
						}
					}

					if(stopPoint-startPoint >= minReadLen)
					{
						memcpy(outReads+currPos, curr, newline-curr);
						memcpy(outQuals+currPos, curr, newline-curr);
						currPos+=(newline-curr);

						outReads[currPos] = '\t'; //CHANGED FROM \n
						outQuals[currPos] = '\t'; //CHANGED FROM \n
						currPos++;

						memcpy(outReads+currPos, (newline+1+startPoint), stopPoint-startPoint);
						memcpy(outQuals+currPos, (newLineQuals+1+startPoint), stopPoint-startPoint);
						currPos+=(stopPoint - startPoint);

						outReads[currPos] = '\n';
						outQuals[currPos] = '\n';
						currPos++;

						memcpy(outReads+currPos, curr, newline-curr);
						memcpy(outQuals+currPos, curr, newline-curr);
						currPos+=(newline-curr);

						outReads[currPos] = '_'; outQuals[currPos] = '_'; currPos++; 	
						outReads[currPos] = 'r'; outQuals[currPos] = 'r'; currPos++; 	
						outReads[currPos] = 'e'; outQuals[currPos] = 'e'; currPos++; 	
						outReads[currPos] = 'v'; outQuals[currPos] = 'v'; currPos++; 	

						outReads[currPos] = '\t'; //CHANGED FROM \n
						outQuals[currPos] = '\t'; //CHANGED FROM \n
						currPos++;

						for(int i = stopPoint; i > startPoint; i--)
						{
							outReads[currPos + (stopPoint-i)] = inverseComplementMapVals[mapVals[newline[i]]]; 
							outQuals[currPos + (stopPoint-i)] = newLineQuals[i]; 
						}	
						currPos+=(stopPoint - startPoint);

						outReads[currPos] = '\n';
						outQuals[currPos] = '\n';
						currPos++;

					}

					curr = next; next = strchr(curr+1, '>');
					currQuals = nextQuals; nextQuals = (strchr(currQuals, '\n'));

					if(nextQuals != '\0')
					{ 
						nextQuals = (strchr(nextQuals+1, '\n')); 
						if(nextQuals != '\0')
						{
							nextQuals = nextQuals+1;
						}
					}

					if(reachedEnd) break;

					if(counts != (11*megabyte) && !reachedEnd && next == '\0')
					{
						next = strchr(curr, '\0'); nextQuals = strchr(currQuals, '\0');
						reachedEnd = true;
					}
				}

				char * end = strchr(curr, '\0');
				fseek(r_inFile, curr-end, SEEK_CUR);

				end = strchr(currQuals, '\0');
				fseek(q_inFile, currQuals-end, SEEK_CUR);


			}while(counts == (11*megabyte)); //Read through until the end of the file

			if(currPos > 0)
			{
				fwrite(outReads, sizeof(char), currPos, r_formatted);
				fwrite(outQuals, sizeof(char), currPos, q_formatted);
				currPos = 0;
			}

			delete [] inReads; delete [] outReads;
			delete [] inQuals; delete [] outQuals;
			fclose(r_formatted); fclose(q_formatted);
		}	
	}

	fclose(r_inFile); fclose(q_inFile);
	cout<<"HERE "<<endl;

	char * cmd = new char[1000];
	sprintf(cmd, "awk '{print length($NF), int(NR), $0}' %s | sort -r -n -k 1,1 -k 2,2 | cut -d\" \" -f3- | awk -F \"\\t\" '{gsub(/^[ \\t]+|[ \\t]+$/, \"\", $1); gsub(/^[ \\t]+|[ \\t]+$/, \"\", $2); print $1\"\\n\"$2}' > %s.tmp" , fmtReadFileName, fmtReadFileName);
	system(cmd);

	sprintf(cmd, "mv %s.tmp %s", fmtReadFileName, fmtReadFileName);
	system(cmd);

	q_inFile = fopen(fmtQualFileName, "r");

	if(q_inFile)
	{
		sprintf(cmd, "awk '{print length($NF), int(NR), $0}' %s | sort -r -n -k 1,1 -k 2,2 | cut -d\" \" -f3- | awk -F \"\\t\" '{gsub(/^[ \\t]+|[ \\t]+$/, \"\", $1); gsub(/^[ \\t]+|[ \\t]+$/, \"\", $2); print $1\"\\n\"$2}' > %s.tmp" , fmtQualFileName, fmtQualFileName);
		system(cmd);

		sprintf(cmd, "mv %s.tmp %s", fmtQualFileName, fmtQualFileName);
		system(cmd);
	}

	fclose(q_inFile);
	delete cmd;

	r_inFile = fopen(fmtReadFileName, "r");
	q_inFile = fopen(fmtQualFileName, "r");

	//Seeing how big the files are
	unsigned long long int datasetSize = 0;
	{
		FILE * pFile = fopen(fmtReadFileName, "r");
		fseek(pFile, 0, SEEK_END);
		datasetSize+=ftell(pFile);
		fclose(pFile);
	}

	//Keep the subsets under 10 megabytes
	unsigned long long int subSetSize = (datasetSize+numJobs-1)/numJobs;
	unsigned int numSubsets = numJobs;
	while(subSetSize > megabyte * 10)
	{
		subSetSize = (subSetSize+1)/2;
		numSubsets = numSubsets * 2;
	}

	//Divide up the files!
	{
		//Tracking the current job that we are on
		unsigned int currJob = 0;
		unsigned int currSet = 0;

		//The input buffers
		char * inReads = new char[11 * megabyte+1];
		char * inQuals = new char[11 * megabyte+1];

		//We are going to tally the number of base pairs and fragments per subset
		unsigned long int numDelims = 0;
		unsigned long long int numCharacters = 0;
		unsigned long long int numBasePairs = 0;
		unsigned long int runningTotalFragments = 0;
		unsigned long long int runningTotalBasePairs = 0;

		if(!q_inFile) //No Quals
		{
			unsigned long long int counts = 0;
			counts = fread(inReads, sizeof(char), 11*megabyte, r_inFile);
			inReads[counts] = '\0';
			char * curr = inReads;

			unsigned int numChunks = min((unsigned long long int) (10*megabyte), counts)/subSetSize;
			if(numChunks == 0) numChunks = 1;

			for(unsigned int i = 0; i < numChunks; i++)
			{
				char outIndexFileName[300];
				char outLabelFileName[300];
				sprintf(outIndexFileName, "%s/Fragments/sequences.job%d.set%d", workDir.c_str(),  currJob, currSet);
				sprintf(outLabelFileName, "%s/Labels/labels.job%d.set%d", workDir.c_str(),  currJob, currSet);
				FILE * oFile = fopen(outIndexFileName, "w");
				FILE * oFile2 = fopen(outLabelFileName, "w");

				unsigned int currentSet = currJob * numSubsets/numJobs + currSet;

				char * placeHolder = curr;
				char * next = strchr(curr+1, '>');

				while(next != '\0' && (numCharacters + currentSet%2 * (next - curr) < subSetSize || (currJob == numJobs-1 && currSet == numSubsets/numJobs-1)))
				{
					numCharacters+=(next-curr);  numDelims++;
					char * newLine = strchr(curr, '\n');
					numBasePairs+=(next - newLine - 2);
					curr = next; next = strchr(curr+1, '>');
				}

				if(next == '\0' && counts != (11*megabyte))
				{
					next = strchr(curr, '\0');
					numCharacters+=(next-curr);  numDelims++;
					char * newLine = strchr(curr, '\n');
					numBasePairs+=(next - newLine - 2);
				}

				Fragment_Index index;
				index.reserve(numDelims, numBasePairs, runningTotalFragments, runningTotalBasePairs);
				curr = placeHolder;

				Labels label;
				label.reserve(seqTagLen, numDelims, runningTotalFragments);

				for(unsigned long int j = 0; j < numDelims; j++)
				{
					label.append(curr+1);
					next = strchr(curr+1, '>');
					if(next == '\0')
						next = strchr(curr, '\0');

					char * newLine = strchr(curr, '\n');
					char * tmpQuals =  new char[next-newLine-2];
					memset(tmpQuals, 73, next-newLine-2);
					index.append(newLine+1, tmpQuals, next-newLine-2);

					delete [] tmpQuals;
					curr = next;
				}

				index.set();
				index.write(oFile); currSet++;
				fclose(oFile);

				label.write(oFile2);
				fclose(oFile2);

				if(currSet == numSubsets/numJobs)
				{
					currSet = 0; currJob++;
				}

				runningTotalFragments+= numDelims;
				runningTotalBasePairs+= numBasePairs;
				numCharacters = 0; numDelims = 0; numBasePairs = 0;
			}

			char * end = strchr(curr, '\0');
			fseek(r_inFile, curr-end, SEEK_CUR);


		}while(counts == (11*megabyte) || currJob < numJobs); //Read through until the end of the file

	}else{
		unsigned long long int counts = 0;
		do{
			counts = fread(inReads, sizeof(char), 11*megabyte, r_inFile);
			inReads[counts] = '\0';
			assert(counts == fread(inQuals, sizeof(char), 11*megabyte, q_inFile));
			inQuals[counts] = '\0';

			char * curr = inReads;
			char * currQuals = inQuals;

			unsigned int numChunks = min((unsigned long long int) (10*megabyte), counts)/subSetSize;
			if(numChunks == 0) numChunks = 1;


			for(unsigned int i = 0; i < numChunks; i++)
			{
				char outIndexFileName[300];
				char outLabelFileName[300];
				sprintf(outIndexFileName, "%s/Fragments/sequences.job%d.set%d", workDir.c_str(),  currJob, currSet);
				sprintf(outLabelFileName, "%s/Labels/labels.job%d.set%d", workDir.c_str(),  currJob, currSet);
				FILE * oFile = fopen(outIndexFileName, "w");
				FILE * oFile2 = fopen(outLabelFileName, "w");

				unsigned int currentSet = currJob * numSubsets/numJobs + currSet;

				char * placeHolder = curr;
				char * placeHolderQuals = currQuals;				

				char * next = strchr(curr+1, '>');
				char * nextQuals = strchr(currQuals+1, '\n');
				nextQuals = (strchr(nextQuals+1, '\n')+1);


				while(next != '\0' && (numCharacters + currentSet%2 * (next - curr) < subSetSize || (currJob == numJobs-1 && currSet == numSubsets/numJobs-1)))
				{
					numCharacters+=(next-curr);  numDelims++;
					char * newLine = strchr(curr, '\n');
					numBasePairs+=(next - newLine - 2);
					curr = next; next = strchr(curr+1, '>');
				}

				if(next == '\0' && counts != (11 * megabyte))
				{
					next = strchr(curr, '\0');
					numCharacters+=(next-curr);  numDelims++;
					char * newLine = strchr(curr, '\n');
					numBasePairs+=(next - newLine - 2);
				}

				Fragment_Index index;
				index.reserve(numDelims, numCharacters, runningTotalFragments, runningTotalBasePairs);
				curr = placeHolder; currQuals = placeHolderQuals;


				Labels label;
				label.reserve(seqTagLen, numDelims, runningTotalFragments);


				for(unsigned long int j = 0; j < numDelims; j++)
				{
					label.append(curr+1);

					char * newLine = strchr(curr, '\n');
					char * newLineQuals = strchr(currQuals, '\n');

					next = strchr(curr+1, '>');
					if(next == '\0')
						next = strchr(curr, '\0');

					nextQuals = strchr(newLineQuals+1, '\n');
					if(nextQuals != '\0')
						nextQuals = nextQuals+1;


					index.append(newLine+1, newLineQuals+1, next-newLine-2);

					curr = next;
					currQuals = nextQuals;
				}


				index.set();
				index.write(oFile); currSet++;
				fclose(oFile);

				label.write(oFile2);
				fclose(oFile2);

				if(currSet == numSubsets/numJobs)
				{
					currSet = 0; currJob++;
				}

				runningTotalFragments+= numDelims;
				runningTotalBasePairs+= numBasePairs;
				numCharacters = 0; numDelims = 0; numBasePairs = 0;

			}

			char * end = strchr(curr, '\0');
			fseek(r_inFile, curr-end, SEEK_CUR);
			fseek(q_inFile, curr-end, SEEK_CUR);


		}while(counts == (11*megabyte) || currJob < numJobs); //Read through until the end of the file

	}

	//clean up
	delete [] inReads;
	delete [] inQuals;
}



//Log the end time
{
	//Get the time
	time_t rawtime;
	struct tm * timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	cerr<<"End time: "<<asctime(timeinfo)<<endl;
	cerr<<endl;
	cerr<<endl;
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

//void makeDir(const char [])
//Description: This function calls the linux mkdir command
//Input: dirName (const char [])
//Output:None
//Return:None
void makeDir(const char dirName []) {

	char * sysCmd = new char [1000];
	sprintf(sysCmd, "mkdir %s", dirName);
	system(sysCmd);
	delete [] sysCmd;
}

//void moveDir(const char [], const char [])
//Description: This function calls the linux mv command
//Input: source (const char []), dest (const char [])
//Output:None
//Return:None
void moveDir(const char source [], const char dest []) {

	char * sysCmd = new char [1000];
	sprintf(sysCmd, "mv %s %s", source, dest);
	system(sysCmd);
	delete [] sysCmd;
}

//void removeDir(const char [])
//Description: This function calls the linux rm -r command
//Input: dirName (const char [])
//Output:None
//Return:None
void removeDir(const char dirName []) {

	char * sysCmd = new char [1000];
	sprintf(sysCmd, "rm -r %s", dirName);
	system(sysCmd);
	delete [] sysCmd;
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


void createDirectories(string & workDir) {

	//Create the workDir
	char * arr = new char [1000];
	sprintf(arr, "mkdir %s", workDir.c_str());
	system(arr);

	//Make the index directory
	sprintf(arr, "mkdir %s/Fragments", workDir.c_str());
	system(arr);

	//Create the labels directory
	sprintf(arr, "mkdir %s/Labels", workDir.c_str());
	system(arr);

	//Make the overlaps directory
	sprintf(arr, "mkdir %s/Overlaps", workDir.c_str());
	system(arr);

	//Create the arrays directory
	sprintf(arr, "mkdir %s/Arrays", workDir.c_str());
	system(arr);

	//Create the clusters directory
	sprintf(arr, "mkdir %s/Clusters", workDir.c_str());
	system(arr);

	//Make the contigs directory
	sprintf(arr, "mkdir %s/Contigs", workDir.c_str());
	system(arr);

	//Make the contigs info directory
	sprintf(arr, "mkdir %s/ContigsInfo", workDir.c_str());
	system(arr);

	//Make the spectrum directory
	sprintf(arr, "mkdir %s/Spectrum", workDir.c_str());
	system(arr);

	//Make the hybrid directory
	sprintf(arr, "mkdir %s/HybSpectrum", workDir.c_str());
	system(arr);

	//Make the overlap graph directory
	sprintf(arr, "mkdir %s/OvlGraph", workDir.c_str());
	system(arr);

	//Make the overlap information directory
	sprintf(arr, "mkdir %s/OvlInfo", workDir.c_str());
	system(arr);
}
