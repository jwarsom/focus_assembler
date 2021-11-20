/*
 * 	ParallelAlign.cpp
 *
 *         Author: Julia
 */

#include<iostream>
#include<vector>
#include <dirent.h>
#include <stdint.h>
#include <limits.h>
using namespace std;
#include "mpi.h"
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
	cerr<<"--highCoverage : High coverage dataset > 50x"<<endl;
	cerr<<"--alignRepsOnly: Align only the representative sequences obtained by clustering when high coverage is set to true"<<endl;
	cerr<<"--time : records run times for the serial and total excution times"<<endl;
	cerr<<"\n\n";

	cerr<<"\n\n"<<endl;

	cerr<<"Exiting program"<<endl;
	MPI_Abort(MPI_COMM_WORLD, help);
}

void recordTime(string &, const int, const bool);
int countFiles(const char []);
unsigned long long int countFragments(string & );
bool fileExists(const char []);


int main ( int argc, char *argv[] )
{
	int rank;
	int nTasks;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nTasks);
	MPI_Status status;

	//Param options and defaults
	string workDir = "NO_WORK_DIR";
	int minIden = 90;
	int minOverlap = 50;
	int kmerSize = 22;
	int stepSize = 11;
	float maxMem = 1024;
	bool verbose = false;
	bool runTime = false;
	bool highCoverage = false;
	bool alignRepsOnly = false;

	//user options
	string workDirS = "--workDir";
	string minIdenS = "--minIden";
	string minOverlapS = "--minOverlap";
	string kmerSizeS =  "--kmerSize";
	string stepSizeS = "--stepSize";
	string maxMemS = "--maxMem";
	string verboseS = "--verbose";
	string timeS = "--time";
	string highCoverageS = "--highCoverage";
	string alignRepsOnlyS = "--alignRepsOnly";

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

		if(argv[i] == verboseS)
		{
			verbose = true;
			i--;
		}

		if(argv[i] == timeS)
		{
			runTime = true;
			i--;
		}

		if(argv[i] == highCoverageS)
		{
			highCoverage = true;
			i--;
		}

		if(argv[i] == alignRepsOnlyS)
		{
			alignRepsOnly = true;
			i--;
		}

	}

	if(nTasks < 2)
	{
		cerr<<"It looks like you have only requested one processor"<<endl;
		cerr<<"Please use the serial version of the iterative program for a single processor task"<<endl;
		cerr<<"The serial program is located at '/home/usr_name/Merge_and_Traverse/bin/SerialAlign"<<endl;
		MPI_Abort(MPI_COMM_WORLD, runSerial);
	}

	if(runTime && rank == 0)
		recordTime(workDir, rank, true);

	//count the subsets
	struct dirent *pDirent;
	DIR *pDir;
	char * cmdStr = new char[1000];
	sprintf(cmdStr, "%s/Fragments", workDir.c_str());

	pDir = opendir(cmdStr);
	int numJobs = 0; int totalSets = 0;
	if(pDir != NULL)
	{
		while ((pDirent = readdir(pDir)) != NULL)
		{
			if(strcmp(pDirent->d_name, ".") != 0 && strcmp(pDirent->d_name, "..") != 0)
			{
				char * pch1 = strstr(pDirent->d_name, "job");
				char * pch2 = strstr(pDirent->d_name, ".set");
				char job [10];
				strncpy(job, pch1+3, (pch2-(pch1+3)));
				job[(pch2-(pch1+3))] = '\0';
				numJobs = max(atoi(job), numJobs);
				totalSets++;
			}
		}

		closedir (pDir);
	}

	numJobs++;
	int numSets = totalSets/numJobs;


	string appLoc = argv[0];
	int endLoc = appLoc.find_last_of('/');
	appLoc.erase(endLoc);

	if(rank == 0)
	{
		for(int i = 0; i < numJobs; i++)
		{ 
			char * ovlInfoDir = new char [1000];
			sprintf(ovlInfoDir, "%s/OvlInfo/", workDir.c_str());

			char * fragmentDir = new char [1000];
			sprintf(fragmentDir, "%s/Fragments/", workDir.c_str());

			char * fileName = new char [1000];

			sprintf(fileName, "%s/sequences.job%d.set0", fragmentDir, i);

			unsigned long long int indexOffset = 0; 

			FILE * pFile = fopen(fileName, "r");
			fread(&indexOffset, sizeof(unsigned long long int), 1, pFile);
			fread(&indexOffset, sizeof(unsigned long long int), 1, pFile);
			fclose(pFile); 

			sprintf(fileName, "%s/sequences.job%d.set%d", fragmentDir, i, numSets-1);

			unsigned long long int indexOffset2 = 0; unsigned long long int totalFragments;

			pFile = fopen(fileName, "r");
			fread(&indexOffset2, sizeof(unsigned long long int), 1, pFile);
			fread(&indexOffset2, sizeof(unsigned long long int), 1, pFile);

			fread(&totalFragments, sizeof(unsigned long long int), 1, pFile);
			fread(&totalFragments, sizeof(unsigned long long int), 1, pFile);
			fclose(pFile);

			int numberFragmentsInSubset = ((indexOffset2+totalFragments) - indexOffset);

			int * cOvl = new int [3*numberFragmentsInSubset+1];

			cOvl[0] = indexOffset;	

			for(int j = 1; j < 3*numberFragmentsInSubset+1; j+=3)
			{
				cOvl[j] = (j-1)/3 + indexOffset; 
				cOvl[j+1] = cOvl[j+2]=0; 
			}

			sprintf(fileName, "%s/sequences.job%d.info", ovlInfoDir, i);

			pFile = fopen(fileName, "w");
			fwrite(cOvl, sizeof(int), 3*numberFragmentsInSubset+1, pFile);
			fclose(pFile);

			delete [] fileName;
			delete [] fragmentDir;
			delete [] ovlInfoDir;		

			delete [] cOvl;
		}

		//Containment Overlapping 
		if(highCoverage == true)
		{
			priority_queue<pair<int, int>, vector<pair<int, int> >, greater<pair<int, int> > > availJobs;
			priority_queue<pair<int, int>, vector<pair<int, int> >, greater<pair<int, int> > > priorityJobs;
			priority_queue<pair<int, bool> > availProcs;

			for(int i = 1; i < nTasks; i++)
				availProcs.push(pair<int, bool>(i, false));


			int queryJob = 0; int refJob = 0; int numSeeded = 0; int ranJobs = 0;
			int processor = availProcs.top().first;

			cout<<"SENT "<<queryJob<<" "<<refJob<<endl;
			MPI_Send(&queryJob, 1, MPI_INT, processor, workTag, MPI_COMM_WORLD);
			MPI_Send(&refJob, 1, MPI_INT, processor, workTag, MPI_COMM_WORLD);

			refJob++; numSeeded++; ranJobs++; availProcs.pop();

			set< pair<int, int> > locks;
			while(ranJobs < ((numJobs * (numJobs+1))/2))
			{
				int fQueryJob = 0; int fRefJob = 0;

				MPI_Recv(&fQueryJob, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&fRefJob, 1, MPI_INT, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				cout<<"RECEIVED "<<fQueryJob<<" "<<fRefJob<<endl;

				int processor = status.MPI_SOURCE;
				availProcs.push(pair<int, bool>(processor, true));

				if(fQueryJob == 0 && fRefJob == 0)
					for(int i = 1; i < numJobs; i++)
						availJobs.push(pair<int, int>(0, i));

				if(fQueryJob + 1 < fRefJob)
				{
					pair<set<pair<int, int> >::iterator, bool > wasInserted;
					wasInserted = locks.insert(pair<int,int>(fQueryJob + 1, fRefJob));

					if(!wasInserted.second)
					{
						availJobs.push(pair<int, int>(fQueryJob + 1, fRefJob));
						locks.erase(wasInserted.first);
					}
				}

				if(fQueryJob + 1 == fRefJob)
				{
					priorityJobs.push(pair<int, int>(fQueryJob + 1, fRefJob));
				}

				if(fQueryJob == fRefJob)
				{
					for(int i = fRefJob+1; i < numJobs; i++)
					{
						pair< set<pair<int, int> >::iterator, bool > wasInserted;
						wasInserted = locks.insert(pair<int,int>(fQueryJob, i));

						if(!wasInserted.second)
						{
							availJobs.push(pair<int, int>(fQueryJob, i));
							locks.erase(wasInserted.first);
						}
					}
				}

				while(!availProcs.empty() && !priorityJobs.empty())
				{
					int processor =  availProcs.top().first;
					queryJob = priorityJobs.top().first;
					refJob = priorityJobs.top().second;

					MPI_Send(&queryJob, 1, MPI_INT, processor, workTag, MPI_COMM_WORLD);
					MPI_Send(&refJob, 1, MPI_INT, processor, workTag, MPI_COMM_WORLD);
					cout<<"SEND "<<queryJob<<" "<<refJob<<endl;

					if(availProcs.top().second == false)
						numSeeded++;

					ranJobs++; availProcs.pop(); priorityJobs.pop();
				}

				while(!availProcs.empty() && !availJobs.empty())
				{
					int processor =  availProcs.top().first;
					queryJob = availJobs.top().first;
					refJob = availJobs.top().second;

					MPI_Send(&queryJob, 1, MPI_INT, processor, workTag, MPI_COMM_WORLD);
					MPI_Send(&refJob, 1, MPI_INT, processor, workTag, MPI_COMM_WORLD);
					cout<<"SEND "<<queryJob<<" "<<refJob<<endl;

					if(availProcs.top().second == false)
						numSeeded++;

					ranJobs++; availProcs.pop(); availJobs.pop();
				}

			}

			int fQueryJob = 0; int fRefJob = 0;
			MPI_Recv(&fQueryJob, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			MPI_Recv(&fRefJob, 1, MPI_INT, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			cout<<"RECIEVED "<<fRefJob<<" "<<fQueryJob<<endl;

			for(int i = 1; i < nTasks; i++)
			{
				MPI_Send(0, 0, MPI_INT, i, dieTag, MPI_COMM_WORLD);

				int completedTask = 0;
				MPI_Recv(&completedTask, 1, MPI_INT, i, dieTag, MPI_COMM_WORLD, &status);
			}

		}else{
			for(int i = 1; i < nTasks; i++)
			{
				MPI_Send(0, 0, MPI_INT, i, dieTag, MPI_COMM_WORLD);

				int completedTask = 0;
				MPI_Recv(&completedTask, 1, MPI_INT, i, dieTag, MPI_COMM_WORLD, &status);
			}
		}

		//Dovetail Overlapping
		{
			priority_queue<pair<int, int>, vector<pair<int, int> >, greater<pair<int, int> > > availJobs;
			priority_queue<pair<int, bool> > availProcs;

			for(int i = 1; i < nTasks; i++)
				availProcs.push(pair<int, bool>(i, false));

			for(int i = 0; i < numJobs; i++)
				availJobs.push(pair<int, int>(0, i)); //refJob + queryJob

			int queryJob = 0; int refJob = 0; int numSeeded = 0; int ranJobs = 0;
			while(!availProcs.empty() && !availJobs.empty())
			{
				int processor = availProcs.top().first;
				queryJob = availJobs.top().first;
				refJob = availJobs.top().second;

				MPI_Send(&queryJob, 1, MPI_INT, processor, workTag, MPI_COMM_WORLD);
				MPI_Send(&refJob, 1, MPI_INT, processor, workTag, MPI_COMM_WORLD);

				refJob++; numSeeded++; ranJobs++; availProcs.pop(); availJobs.pop();
			}

			while(ranJobs < numJobs * numJobs)
			{
				int fQueryJob = 0; int fRefJob = 0;

				MPI_Recv(&fQueryJob, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&fRefJob, 1, MPI_INT, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

				int processor = status.MPI_SOURCE;
				availProcs.push(pair<int, bool>(processor, true));


				if(fQueryJob == 0)
					for(int i = 1; i < numJobs; i++)
						availJobs.push(pair<int, int>(i, fRefJob));

				while(!availProcs.empty() && !availJobs.empty())
				{
					int processor =  availProcs.top().first;
					queryJob = availJobs.top().first;
					refJob = availJobs.top().second;

					MPI_Send(&queryJob, 1, MPI_INT, processor, workTag, MPI_COMM_WORLD);
					MPI_Send(&refJob, 1, MPI_INT, processor, workTag, MPI_COMM_WORLD);

					if(availProcs.top().second == false)
						numSeeded++;

					ranJobs++; availProcs.pop(); availJobs.pop();
				}
			}

			for(int i = 0; i < numSeeded; i++)
			{
				int fQueryJob = 0; int fRefJob = 0;
				MPI_Recv(&fQueryJob, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&fRefJob, 1, MPI_INT, status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				MPI_Send(0, 0, MPI_INT, status.MPI_SOURCE, dieTag, MPI_COMM_WORLD);
			}

			delete [] cmdStr;
		}

	}else{

		while(true)
		{
			int queryJob = 0; int refJob = 0;
			MPI_Recv(&queryJob, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			if(status.MPI_TAG == dieTag)
			{
				MPI_Send(0, 0, MPI_INT, 0, dieTag, MPI_COMM_WORLD);
				break;
			}

			MPI_Recv(&refJob, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			char * cmdStr = new char[1000];
			sprintf(cmdStr, "%s/align --workDir %s --queryJob %d --numSets %d --refJob %d --containment --minIden %d --minOverlap %d --kmerSize %d --stepSize %d --maxMem %f", appLoc.c_str(), workDir.c_str(), queryJob, numSets, refJob, minIden, minOverlap, kmerSize, stepSize, maxMem);

			if(runTime)
				recordTime(workDir, rank, true);

			system(cmdStr);  delete [] cmdStr;

			if(runTime)
				recordTime(workDir, rank, false);

			MPI_Send(&queryJob, 1, MPI_INT, 0, workTag, MPI_COMM_WORLD);
			MPI_Send(&refJob, 1, MPI_INT, 0, workTag, MPI_COMM_WORLD);
		}

		while(true)
		{
			int queryJob = 0; int refJob = 0;
			MPI_Recv(&queryJob, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			if(status.MPI_TAG == dieTag)
				break;

			MPI_Recv(&refJob, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			char * cmdStr = new char[1000];

			if(alignRepsOnly)
				sprintf(cmdStr, "%s/align --workDir %s --queryJob %d --numSets %d --refJob %d --minIden %d --minOverlap %d --kmerSize %d --stepSize %d --maxMem %f --alignRepsOnly", appLoc.c_str(), workDir.c_str(), queryJob, numSets, refJob, minIden, minOverlap, kmerSize, stepSize, maxMem);
			else
				sprintf(cmdStr, "%s/align --workDir %s --queryJob %d --numSets %d --refJob %d --minIden %d --minOverlap %d --kmerSize %d --stepSize %d --maxMem %f", appLoc.c_str(), workDir.c_str(), queryJob, numSets, refJob, minIden, minOverlap, kmerSize, stepSize, maxMem);

			if(runTime)
				recordTime(workDir, rank, true);

			system(cmdStr);  delete [] cmdStr;

			if(runTime)
				recordTime(workDir, rank, false);

			MPI_Send(&queryJob, 1, MPI_INT, 0, workTag, MPI_COMM_WORLD);
			MPI_Send(&refJob, 1, MPI_INT, 0, workTag, MPI_COMM_WORLD);
		}
	}

	if(runTime && rank == 0)
		recordTime(workDir, rank, false);


	MPI_Finalize();
	return 0;
}

//void countFragments(const char [])
//Description: This function counts the number of fragments in the fragment set
//Input: workDir (string &) : The working directory
//Output:None
//Return:None
unsigned long long int countFragments(string & workDir)
{
	char * fileName = new char [1000];
	sprintf(fileName, "%s/Fragments", workDir.c_str());

	int numFiles = countFiles(fileName);
	int numJobs = 0; int numSets = 0;

	for(int i = 0; i < numFiles; i++)
	{
		sprintf(fileName, "%s/Fragments/sequences.job0.set%d", workDir.c_str(), i);
		if(!fileExists(fileName))
			break;

		numSets++;
	}

	numJobs = numFiles/numSets;

	sprintf(fileName, "%s/Fragments/sequences.job%d.set%d", workDir.c_str(), numJobs-1, numSets-1);

	unsigned long long int indexOffset = 0; unsigned long long int totalFragments;

	FILE * pFile = fopen(fileName, "r");
	fread(&indexOffset, sizeof(unsigned long long int), 1, pFile);
	fread(&indexOffset, sizeof(unsigned long long int), 1, pFile);

	fread(&totalFragments, sizeof(unsigned long long int), 1, pFile);
	fread(&totalFragments, sizeof(unsigned long long int), 1, pFile);

	delete [] fileName;

	return (indexOffset + totalFragments);
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


//void recordTime(string &, const int, const bool)
//Description: This function records execution times
//Input: workDir (string &): The working directory, rank (const int):
//The nodes rank, begin (bool): execution start or finish
//Output:None
//Return:None
void recordTime(string & workDir, const int rank, const bool begin) {

	time_t rawtime;
	struct tm * timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);

	char * fileName = new char [1000];
	if(rank == 0)
		sprintf(fileName, "%s/timeAll", workDir.c_str());
	else
		sprintf(fileName, "%s/timeNode%d", workDir.c_str(), rank);

	ofstream myfile;
	myfile.open (fileName, ios::app);

	if(begin)
		myfile<<"start time: "<<asctime(timeinfo)<<endl;
	else
		myfile<<"end time: "<<asctime(timeinfo)<<endl;

	myfile.close();

	delete [] fileName;
}
