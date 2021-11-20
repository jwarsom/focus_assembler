/*
 * 	ParallelAlign.cpp
 *
 *    Created on: Jan 12, 2013
 *         Author: Julia
 */

#include<vector>
#include <dirent.h>
#include <stdint.h>
#include <limits.h>
#include <cstdlib>
#include <math.h>
#include <fstream>
#include <map>
#include <algorithm>
using namespace std;
#include "mpi.h"
#include "MappingValues.h"
#include "Dictionary.h"
#include "Labels.h"
#include "Fragment_Index.h"
#include "Graph.h"


//The help message
void help_message(){

	cerr<<"\n\n";
	cerr<<"\t\t\t Parallel Cluster Analysis"<<endl;
	cerr<<"\t\t\tWritten by Julia Warnke\n"<<endl;
	cerr<<"Command line arguments\n"<<endl;
	cerr<<"--workDir :working directory :default dat"<<endl;
	cerr<<"--graphLevel : current graph level : default 0"<<endl;
	cerr<<"--GCcount :analyze for GC counts"<<endl;
	cerr<<"--tmerFreq :analyze for tetramer frequencies"<<endl;
	cerr<<"--minClustSizePercentile : return clusters with size percentile greater than a minimum (int) [0-100] "<<endl;
	cerr<<"--minClustDenPercentile :return clusters with densities greater than a minimum (int) [0-100]"<<endl;
	cerr<<"--maxClustSizePercentile : return clusters with size percentile less than a maximum (int) [0-100] "<<endl;
	cerr<<"--maxClustDenPercentile :return clusters with densities less than a maximum (int) [0-100]"<<endl;
	cerr<<"\n\n";

	cerr<<"\n\n"<<endl;

	cerr<<"Exiting program"<<endl;
	MPI_Abort(MPI_COMM_WORLD, help);
}

bool fileExists(const char []);
void recordTime(string &, const int, const bool);
int recoverClusters(string &, int * &, int);
void calculateBins(string &, multimap<int, int> &, int * &, int, int);
int countFiles(const char []);
void createFragmentIndex(string &, Fragment_Index &);
void createLabels(string &, Labels &); 
void gcCount(Fragment_Index &, int * &, int, int, ofstream &, int & numElements);
void tetraCount(Fragment_Index &, int * &, int, int, ofstream &);
void outputCluster(Fragment_Index &, Labels & label, int * &, int, int, ofstream &);

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
	int graphLevel = 0;
	bool GCcount = false;
	bool tmerFreq = false;
	float minClustSizePercentile = 0;
	float minClustDenPercentile = 0;
	float maxClustSizePercentile = 100;
	float maxClustDenPercentile = 100;
	bool verbose = false;
	bool runTime = false; 

	//user options
	string workDirS = "--workDir";
	string timeS = "--time";
	string graphLevelS = "--graphLevel";
	string GCcountS = "--GCcount";
	string tmerFreqS = "--tmerFreq";
	string minClustSizePercentileS = "--minClustSizePercentile";
	string minClustDenPercentileS = "--minClustDenPercentile";
	string maxClustSizePercentileS = "--maxClustSizePercentile";
	string maxClustDenPercentileS = "--maxClustDenPercentile";
	string verboseS = "--verbose";

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

		//The graph level
		if(argv[i] == graphLevelS)
		{
			graphLevel = atoi(argv[i+1]);
		}

		//Count GC content
		if(argv[i] == GCcountS){

			GCcount = true;
			i--;
		}

		//Count tetramers
		if(argv[i] == tmerFreqS)
		{
			tmerFreq = true;
			i--;
		}

		//The minimum identity
		if(argv[i] == graphLevelS)
		{
			graphLevel = atoi(argv[i+1]);
		}

		//The min cluster percentile size
		if(argv[i] == minClustSizePercentileS)
		{
			minClustSizePercentile = atof(argv[i+1]);
		}

		//The min cluster percentile density
		if(argv[i] == minClustDenPercentileS)
		{
			minClustDenPercentile = atof(argv[i+1]);
		}

		//The max cluster percentile size
		if(argv[i] == maxClustSizePercentileS)
		{
			maxClustSizePercentile = atof(argv[i+1]);
		}

		//The min cluster density
		if(argv[i] == maxClustDenPercentileS)
		{
			maxClustDenPercentile = atof(argv[i+1]);
		}

		//If verbose
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

	if(rank == 0)
	{
		multimap<int, int> binAssignments;

		int * binWeights;
		calculateBins(workDir, binAssignments, binWeights, graphLevel, nTasks);

		int * buff = new int [binAssignments.size() + nTasks];

		int count = 0;
		for(int i = 0; i < nTasks; i++)
		{
			buff[count++] =  binAssignments.count(i);

			pair<multimap<int,int>::iterator, multimap<int,int>::iterator> ret;
			ret = binAssignments.equal_range(i);

			for(multimap<int, int>::iterator it = ret.first; it != ret.second; it++)
				buff[count++] = (*it).second;
		}

		{
			char * arr = new char[200];
			sprintf(arr, "%s/Analysis/bins.dat", workDir.c_str());

			FILE * pFile = fopen(arr, "w");
			fwrite(buff, sizeof(int), count, pFile);
			fclose(pFile);
			delete [] arr;
		}

		delete [] buff;


		if(verbose)
		{
			ofstream myFile;

			char * arr = new char [200];
			sprintf(arr, "%s/Analysis/stats.txt", workDir.c_str());			
			myFile.open(arr);
			delete [] arr;

			float numElements = 0;
			for(int i = 0; i < nTasks; i++) numElements+=binWeights[i];

			float elementsPerBin = numElements/nTasks;

			float sumOfSquares = 0;

			myFile<<"Expected Elements per Bin"<<endl;

			myFile<<"Elements per Bin "<<endl;
			for(int i = 0; i < nTasks; i++)
			{
				myFile<<binWeights[i]<<endl;
				sumOfSquares+=pow(binWeights[i] - elementsPerBin, 2);
			}

			float variance = sumOfSquares/nTasks;

			myFile<<"Binning Variance : "<<variance<<endl;

			myFile.close();

		}

		delete [] binWeights;

		int work = 0;
		for(int i = 1; i < nTasks; i++)
			MPI_Send(&work, 1, MPI_INT, i, workTag, MPI_COMM_WORLD);

		char * arr = new char [200];
		sprintf(arr, "%s/Analysis/bins.dat", workDir.c_str());

		FILE * pFile = fopen(arr, "r");
		fseek(pFile, 0, SEEK_END);
		delete [] arr;

		count = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);

		buff = new int [count];
		fread(buff, count, sizeof(int), pFile);

		fclose(pFile);

		int startPos = 0;
		for(int i = 0; i < rank; i++)
		{
			int numIter = buff[startPos];
			for(int j = 0; j < numIter+1; j++)
				startPos++;
		}

		int * clusterSubset = new int [buff[startPos]];

		int numInSubset = buff[startPos]; startPos++;
		for(int i = 0; i < numInSubset; i++)
		{
			clusterSubset[i] = buff[startPos++];
		}

		delete [] buff;

		sort(clusterSubset, clusterSubset + numInSubset);

		arr = new char [200];
		sprintf(arr, "%s/Analysis/gcCount%d.txt", workDir.c_str(), rank);
		ofstream gcFile(arr);

		sprintf(arr, "%s/Analysis/tFreq%d.txt", workDir.c_str(), rank);
		ofstream tetraFile(arr);

		sprintf(arr, "%s/Analysis/clusterDen%d.txt", workDir.c_str(), rank);
		ofstream cDenFile(arr);

		sprintf(arr, "%s/Edges/Graph%d/nDen", workDir.c_str(), graphLevel);
		pFile = fopen(arr, "r");

		fseek(pFile, 0, SEEK_END);
		int numClusters = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);

		int * nodeDen = new int [numClusters];

		fread(nodeDen, sizeof(int), numClusters, pFile);
		fclose(pFile);

		sprintf(arr, "%s/Edges/Graph%d/eDen", workDir.c_str(), graphLevel);
		pFile = fopen(arr, "r");

		int * edgeDen = new int[numClusters];

		fread(edgeDen, sizeof(int), numClusters, pFile);
		fclose(pFile);
		delete [] arr;

		multimap<float, int> clusterSizePercentile;
		multimap<float, int> clusterDenPercentile;

		for(int i = 0; i < numClusters; i++)
		{
			if(nodeDen[i] > 1)
			{

				float size = nodeDen[i];
				float den = edgeDen[i]/(nodeDen[i] * (nodeDen[i] - 1));

				clusterSizePercentile.insert(pair<float, int>(size, i ));
				clusterDenPercentile.insert(pair<float, int>(den, i));

			}else{
				clusterSizePercentile.insert(pair<float, int>(1, i ));
				clusterDenPercentile.insert(pair<float, int>(0, i));
			}
		}

		{
			int i = 0;
			for(multimap<float, int>::iterator it = clusterSizePercentile.begin(); it != clusterSizePercentile.end(); it++)
			{
				nodeDen[(*it).second] = i++;
			}

			i = 0;
			for(multimap<float, int>::iterator it = clusterDenPercentile.begin(); it != clusterDenPercentile.end(); it++)
			{
				edgeDen[(*it).second] = i++;
			}
		}

		int minSizeIndex = clusterSizePercentile.size() * (minClustSizePercentile/100);
		int maxSizeIndex = clusterSizePercentile.size() * (maxClustSizePercentile/100);
		int minDenIndex = clusterDenPercentile.size() * (minClustDenPercentile/100);
		int maxDenIndex = clusterDenPercentile.size() * (maxClustDenPercentile/100);

		clusterSizePercentile.clear();
		clusterDenPercentile.clear();

		Fragment_Index index;
		createFragmentIndex(workDir, index);

		Labels label;
		createLabels(workDir, label); 

		int * clusters;
		count = recoverClusters(workDir, clusters, graphLevel);

		int currentCluster = 0; startPos = 0;
		int numElements = 0;
		for(int i = 0; i < count; i++)
		{
			if(startPos < numInSubset && currentCluster == clusterSubset[startPos])
			{
				//Indexes of clusters
				gcFile<<currentCluster<<" ";
				tetraFile<<currentCluster<<" ";
				cDenFile<<currentCluster<<" ";

				gcCount(index, clusters, i, count, gcFile, numElements);
				tetraCount(index, clusters, i, count, tetraFile);

				if(nodeDen[currentCluster] >= minSizeIndex && nodeDen[currentCluster] <= maxSizeIndex)
					if((edgeDen[currentCluster] >= minDenIndex && edgeDen[currentCluster] <= maxDenIndex))
						outputCluster(index, label, clusters, i, count, cDenFile);

				startPos++;
			}

			if(clusters[i] == -1)
			{
				currentCluster++;
			}
		}


		gcFile.close();
		tetraFile.close();
		cDenFile.close();

		delete [] edgeDen;
		delete [] nodeDen;
		delete [] clusters;
		delete [] clusterSubset;

		for(int i = 1; i < nTasks; i++)
			MPI_Recv(&work, 1, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

	}else{

		int work = 0;
		MPI_Recv(&work, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		char * arr = new char [200];
		sprintf(arr, "%s/Analysis/bins.dat", workDir.c_str());

		FILE * pFile = fopen(arr, "r");
		fseek(pFile, 0, SEEK_END);
		delete [] arr;

		int count = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);

		int * buff = new int [count];
		fread(buff, count, sizeof(int), pFile);

		fclose(pFile);

		int startPos = 0;
		for(int i = 0; i < rank; i++)
		{
			int numIter = buff[startPos];
			for(int j = 0; j < numIter+1; j++)
				startPos++;
		}

		int * clusterSubset = new int [buff[startPos]];

		int numInSubset = buff[startPos]; startPos++;
		for(int i = 0; i < numInSubset; i++)
		{
			clusterSubset[i] = buff[startPos++];
		}

		sort(clusterSubset, clusterSubset + numInSubset);

		arr = new char [200];
		sprintf(arr, "%s/Analysis/gcCount%d.txt", workDir.c_str(), rank);
		ofstream gcFile(arr);

		sprintf(arr, "%s/Analysis/tFreq%d.txt", workDir.c_str(), rank);
		ofstream tetraFile(arr);

		sprintf(arr, "%s/Analysis/clusterDen%d.txt", workDir.c_str(), rank);
		ofstream cDenFile(arr);

		sprintf(arr, "%s/Edges/Graph%d/nDen", workDir.c_str(), graphLevel);
		pFile = fopen(arr, "r");

		fseek(pFile, 0, SEEK_END);
		int numClusters = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);

		int * nodeDen = new int [numClusters];

		fread(nodeDen, sizeof(int), numClusters, pFile);
		fclose(pFile);

		sprintf(arr, "%s/Edges/Graph%d/eDen", workDir.c_str(), graphLevel);
		pFile = fopen(arr, "r");

		int * edgeDen = new int[numClusters];

		fread(edgeDen, sizeof(int), numClusters, pFile);
		fclose(pFile);

		multimap<float, int> clusterSizePercentile;
		multimap<float, int> clusterDenPercentile;

		for(int i = 0; i < numClusters; i++)
		{
			if(nodeDen[i] > 1)
			{
				float size = nodeDen[i];
				float den = edgeDen[i]/(nodeDen[i] * (nodeDen[i] - 1));

				clusterSizePercentile.insert(pair<float, int>(size, i ));
				clusterDenPercentile.insert(pair<float, int>(den, i));

			}else{

				clusterSizePercentile.insert(pair<float, int>(1, i));
				clusterDenPercentile.insert(pair<float, int>(0, i));
			}
		}

		{
			int i = 0;
			for(multimap<float, int>::iterator it = clusterSizePercentile.begin(); it != clusterSizePercentile.end(); it++)
			{
				nodeDen[(*it).second] = i++;
			}

			i = 0;
			for(multimap<float, int>::iterator it = clusterDenPercentile.begin(); it != clusterDenPercentile.end(); it++)
			{
				edgeDen[(*it).second] = i++;
			}
		}

		int minSizeIndex = clusterSizePercentile.size() * (minClustSizePercentile/100);
		int maxSizeIndex = clusterSizePercentile.size() * (maxClustSizePercentile/100);
		int minDenIndex = clusterDenPercentile.size() * (minClustDenPercentile/100);
		int maxDenIndex = clusterDenPercentile.size() * (maxClustDenPercentile/100);

		clusterSizePercentile.clear();
		clusterDenPercentile.clear();

		Fragment_Index index;
		createFragmentIndex(workDir, index);

		Labels label;
		createLabels(workDir, label); 

		int * clusters;
		count = recoverClusters(workDir, clusters, graphLevel);

		int currentCluster = 0; startPos = 0; int numElements = 0;
		for(int i = 0; i < count; i++)
		{
			if(startPos < numInSubset && currentCluster == clusterSubset[startPos])
			{
				//Indexes of clusters
				gcFile<<currentCluster<<" ";
				tetraFile<<currentCluster<<" ";
				cDenFile<<currentCluster<<" ";

				gcCount(index, clusters, i, count, gcFile, numElements);
				tetraCount(index, clusters, i, count, tetraFile);

				if(nodeDen[currentCluster] >= minSizeIndex && nodeDen[currentCluster] <= maxSizeIndex)
					if(edgeDen[currentCluster] >= minDenIndex && edgeDen[currentCluster] <= maxDenIndex)
						outputCluster(index, label, clusters, i, count, cDenFile);

				startPos++;
			}

			if(clusters[i] == -1)
			{
				currentCluster++;
			}
		}

		gcFile.close();
		tetraFile.close();
		cDenFile.close();

		delete [] clusters;
		delete [] clusterSubset;
		delete [] edgeDen;
		delete [] nodeDen;

		MPI_Send(&work, 1, MPI_INT, 0, workTag, MPI_COMM_WORLD);
	}

	if(runTime && rank == 0)
		recordTime(workDir, rank, false);

	MPI_Finalize();
	return 0;
}

//void gcCount(int *, const int)
//Description: This function counts GC content.
//Input: workDir (string &): The working directory, rank (const int):
//The nodes rank, begin (bool): execution start or finish
//Output:None
//Return:None
void gcCount(Fragment_Index & index, int * & clusters, int startPos, int count, ofstream & gcFile, int & numElements)
{
	float A = 0; float T = 0; float C = 0; float G = 0;

	while(clusters[startPos]  != -1 && startPos < count)
	{

		int fragment = clusters[startPos];

		numElements++;
		pair<long long int, long long int> myRange = index.indexBounds(fragment);


		for(int i = myRange.first; i != myRange.second; i++)
		{
			if(index.at(i) == 'A') A++;
			if(index.at(i) == 'T') T++;
			if(index.at(i) == 'C') C++;
			if(index.at(i) == 'G') G++;
		}


		startPos++;
	}

	gcFile<<((G+C)/(G+C+A+T))<<endl;
}

//void tetraCount(int *, const int)
//Description: This function counts GC content.
//Input: workDir (string &): The working directory, rank (const int):
//The nodes rank, begin (bool): execution start or finish
//Output:None
//Return:None
void tetraCount(Fragment_Index & index, int * & clusters, int startPos, int count, ofstream & tetraFile)
{
	float totalTetramers = 0; float tCounts[256];
	for(int i = 0; i < 256; i++) tCounts[i] = 0;

	while(clusters[startPos]  != -1 && startPos < count)
	{
		int fragment = clusters[startPos];
		pair<long long int, long long int> myRange = index.indexBounds(fragment);

		for(int i = myRange.first; i < myRange.second-4; i++)
		{
			int tHash = 0;
			for(int j = 0; j < 4; j++)
			{
				if(index.at(i+j) == 'C')
				{
					tHash+=(pow(2, 6-j*2));
				}

				if(index.at(i+j) == 'G')
				{
					tHash+=(2*pow(2, 6-j*2));
				}

				if(index.at(i+j) == 'T')
				{
					tHash+=(3*pow(2, 6-j*2));
				}
			}
			tCounts[tHash]++;
			totalTetramers++;
		}
		startPos++;
	}

	for(int i = 0; i < 256; i++)
		tetraFile<<tCounts[i]/totalTetramers<<" ";

	tetraFile<<endl;
}

//void gcCount(int *, const int)
//Description: This function counts GC content.
//Input: workDir (string &): The working directory, rank (const int):
//The nodes rank, begin (bool): execution start or finish
//Output:None
//Return:None
void outputCluster(Fragment_Index & index, Labels & label, int * & clusters, int startPos, int count, ofstream & myFile)
{
	myFile<<"*"<<endl;
	while(clusters[startPos] != -1 && startPos < count)
	{
		int fragment = clusters[startPos];
		pair<long long int, long long int> myRange = index.indexBounds(fragment);

		myFile<<label.labelAt(fragment)<<endl;

		for(int i = myRange.first; i != myRange.second; i++)
		{
			myFile<<index.at(i);
		}
		myFile<<endl;
		startPos++;
	}
}

//void createFragmentIndex()
//Description: This fills the fragment index
//Input: fIndex, (Fragment_Index)
//Output: None
//Return: None
void createFragmentIndex(string & workDir, Fragment_Index & index)
{
	char * fragmentDir = new char [1000];
	sprintf(fragmentDir, "%s/Fragments", workDir.c_str());

	int numFiles = countFiles(fragmentDir);

	int numSets = 0;
	char * fileName = new char [1000];
	sprintf(fileName, "%s/sequences.job0.set%d", fragmentDir, numSets++);

	while(fileExists(fileName))
	{
		sprintf(fileName, "%s/sequences.job0.set%d", fragmentDir, numSets++);
	}
	numSets-=1;

	int numJobs = numFiles/numSets;

	sprintf(fileName, "%s/sequences.job%d.set%d", fragmentDir, numJobs-1, numSets-1);

	FILE * pFile  = fopen(fileName, "r");
	unsigned long long int fragmentOffset;
	fread(&fragmentOffset, sizeof(unsigned long long int), 1, pFile);

	unsigned long long int indexOffset;
	fread(&indexOffset, sizeof(unsigned long long int), 1, pFile);

	unsigned long long int intervalLength;
	fread(&intervalLength, sizeof(unsigned long long int), 1, pFile);

	unsigned long long int numIntervalFragments;
	fread(&numIntervalFragments, sizeof(unsigned long long int), 1, pFile);
	fclose(pFile);

	int numFragments = numIntervalFragments + indexOffset;
	int length = intervalLength + fragmentOffset;

	index.reserve(numFragments, length, 0, 0);

	for(int i = 0; i < numJobs; i++)
	{
		for(int j = 0; j < numSets; j++)
		{
			sprintf(fileName, "%s/sequences.job%d.set%d", fragmentDir, i, j);
			pFile = fopen(fileName, "r");
			index.read(pFile);
			fclose(pFile);
		}
	}

	index.set();

	delete [] fragmentDir;
	delete [] fileName;
}

//void createLabels()
//Description: This fills the labels
//Input: label, (Labels)
//Output: None
//Return: None
void createLabels(string & workDir, Labels & label)
{
	char * fragmentDir = new char [1000];
	sprintf(fragmentDir, "%s/Fragments", workDir.c_str());

	int numFiles = countFiles(fragmentDir);

	int numSets = 0;
	char * fileName = new char [1000];
	sprintf(fileName, "%s/sequences.job0.set%d", fragmentDir, numSets++);

	while(fileExists(fileName))
	{
		sprintf(fileName, "%s/sequences.job0.set%d", fragmentDir, numSets++);
	}
	numSets-=1;

	int numJobs = numFiles/numSets;

	sprintf(fileName, "%s/sequences.job%d.set%d", fragmentDir, numJobs-1, numSets-1);

	FILE * pFile  = fopen(fileName, "r");
	unsigned long long int fragmentOffset;
	fread(&fragmentOffset, sizeof(unsigned long long int), 1, pFile);

	unsigned long long int indexOffset;
	fread(&indexOffset, sizeof(unsigned long long int), 1, pFile);

	unsigned long long int intervalLength;
	fread(&intervalLength, sizeof(unsigned long long int), 1, pFile);

	unsigned long long int numIntervalFragments;
	fread(&numIntervalFragments, sizeof(unsigned long long int), 1, pFile);
	fclose(pFile);

	int numFragments = numIntervalFragments + indexOffset;

	int tagLen;
	sprintf(fileName, "%s/Labels/labels.job%d.set%d", workDir.c_str(), numJobs-1, numSets-1);
	pFile = fopen(fileName, "r");
	fread(&tagLen, sizeof(int), 1, pFile);
	fclose(pFile);

	label.reserve(tagLen, numFragments, 0);

	for(int i = 0; i < numJobs; i++)
	{
		for(int j = 0; j < numSets; j++)
		{
			sprintf(fileName, "%s/Labels/labels.job%d.set%d", workDir.c_str(), i, j);
			pFile = fopen(fileName, "r");
			label.read(pFile);
			fclose(pFile);
		}
	}

	delete [] fragmentDir;
	delete [] fileName;
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


//void calculateBins()
//Description: This calculates bins for the parallel algorithm
//Input: workDir (string &): the working directory, minMerge (int):
//the minimum overlap length for merging, minDensity (double): The minimum inter-node
//density
//Output:None
//Return:None
void calculateBins(string & workDir, multimap<int, int> & binAssignments, int * & binWeights, int graphLevel, int nTasks)
{
	int * clusters;
	int count = recoverClusters(workDir, clusters, graphLevel);

	binWeights = new int [nTasks];
	for(int i = 0; i < nTasks; i++) binWeights[i] = 0;
	multimap<int, int> clusterWeights;

	int tally = 0; int numClusters = 0; float numElements = 0;
	for(int i = 0; i < count; i++)
	{
		if(clusters[i] != -1)
		{
			tally++; numElements++;
		}
		else
		{
			clusterWeights.insert(pair<int, int>(-(tally), numClusters));
			numClusters++; tally = 0;
		}
	}

	if(tally != 0)
	{
		clusterWeights.insert(pair<int, int>(-(tally), numClusters));
		numClusters++; tally = 0;
	}

	delete [] clusters;

	float elementsPerBin = numElements/nTasks;

	int * randomArr = new int [nTasks];
	for(int i = 0; i < nTasks; i++)
		randomArr[i] = i;


	for(multimap<int, int>::iterator it = clusterWeights.begin(); it != clusterWeights.end(); it++)
	{
		bool placed = false;
		random_shuffle(randomArr, randomArr+nTasks);

		for(int i = 0; i < nTasks; i++)
		{
			int clusterSize = -(*it).first;
			if(binWeights[randomArr[i]]+clusterSize < elementsPerBin)
			{
				binAssignments.insert(pair<int, int>(randomArr[i], (*it).second));
				binWeights[randomArr[i]]+=clusterSize;
				placed = true; break;
			}
		}

		if(!placed)
		{
			int min = INT_MAX; int hold = 0;
			for(int i = 0; i < nTasks; i++)
			{
				int clusterSize = -(*it).first;
				if(binWeights[randomArr[i]]+clusterSize < min)
				{
					min = binWeights[randomArr[i]]+clusterSize; hold = randomArr[i];
				}
			}
			binAssignments.insert(pair<int, int>(hold, (*it).second));
			binWeights[hold] = min;
		}
	}
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


//void recoverClusters()
//Description: This runs the the parent thread in a parallel program
//Input: workDir (string &): the working directory, minMerge (int):
//the minimum overlap length for merging, minDensity (double): The minimum inter-node
//density
//Output:None
//Return:None
int recoverClusters(string & workDir, int * & clusters, int graphLevel) {

	char * graphDir = new char [1000];
	sprintf(graphDir, "%s/Edges/Graph0", workDir.c_str());

	int numFiles = countFiles(graphDir) - 2;
	long int numNodes = 0; long int numEdges = 0;

	cout<<numFiles-1<<endl;
	for(int i = 0; i < 1; i++)
	{
		char * currentGraphDir = new char [1000];
		sprintf(currentGraphDir, "%s/Edges/Graph%d", workDir.c_str(), graphLevel);

		Graph graph;
		char * fileName = new char [1000];
		sprintf(fileName, "%s/graph%d", currentGraphDir, i);
		cout<<fileName<<endl;
		FILE * pFile = fopen(fileName, "r");
		delete [] fileName;

		graph.read(pFile);
		fclose(pFile);
		numEdges+= graph.getNumEdges();
	
		delete [] currentGraphDir;

	}

	{
		Graph graph;
		char * fileName = new char [1000];
		sprintf(fileName, "%s/graph%d", graphDir, numFiles-1);
		FILE * pFile = fopen(fileName, "r");
		delete [] fileName;

		graph.read(pFile);
		fclose(pFile);
		numNodes = graph.getOffset() + graph.getNumNodes();
	}

	sprintf(graphDir, "%s/Edges", workDir.c_str());
	cout<<"NumNodes :"<<numNodes<<endl;
	cout<<"NumEdges :"<<numEdges<<endl;

	clusters = new int [2 * numNodes];

	int count = 0;
	for(int i = graphLevel; i > 0; i--)
	{
		char * fileName = new char [1000];
		sprintf(fileName, "%s/Graph%d/map", graphDir, i);

		FILE * pFile = fopen(fileName, "r");
		fseek(pFile, 0, SEEK_END);

		int buffSize = ftell(pFile)/sizeof(int);
		fseek(pFile, 0, SEEK_SET);

		int * buff = new int [buffSize];
		fread(buff, sizeof(int), buffSize, pFile);
		fclose(pFile);

		if(i == graphLevel)
		{
			for(int j = 0; j < buffSize; j+=2)
			{
				clusters[j] = j/2;
				clusters[j+1] = -1;
			}
			count = buffSize;
		}

		int * tmp = new int [2 * numNodes];

		int cHold = count; count = 0;
		for(int j = 0; j < cHold; j++)
		{
			if(clusters[j]  != -1)
			{
				tmp[count++] = buff[clusters[j]*2];
				if(buff[clusters[j]*2+1] != -1)
					tmp[count++] = buff[clusters[j]*2+1];
			}else{

				tmp[count++] = -1;
			}
		}


		delete [] clusters;

		clusters = tmp;

		delete [] buff;
		delete [] fileName;
	}

	delete [] graphDir;

	return count;
}
