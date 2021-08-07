#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <csignal>
#include <unistd.h>
#include "TSP.h"

using namespace std;

bool die = false;
const int TIME_LIMIT = 180;

void sigHandler(int p){
	die = true;
}	

int main(int argc, char **argv){
	time_t tempTime;
	char *filename = NULL, *algorithm = NULL/*, *localOpt = NULL*/;
	unsigned int seed = 0;
	int cost, numIter = 100, k = 15;
	double T = 100.0, t = 1.0, decreaseFactor = 0.95;
	bool firstImprovement = false;
	
	signal(SIGALRM, sigHandler);
	//
	alarm(TIME_LIMIT);
	// Getting info from arguments
	
	if(argc < 5 || !(argc%2)){
		cerr << "Please be more careful with arguments...\n";
		return -1;
	}
	for(int i = 1; i < argc; i += 2){
		if(!strcmp(argv[i], "-f")){
			filename = new char[strlen(argv[i+1]) + 1];
			strcpy(filename, argv[i+1]);
		}
		else if(!strcmp(argv[i], "-a")){
			algorithm = new char[strlen(argv[i+1]) + 1];
			strcpy(algorithm, argv[i+1]);
		}
//		else if(!strcmp(argv[i], "-l")){
//			localOpt = new char[strlen(argv[i+1]) + 1];
//			strcpy(localOpt, argv[i+1]);
//		}
		else if(!strcmp(argv[i], "-s")){
			sscanf(argv[i+1], "%u", &seed);
		}
		else if(!strcmp(argv[i], "-T")){
			T = atof(argv[i+1]);
		}
		else if(!strcmp(argv[i], "-t")){
			t = atof(argv[i+1]);
		}
		else if(!strcmp(argv[i], "-d")){
			decreaseFactor = atof(argv[i+1]);
		}
		else if(!strcmp(argv[i], "-n")){
			sscanf(argv[i+1], "%d", &numIter);
		}
		else if(!strcmp(argv[i], "-k")){
			sscanf(argv[i+1], "%d", &k);
		}	
		else if(!strcmp(argv[i], "-i")){
			firstImprovement = (!strcmp(argv[i+1], "true")) ? true : false;
		}
		else{
			cerr << "Wrong arguments were given...\nPlease retry...\n";
			return -1;
		}
	}		
	if(filename == NULL){
		cerr << "No file was given...\n";
		return -1;
	}	
	if(algorithm == NULL){
		cerr << "No algorithm was provided...\n";
		return -1;
	}
//	if(localOpt == NULL){
//		localOpt = new char[4];
//		strcpy(localOpt, "2opt");
//	}	
	
	//	Initiating program
	cout << "Initiating TSP program...\n" << endl;
	cout << "(Keep in mind that any algorithm will return the so far solution after approximately a 3-minute time period)" << endl;
	
	// Creating instance from file
	Instance *instance = new Instance(filename);
	cout << "\nNew instance was created with the following description: " << endl;
	cout << "Filename: " << instance->getName() << endl;
	cout << "Comment: " << instance->getComment() << endl;
	cout << "Type: " << instance->getType() << endl;
	cout << "Dimension: " << instance->dim() << endl;
	cout << "EdgeWeightType: " << instance->getEdgeWeightType() << endl;
	cout << "BestKnown: " << instance->getBestKnown() << endl;
	cout << endl;
	// Seed for random numbers
	if(!seed)
		seed = (unsigned int) time(NULL);
	srand(seed);
	cout << "Random seed to be used is " << seed << endl << endl;
	
	cout << "Waiting for algorithm to converge or time limit to be reached...\n";
	
	if(!strcmp(algorithm, "rand")){
		tempTime = time(NULL);
		instance->initialSolution();
		tempTime = time(NULL) - tempTime;
		cout << "------------------------------------------------------------------------------------------------" << endl;
		cout << "Cost of (Random) Initial Solution: " << (cost = instance->evaluateSolution()) << endl;
		cout << "Relative error: " << (cost - instance->getBestKnown())/((double) instance->getBestKnown())*100 << "%" << endl;
		cout << "Running time: " << tempTime << " seconds" << endl;
		cout << "------------------------------------------------------------------------------------------------" << endl;
	}
	else if(!strcmp(algorithm, "nn")){
		tempTime = time(NULL);
		instance->nearestNeighbour();
		tempTime = time(NULL) - tempTime;
		cout << "------------------------------------------------------------------------------------------------" << endl;
		cout << "Cost of Nearest Neighbour: " << (cost = instance->evaluateSolution()) << endl;
		cout << "Relative error: " << (cost - instance->getBestKnown())/((double) instance->getBestKnown())*100 << "%" << endl;
		cout << "Running time: " << tempTime << " seconds" << endl;
		cout << "------------------------------------------------------------------------------------------------" << endl;
	}
	else if(!strcmp(algorithm, "ni")){
		tempTime = time(NULL);
		instance->nearestInsertion();
		tempTime = time(NULL) - tempTime;
		cout << "------------------------------------------------------------------------------------------------" << endl;
		cout << "Cost of Nearest Insertion: " << (cost = instance->evaluateSolution()) << endl;
		cout << "Relative error: " << (cost - instance->getBestKnown())/((double) instance->getBestKnown())*100 << "%" << endl;
		cout << "Running time: " << tempTime << " seconds" << endl;
		cout << "------------------------------------------------------------------------------------------------" << endl;
	}	
	else if(!strcmp(algorithm, "fi")){
		tempTime = time(NULL);
		instance->farthestInsertion();
		tempTime = time(NULL) - tempTime;
		cout << "------------------------------------------------------------------------------------------------" << endl;
		cout << "Cost of Farthest Insertion: " << (cost = instance->evaluateSolution()) << endl;
		cout << "Relative error: " << (cost - instance->getBestKnown())/((double) instance->getBestKnown())*100 << "%" << endl;
		cout << "Running time: " << tempTime << " seconds" << endl;
		cout << "------------------------------------------------------------------------------------------------" << endl;
	}
	else if(!strcmp(algorithm, "2opt")){
		tempTime = time(NULL);
		instance->neighborList(k);
		instance->initialSolution();
		instance->opt2(firstImprovement, die);
		tempTime = time(NULL) - tempTime;
		cout << "------------------------------------------------------------------------------------------------" << endl;
		cout << "Cost of Local Search(2-opt): " << (cost = instance->evaluateSolution()) << endl;
		cout << "Relative error: " << (cost - instance->getBestKnown())/((double) instance->getBestKnown())*100 << "%" << endl;
		cout << "Running time: " << tempTime << " seconds" << endl;
		cout << "------------------------------------------------------------------------------------------------" << endl;
	}
	else if(!strcmp(algorithm, "3opt")){
		tempTime = time(NULL);
		instance->initialSolution();
		instance->opt3(firstImprovement, die);
		tempTime = time(NULL) - tempTime;
		cout << "------------------------------------------------------------------------------------------------" << endl;
		cout << "Cost of Local Search(3-opt): " << (cost = instance->evaluateSolution()) << endl;
		cout << "Relative error: " << (cost - instance->getBestKnown())/((double) instance->getBestKnown())*100 << "%" << endl;
		cout << "Running time: " << tempTime << " seconds" << endl;
		cout << "------------------------------------------------------------------------------------------------" << endl;
	}
	else if(!strcmp(algorithm, "rr")){
		tempTime = time(NULL);
//		if(!strcmp(localOpt, "2opt"))
			instance->neighborList(k);
			instance->randomRestart(/*&opt2,*/ firstImprovement, numIter, die);
//		else
//			instance->randomRestart(&opt3, firstImprovement, numIter);
		tempTime = time(NULL) - tempTime;
		cout << "------------------------------------------------------------------------------------------------" << endl;
		cout << "Cost of Random Restart(2opt): " << (cost = instance->evaluateSolution()) << endl;
		cout << "Relative error: " << (cost - instance->getBestKnown())/((double) instance->getBestKnown())*100 << "%" << endl;
		cout << "Running time: " << tempTime << " seconds" << endl;
		cout << "------------------------------------------------------------------------------------------------" << endl;
	}
	else if(!strcmp(algorithm, "ils")){
		cout << "Starting an iterated local search with initial temperature " << T << ", threshold " << t << " and decrease factor " << decreaseFactor << endl;
		cout << "Number of iterations used: " << numIter << endl;
		tempTime = time(NULL);
//		if(!strcmp(localOpt, "2opt"))
			instance->neighborList(k);
			instance->iteratedLocalSearch(/*&opt2,*/ firstImprovement, T, t, decreaseFactor, numIter, die);
//		else
//			instance->iteratedLocalSearch(&opt3, firstImprovement, T, t, decreaseFactor, numIter);
		tempTime = time(NULL) - tempTime;
		cout << "------------------------------------------------------------------------------------------------" << endl;
		cout << "Cost of Iterated Local Search(double bridge, 2opt): " << (cost = instance->evaluateSolution()) << endl;
		cout << "Relative error: " << (cost - instance->getBestKnown())/((double) instance->getBestKnown())*100 << "%" << endl;
		cout << "Running time: " << tempTime << " seconds" << endl;
		cout << "------------------------------------------------------------------------------------------------" << endl;
	}
	else{
		cerr << "None of the known algorithms was given..." << endl;
		delete []filename;
		delete []algorithm;
//		delete []localOpt;
		delete instance;
		return -1;
	}	
	
	delete []filename;
	delete []algorithm;
//	delete []localOpt;
	
	delete instance;
	
	return 0;
}
