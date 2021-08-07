#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <climits>
#include <fstream>
#include <iostream>
#include <algorithm>
#include "TSP.h"

using namespace std;

static const int BUF_SIZE = 128;

//********************* Instance class functions ***********************************************************//

Instance::Instance(char *fileName){	// Constructor - reads file, etc.                
	char *tempBuf;
	// Read file and keep information in storage
	ifstream inputFile(fileName);
	if(inputFile.fail()){
	      cerr << "Failed to open " << fileName << endl;
	      exit(EXIT_FAILURE);
	}
	// Reading initial info
	for(int i = 0; i < 6; ++i){
		inputFile.ignore(BUF_SIZE,':');
		inputFile.ignore(1,' ');
		tempBuf = new char[BUF_SIZE];
		inputFile.getline(tempBuf, BUF_SIZE, '\n');
		switch(i){
			case 0:	name = new char[inputFile.gcount() + 1];
					strncpy(name, tempBuf, strlen(tempBuf) + 1);
					break;
			case 1:	comment = new char[inputFile.gcount() + 1];
					strncpy(comment, tempBuf, strlen(tempBuf) + 1);
					break;
			case 2:	type = new char[inputFile.gcount() + 1]; 
					strncpy(type, tempBuf, strlen(tempBuf) + 1);
					break;
			case 3:	dimension = atoi(tempBuf); 
					break;
			case 4: edgeWeightType = new char[inputFile.gcount() + 1]; 
					strncpy(edgeWeightType, tempBuf, strlen(tempBuf) + 1);
					break;
			case 5:	bestKnown = atoi(tempBuf); 
					break;		
		}	
		delete []tempBuf;
	}
	// Reading nodes' coordinates
	inputFile.ignore(BUF_SIZE, '\n'); // ignore NODE_COORD_SECTION line
	for(int i = 0; i < dimension; ++i){
		inputFile.ignore(BUF_SIZE, ' ');
		pair<double, double> tempPair;
		tempBuf = new char[BUF_SIZE];
		inputFile.getline(tempBuf, BUF_SIZE, ' ');
		tempPair.first = strtod(tempBuf, NULL);
		delete []tempBuf;
		tempBuf = new char[BUF_SIZE];
		inputFile.getline(tempBuf, BUF_SIZE, '\n');
		tempPair.second = strtod(tempBuf, NULL);
		delete []tempBuf;
		coordinates.push_back(tempPair);
	}
	inputFile.close();
	if(inputFile.fail()){
		cerr << "Failed to close " << fileName << endl;
		exit(EXIT_FAILURE);
	}
	cout << "Successfully read information from " << fileName << endl;
	// Construct distances' matrix
	distances.resize(dimension);
	for(int i = 0; i < dimension; ++i){
		distances[i].resize(dimension);
	}
	for(int i = 0; i < dimension; ++i){
		for(int j = 0; j < i; ++j){
			distances[i][j] = nearInt(euclideanDist(coordinates[i], coordinates[j]));
			distances[j][i] = distances[i][j];
		}
		distances[i][i] = 0;
	}
	solution.resize(dimension);
	position.resize(dimension);
	///////////////////////////////////////////////////////////
	// Construct neighbor lists of fixed size k = 15
/*	this->k = 15;
	neighbor.resize(dimension);
	for(int i = 0; i < dimension; ++i){
		neighbor[i].resize(k);
		vector< pair<int, int> > neigh;
		for(int o = 0; o < dim(); ++o){
			if(o != i){
				pair<int, int> tempPair;
				tempPair.first = o;
				tempPair.second = distances[i][o];
				neigh.push_back(tempPair);
			}
		}
		for(int o = 0; o < this->k; ++o)
			neighbor[i][o] = specialMin(neigh);
	}
	*/		
} 

Instance::~Instance(){		// Destructor			
	delete []name;
	delete []comment;
	delete []type;
	delete []edgeWeightType;
	coordinates.clear();
	solution.clear();
	position.clear();
	for(int i = 0; i < dimension; ++i)
		distances[i].clear();
	distances.clear();
	if(!neighbor.empty()){
		for(int i = 0; i < dimension; ++i)
			neighbor[i].clear();
		neighbor.clear();
	}		
	
}

void Instance::neighborList(int num){
	this->k = num;
	neighbor.resize(dimension);
	for(int i = 0; i < dimension; ++i){
		neighbor[i].resize(k);
		vector< pair<int, int> > neigh;
		for(int o = 0; o < dim(); ++o){
			if(o != i){
				pair<int, int> tempPair;
				tempPair.first = o;
				tempPair.second = distances[i][o];
				neigh.push_back(tempPair);
			}
		}
		for(int o = 0; o < this->k; ++o)
			neighbor[i][o] = specialMin(neigh);
	}
}	

int Instance::specialMin(vector< pair<int,int> > &vec){
	int min = INT_MAX, argMin, toReturn;
	
	for(int i = 0; i < vec.size(); ++i){
		if(vec[i].second <= min){
			min = vec[i].second;
			argMin = i;
			toReturn = vec[i].first;
		}
	}
	vec.erase(vec.begin() + argMin);
	return toReturn;
}			

void Instance::initialSolution(){	// Return a random initial solution 
	for(int i = 0; i < dimension; ++i){
		solution[i] = i;
		position[i] = i;
	}	
	this->randPerm();
}

double Instance::euclideanDist(pair<double, double> pair1, pair<double, double> pair2){	
	return	sqrt((pair1.first - pair2.first)*(pair1.first - pair2.first) + (pair1.second - pair2.second)*(pair1.second - pair2.second));
}

int Instance::nearInt(double doub){
	double intPart, fractPart;
	
	fractPart = modf(doub, &intPart);
	if(fractPart < 0.5)
		return (int) intPart;
	else
		return (int) intPart + 1.0;
}		

int Instance::dim(){
	return this->dimension;
}	

int Instance::d(int i, int j){
	return this->distances[i][j];
}	

void Instance::randPerm(){		
	// Knuth shuffle method
	for(int i = 0; i < dimension; ++i){
		int random = (rand() % (dimension - i)) + i;	// produce a random position in [i,dimension-1]
		swap(i, random);
	}
}	

void Instance::swap(int i, int j){	// swap vertices in positions i and j	
	int temp;
	
	temp = solution[j];
	solution[j] = solution[i];
	solution[i] = temp;
	
	position[solution[j]] = j;
	position[solution[i]] = i;
}

int Instance::evaluateSolution(){
	int sum = 0;
	
	for(int i = 0; i < dimension; ++i)
		sum += d(solution[i],solution[(i+1)%dimension]);
	return sum;	
}

int Instance::evaluateSolution(vector<int> v){ 	
	int sum = 0;
	
	for(int i = 0; i < dimension; ++i)
		sum += d(v[i],v[(i+1)%dimension]);
	return sum;
}	

int Instance::next(int i){	
	return solution[(position[i] + 1)%dimension];
}

int Instance::prev(int i){
	return solution[(position[i] - 1 + dimension)%dimension];
}		

bool Instance::between(int i, int j, int k){
	return (position[i] < position[k] && position[k] < position[j]);//
}

void Instance::shift(int i, int j){	// Delete the vertex i from its current position and insert it after vertex j	
	int p, nextP;
	
	p = i;
	do{
		nextP = next(p);
		swap(position[p], position[nextP]);
	}while(nextP != j);
}
	
void Instance::invert(int i, int j){	// Invert the segment starting at vertex i and ending at vertex j	
	int dist = (position[j] < position[i]) ? (position[i] - position[j]) : (position[j] - position[i]);
	int plusOne = 0;
	int start = position[i];
	int end = position[j];
	
	if(dist%2 == 1)
		plusOne = 1;
	for(int k = 0; k < (dist/2 + plusOne); k++)
		swap((start + k)%dimension, (end - k + dimension)%dimension);
}	

void Instance::nearestNeighbour(){	
	bool *added = new bool[dim()];
	
	for(int i = 0; i < dim(); ++i)
		added[i] = false;
	solution[0] = rand() % dim();
	position[solution[0]] = 0;
	added[solution[0]] = true;
	for(int i = 1; i < dim(); ++i){
		solution[i] = closestNotAddedVertex(solution[i-1], added);
		position[solution[i]] = i;
		added[solution[i]] = true;
	}
	delete []added;
}

int Instance::closestNotAddedVertex(int v, bool *added){
	int arg = -1;
	int min = INT_MAX;
	
	for(int i = 0; i < dim(); ++i){
		if(i != v){	// do not look on diagonal
			if(d(v,i) <= min && !added[i]){
				min = d(v,i);
				arg = i;
			}	
		}
	}	
	return arg;
}

void Instance::farthestInsertion(){	
	int *distance = new int[dim()];
	int pos;
	
	for(int i = 0; i < dim(); ++i)
		distance[i] = INT_MAX;
	
	for(int i = 0; i < dim(); ++i){
		if(!i){
			pos = rand() % dim();
			solution[i] = pos;
			position[pos] = i;
			distance[pos] = -1;
		}	
		else{
			pos = argMax(distance); // pick farthest	
			// find where to add node
			int min = INT_MAX;
			int aMin;
			for(int j = 0; (j+1) < i; ++j){
				if((d(solution[j],pos) + d(pos, solution[j+1]) - d(solution[j], solution[j+1])) < min){
					min = d(solution[j],pos) + d(pos, solution[j+1]) - d(solution[j], solution[j+1]);
					aMin = j;
				}
			}
			if(d(solution[i-1], pos) < min){
				aMin = i-1;
			}
			for(int j = i-1 ; j > aMin; --j){
				solution[j+1] = solution[j];
				position[solution[j+1]] = j+1;
			}
			solution[aMin+1] = pos;
			position[pos] = aMin + 1;
			distance[pos] = -1;	 
		}				
		// update distance table
		for(int j = 0; j < dim(); ++j){
			distance[j] = (distance[j] < d(pos, j)) ? distance[j] : d(pos, j);	
		}
	}
	delete []distance;
}

void Instance::nearestInsertion(){	
	int *distance = new int[dim()];
	int pos;
	
	for(int i = 0; i < dim(); ++i)
		distance[i] = INT_MAX;
	
	pos = rand() % dim();
	solution[0] = pos;
	position[pos] = 0;
	distance[pos] = -1;
	// update distance table
	for(int j = 0; j < dim(); ++j){
		distance[j] = (distance[j] < d(pos, j)) ? distance[j] : d(pos, j);	
	}
	for(int i = 1; i < dim(); ++i){	// main loop
		pos = argMin(distance); // pick nearest	
		// find where to add node
		int min = INT_MAX;
		int aMin;
		for(int j = 0; j < i-1; ++j){
			if((d(solution[j],pos) + d(pos, solution[j+1]) - d(solution[j], solution[j+1])) < min){
				min = d(solution[j],pos) + d(pos, solution[j+1]) - d(solution[j], solution[j+1]);
				aMin = j;
			}
		}
		if(d(solution[i-1], pos) < min){
			aMin = i-1;
		}
		for(int j = i-1 ; j > aMin; --j){
			solution[j+1] = solution[j];
			position[solution[j+1]] = j+1;
		}
		solution[aMin + 1] = pos;
		position[pos] = aMin + 1;
		distance[pos] = -1;
		// update distance table
		for(int j = 0; j < dim(); ++j){
			distance[j] = (distance[j] < d(pos, j)) ? distance[j] : d(pos, j);	
		}	 
	}				
	
	delete []distance;
}

int Instance::argMax(int *a){	
	int max = INT_MIN;
	int arg;
	
	for(int i = 0; i < dim(); ++i){
		if(max <= a[i] && a[i] != -1){
			max = a[i];
			arg = i;
		}
	}
	return arg;
}

int Instance::argMin(int *a){	
	int min = INT_MAX;
	int arg;
	
	for(int i = 0; i < dim(); ++i){
		if(a[i] <= min && a[i] != -1){
			min = a[i];
			arg = i;
		}
	}
	return arg;
}	


void Instance::opt2(bool firstImprovement, bool &die){
	int bestGain;
	int gain, best_i, best_j;
/*	
	do{
		bestGain = 0;
		for(int i = 0; i < dim() && !die; ++i){
			for(int j = 0; j < dim() && !die; ++j){
				if(distanceConstraints(i, j)){
					gain = d(solution[i], solution[j]) + d(solution[(i+1)%dim()], solution[(j+1)%dim()]) - d(solution[(i+1)%dim()], solution[i]) - d(solution[(j+1)%dim()], solution[j]);
					if(gain < bestGain){
						bestGain = gain;
						best_i = (i<j) ? solution[i] : solution[j];
						best_j = (i<j) ? solution[j] : solution[i];
						if(firstImprovement) break;
					}
				}
			}
			if(bestGain < 0 && firstImprovement) break;
		}
		if(bestGain)
			invert(next(best_i), best_j);
	}while(bestGain && !die);
*/
	
	do{
		bestGain = 0;
		for(int i = 0; i < dim() && !die; ++i){
			for(int j = 0; j < this->k && !die; ++j){
				if(d(solution[(i+1)%dim()], neighbor[solution[(i+1)%dim()]][j]) >= d(solution[i], solution[(i+1)%dim()]))
					break;
				int m = (position[neighbor[solution[(i+1)%dim()]][j]] - 1 + dim())%dim();	
				gain = d(solution[i], solution[m]) + d(solution[(i+1)%dim()], solution[(m+1)%dim()]) - d(solution[(i+1)%dim()], solution[i]) - d(solution[(m+1)%dim()], solution[m]);
				if(gain < bestGain){
					bestGain = gain;
					best_i = (i < m) ? solution[i] : solution[m];
					best_j = (i < m) ? solution[m] : solution[i];
					if(firstImprovement) break;
				}
			}
			if(bestGain < 0 && firstImprovement) break;
		}
		if(bestGain)
			invert(next(best_i), best_j);
	}while(bestGain && !die);			
}

void Instance::opt3(bool firstImprovement, bool &die){
	int bestGain, gain1, gain2, minGain, argMinGain, best_i, best_j, best_k;
		
	do{
		bestGain = 0;
		for(int i = 0; i < dim() && !die; ++i){
			for(int j = 0; j < dim() && !die; ++j){
				for(int k = 0; k < dim() && !die; ++k){
					if(distanceConstraints(i, j, k)){
						gain1 = d(solution[i], solution[(j+1)%dim()]) + d(solution[j], solution[(k+1)%dim()]) + d(solution[k], solution[(i+1)%dim()]) - d(solution[(i+1)%dim()], solution[i]) - d(solution[(j+1)%dim()], solution[j]) - d(solution[(k+1)%dim()], solution[k]);
						gain2 = d(solution[i], solution[(j+1)%dim()]) + d(solution[j], solution[k]) + d(solution[(k+1)%dim()], solution[(i+1)%dim()]) - d(solution[(i+1)%dim()], solution[i]) - d(solution[(j+1)%dim()], solution[j]) - d(solution[(k+1)%dim()], solution[k]);
						minGain = (gain1 <= gain2) ? gain1 : gain2;
						if(minGain < bestGain){
							bestGain = minGain;
							argMinGain = (gain1 <= gain2) ? 1 : 2;
							vector<int> toSort;
							toSort.push_back(i);toSort.push_back(j);toSort.push_back(k);
							sort(toSort.begin(), toSort.end());
							best_i = toSort[0];
							best_j = toSort[1];
							best_k = toSort[2];
							toSort.clear();
							if(firstImprovement) break;
						}
					}
				}
				if(bestGain < 0 && firstImprovement) break;
			}
			if(bestGain < 0 && firstImprovement) break;
		}	
		if(bestGain){
			vector<int>	iToj(solution.begin() + best_i + 1, solution.begin() + best_j + 1);
			vector<int> jTok(solution.begin() + best_j + 1, solution.begin() + best_k + 1);
			
			for(int o = 0; o < jTok.size(); ++o){
				solution[best_i + 1 + o] = jTok[o];
				position[solution[best_i + 1 + o]] = best_i + 1 + o;
			}
			for(int o = 0; o < iToj.size(); ++o){
				solution[best_i + 1 + jTok.size() + o] = iToj[o];
				position[solution[best_i + 1 + jTok.size() + o]] = best_i + 1 + jTok.size() + o;
			}
			if(argMinGain == 2)
				invert(solution[best_i + 1 +jTok.size()], solution[best_i + jTok.size() + iToj.size()]);
		}
	}while(bestGain && !die);
}	

void Instance::doubleBridgeMove(bool firstImprovement, bool &die){
	int bestGain = 0, gain, count = 0;
	int best_i, best_j, best_k, best_m;
/*	
	for(int i = 0; i < dim() && !die; ++i){
		for(int j = 0; j < dim() && !die; ++j){
			for(int k = 0; k < dim() && !die; ++k){
				for(int m = 0; m < dim() && !die; ++m){
					if(distanceConstraints(i, j, k, m)){
						gain = d(solution[i], solution[(k+1)%dim()]) + d(solution[(i+1)%dim()], solution[k]) + d(solution[(j+1)%dim()], solution[m]) + d(solution[(m+1)%dim()], solution[j]) - d(solution[(i+1)%dim()], solution[i]) - d(solution[(j+1)%dim()], solution[j]) - d(solution[(k+1)%dim()], solution[k]) - d(solution[(m+1)%dim()], solution[m]);
						if(gain < bestGain){
							bestGain = gain;
							vector<int> toSort;
							toSort.push_back(i);toSort.push_back(j);toSort.push_back(k);toSort.push_back(m);
							sort(toSort.begin(), toSort.end());
							best_i = toSort[0];
							best_j = toSort[1];
							best_k = toSort[2];
							best_m = toSort[3];
							toSort.clear();
							if(firstImprovement) break;
						}
					}
				}
				if(bestGain < 0 && firstImprovement) break;
			}
			if(bestGain < 0 && firstImprovement) break;
		}
		if(bestGain < 0 && firstImprovement) break;
	}
*/
	int i = 0, j = 0, k = 0, m = 0;
	while(!distanceConstraints(i, j, k, m)){
		i = rand()%dim();
		j = rand()%dim();
		k = rand()%dim();
		m = rand()%dim();
	}
	vector<int> toSort;
	toSort.push_back(i);toSort.push_back(j);toSort.push_back(k);toSort.push_back(m);
	sort(toSort.begin(), toSort.end());
	best_i = toSort[0];
	best_j = toSort[1];
	best_k = toSort[2];
	best_m = toSort[3];
	toSort.clear();
	
	bestGain = -1;//
	if(bestGain){//
//		cout << "Ready to perform double bridge\n";
		vector<int>	iToj(solution.begin() + best_i + 1, solution.begin() + best_j + 1);
		vector<int> jTok(solution.begin() + best_j + 1, solution.begin() + best_k + 1);
		vector<int> kTom(solution.begin() + best_k + 1, solution.begin() + best_m + 1);
		
		for(int o = 0; o < kTom.size(); ++o){
			solution[best_i + 1 + o] = kTom[o];
			position[solution[best_i + 1 + o]] = best_i + 1 + o;
		}		
		for(int o = 0; o < jTok.size(); ++o){
			solution[best_i + 1 + kTom.size() + o] = jTok[o];
			position[solution[best_i + 1 + kTom.size() + o]] = best_i + 1 + kTom.size() + o;
		}
		for(int o = 0; o < iToj.size(); ++o){
			solution[best_i + 1 + kTom.size() + jTok.size() + o] = iToj[o];
			position[solution[best_i + 1 + kTom.size() + jTok.size() + o]] = best_i + 1 + kTom.size() + jTok.size() + o;
		}
//		cout << " Double Bridge done...\n";	
	}
}

void Instance::randomRestart(/*void (Instance::*ls)(bool),*/ bool firstImprovement, int numIter, bool &die){
	int min = INT_MAX;
	vector<int> bestSol(dim());
	
	for(int i = 0; i < numIter && !die; ++i){
		initialSolution();
		//(*ls)(firstImprovement);
		opt2(firstImprovement, die);
		if(evaluateSolution() < min){
			min = evaluateSolution();
			bestSol = getSolution();
		}	
	}
	solution = bestSol;

}

void Instance::iteratedLocalSearch(/*void (Instance::*ls)(bool),*/ bool fi, double T, double frozenThreshold, double decreaseFactor, int numIterations, bool &die){ 
	int delta;
	time_t tempTime;
	
	tempTime = time(NULL);
	initialSolution();
	vector<int> sol = getSolution();
	vector<int> pos = getPosition();
	vector<int> bestSol = getSolution();
	vector<int> bestPos = getPosition();
	while(T > frozenThreshold && !die){
		for(int i = 0; i < numIterations && !die ; ++i){
			doubleBridgeMove(fi, die);
			//(*this).*ls(fi);
			opt2(fi, die);
			//opt3(fi, die);
			delta = evaluateSolution() - evaluateSolution(sol);
			if((delta <= 0)){
				sol = solution;
				pos = position;
				bestSol = solution;
				bestPos = position;
			}	
			else if(delta > 0 && (rand()/(double) RAND_MAX) < exp(-delta/T)){
				sol = solution;
				pos = position;
			}
			else{
				solution = sol;
				position = pos;
			}	
		}
		T *= decreaseFactor;
		//numIterations = (int) ceil(numIterations * 1.1);//
	}
	solution = bestSol;
	position = bestPos;
	cout << "Temperature when exiting is " << T << endl;
}

void Instance::antColonization(int m, bool &die){
	int *position = new int[m];
	

	while(1){
		// Randomly position m agents on n cities
		for(int ant = 0; ant < m; ++ant)
			position[ant] = rand()%dim();
		for(int i = 0; i < dim() && !die; ++i){
			for(int ant = 0; ant < m && !die; ++ant){
				// Apply state transition rule
				// Apply local trail updating rule
			}
		}
		// local search - not necessary
		// Apply global trail updating rule
	}
	delete []position;
}

bool Instance::distanceConstraints(int a, int b){
	if(a == b) return false;
	if(fabs((double)(a-b)) < 2) return false;
	return true;
}	

bool Instance::distanceConstraints(int a, int b, int c){
	if(a == b) return false;
	if(fabs((double)(a-b)) < 2) return false;
	if(a == c) return false;
	if(fabs((double)(a-c)) < 2) return false;
	if(b == c) return false;
	if(fabs((double)(c-b)) < 2) return false;
	return true;
}	

bool Instance::distanceConstraints(int a, int b, int c, int d){
	if(a == b) return false;
	if(fabs((double)(a-b)) < 2) return false;
	if(a == c) return false;
	if(fabs((double)(a-c)) < 2) return false;
	if(a == d) return false;
	if(fabs((double)(a-d)) < 2) return false;
	if(b == c) return false;
	if(fabs((double)(c-b)) < 2) return false;
	if(b == d) return false;
	if(fabs((double)(d-b)) < 2) return false;
	if(c == d) return false;
	if(fabs((double)(c-d)) < 2) return false;
	return true;
}

// accessors
char* Instance::getName(){ return this->name;}
char* Instance::getComment(){ return this->comment;}
char* Instance::getType(){ return this->type;}
char* Instance::getEdgeWeightType(){ return this->edgeWeightType;}
int Instance::getBestKnown(){ return this->bestKnown;}
vector< pair<double, double> > Instance::getCoordinates(){ return this->coordinates;}
vector< vector<int> > Instance::getDistances(){ return this->distances;}
vector<int> Instance::getSolution(){ return this->solution;}
vector<int> Instance::getPosition(){ return this->position;}
