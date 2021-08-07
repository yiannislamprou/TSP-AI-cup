#ifndef TSP_H_
#define TSP_H_

#include <vector>
#include <utility>

using namespace std;

class Instance{
	
	private:
		// info from file
		char *name;
		char *comment;
		char *type;
		int dimension;
		char *edgeWeightType;
		int bestKnown;
		vector< pair<double, double> > coordinates;
		
		// instance representation
		vector< vector<int> > distances;
		vector<int> solution; // for simple or array representation
		vector<int> position; // for array representation
		
		
		int k; // keeping k-nearest neighbors for each city
		vector< vector<int> > neighbor;
	
	public:
		Instance(char*);
		~Instance();
		void initialSolution();
		
		// auxiliary functions
		double euclideanDist(pair<double, double>, pair<double, double>);
		int nearInt(double);
		void randPerm();
		void swap(int, int);
		int closestNotAddedVertex(int, bool*);
		int argMax(int*);
		int argMin(int*);
		bool distanceConstraints(int, int);
		bool distanceConstraints(int, int, int);
		bool distanceConstraints(int, int, int, int);
		void neighborList(int);
		int specialMin(vector< pair<int, int> >&); 
		
		// solution evaluation
		int evaluateSolution();
		int evaluateSolution(vector<int>);
		
		// operations
		int next(int);
		int prev(int);
		bool between(int, int, int);
		void shift(int, int);
		void invert(int, int);
		
		// heuristics
		void nearestNeighbour();
		void farthestInsertion();
		void nearestInsertion();
		
		// local search
		void opt2(bool, bool&);
		void opt3(bool, bool&);
		void doubleBridgeMove(bool, bool&);
		
		// Metaheuristics
		void randomRestart(/*void (Instance::*ls)(bool),*/ bool firstImprovement, int numIter, bool &die);
		void iteratedLocalSearch(/*void (Instance::*ls)(bool),*/ bool firstImprovement, double T, double t, double ft, int numIter, bool &die);
		void antColonization(int numAnts, bool &die);
		
		// accessors
		char* getName();
		char* getComment();
		char* getType();
		char* getEdgeWeightType();
		int getBestKnown();
		int dim();
		int d(int, int);
		vector< pair<double, double> > getCoordinates();
		vector< vector<int> > getDistances();
		vector<int> getSolution();
		vector<int> getPosition();
		
};	

#endif // TSP_H_
