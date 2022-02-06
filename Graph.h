#pragma once
#include <iostream>
#include <vector>
#include <math.h>
#include <list>
#include <iterator>
#include <random>

using namespace std;
class Graph
{
private:
	 vector<vector<int>> adjencyList; //adjency list matrix
	 vector<int> indexes;//will help us choose vertex which is not neccesserily an adjent to another vertex
	 int n;

public:
	 //Constructors:
	 //normal constructor:
	 Graph(int dim);
	 //Generate graph for first family to discuss (Q6 and on)
	 Graph(int dim, double p);

	 double CalculateAverageRank();
	 int RandomWalk(int& N, double& epsilon);
	 vector<double> t_times_RandomWalk(int N, int t, double epsilon);
	 vector<double> divide_vector_by_t(vector<double>& vector, int t);
	 int ChooseNeighborVertex(int v);
	 int ChooseRandomVertex();
	 double GenerateFraction();
	 vector<double> quest_1(int N, int t, double epsilon);
	 void instert_neigbor(int add_to, int add);
	 //quest_1 related functions:
	 void InitLists();
	 void check_if_probability_vector(vector<double>& d_array);

	 //First Family related functions:
	 void InitFirstFamily(double p);
	 vector<double> Quest_7(int N, double epsilon = 0.0, bool quest_8 = false);
	 double Calculate_Distance(vector<double>& current, vector<double>& previous);

	 //Second Family related functions:
	 void JoinRingGraph(int connect_number);
	 int ChooseRandomVertexIndex(int min, int max);
	 void check_average_of_two_groups(int rings_end, vector<double> average_vector);
};