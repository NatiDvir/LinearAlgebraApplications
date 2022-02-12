#include "Graph.h"
using namespace std;
//Generate manual graph:
Graph::Graph(int dim)
{
	 this->n = dim;
	 InitLists();
}

//Generate graph for first family to discuss (Q6 and on)
Graph::Graph(int dim, double p)
{
	 this->n = dim;
	 adjencyList = vector<vector<int>>(n, vector<int>());
	 indexes = vector<int>(n, 0);
	 InitFirstFamily(dim, p);
}

Graph::Graph(int dim, double p, int non_ring_index)
{
	 this->n = dim;
	 adjencyList = vector<vector<int>>(n, vector<int>());
	 indexes = vector<int>(n, 0);
	 InitFirstFamily(non_ring_index, p);
	 JoinRingGraph(non_ring_index);
	 int outside_ring = ChooseRandomVertexIndex(0, non_ring_index);
	 int from_ring = ChooseRandomVertexIndex(non_ring_index + 1, n);
	 adjencyList[outside_ring].push_back(from_ring);
}

Graph::Graph(int dim, bool nine_quest)
{
	 this->n = dim;
	 adjencyList = vector<vector<int>>(n, vector<int>());
	 indexes = vector<int>(n, 0);
	 InitThirdFamily(dim);
}

vector<double> Graph::quest_1(int N, int t, double epsilon)
{
	 vector<double> d_array = t_times_RandomWalk(N, t, epsilon);
	 vector<double> divided_vector = divide_vector_by_t(d_array, t);

	 check_if_probability_vector(divided_vector);
	 return divided_vector;
}

vector<double> Graph::t_times_RandomWalk(int N, int t, double epsilon)
{
	 vector<double> d_array = vector<double>(n, 0);
	 int random_walk_result;
	 for (int i = 0; i < t; i++)
	 {
		  random_walk_result = RandomWalk(N, epsilon);
		  d_array[random_walk_result] += 1.0;
	 }
	 return d_array;
}

vector<double> Graph::divide_vector_by_t(vector<double>& v, int t)
{
	 vector<double> d_array = vector<double>(n, 0);
	 for (int i = 0; i < n; i++) //Divide all vectors values by t
	 {
		  d_array[i] = v[i] / (double)t;
	 }
	 return d_array;
}

int Graph::RandomWalk(int& N, double& epsilon)
{
	 double fraction;
	 int current_vertex = ChooseRandomVertex();
	 for (int step = 0; step < N; step++)
	 {
		  fraction = GenerateFraction();
		  if (fraction >= epsilon && adjencyList[current_vertex].size() > 0)
		  {
			   current_vertex = ChooseNeighborVertex(current_vertex);
		  }
		  else// 1-epsilon change or current_vertex doesn't have neighbors
		  {
			   current_vertex = ChooseRandomVertex();
		  }
	 }
	 return current_vertex;
}

int Graph::ChooseNeighborVertex(int v)
{
	 vector<int>neighbors_of_v = vector<int>(adjencyList[v]);
	 random_shuffle(neighbors_of_v.begin(), neighbors_of_v.end());
	 return neighbors_of_v[0];
}

int Graph::ChooseRandomVertex()
{
	 random_shuffle(indexes.begin(), indexes.end());
	 return indexes[0];
}

void Graph::check_if_probability_vector(vector<double>& d_array)
{
	 double prob_sum = 0;

	 for (int i = 0; i < n - 1; i++)
	 {
		  prob_sum += d_array[i];
		  cout << d_array[i] << " + ";
	 }
	 prob_sum += d_array[n - 1];
	 cout << d_array[n - 1] << " = " << prob_sum << endl;
}

void Graph::InitLists()
{
	 vector<vector<int>> adjencyMatrix = vector<vector<int>>(n, vector<int>(n, 0));
	 adjencyList = vector<vector<int>>(n, vector<int>());
	 indexes = vector<int>(n, 0);
	 int row, col;
	 cin >> row >> col;
	 while (row >= 0 && col >= 0)
	 {
		  if (row >= n || col >= n)
		  {
			   cout << "Out of verteces range.";
			   exit(1);
		  }
		  adjencyMatrix[row][col] = 1;
		  cin >> row >> col;
	 }
	 for (int i = 0; i < n; i++)
	 {
		  for (int j = 0; j < n; j++)
		  {
			   if (adjencyMatrix[i][j] == 1)
					adjencyList[i].push_back(j);
		  }
		  indexes[i] = i;
	 }
}

void Graph::InitFirstFamily(int dim, double p)
{
	 for (int from_vertex = 0; from_vertex < dim; from_vertex++)
	 {
		  for (int to_vertex = 0; to_vertex < dim; to_vertex++)
		  {
			   double current_probabilty = GenerateFraction();
			   if (current_probabilty <= p)
			   {
					adjencyList[from_vertex].push_back(to_vertex);
			   }
			   /*Else -> Does nothing...*/
		  }
	 }
}

vector<double> Graph::Quest_7(int N, double epsilon, bool quest_8)
{
	 bool quest_7 = !quest_8;
	 int t = 1, cntr = 0;
	 int itr = 1;
	 //for t = 2:
	 vector<double> d_previous = vector<double>(),
		  d_previous_divided = vector<double>();
	 //for t = 4:
	 vector<double> d_current = t_times_RandomWalk(N, 2, epsilon);
	 vector<double> d_current_divided = divide_vector_by_t(d_current, 2);

	 double distance = 1;
	 while (distance >= (1.0 / 256.0))
	 {
		  if (quest_8 == true && t > pow(2, 14))
			   break; //In case we're in quest 8 and the loop didn't stop - break out.
		  t *= 2;
		  d_previous.swap(d_current); d_previous_divided.swap(d_current_divided);
		  d_current.clear(); d_current_divided.clear();
		  d_current = t_times_RandomWalk(N, t, epsilon);
		  d_current_divided = divide_vector_by_t(d_current, t);
		  distance = Calculate_Distance(d_current_divided, d_previous_divided);
		  cout << "t = " << t << ";		iteration = " << itr << ";	  Distance = " << distance << endl;
		  itr++;
	 }
	 //FOR QUEST 7: if out -> distance < 1.0 / 256.0
	 double avg = 0.0;
	 for (int i = 0; i < n; i++)
	 {
		  avg += d_current[i];
	 }
	 avg = avg / (double)n;
	 cout << "average d is: " << avg << endl;

	 //FOR QUEST 8: if out -> Check distance.
	 if (quest_8)
	 {
		  cout << "As requested for Task #8: Distance noted in last iteration of the loop." << endl;
	 }
	 return d_current_divided;
}

double Graph::Calculate_Distance(vector<double>& current, vector<double>& previous)
{
	 double distance = 0;
	 for (int i = 0; i < n; i++)
	 {
		  distance += pow(current[i] - previous[i], 2.0);
	 }
	 return sqrt(distance);
}

void Graph::JoinRingGraph(int ring_index)
{
	 for (int i = ring_index; i < n - 1; i++)
	 {
		  adjencyList[i].push_back(i + 1);
	 }
	 adjencyList[n - 1].push_back(ring_index);
}

int Graph::ChooseRandomVertexIndex(int min, int max)
{
	 std::random_device rd;     // only used once to initialise (seed) engine
	 std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
	 std::uniform_int_distribution<int> uni(min, max); // guaranteed unbiased
	 auto random_integer = uni(rng);
	 return random_integer;
}

void Graph::check_average_of_two_groups(int rings_start, vector<double> average_vector)
{
	 cout << "Ring's group average rank: ";
	 double sum = 0.0, ring_avg, non_ring_avg;
	 //Calculate average rank of ring's verteces:
	 cout << rings_start << endl;
	 for (int i = rings_start; i < n; i++)
	 {
		  sum += average_vector[i];
	 }
	 non_ring_avg = sum / (double)(n - rings_start);
	 cout << non_ring_avg << endl;
	 sum = 0.0;
	 //Calculate average rank of non-ring's verteces:
	 for (int i = 0; i < rings_start; i++)
	 {
		  sum += average_vector[i];
	 }
	 ring_avg = sum / (double)(rings_start);
	 cout << "Non-Ring's group average rank: " << non_ring_avg << endl;
}

void Graph::InitThirdFamily(int dim)
{
	 for (int i = 0; i < n; i++)
	 {
		  for (int j = 0; j < dim; j++)
		  {
			   double p = 1 / (log2(j + 1));
			   double current_probabilty = GenerateFraction();
			   if (current_probabilty <= p)
			   {
					adjencyList[i].push_back(j);
			   }
			   /*Else -> Does nothing...*/
		  }
	 }
}

double Graph::CalculateAverageRank()
{
	 double sum = 0;
	 for (int i = 0; i < n; i++)
	 {
		  //Each "out" edge of a vertex enters to another vertex,
		  //therefore increase their rank by one.
		  sum += 2.0 * adjencyList[i].size();
	 }
	 return sum / (double)n;
}

double Graph::GenerateFraction()
{
	 double mone = (double)rand();
	 return mone / (double)RAND_MAX; //generate random fraction
}