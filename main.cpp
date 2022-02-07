#include <iostream>
#include <algorithm>
#include <math.h>
#include <vector>
#include "Graph.h"
using namespace std;

void printVector(vector<double> v)
{
	 for (int i = 0; i < v.size(); i++)
	 {
		  cout << "d[" << i << "] = " << v[i] << endl;
	 }
	 cout << endl;
}

void first_part()
{
	 int vnum = 6;
	 Graph g(vnum);
	 vector<double> dest_vector1, dest_vector2;
	 vector<double> avg_vector1 = vector<double>(vnum, 0),
		  avg_vector2 = vector<double>(vnum, 0);
	 int N1 = pow(2, 17), t = 1024, exp_itr = 1;
	 int N2 = 10;
	 double frac1 = 23.0 / 100.0;
	 double frac2 = 400.0 / 2000.0;
	 cout << "N = " << N1 << endl;
	 cout << "t = " << t << endl;
	 cout << "Epsilon = " << frac1 << endl;
	 for (int i = 0; i < exp_itr; i++)
	 {
		  //t = (i + 1) * 5;
		  cout << "t is: " << t << endl;
		  dest_vector1 = g.quest_1(N1, t, frac1);
		  printVector(dest_vector1);
		  cout << endl;
		  //dest_vector2 = g.quest_1(N2, t, frac2);
		  //printVector(dest_vector2);
		  for (int j = 0; j < vnum; j++)
		  {
			   avg_vector1[j] += dest_vector1[j];
			   //avg_vector2[j] += dest_vector2[j];
		  }
	 }

	 for (int j = 0; j < vnum; j++)
	 {
		  avg_vector1[j] /= (double)exp_itr;
		  //avg_vector2[j] /= (double)exp_itr;
	 }
	 cout << "**********************************************************************************\n";
	 cout << "AVG VECTOR1" << endl;
	 printVector(avg_vector1);
	 cout << "**********************************************************************************\n";
	 //cout << "AVG VECTOR2" << endl;
	 //printVector(avg_vector2);
}

void Quest_7()
{
	 int vnum = 1024;
	 double p = 1.0 / 64.0;
	 Graph g(vnum, p);
	 cout << "Average vertex rank: " << g.CalculateAverageRank() << endl;
	 /***********************************************************************
	 The averge is almost 32, which is 2^5.
	 The probability of edge to be added to the graph is 1/64,
	 which equals to 1/(2^6).
	 ************************************************************************/
	 int N = 100;
	 vector<double> dest_vector1 = g.Quest_7(N);
}

void Quest_8_TemplatedFunction(int N)
{
	 int vnum = 1024;
	 double p = 1.0 / 64.0;
	 Graph g(vnum, p);
	 vector<double> vec = g.Quest_7(N, 0.0, false);
	 printVector(vec);
}

void Quest_8()
{
	 cout << "Quest 8 Part 1:" << endl;
	 Quest_8_TemplatedFunction(/*N = */64);
	 cout << "Quest 8 Part 2:" << endl;
	 Quest_8_TemplatedFunction(/*N = */128);
}

void Quest_9()
{
	 //Choose random vertex outside ring and one inside ring,
	 //and add the inside-one to be neighbor of the outside one (means there's an edge outside->inside)
	 int vnum = 1024, N = 650;
	 double p = 1.0 / 64.0;
	 Graph g(vnum, p, 1024 - 64);
	 vector<double> vec = g.Quest_7(N, 0.0);
	 printVector(vec);
	 g.check_average_of_two_groups(1024 - 64, vec);
}

int main()
{
	 Quest_9();
	 return 0;
}