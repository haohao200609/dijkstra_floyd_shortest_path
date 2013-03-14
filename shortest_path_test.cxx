/**
 * Shortest Path implementation using Dijkastra & Floyd Algorithm
 * This is the Unit Test for the implementation
 * All the test cases would be run in this Unit Test file
 *
 * File: shortest_path_test.cxx
 * Author: iJab(Zhan Caibao) zhancaibaoATgmail.com
 * Date: 2013/02/20
 *
 */

#include <UnitTest++.h>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <limits>

#include "graph.hxx"

using namespace std;

double d_rand(double d_min, double d_max)
{
    double f = (double)rand() / RAND_MAX;
    return d_min + f * (d_max - d_min);
}

void add_edges(graph *g, int n)
{
	int r_edges = n - 1; //rand() % n;  // random number 0 to n-1, try to add edges to this number

	// Random edges to add
	/*
	for(int i = 0; i <= r_edges; ++i)
	{
		int src = rand() % n;
		int dst = rand() % (src+1);
		double w = d_rand(1.0, 100.0);
		g->add_edge(src, dst, w);
		g->add_edge(dst, src, w);
	}
	*/
	for(int i = 0; i < r_edges; ++i)
	{
		for(int j = i+1; j <= r_edges; ++j)
		{
			if(i != j && i < j)
			{
				double cost = d_rand(1.0, 100.0);
				g->add_edge(i, j, cost);
				g->add_edge(j, i, cost);
			}
		}
	}
}


// Test Case 1: 1000 nodes for two-dimension array based dijkarstra algorithm
// 1000
TEST(TwoArrayDij1000)
{
	cout << "Init Two-dimension array based graph ...";
	print_current_time();
	cout << endl;

	int n = 1000;
	twodarray_graph _aTest(n);
	add_edges(&_aTest, n);

	cout << "Init  Two-dimension array based graph completed.";
	print_current_time();
	cout << endl;

	int path_ct = _aTest.sp_dij_all_paths();
	CHECK(path_ct + 1);

	cout << "path number:" << path_ct << endl;
	//_aTest.print_result();
}

// Test Case 1: 100 nodes for two-dimension array based Floyd algorithm
// 1000 
TEST(TwoArrayFloyd1000)
{
	int n = 2000;
	twodarray_graph _aTest(n);
	add_edges(&_aTest, n);
	int path_ct = _aTest.sp_floyd();
	CHECK(path_ct + 1);

	cout << "path number:" << path_ct << endl;
	//_aTest.print_result();
}

// Test Case 1: 1000 nodes for Linked List based dijkarstra algorithm
// 1000
TEST(ListDij1000)
{
	cout << "Init Linked List based graph ...";
	print_current_time();
	cout << endl;

	int n = 200;
	list_graph _aTest(n);
	add_edges(&_aTest, n);

	cout << "Init Linked List based graph comleted.";
	print_current_time();
	cout << endl;

	int path_ct = _aTest.sp_dij_all_paths();
	CHECK(path_ct + 1);

	cout << "path number:" << path_ct << endl;
	//_aTest.print_result();
}

// Test Case 1: 1000 nodes for Linked List based Floyd algorithm
// 1000 
TEST(ListFloyd1000)
{
	int n = 2000;
	list_graph _aTest(n);
	add_edges(&_aTest, n);
	int path_ct = _aTest.sp_floyd();
	CHECK(path_ct + 1);

	cout << "path number:" << path_ct << endl;
	//_aTest.print_result();
}

// Test Case 1: 1000 nodes for One-d array based dijkstra algorithm
// 1000 
TEST(OneDArrayDij1000)
{
	cout << "Init One-dimension array based graph ...";
	print_current_time();
	cout << endl;

	int n = 1000;
	onedarray_graph _aTest(n);
	add_edges(&_aTest, n);

	cout << "Init One-dimension array based graph completed.";
	print_current_time();
	cout << endl;

	int path_ct = _aTest.sp_dij_all_paths();
	CHECK(path_ct + 1);

	cout << "path number:" << path_ct << endl;
	//_aTest.print_result();
}


// Test Case 1: 1000 nodes for Linked List based Floyd algorithm
// 1000 
TEST(OneDArrayFloyd1000)
{
	int n = 2000;
	onedarray_graph _aTest(n);
	add_edges(&_aTest, n);
	int path_ct = _aTest.sp_floyd();
	CHECK(path_ct + 1);

	cout << "path number:" << path_ct << endl;
	//_aTest.print_result();
}

int main(int argc, char* argv[])
{
	int rv = 0;
	rv = UnitTest::RunAllTests();

	getchar();

	return rv;
}