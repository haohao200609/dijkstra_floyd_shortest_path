/**
 * Dijkastra & Floyd Algorithm implementation
 * File: graph.hxx
 * Author: iJab(Zhan Caibao) zhancaibaoATgmail.com
 * Date: 2013/02/20
 *
 */

#ifndef __DIJ_FLOYD_SHORTES_PATH_GRAPH__
#define __DIJ_FLOYD_SHORTES_PATH_GRAPH__

#include <map>
#include <set>
#include <list>
#include <limits>
#include <ctime>

using namespace std;

const bool VERBOSE_DEBUG = false;
const double INF_WEIGHT = numeric_limits<double>::infinity();

// Helper function to print time
void print_current_time(); 


/**
 * Graph class
 *
 */
class graph
{
public:
	graph(int num){this->size = num; this->f_path = NULL;};

	virtual ~graph();

	virtual void add_edge(int source, int target, double w) = 0;
	virtual void remove_edge(int source, int target) = 0;
	int get_size() {return this->size;};

	virtual int sp_dijkstra(int source);
	virtual int sp_dij_all_paths();
	virtual int sp_floyd();
	virtual void get_nearest_vertex(const set<int> &visiting_v, int source, int current_v, double current_min_d) = 0;
	virtual void init_result_matrix(set<int> &visiting_v) = 0;
	
	virtual bool is_empty_graph() = 0;
	virtual void print_graph() = 0;
	void print_result();
	void print_path(int source, int target);

protected:
	int size;
	double **f_path;							// Store result of cost results
	int **prev_v;								// To track path
};

/**
 * Two dimensions based Graph implementation
 */
class twodarray_graph : public graph
{
public:
	twodarray_graph(int num);
	
	virtual ~twodarray_graph();

public:
	virtual void add_edge(int source, int target, double w);
	virtual void remove_edge(int source, int target);

	virtual void get_nearest_vertex(const set<int> &visiting_v, int source, int current_v, double current_min_d);
	virtual void init_result_matrix(set<int> &visiting_v);

	virtual bool is_empty_graph() {return this->adj_matrix == NULL;};
	virtual void print_graph();
private:
	double **adj_matrix;
};

/**
 * Linked list based Graph implementation
 */

/**
 * Node structure of linked list
 */
typedef struct node
{
	int v_id;
	double w;
	struct node *next;
} node;

class list_graph : public graph
{
	public:
	list_graph(int num);
	
	virtual ~list_graph();

public:
	virtual void add_edge(int source, int target, double w);
	virtual void remove_edge(int source, int target);

	virtual void get_nearest_vertex(const set<int> &visiting_v, int source, int current_v, double current_min_d);
	virtual void init_result_matrix(set<int> &visiting_v);

	virtual bool is_empty_graph() {return this->adj_matrix == NULL;};
	virtual void print_graph();

private:
	node **adj_matrix;
};

/**
 * One dimensions based Graph implementation
 */
class onedarray_graph : public graph
{
public:
	onedarray_graph(int num);
	
	virtual ~onedarray_graph();

public:
	virtual void add_edge(int source, int target, double w);
	virtual void remove_edge(int source, int target);

	virtual void get_nearest_vertex(const set<int> &visiting_v, int source, int current_v, double current_min_d);
	virtual void init_result_matrix(set<int> &visiting_v);

	virtual bool is_empty_graph() {return this->adj_matrix == NULL;};
	virtual void print_graph();

private:
	double *adj_matrix;
};

#endif