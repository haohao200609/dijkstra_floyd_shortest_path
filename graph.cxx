/**
 * Dijkastra & Floyd Algorithm implementation
 * File: graph.cxx
 * Author: iJab(Zhan Caibao) zhancaibaoATgmail.com
 * Date: 2013/02/20
 *
 */

#include <iostream>
#include <stack>
#include <list>
#include <set>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <limits>

#include "graph.hxx"

using namespace std;

// Helper function to print time
void print_current_time()
{
		time_t rawtime;
		struct tm * timeinfo;

		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		printf ( "The current date/time is: %s", asctime (timeinfo) );
};

///////////////////////////////////////////////////////////////
/**
 * Graph base class implementation
 *
 */
graph::~graph()
{
	// Release resource allocated
	if(this->f_path != NULL)
	{
		int v_size = this->get_size();
		for(int i = 0; i < v_size; ++i)
		{
			if(this->f_path[i] != NULL)
				delete []this->f_path[i];
		}
		delete []this->f_path;
	}

	if(this->prev_v != NULL)
	{
		int v_size = this->get_size();
		for(int i = 0; i < v_size; ++i)
		{
			if(this->prev_v[i] != NULL)
				delete []this->prev_v[i];
		}
		delete []this->prev_v;
	}
}

/**
 * Calculate shortest path using Dijkastra
 * @param int source : source vertex id
 * @return int : number of path from source to other nodes
 */
int graph::sp_dijkstra(int source)
{
	//cout << "Dijkastra Try Calculating shortest path from vertex " << source << endl;
	//print_current_time();
	//cout << endl;

	int _v_size = this->get_size();
	int path_ct = 0;

	if(source < 0 || source >= _v_size) return 0;

	set<int> visiting_v;		// vertices which are waiting to visited
	this->init_result_matrix(visiting_v);	// Init path cost results and path track matrix
		
	
	// Size of vertices
	int n_visiting = visiting_v.size();

	// Loop from source vertex
	int current_v = source;
	
	while(n_visiting > 0)
    {
		// Find minimum distance in handling data
		double current_min_d = INF_WEIGHT;
		set<int>::iterator cur_v_it = visiting_v.end();
		for(set<int>::iterator v_it = visiting_v.begin(); v_it != visiting_v.end(); ++v_it)
		{
			if(this->f_path[source][*v_it] < current_min_d)
			{
				current_min_d = this->f_path[source][*v_it];
				current_v = *v_it;
				cur_v_it = v_it;
			}
		}
		
		if(current_min_d == INF_WEIGHT)			// No path from source to current_v
			break;

		// Remove current_v from visiting vertex list
		if(cur_v_it != visiting_v.end())
			visiting_v.erase(cur_v_it);
		else									// No path from source to current vertex
			break;
		n_visiting = visiting_v.size();
	

		// Get the nearest neighbourhood of current vertex: current_v
		this->get_nearest_vertex(visiting_v, source, current_v, current_min_d);
		
		// Increase path count
		path_ct++;
    }

	//cout << "Dijkastra Complete Calculating shortest path from vertex " << source << endl;
	//print_current_time();
	//cout << endl;

	return path_ct;
}

/**
 * Calculate all shortest path using Dijkastra
 * @param void
 * @return int : all paths
 */
int graph::sp_dij_all_paths()
{
	int _v_size = this->get_size();
	int all_path_ct = 0;

	for(int i = 0; i < _v_size; ++i)
	{
		all_path_ct += this->sp_dijkstra(i);
	}

	return all_path_ct;
}

/**
 * Shortest path through Floyd algorithm
 * @param void
 * @return int : all paths
 */
int graph::sp_floyd()
{
	if(this->is_empty_graph()) return 0;

	cout << "Floyd try Calculating all shortest path ... " << endl;
	print_current_time();
	cout << endl;

	int _v_size = this->get_size();
	int path_ct = 0;
	
	// Init all other path as unknown cost
	set<int> visiting_v;	// Dummy for Floyd algorithm
	this->init_result_matrix(visiting_v);

	if(this->f_path == NULL) return 0;

	// Three nested loop to calculate all possible shortest path
	for(int k = 0; k < _v_size; ++k)
	{
		for(int i = 0; i < _v_size; ++i)
		{
			for(int j = 0; j < _v_size; ++j)
			{
				// Find least cost from i to j
				if((this->f_path[i][k] + this->f_path[k][j]) < this->f_path[i][j])
				{
					this->f_path[i][j] = this->f_path[i][k] + this->f_path[k][j];
					this->prev_v[i][j] = this->prev_v[k][j];
				}
			}
			path_ct++;
		}
	}

	cout << "Floyd Completed Calculating all shortest path." << endl;
	print_current_time();
	cout << endl;

	return path_ct;
}

/**
 * Print the shortest path cost result
 */
void graph::print_result()
{
	if(this->f_path != NULL)
	{
		int v_size = this->get_size();
		for(int i = 0; i < v_size; ++i)
		{
			for(int j = 0; j < v_size; ++j)
			{
				if(this->f_path[i][j] == INF_WEIGHT) continue;

				cout << "Shortest path from " << i << " to " << j;
				cout << " Cost(" << this->f_path[i][j] << "):";
				this->print_path(i, j);
				cout << endl;
			}
			cout << endl;
		}
	}
}

/**
 * Print the shortest path from source vertex to target vertex
 * @param int source : source vertex
 * @param int target : target vertex
 */
void graph::print_path(int source, int target)
{
	// Use stack to track the path
	stack<int> stack_p;
	int _p_target = target;
	stack_p.push(target);

	if(this->prev_v != NULL)
	{
		while(true)
		{
			int _p_vertex = this->prev_v[source][_p_target];
			if(_p_vertex != -1)
			{
				stack_p.push(_p_vertex);

				if(_p_vertex == source)		// Get the path
					break;
			}
			else
				break;

			// Set target to new previous vertex
			_p_target = _p_vertex;
		}
	}

	// Try to print it
	int node_size = stack_p.size();
	while(node_size > 0)
	{
		cout << stack_p.top();
		stack_p.pop();

		node_size = stack_p.size();
		if(node_size > 0)
			cout << ">";
	}
}

///////////////////////////////////////////////////////////////////////////////////
/**
 * Implementation of 2-d array based graph
 *
 */

/**
 * Constructor of 2-d array based graph
 * @param num int: size of vertices in the graph
 * @return void
 */
twodarray_graph::twodarray_graph(int num) : graph(num)
{
	if(num < 1)
	{
		this->adj_matrix = NULL;
	}
	else
	{
		// Init and new memory for graph
		this->adj_matrix = new double*[num];

		for(int i = 0; i < num; ++i)
		{
			this->adj_matrix[i] = new double[num];
			// Set init value to INF_WEIGHT which means no path between vertices
			for( int j = 0; j < num; ++j)
			{
				this->adj_matrix[i][j] = INF_WEIGHT;
			}
		}
	}
}

/**
 * Desstructor of 2-d array based graph
 * @param num int: size of vertices in the graph
 * @return void
 */
twodarray_graph::~twodarray_graph()
{
	// Release resource
	if(this->adj_matrix != NULL)
	{
		int _size = this->get_size();
		for(int i = 0; i < _size; ++i)
		{
			delete []this->adj_matrix[i];
			this->adj_matrix[i] = NULL;
		}
		delete []this->adj_matrix;
		this->adj_matrix = NULL;
	}
}

/**
 * Add connected edge to graph
 * @param int source : source vertex id
 * @param int target : target vertex id
 * @param double w : weight between these two vertices
 * @return void
 */
void twodarray_graph::add_edge(int source, int target, double w)
{
	// Not allocated memory for graph
	if(this->adj_matrix == NULL)
		return;

	// Negative source or target
	if(source < 0 || target < 0)
		return;

	// source or target is greater than number of vertices predefined
	int _size = this->get_size();
	if(source >= _size || target >= _size)
		return;

	if(source == target)
		this->adj_matrix[source][target] = 0;
	else
		this->adj_matrix[source][target] = w < 1 ? INF_WEIGHT : w;
}

/**
 * Remove connected edge between vertices from graph
 * @param int source : source vertex id
 * @param int target : target vertex id
 * @return void
 */
void twodarray_graph::remove_edge(int source, int target)
{
	// Make the weight between vertices to be -1 to disconnect
	this->add_edge(source, target, INF_WEIGHT);
}

/**
 * Get current node's nearest neighbour
 * @param const set<int> &visiting_v : Still not visited vertices
 * @param int source : original source vertex
 * @param int current_v : Current vertex id
 * @param double current_min_d : Current minimum cost
 * @return void
 *
 */
void twodarray_graph::get_nearest_vertex(const set<int> &visiting_v, int source, int current_v, double current_min_d)
{
	// Scan all neighbourhood of current vertex: current_v
	double alt_d = INF_WEIGHT;
	
	for(set<int>::const_iterator v_it = visiting_v.begin(); v_it != visiting_v.end(); ++v_it)
	{
		int j = *v_it;
		if(j == current_v) continue;

		// If not connected to current vertex
		if(this->adj_matrix[current_v][j] == INF_WEIGHT) 
			continue;

		// If cost of connected to this node is smaller
		alt_d = current_min_d + this->adj_matrix[current_v][j];

		if(alt_d < this->f_path[source][j])
		{
			// update result path cost and track matrix
			this->f_path[source][j] = alt_d;
			this->prev_v[source][j] = current_v;
		}
	}
}


/**
 * Init Floyd path matrix
 * 
 */
void twodarray_graph::init_result_matrix(set<int> &visiting_v)
{
	if(this->adj_matrix == NULL) return;

	int _v_size = this->get_size();

	this->f_path = new double*[_v_size];
	this->prev_v = new int*[_v_size];
	
	// Init all other path as unknown cost
	for(int i = 0; i < _v_size; ++i)
	{	
		// Add all vertices into undetermined list
		visiting_v.insert(i);

		// Init path and path track matrix
		this->f_path[i] = new double[_v_size];
		this->prev_v[i] = new int[_v_size];
		for(int j = 0; j < _v_size; ++j)
		{
			// Init cost and path
			double cost = INF_WEIGHT;
			if(i == j)
			{
				cost = 0;
				this->prev_v[i][j] = -1;
			}
			else
			{
				cost = this->adj_matrix[i][j];
				if(cost == INF_WEIGHT)
					this->prev_v[i][j] = -1;
				else
					this->prev_v[i][j] = i;
			}
			this->f_path[i][j] = cost;
		}
	}
}

/**
 * Print graph
 * @param void
 * @return void
 */
void twodarray_graph::print_graph()
{
	int _v_size = this->get_size();
	
	for(int i = 0; i < _v_size; ++i)
	{
		for(int j = 0; j < _v_size; ++j)
		{
			cout << this->adj_matrix[i][j] << "\t\t";
		}
		cout << endl;
	}
}

///////////////////////////////////////////////////////////////////////////////////
/**
 * Implementation of Linked List based graph
 *
 */

/**
 * Constructor of Linked List based graph
 * @param num int: size of vertices in the graph
 * @return void
 */
list_graph::list_graph(int num) : graph(num)
{
	if(num < 1)
	{
		this->adj_matrix = NULL;
	}
	else
	{
		// Init and new memory for graph
		this->adj_matrix = new node*[num];

		for(int i = 0; i < num; ++i)
		{
			this->adj_matrix[i] = NULL;
		}
	}
}

/**
 * Desstructor of Linked List based graph
 * @param num int: size of vertices in the graph
 * @return void
 */
list_graph::~list_graph()
{
	// Release resource
	if(this->adj_matrix != NULL)
	{
		int _size = this->get_size();
		for(int i = 0; i < _size; ++i)
		{
			node *tmp;
			node *cur = this->adj_matrix[i];
			while(cur != NULL)
			{
				tmp = cur->next;
				delete cur;
				cur = NULL;
				cur = tmp;
			}
		}
		delete []this->adj_matrix;
		this->adj_matrix = NULL;
	}
}

/**
 * Add connected edge to graph
 * @param int source : source vertex id
 * @param int target : target vertex id
 * @param double w : weight between these two vertices
 * @return void
 */
void list_graph::add_edge(int source, int target, double w)
{
	// Not allocated memory for graph
	if(this->adj_matrix == NULL)
		return;

	// Self assignment
	if(source == target)
		return;

	// Negative source or target
	if(source < 0 || target < 0)
		return;

	// Negative or zero w means no connection, just ignore it
	if(w <= 0)
		return;

	// source or target is greater than number of vertices predefined
	int _size = this->get_size();
	if(source >= _size || target >= _size)
		return;

	// Try to add new target to source's list
	node *new_node = new node;
	new_node->v_id = target;
	new_node->w = w;
	new_node->next = NULL;

	if(this->adj_matrix[source] == NULL)	// List is empty
		this->adj_matrix[source] = new_node;
	else
	{
		node *head = this->adj_matrix[source];
		while(head->next != NULL)	// Move to the list's tail
			head = head->next;

		head->next = new_node;
	}
}

/**
 * Remove connected edge between vertices from graph
 * @param int source : source vertex id
 * @param int target : target vertex id
 * @return void
 */
void list_graph::remove_edge(int source, int target)
{
	int _v_size = this->get_size();

	// check out of index range or not
	if(source < 0 || target < 0)
		return;

	if(source >= _v_size || target >= _v_size)
		return;

	// Check NULL pointer
	if(this->adj_matrix == NULL || this->adj_matrix[source] == NULL)
		return;

	// Remove from soure's edge list
	node *head = this->adj_matrix[source];
	node *prev = NULL;
	while(head != NULL)
	{
		if(head->v_id == target)
		{
			if(prev != NULL)
			{
				// Set previous' next to head->next
				prev->next = head->next;
			}
			delete head;	// Free memory
		}

		// Set prev
		prev = head;
		head = head->next;
	}
}

/**
 * Get current node's nearest neighbour
 * @param const set<int> &visiting_v : Still not visited vertices
 * @param int source : orginal source vertex
 * @param int current_v : Current vertex id
 * @param double current_min_d : Current minimum cost
 * @return void
 *
 */
void list_graph::get_nearest_vertex(const set<int> &visiting_v, int source, int current_v, double current_min_d)
{
	// Scan all neighbourhood of current vertex: current_v
	double alt_d = INF_WEIGHT;
	int alt_current_v = current_v;
	double alt_cur_min_d = INF_WEIGHT;

	node *head = this->adj_matrix[current_v];
	
	while(head != NULL)
	{
		int j = head->v_id;

		// If it's not in waiting visiting vertices list, it should have been visited
		if(visiting_v.find(j) != visiting_v.end() && j != current_v)
		{
			// If cost of connected to this node is smaller
			alt_d = current_min_d + head->w;

			if(alt_d < this->f_path[source][j])
			{
				// Update path cost results and track matrix
				this->f_path[source][j] = alt_d;
				this->prev_v[source][j] = current_v;
			}
		}
		head = head->next;
	}
}


/**
 * Init Floyd path matrix
 * 
 */
void list_graph::init_result_matrix(set<int> &visiting_v)
{
	if(this->adj_matrix == NULL) return;

	int _v_size = this->get_size();

	this->f_path = new double*[_v_size];
	this->prev_v = new int*[_v_size];
	
	// Init all other path as unknown cost
	for(int i = 0; i < _v_size; ++i)
	{	
		// Add all vertices into undetermined list
		visiting_v.insert(i);

		// Init path and path track matrix
		this->f_path[i] = new double[_v_size];
		this->prev_v[i] = new int[_v_size];
		for(int j = 0; j < _v_size; ++j)
		{
			// Init cost and path
			double cost = INF_WEIGHT;
			if(i == j)
			{
				cost = 0;
				this->prev_v[i][j] = -1;
			}
			else
			{
				node *head = this->adj_matrix[i];
				while(head != NULL)
				{
					if(head->v_id == j)
					{
						cost = head->w;
						if(cost == INF_WEIGHT)
							this->prev_v[i][j] = -1;
						else
							this->prev_v[i][j] = i;
						break;
					}
					head = head->next;
				}
			}

			this->f_path[i][j] = cost;
		}
	}
}

/**
 * Print graph
 * @param void
 * @return void
 */
void list_graph::print_graph()
{
	int _v_size = this->get_size();
	
	for(int i = 0; i < _v_size; ++i)
	{
		node *head = this->adj_matrix[i];
		while(head != NULL)
		{
			cout << head->w << "\t\t";
		}
		cout << endl;
	}
}


///////////////////////////////////////////////////////////////////////////////////
/**
 * Implementation of One-d array based graph
 *
 */

/**
 * Constructor of One-d array based graph
 * @param num int: size of vertices in the graph
 * @return void
 */
onedarray_graph::onedarray_graph(int num) : graph(num)
{
	if(num < 1)
	{
		this->adj_matrix = NULL;
	}
	else
	{
		int _v_size = (num*num-num)/2;
		// Init and new memory for graph
		this->adj_matrix = new double[_v_size];

		for(int i = 0; i < _v_size; ++i)
		{
			this->adj_matrix[i] = INF_WEIGHT;
		}
	}
}

/**
 * Desstructor of One-d array based graph
 * @param num int: size of vertices in the graph
 * @return void
 */
onedarray_graph::~onedarray_graph()
{
	// Release resource
	if(this->adj_matrix != NULL)
	{
		delete []this->adj_matrix;
		this->adj_matrix = NULL;
	}
}

/**
 * Add connected edge to graph
 * @param int source : source vertex id
 * @param int target : target vertex id
 * @param double w : weight between these two vertices
 * @return void
 */
void onedarray_graph::add_edge(int source, int target, double w)
{
	// Not allocated memory for graph
	if(this->adj_matrix == NULL)
		return;

	// Self assignment
	if(source == target)
		return;

	// Negative source or target
	if(source < 0 || target < 0)
		return;

	// source or target is greater than number of vertices predefined
	int _size = this->get_size();
	if(source >= _size || target >= _size)
		return;

	// source should less than target for we only store the left bottom of the matrix
	if(source < target)
	{
		int tmp = source;
		source = target;
		target = tmp;
	}

	int _v_size = (_size*_size - _size)/2;
	int _v_ix = source*(source - 1)/2 + target;
	if(_v_ix >= _v_size) return;

	this->adj_matrix[_v_ix] = w < 1 ? INF_WEIGHT : w;
}

/**
 * Remove connected edge between vertices from graph
 * @param int source : source vertex id
 * @param int target : target vertex id
 * @return void
 */
void onedarray_graph::remove_edge(int source, int target)
{
	// Just set the weight to negative to make them unconnected
	this->add_edge(source, target, -1);
}

/**
 * Get current node's nearest neighbour
 * @param const set<int> &visiting_v : Still not visited vertices
 * @param int source : orginal source vertex
 * @param int current_v : Current vertex id
 * @param double current_min_d : Current minimum cost
 * @return pair<int, double> : Nearest neighbour and the new minimum cost to this neighbour node
 *
 */
void onedarray_graph::get_nearest_vertex(const set<int> &visiting_v, int source, int current_v, double current_min_d)
{
	// Scan all neighbourhood of current vertex: current_v
	double alt_d = INF_WEIGHT;

	for(set<int>::const_iterator v_it = visiting_v.begin(); v_it != visiting_v.end(); ++v_it)
	{
		int j = *v_it;
		if(j == current_v) continue;

		// Get index in one-dimension array
		int ix = current_v > j ? current_v*(current_v-1)/2 + j : j*(j-1)/2 + current_v;

		// If not connected to current vertex
		if(this->adj_matrix[ix] == INF_WEIGHT) 
			continue;

		// If cost of connected to this node is smaller
		alt_d = current_min_d + this->adj_matrix[ix];

		if(alt_d < this->f_path[source][j])
		{
			// update shortest path result and track matrix
			this->f_path[source][j] = alt_d;
			this->prev_v[source][j] = current_v;
		}
	}
}


/**
 * Init Floyd path matrix
 * 
 */
void onedarray_graph::init_result_matrix(set<int> &visiting_v)
{
	if(this->adj_matrix == NULL) return;

	int _v_size = this->get_size();

	this->f_path = new double*[_v_size];
	this->prev_v = new int*[_v_size];
	
	// Init all other path as unknown cost
	for(int i = 0; i < _v_size; ++i)
	{	
		// Add all vertices into undetermined list
		visiting_v.insert(i);

		// Init path and path track matrix
		this->f_path[i] = new double[_v_size];
		this->prev_v[i] = new int[_v_size];
		for(int j = 0; j < _v_size; ++j)
		{
			// Init cost and path
			double cost = INF_WEIGHT;
			if(i == j)
			{
				cost = 0;
				this->prev_v[i][j] = -1;
			}
			else
			{
				// Get index in one-dimension array
				int ix = i > j ? i*(i-1)/2 + j : j*(j-1)/2 + i;
				cost = this->adj_matrix[ix];
				if(cost == INF_WEIGHT)
					this->prev_v[i][j] = -1;
				else
					this->prev_v[i][j] = i;
			}

			this->f_path[i][j] = cost;
		}
	}
}

/**
 * Print graph
 * @param void
 * @return void
 */
void onedarray_graph::print_graph()
{
	int _v_size = this->get_size();
	
	for(int i = 0; i < _v_size; ++i)
	{
		for(int j = 0; j < _v_size; ++j)
		{
			// Get index in one-dimension array
			int ix = i > j ? i*(i-1)/2 + j : j*(j-1)/2 + i;
			cout << this->adj_matrix[ix] << "\t\t";
		}
		cout << endl;
	}
}