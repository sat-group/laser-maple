/*
 * proof_analyzer.cpp
 *
 *  Created on: May 23, 2017
 *      Author: ezulkosk
 */

using namespace std;
#include<iostream>
#include<fstream>
#include<vector>
#include <utility>
#include<map>
#include<queue>
#include <algorithm>
#include<set>


typedef map<int,vector<int> > Graph;


//================================= FILE PARSING ================================================

void process_graph_file_p1(const char* proof_file_name, Graph& graph, int& total_clauses){
	/*
	 * Parse the graph file to get the initial graph abstraction (ignore actual clauses for now).
	 *
	 * proof_file_name: the infile
	 * graph:           adjacency-list based DAG with edges pointing to clauses used to derive each clause
	 * total_clauses:   the number of clauses found by the solver, where the last clause should be bot
	 */
	fstream pf_file(proof_file_name, std::ios_base::in);
	int t;

	enum Mode {cid, lits, deps, min};
	Mode m = cid;

	int curr_cid = -1;
	total_clauses = -1;
	vector<int>* curr_adj_list;
	while(pf_file >> t){
		switch(m){
		case cid:
			curr_cid = t;
			if(curr_cid > total_clauses){
				total_clauses = curr_cid;
			}
			m = lits;
			break;
		case lits:
			if(t == 0){
				curr_adj_list = new vector<int>;
				m = deps;
			}
			break;
		case deps:
			if(t == 0){
				graph[curr_cid] = *curr_adj_list;
				m = min;
			}
			else{
				curr_adj_list->push_back(t);
			}
			break;
		case min:
			if(t == 0)
				m = cid;
			break;
		}
	}
}


void process_graph_file_p2(const char* proof_file_name,
		                   set<int>& useful_clauses,
						   map<int, vector<int> >& clause_literals,
						   int& proof_width){
	/*
	 * Parse the graph file to get the literals in each useful clause.
	 *
	 * proof_file_name: the infile
	 * useful_clauses:  the set of IDs of clauses used to derive unsat
	 * clause_literals: map from clause IDs to the literals in the clause
	 * proof_width:     the number of literals in the largest clause
	 */
	fstream pf_file(proof_file_name, std::ios_base::in);
	int t;

	enum Mode {cid, lits, deps, min};
	Mode m = cid;

	int curr_cid;
	int curr_width;
	vector<int>* curr_lits;
	bool curr_useful;
	while(pf_file >> t){
		switch(m){
		case cid:
			curr_cid = t;
			if(useful_clauses.find(curr_cid) != useful_clauses.end()){
				curr_lits = new vector<int>;
				curr_width = 0;
				curr_useful = true;
			}
			else
				curr_useful = false;
			m = lits;
			break;
		case lits:
			if(t == 0){
				if(curr_useful){
					clause_literals[curr_cid] = *curr_lits;
					if(curr_width > proof_width)
						proof_width = curr_width;
				}
				m = deps;
			}
			else if(curr_useful){
				curr_lits->push_back(t);
				curr_width++;
			}
			break;
		case deps:
			if(t == 0){
				m = min;
			}
			break;
		case min:
			if(t == 0)
				m = cid;
			break;
		}
	}
}


//============================= END FILE PARSING ================================================

//============================== CLAUSE PRUNING =================================================

void find_useful_clauses(Graph& graph,
		                 int& bot_clause,
						 set<int>& useful_clauses,
						 vector<int>& ordered_useful_clauses,
						 int& proof_length){
	/*
	 * Finds the set of clauses are actually used to derive unsat.
	 *
	 * graph:                  adjacency-list based DAG with edges pointing to clauses used to derive each clause
	 * bot_clause:             index of the final empty clause
	 * useful_clauses:         the set of IDs of clauses used to derive unsat
	 * ordered_useful_clauses: a vector abstraction of the useful clauses, ordered by index
	 * proof_length:           the number of useful clauses
	 */
	int c;
	queue<int> q;
	q.push(bot_clause);
	while(!q.empty()){
		c = q.front();
		q.pop();
		proof_length++;
		useful_clauses.insert(c);
		vector<int> adj = graph[c];
		for(int i : adj){
			if(useful_clauses.find(i) == useful_clauses.end())
				q.push(i);
		}
	}
	for(int i: useful_clauses)
		ordered_useful_clauses.push_back(i);
	std::sort(ordered_useful_clauses.begin(), ordered_useful_clauses.end());
}

void prune_useless_clauses(Graph& graph, set<int>& useful_clauses){
	/*
	 * Delete any useless clauses from the graph.
	 *
	 * graph:          adjacency-list based DAG with edges pointing to clauses used to derive each clause
	 * useful_clauses: the set of IDs of clauses used to derive unsat
	 */
	for(map<int, vector<int> >::iterator git = graph.begin(); git != graph.end();){
		int f = (*git).first;
		// delete any entries in the map that are not useful
		if(useful_clauses.find(f) == useful_clauses.end()){
			(*git).second.clear();
			git = graph.erase(git);
		}
		else{
			vector<int>& v = (*git).second;
			vector<int>::iterator vit = v.begin();
			 while (vit != v.end()) {
				 if(useful_clauses.find(*vit) == useful_clauses.end())
					 vit = v.erase(vit);
				 else
					 ++vit;
			 }
			 ++git;
		}
	}
}

//========================== END CLAUSE PRUNING =================================================

//============================ PROOF ANALYSES ===================================================

void ordered_proof_space(Graph& graph, vector<int>& ordered_useful_clauses, int& proof_ordered_space){
	/*
	 * Computes the space of the proof according to the order derived by the solver.
	 * For each node at index t, compute the number of nodes at index <= t that are used at index >= t.
	 * The space of the proof is the max over all nodes.
	 *
	 * graph:               adjacency-list based DAG with edges pointing to clauses used to derive each clause
	 *	 	 	 	        assumed to be pruned at this point
	 * bot_clause:          index of the final empty clause
	 * proof_ordered_space: the space of the proof
	 */
	set<int> live_clauses;
	proof_ordered_space = 0;
	for (auto it = ordered_useful_clauses.rbegin(); it != ordered_useful_clauses.rend(); ++it)
	{
		live_clauses.insert(*it);
		vector<int> deps = graph[*it];
		for(int i: deps)
			live_clauses.insert(i);

		if(live_clauses.size() > proof_ordered_space){
			proof_ordered_space = live_clauses.size();
		}
		live_clauses.erase(*it);
	}
}


//========================== END PROOF ANALYSES =================================================


int main(int argc, char** argv){

	const char* proof_file_name = argv[1];
	Graph graph;
	set<int> useful_clauses; //the set of IDs of clauses used to derive unsat
	vector<int> ordered_useful_clauses;
	map<int, vector<int> > clause_literals;
	int total_clauses; // should index the empty clause
	int proof_length = 0;
	int proof_width = -1;
	int proof_ordered_space = -1; // space of proof w.r.t. the order the solver found the clauses
	int unused_clauses = -1;

	process_graph_file_p1(proof_file_name, graph, total_clauses);
	find_useful_clauses(graph, total_clauses, useful_clauses, ordered_useful_clauses, proof_length);
	unused_clauses = total_clauses - proof_length;
	prune_useless_clauses(graph, useful_clauses);
	process_graph_file_p2(proof_file_name, useful_clauses, clause_literals, proof_width);
	ordered_proof_space(graph, ordered_useful_clauses, proof_ordered_space);

	//cout <<graph.size()<<" "<<useful_clauses.size()<<endl;
	//for(auto elem: graph)
	//	cout<<elem.first<<" "<<elem.second.size()<<endl;

	printf("%-15s: %d\n", "Length", proof_length);
	printf("%-15s: %d\n", "Unused", unused_clauses);
	printf("%-15s: %d\n", "Width", proof_width);
	printf("%-15s: %d\n", "OrderedSpace", proof_ordered_space);
}


