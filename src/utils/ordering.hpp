#ifndef ORDERING_HPP
#define ORDERING_HPP

#pragma once

#include "linearheap.hpp"
#include "graph.hpp"
#include "bigraph.hpp"

#include <vector>

class Ordering {
	LinearHeap vHeap;
public:
	std::vector<int> order, ordered, value;
	int numOrdered, capacity;
	Ordering(int size=0):capacity(0) {
		resize(size);
		numOrdered = 0;
	}
	void push(int u, int val) {
		ordered[numOrdered] = u;
		order[u] = numOrdered++;
		value[u] = val;
	}
	void resize(int size) {
		capacity = size;
		order.resize(size);
		ordered.resize(size);
		value.resize(size);
	}
	void degeneracyOrdering(Graph &G) {
		if (capacity < G.n) resize(G.n);
		vHeap.build(G.V, [&](int v) { return (int)G.nbr[v].size(); });
		numOrdered = 0;
		while (!vHeap.empty()) {
			int u = vHeap.min();
			vHeap.popMin();
			for (int v : G.nbr[u]) {
				if (vHeap.value(v) <= vHeap.value(u)) continue;
				vHeap.dec(v);
			}
			push(u, vHeap.value(u));
		}
	}
	void degreeOrdering(Graph &G) {
		if (capacity < G.n) resize(G.n);
		vHeap.build(G.V, [&](int v) { return (int)G.nbr[v].size(); });
		numOrdered = 0;
		while (!vHeap.empty()) {
			int u = vHeap.min();
			vHeap.popMin();
			push(u, vHeap.value(u));
		}
	}
	void degeneracyOrdering(BiGraph &G) {
		Graph H;
		for (int u : G.V[0])
			for (int v : G.nbr[0][u])
				H.addEdge(u, v+G.n[0]);
		degeneracyOrdering(H);
	}
	void degreeOrdering(BiGraph &G) {
		Graph H;
		for (int u : G.V[0])
			for (int v : G.nbr[0][u])
				H.addEdge(u, v+G.n[0]);
		degreeOrdering(H);
	}
	static Ordering DegeneracyOrdering(Graph &G) {
		Ordering o;
		o.degeneracyOrdering(G);
		return o;
	}
	static Ordering DegreeOrdering(Graph &G) {
		Ordering o;
		o.degreeOrdering(G);
		return o;
	}

	std::vector<int>::iterator begin() {
		return ordered.begin();
	}
	std::vector<int>::iterator end() {
		return ordered.begin() + numOrdered;
	}
	std::vector<int>::reverse_iterator rbegin() {
		return ordered.rend() - numOrdered;
	}
	std::vector<int>::reverse_iterator rend() {
		return ordered.rend();
	}

};

#endif
