#ifndef BIGRAPH_HPP
#define BIGRAPH_HPP

#include <cstdint>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <memory>

#include "../utils/hash.hpp"
#include "../utils/fastio.hpp"
#include "../utils/vertexset.hpp"

struct BiGraph {
	VertexSet V[2];
	int n[2], m, maxDeg[2];
	std::vector<std::vector<int>> nbr[2];
	std::shared_ptr<std::vector<CuckooHash>> nbrMap[2];

	BiGraph(): n{0}, m(0), maxDeg{0} {
		nbrMap[0] = std::make_shared<std::vector<CuckooHash>>();
		nbrMap[1] = std::make_shared<std::vector<CuckooHash>>();
	}

	BiGraph(const std::string &filename): n{0}, m(0), maxDeg{0} {
		nbrMap[0] = std::make_shared<std::vector<CuckooHash>>();
		nbrMap[1] = std::make_shared<std::vector<CuckooHash>>();
		loadFromFile(filename);
	}

	int capacity(int side) {
		return nbr[side].size();
	}

	inline int degree(int side, int v) {
		return nbr[side][v].size();
	}

	inline int numVertices(int side) {
		return V[side].size();
	}

	inline int numEdges() {
		return m;
	}

	void clear() {
		m = 0;
		for (int s = 0; s <= 1; ++s) {
			n[s] = maxDeg[s] = 0;
			V[s].clear();
			nbr[s].clear();
			nbrMap[s] = std::make_shared<std::vector<CuckooHash>>();
		}

	}

	void resize(int side, int size) {

		nbr[side].resize(size);
		nbrMap[side]->resize(size);
		V[side].reserve(size);
	}

	void addEdge(int u, int v) {
		if (u >= capacity(0)) resize(0, u+1);
		if (v >= capacity(1)) resize(1, v+1);
		if ((*nbrMap[0])[u].find(v)) return;
		++m;
		n[0] = std::max(n[0], u+1);
		n[1] = std::max(n[1], v+1);
		V[0].push(u);
		V[1].push(v);
		nbr[0][u].push_back(v);
		nbr[1][v].push_back(u);
		(*nbrMap[0])[u].insert(v);
		(*nbrMap[1])[v].insert(u);
		maxDeg[0] = std::max(maxDeg[0], degree(0, u));
		maxDeg[1] = std::max(maxDeg[1], degree(1, v));
	}

	bool connect(int u, int v) {
		return (*nbrMap[0])[u].find(v);
	}

	bool connect(int uSide, int u, int v) {
		return (*nbrMap[uSide])[u].find(v);
	}

	void loadFromFile(const std::string &fname) {
		clear();
		FastIO fio(fname, "r");
		int numEdges = fio.getInt();
		int numVertsL = fio.getInt();
		int numVertsR = fio.getInt();
		resize(0, numVertsL);
		resize(1, numVertsR);
		for (int i = 0; i < numEdges; ++i) {
			int u = fio.getInt(), v = fio.getInt();
			addEdge(u, v);
		}
	}

	void addVertices(VertexSet S[2]) {
		for (int s = 0; s <= 1; ++s)
			for (int v : S[s])
				V[s].push(v);
	}

	template<typename... Args>
	void addVertices(VertexSet S[2], Args... Ss) {
		addVertices(S);
		addVertices(Ss...);
	}

	template <typename... Args>
	void subGraph(BiGraph &G, Args... Ss) {
		auto add_edge = [this](int u, int v) {
			++m;
			nbr[0][u].push_back(v);
			nbr[1][v].push_back(u);
		};

		V[0].clear(); V[1].clear();

		addVertices(Ss...);

		for (int s = 0; s <= 1; ++s) {
			n[s] = 0;
			for (int v : V[s]) {
				n[s] = std::max(n[s], v+1);
				if (v < capacity(s)) nbr[s][v].clear();
			}
			if (n[s] > capacity(s)) {
				nbr[s].resize(n[s]);
			}
		}

		for (int u : V[0]) {
			if (V[1].size() < G.degree(0, u)) {
				for (int v : V[1]) if (G.connect(u, v))
					add_edge(u, v);
			}
			else {
				for (int v : G.nbr[0][u]) if (V[1].inside(v))
					add_edge(u, v);
			}
		}
		for (int s = 0; s <= 1; ++s) {
			maxDeg[s] = 0;
			for (int u : V[s])
				maxDeg[s] = std::max(maxDeg[s], (int)nbr[s][u].size());
		}
	}

};

#endif
