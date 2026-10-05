#pragma once
#ifndef MDB_H
#define MDB_H

#include "../utils/vertexset.hpp"
#include "../utils/hash.hpp"
#include "../utils/bigraph.hpp"
#include "result.h"

#define numEdgesSs (Ss[0].size() * Ss[1].size() - numNnbSs)
#define numEdgesS (S[0].size() * S[1].size() - numNnbS)
#define nnbS(s, v) (S[s ^ 1].size() - degS[s][v])
#define nnbC(s, v) (C[s ^ 1].size() - degC[s][v])
#define degSub(s, v) (degS[s][v] + degC[s][v])
#define nnbSub(s, v) (nnbS(s, v) + nnbC(s, v))

#define COMM_ROUNDS 1

#define FLAG_DEBUG (1 << 0)
#define FLAG_UB (1 << 1)
#define FLAG_CORE (1 << 2)
#define FLAG_CN (1 << 3)
#define FLAG_1NN (1 << 4)
#define FLAG_QUEUE (1 << 5)
#define FLAG_ORDER (1 << 6)
#define FLAG_PB (1 << 7)
#define FLAG_HEU (1 << 8)
#define FLAG_BR (1 << 9)

class MDB
{
public:
	void findMDB(const std::string &dataPath, int q[2], int k, int flags = 7, int numThreads = 1);
	SolverResult getResult() const;
	VertexSet Ss[2];
	int numNnbSs, k, lb[2], flags, numThreads;
	thread_local static VertexSet S[2], C[2], X[2];
	thread_local static int numNnbS;
	thread_local static std::vector<int> degS[2], degC[2];
	thread_local static BiGraph G;
	int branchTime, reductionTime, numBranches, numUbPruned, numPivoting, numBipartite;

	struct BakPos
	{
		int p[2][2];
		inline void backup(const VertexSet *V)
		{
			p[0][0] = V[0].frontPos();
			p[0][1] = V[1].frontPos();
			p[1][0] = V[0].backPos();
			p[1][1] = V[1].backPos();
		}
	};

	virtual void branch(int dep) = 0;

	BiGraph core(BiGraph &G, int a, int b);
	BiGraph comm(BiGraph &G, int q[2], int k);

	void heuristic(BiGraph &G);

	bool upperbound();

	void russianDoll();
	void searchAll();

	BakPos backup(const VertexSet V[2]);
	BakPos update(int uSide, int u);
	BakPos minus(int uSide, int u);
	void restore(BakPos &pos);

	inline void addS(int s, int u)
	{
		S[s].push(u);
		for (int v : G.nbr[s][u])
			++degS[s ^ 1][v];
	}
	inline void addC(int s, int u)
	{
		C[s].push(u);
		for (int v : G.nbr[s][u])
			++degC[s ^ 1][v];
	}
	inline void subS(int s, int u)
	{
		S[s].pop(u);
		for (int v : G.nbr[s][u])
			--degS[s ^ 1][v];
	}
	inline void subC(int s, int u)
	{
		C[s].pop(u);
		for (int v : G.nbr[s][u])
			--degC[s ^ 1][v];
	}

	inline void moveC2S(int s, int u)
	{
		numNnbS += nnbS(s, u);
		C[s].popBack(u);
		S[s].push(u);
		for (int v : G.nbr[s][u])
		{
			--degC[s ^ 1][v];
			++degS[s ^ 1][v];
		}
	}

	inline void moveS2C(int s, int u)
	{
		S[s].pop(u);
		C[s].push(u);
		for (int v : G.nbr[s][u])
		{
			--degS[s ^ 1][v];
			++degC[s ^ 1][v];
		}
		numNnbS -= nnbS(s, u);
	}
};

#endif
