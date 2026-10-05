#include "mdb.h"
#include "../utils/vertexset.hpp"
#include "../utils/bigraph.hpp"
#include "../utils/ordering.hpp"
#include "../utils/log.hpp"
#include "../utils/traversal.hpp"
#include <chrono>
#include <queue>
#include <algorithm>
#include <cassert>
#include <utility>

using namespace logging;
namespace jsbcrit
{
	thread_local bool witnessValid = false;
	thread_local int witnessI = 0;
	thread_local int witnessJ = 0;
	thread_local unsigned char activeMask = 0;
	thread_local int attackH[2] = {-1, -1};
	thread_local long long attackDrop[2] = {0, 0};
}

namespace
{
	const long long SHORTAGE_EXACT_PAIR_CAP = 2000000LL;

	int maxPrefixCount(const std::vector<long long> &pref, long long budget)
	{
		int lo = 0, hi = (int)pref.size() - 1;
		while (lo < hi)
		{
			int mid = (lo + hi + 1) >> 1;
			if (pref[mid] <= budget)
				lo = mid;
			else
				hi = mid - 1;
		}
		return lo;
	}

	void fillPrefixFromDSCount(const std::vector<int> &cnt, std::vector<long long> &pref)
	{
		std::size_t total = 0;
		for (int x : cnt)
			total += (std::size_t)x;

		pref.clear();
		if (pref.capacity() < total + 1)
			pref.reserve(total + 1);

		long long sum = 0;
		pref.push_back(0);
		for (int ds = 0; ds < (int)cnt.size(); ++ds)
		{
			for (int t = 0; t < cnt[ds]; ++t)
			{
				sum += ds;
				pref.push_back(sum);
			}
		}
	}

	std::vector<long long> buildTopMinRowFromDegreeHist(
		const std::vector<int> &degHist, int maxCount, int cap)
	{
		maxCount = std::max(0, maxCount);
		int maxDeg = (int)degHist.size() - 1;
		cap = std::max(0, std::min(cap, maxDeg));

		std::vector<long long> pref(maxCount + 1, 0);
		if (maxCount == 0 || cap == 0)
			return pref;

		long long high = 0;
		for (int d = cap; d <= maxDeg; ++d)
			high += degHist[d];

		int curCount = 0;
		long long curSum = 0;

		int take = (int)std::min<long long>(high, maxCount);
		for (int t = 1; t <= take; ++t)
			pref[t] = (long long)t * cap;
		curCount = take;
		curSum = (long long)take * cap;

		for (int d = cap - 1; d >= 0 && curCount < maxCount; --d)
		{
			int cnt = std::min<int>(degHist[d], maxCount - curCount);
			for (int t = 1; t <= cnt; ++t)
				pref[curCount + t] = curSum + (long long)t * d;
			curSum += (long long)cnt * d;
			curCount += cnt;
		}

		return pref;
	}

	std::vector<std::vector<long long>> buildJointSuffix(const std::vector<std::vector<int>> &hist,
														 int maxDS, int maxCap)
	{
		std::vector<std::vector<long long>> suffix(maxDS + 1, std::vector<long long>(maxCap + 2, 0));
		for (int ds = 0; ds <= maxDS; ++ds)
			for (int dc = maxCap; dc >= 0; --dc)
				suffix[ds][dc] = suffix[ds][dc + 1] + hist[ds][dc];
		return suffix;
	}

	std::vector<long long> buildJointRowFromBucket(
		const std::vector<std::vector<int>> &hist,
		const std::vector<std::vector<long long>> &suffix,
		int maxDS, int maxCount, int maxCap, int cap)
	{
		maxCount = std::max(0, maxCount);
		maxCap = std::max(0, maxCap);
		cap = std::max(0, std::min(cap, maxCap));

		std::vector<long long> pref(maxCount + 1, 0);
		if (maxCount == 0)
			return pref;

		const int offset = maxDS;
		std::vector<long long> valueCnt(maxCap + maxDS + 1, 0);

		for (int ds = 0; ds <= maxDS; ++ds)
		{

			for (int dc = 0; dc < cap; ++dc)
				if (hist[ds][dc])
					valueCnt[dc - ds + offset] += hist[ds][dc];

			long long high = suffix[ds][cap];
			if (high)
				valueCnt[cap - ds + offset] += high;
		}

		int curCount = 0;
		long long curSum = 0;
		for (int value = cap; value >= -maxDS && curCount < maxCount; --value)
		{
			long long available = valueCnt[value + offset];
			if (available <= 0)
				continue;
			int take = (int)std::min<long long>(available, maxCount - curCount);
			for (int t = 1; t <= take; ++t)
				pref[curCount + t] = curSum + (long long)t * value;
			curSum += (long long)take * value;
			curCount += take;
		}

		return pref;
	}

	long long sumTopMinFromDegreeHist(const std::vector<int> &degHist, int count, int cap)
	{
		if (count <= 0 || cap <= 0)
			return 0;
		int maxCap = (int)degHist.size() - 1;
		cap = std::min(cap, maxCap);

		int curCount = 0;
		long long curSum = 0;
		long long high = 0;
		for (int d = cap; d <= maxCap; ++d)
			high += degHist[d];

		long long take = std::min<long long>(high, count - curCount);
		curSum += take * cap;
		curCount += (int)take;

		for (int d = cap - 1; d >= 0 && curCount < count; --d)
		{
			take = std::min<long long>(degHist[d], count - curCount);
			curSum += take * d;
			curCount += (int)take;
		}
		return curSum;
	}
}

thread_local VertexSet MDB::S[2], MDB::C[2], MDB::X[2];
thread_local std::vector<int> MDB::degS[2], MDB::degC[2];
thread_local BiGraph MDB::G;
thread_local int MDB::numNnbS;

SolverResult MDB::getResult() const
{
	return SolverResult{numEdgesSs, numBranches, numUbPruned};
}

void MDB::findMDB(const std::string &dataPath, int q[2], int k, int flags, int numThreads)
{
	BiGraph G(dataPath);
	logBiGraph(G);

	MDB::flags = flags;
	MDB::numThreads = numThreads;
	MDB::k = k;
	MDB::lb[0] = q[0];
	MDB::lb[1] = q[1];

	numNnbS = numNnbSs =
		numBranches = numPivoting = numBipartite =
			numUbPruned = branchTime = reductionTime = 0;

	for (int s = 0; s <= 1; ++s)
	{
		Ss[s].reserve(G.n[s]);
#pragma omp parallel num_threads(numThreads)
		{
			S[s].reserve(G.n[s]);
			C[s].reserve(G.n[s]);

			degS[s].resize(G.n[s]);
			degC[s].resize(G.n[s]);
		}

	}

	auto searchStartTime = std::chrono::steady_clock::now();

	log("Start reductions...");
	BiGraph g = flags & FLAG_CORE ? core(G, q[1] - k, q[0] - k) : G;
	logBiGraph(g);
	if (flags & FLAG_CN)
		g = comm(g, q, k);
	logBiGraph(g);

	if (flags & FLAG_HEU)
	{
		log("Start heuristic...");
		auto heuristicStartTime = std::chrono::steady_clock::now();
		heuristic(g);
		int heuristicTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - heuristicStartTime).count();
		log("Heuristic done! Time: %d ms", heuristicTime);
		log("G[S*] info: |U|=%d, |V|=%d, |E|=%d", Ss[0].size(), Ss[1].size(), numEdgesSs);
		logSetInfo(Ss[0], "U*");
		logSetInfo(Ss[1], "V*");
	}

	log("Start searching...");
	BiGraph &Sub = this->G;

	if (flags & FLAG_PB)
	{
		int s = (int)(Ss[0].size() < Ss[1].size());
		lb[s] = g.maxDeg[s ^ 1] + k;
		while (lb[s] > q[s])
		{
			MDB::lb[s ^ 1] = std::max(q[s ^ 1], std::max(numEdgesSs, q[0] * q[1] - k) / lb[s]);
			MDB::lb[s] = std::max(q[s], lb[s] >> 1);
			log("*** lb = (%d, %d) ***", lb[0], lb[1]);

			Sub = (flags & FLAG_CORE) ? core(g, lb[1] - k, lb[0] - k) : g;
			logBiGraph(Sub);

			if (Sub.numVertices(0) < lb[0] || Sub.numVertices(1) < lb[1])
				continue;

			if (flags & FLAG_CN)
				Sub = comm(Sub, lb, k);
			logBiGraph(Sub);

			if (Sub.numVertices(0) < lb[0] || Sub.numVertices(1) < lb[1])
				continue;

			auto branchStartTime = std::chrono::steady_clock::now();
			if ((flags & FLAG_ORDER) && lb[0] > k && lb[1] > k)
			{
				log("Search scheme: RussianDoll");
				russianDoll();
			}
			else
			{
				log("Search scheme: All");
				searchAll();
			}
			branchTime += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - branchStartTime).count();
		}
	}
	else
	{
		lb[0] = q[0];
		lb[1] = q[1];
		log("*** lb = (%d, %d) ***", lb[0], lb[1]);

		Sub = (flags & FLAG_CORE) ? core(g, lb[1] - k, lb[0] - k) : g;
		logBiGraph(Sub);

		if (flags & FLAG_CN)
			Sub = comm(Sub, lb, k);
		logBiGraph(Sub);

		auto branchStartTime = std::chrono::steady_clock::now();
		if ((flags & FLAG_ORDER) && lb[0] > k && lb[1] > k)
		{
			log("Search scheme: RussianDoll");
			russianDoll();
		}
		else
		{
			log("Search scheme: All");
			searchAll();
		}
		branchTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - branchStartTime).count();
	}

	int searchTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - searchStartTime).count();
	log("Searching done! Total time: %d ms, Branch time: %d ms, branch num: %d (pivoting: %d, bipartite: %d), ub pruned num: %d",
		searchTime, branchTime, numBranches, numPivoting, numBipartite, numUbPruned);
	log("G[S*] info: |U|=%d, |V|=%d, |E|=%d", Ss[0].size(), Ss[1].size(), numEdgesSs);

	logSetInfo(Ss[0], "U*");
	logSetInfo(Ss[1], "V*");

	log("Missing edges:");
	for (int u : Ss[0])
		for (int v : Ss[1])
			if (!G.connect(u, v))
				log("(%d, %d)", u, v);
}

BiGraph MDB::core(BiGraph &G, int a, int b)
{
	std::queue<int> q[2];
	std::vector<bool> vis[2] = {std::vector<bool>(G.n[0]), std::vector<bool>(G.n[1])};
	std::vector<int> deg[2] = {std::vector<int>(G.n[0]), std::vector<int>(G.n[1])};
	int c[2] = {std::max(a, 0), std::max(b, 0)};

	for (int s = 0; s <= 1; ++s)
	{
		for (int v : G.V[s])
		{
			deg[s][v] = G.degree(s, v);
			if (deg[s][v] < c[s])
			{
				vis[s][v] = true;
				q[s].push(v);
			}
		}
	}

	while (!q[0].empty() || !q[1].empty())
	{
		for (int s = 0; s <= 1; ++s)
		{
			while (!q[s].empty())
			{
				int u = q[s].front();
				q[s].pop();
				for (int v : G.nbr[s][u])
					if (!vis[s ^ 1][v] && --deg[s ^ 1][v] < c[s ^ 1])
					{
						vis[s ^ 1][v] = true;
						q[s ^ 1].push(v);
					}
			}
		}
	}

	BiGraph Core;

	for (int u : G.V[0])
		if (!vis[0][u])
			for (int v : G.nbr[0][u])
				if (!vis[1][v])
					Core.addEdge(u, v);

	return Core;
}

BiGraph MDB::comm(BiGraph &G, int q[2], int k)
{
	std::vector<bool> visv[2] = {std::vector<bool>(G.n[0]), std::vector<bool>(G.n[1])};
	std::vector<int> deg[2] = {std::vector<int>(G.n[0]), std::vector<int>(G.n[1])};

	thread_local static std::vector<int> cn(std::max(G.n[0], G.n[1]));

#pragma omp parallel num_threads(numThreads)
	{
		cn.resize(std::max(G.n[0], G.n[1]));
	}

	for (int s = 0; s <= 1; ++s)
	{

		for (int u : G.V[s])
		{

			deg[s][u] = G.degree(s, u);

		}
	}

	for (int s = 0; s <= 1; ++s)
	{
#pragma omp parallel for num_threads(numThreads) schedule(dynamic)

		for (int u : G.V[s])
		{

			if (visv[s][u])
				continue;
			for (int v : G.nbr[s][u])
				if (!visv[s ^ 1][v]  )
					for (int w : G.nbr[s ^ 1][v])
						if (!visv[s][w]  )
							cn[w] = 0;
			for (int v : G.nbr[s][u])
				if (!visv[s ^ 1][v]  )
					for (int w : G.nbr[s ^ 1][v])
						if (!visv[s][w]  )
							++cn[w];

			int cntv = 0;
			for (int v : G.nbr[s][u])
			{
				int cntw = 0;
				if (visv[s ^ 1][v]  )
					continue;
				for (int w : G.nbr[s ^ 1][v])
					if (!visv[s][w] && cn[w] >= q[s ^ 1] - k && ++cntw >= q[s] - k)
						break;
				if (cntw < q[s] - k)
				{

				}
				else if (++cntv >= q[s ^ 1] - k)
					break;
			}

			if (cntv < q[s ^ 1] - k)
			{
				visv[s][u] = true;
				for (int v : G.nbr[s][u])
				{
#pragma omp atomic
					--deg[s ^ 1][v];
					if (deg[s ^ 1][v] < q[s] - k)
						visv[s ^ 1][v] = true;
				}
			}
		}
	}

	BiGraph CN;
	for (int u : G.V[0])
		if (!visv[0][u])
			for (int v : G.nbr[0][u])
				if (!visv[1][v]  )
					CN.addEdge(u, v);
	return CN;
}

void MDB::heuristic(BiGraph &G)
{

	for (int s = 0; s <= 1; ++s)
	{
		S[s].clear();
		C[s].clear();
		degS[s].assign(G.n[s], 0);
	}

	Ordering o;
	o.degeneracyOrdering(G);

	for (int i = o.numOrdered - 1; i >= 0; --i)
	{
		int s = o.ordered[i] >= G.n[0];
		int u = o.ordered[i] - s * G.n[0];
		if (numNnbS + nnbS(s, u) <= k)
		{
			S[s].push(u);
			for (int v : G.nbr[s][u])
				++degS[s ^ 1][v];
			numNnbS += nnbS(s, u);
		}
	}

	if (S[0].size() >= lb[0] && S[1].size() >= lb[1] && numEdgesS > numEdgesSs)
	{
		numNnbSs = numNnbS;
		for (int s = 0; s <= 1; ++s)
		{
			Ss[s].clear();
			for (int v : S[s])
				Ss[s].push(v);
		}
	}

	for (int s = 0; s <= 1; ++s)
	{

		S[s].clear();

		for (int v : G.V[s ^ 1])
		{
			S[s ^ 1].push(v);
			degS[s ^ 1][v] = 0;
		}

		for (int v : G.V[s])
		{
			C[s].push(v);
			degS[s][v] = G.degree(s, v);
		}

		LinearHeap U(C[s], [&](int v)
					 { return degS[s][v] + G.maxDeg[s] * o.value[v + s * G.n[0]]; });
		LinearHeap V(S[s ^ 1], [&](int v)
					 { return degS[s ^ 1][v]; });

		numNnbS = 0;

		std::vector<std::vector<int>> bin(k + 2);

		while (S[s].size() < lb[s] && C[s].size() > 0 && !U.empty())
		{
			int u = U.max();
			U.popMax();

			if (degS[s][u] < lb[s ^ 1] - k)
				break;

			S[s].push(u);
			numNnbS += nnbS(s, u);
			for (int v : G.nbr[s][u])
			{
				V.inc(v);
				++degS[s ^ 1][v];
			}

			int nd = 0;

			for (int i = 0; i < lb[s ^ 1]; ++i)
				nd += nnbS(s ^ 1, V.max(i));

			if (nd <= k)
			{
				while (numNnbS > k)
				{
					int v = V.min();
					V.popMin();
					numNnbS -= nnbS(s ^ 1, v);
					S[s ^ 1].pop(v);
					for (int w : G.nbr[s ^ 1][v])
					{
						U.dec(w);
						--degS[s][w];
					}
				}
			}
			else
			{
				S[s].pop(u);
				numNnbS -= nnbS(s, u);
				for (int v : G.nbr[s][u])
				{
					V.dec(v);
					--degS[s ^ 1][v];
				}
			}
		}

		if (S[0].size() >= lb[0] && S[1].size() >= lb[1] && numEdgesS > numEdgesSs)
		{
			numNnbSs = numNnbS;
			for (int s = 0; s <= 1; ++s)
			{
				Ss[s].clear();
				for (int v : S[s])
					Ss[s].push(v);
			}
		}
	}
}

void MDB::searchAll()
{
	numNnbS = 0;
	for (int s = 0; s <= 1; ++s)
	{
		S[s].clear();
		C[s].clear();
		for (int u : G.V[s])
		{
			degS[s][u] = 0;
			degC[s][u] = G.degree(s, u);
			C[s].push(u);
		}
	}

	branch(0);
}

void MDB::russianDoll()
{
	static Ordering o;
	thread_local static std::vector<int> q[2];

	int tSide = G.maxDeg[1] > G.maxDeg[0];

	BiGraph G0 = std::move(G);

	auto prune = [&]()
	{
		for (int s = 0; s <= 1; ++s)
		{
			q[s].clear();
			for (int v : S[s])
				if (degC[s][v] + degS[s][v] < lb[s ^ 1] - k)
					return false;
			for (int v : C[s])
				if (degC[s][v] + degS[s][v] < lb[s ^ 1] - k)
				{
					q[s].push_back(v);
					C[s].pop(v);
				}
		}

		while (!q[0].empty() || !q[1].empty())
		{
			for (int s = 0; s <= 1; ++s)
			{
				while (!q[s].empty())
				{
					int v = q[s].back();
					q[s].pop_back();

					for (int w : G0.nbr[s][v])
						if (S[s ^ 1].inside(w) || C[s ^ 1].inside(w))
						{
							if (--degC[s ^ 1][w] + degS[s ^ 1][w] < lb[s] - k)
							{
								if (C[s ^ 1].inside(w))
								{
									q[s ^ 1].push_back(w);
									C[s ^ 1].pop(w);
								}
								else
									return false;
							}
						}
				}
			}
		}
		return true;
	};

	o.degeneracyOrdering(G0);

#pragma omp parallel num_threads(numThreads)
	{
		for (int s = 0; s <= 1; ++s)
		{
			degS[s].assign(G0.n[s], 0);
			degC[s].assign(G0.n[s], 0);
		}
		G.clear();
		for (int s = 0; s <= 1; ++s)
			G.nbrMap[s] = G0.nbrMap[s];
	}

#pragma omp parallel for num_threads(numThreads) schedule(dynamic)
	for (int i = o.numOrdered - 1; i >= 0; --i)
	{
		int uSide = o.ordered[i] >= G0.n[0];
		int u = o.ordered[i] - uSide * G0.n[0];

		if (uSide != tSide)
			continue;

		numNnbS = 0;
		for (int s = 0; s <= 1; ++s)
		{
			S[s].clear();
			C[s].clear();
		}

		S[uSide].push(u);
		degS[uSide][u] = degC[uSide][u] = 0;

		for (int v : G0.nbr[uSide][u])
		{
			if (!C[uSide ^ 1].inside(v))
			{
				degS[uSide ^ 1][v] = degC[uSide ^ 1][v] = 0;
				C[uSide ^ 1].push(v);
			}
			++degC[uSide][u];
			++degS[uSide ^ 1][v];
			for (int w : G0.nbr[uSide ^ 1][v])
				if (w != u && o.order[w + uSide * G0.n[0]] > i)
				{
					if (!C[uSide].inside(w))
					{
						degS[uSide][w] = degC[uSide][w] = 0;
						C[uSide].push(w);
					}
					++degC[uSide ^ 1][v];
					++degC[uSide][w];
				}
		}

		if (!prune())
			continue;

		if (k > 0)
		{

			for (int w : C[uSide])
			{
				for (int x : G0.nbr[uSide][w])
					if (!G0.connect(uSide, u, x))
					{
						if (!C[uSide ^ 1].inside(x))
						{
							degS[uSide ^ 1][x] = degC[uSide ^ 1][x] = 0;
							C[uSide ^ 1].push(x);
						}

						++degC[uSide][w];
						++degC[uSide ^ 1][x];
					}
			}

			if (!prune())
				continue;
		}

		G.subGraph(G0, S, C);

		branch(0);

	}
}

bool MDB::upperbound()
{
	jsbcrit::witnessValid = false;
	jsbcrit::witnessI = 0;
	jsbcrit::witnessJ = 0;
	jsbcrit::activeMask = 0;
	jsbcrit::attackH[0] = jsbcrit::attackH[1] = -1;
	jsbcrit::attackDrop[0] = jsbcrit::attackDrop[1] = 0;

	auto oldUpperboundOnly = [&]() -> bool
	{
		for (int s = 0; s <= 1; ++s)
		{
			for (int u : S[s])
				if (degSub(s, u) < lb[s ^ 1] - k)
					return false;
		}

		int rem = k - numNnbS;
		if (rem < 0)
			return false;

		thread_local static std::vector<int> bin[2];
		for (int s = 0; s <= 1; ++s)
		{
			if ((int)bin[s].size() < k + 1)
				bin[s].resize(k + 1, 0);
		}

		int idx[2] = {rem, rem}, sum[2] = {0, 0}, num[2] = {0, 0};

		for (int s = 0; s <= 1; ++s)
		{
			for (int i = 0; i <= rem; ++i)
				bin[s][i] = 0;
			for (int v : C[s])
				++bin[s][nnbS(s, v)];
			for (int i = 0; i <= rem; ++i)
			{
				if (sum[s] + bin[s][i] * i > rem)
				{
					bin[s][i] = (k - sum[s] - numNnbS) / i;
					sum[s] += bin[s][i] * i;
					num[s] += bin[s][i];
					idx[s] = i;
					break;
				}
				sum[s] += bin[s][i] * i;
				num[s] += bin[s][i];
			}
		}

		if (num[0] + S[0].size() < lb[0] || num[1] + S[1].size() < lb[1])
			return false;

		int s = (int)(num[0] > num[1]);
		long long numEdges = (long long)S[s].size() * (num[s ^ 1] + S[s ^ 1].size()) - numNnbS - sum[s ^ 1];

		sum[s] = num[s] = 0;

		for (int i = 0, j = idx[s ^ 1]; i <= idx[s] && j >= 0; ++i)
		{
			while (bin[s][i]--)
			{
				while (j > 0 && sum[0] + sum[1] + i > rem)
				{
					int x = std::min(bin[s ^ 1][j], (sum[0] + sum[1] + i - k + numNnbS - 1) / j + 1);
					sum[s ^ 1] -= x * j;
					num[s ^ 1] -= x;
					bin[s ^ 1][j] -= x;
					j -= (int)(!bin[s ^ 1][j]);
				}
				if (sum[0] + i + sum[1] > rem)
					return numEdges > numEdgesSs;
				sum[s] += i;
				++num[s];
				numEdges = std::max(numEdges, (long long)(num[0] + S[0].size()) * (num[1] + S[1].size()) - numNnbS - sum[0] - sum[1]);
				if (numEdges > numEdgesSs)
					return true;
			}
		}

		return numEdges > numEdgesSs;
	};

	if (!oldUpperboundOnly())
		return false;

	long long sU = S[0].size(), sV = S[1].size();
	long long edgesS = sU * sV - numNnbS;
	long long kRemainLL = (long long)k - numNnbS;
	if (kRemainLL < 0)
		return false;
	int kRemain = (int)kRemainLL;

	thread_local static std::vector<int> dCnt[2];
	thread_local static std::vector<std::pair<int, int>> cand[2];
	thread_local static std::vector<long long> prefixD[2];
	thread_local static std::vector<int> degHist[2];

	for (int s = 0; s <= 1; ++s)
	{
		dCnt[s].assign(kRemain + 1, 0);
		cand[s].clear();
		if (cand[s].capacity() < C[s].size())
			cand[s].reserve(C[s].size());

		for (int v : C[s])
		{
			int ds = (int)S[s ^ 1].size() - degS[s][v];
			if (ds < 0)
				ds = 0;
			if (ds > kRemain)
				continue;

			int dc = degC[s][v];
			cand[s].push_back(std::make_pair(ds, dc));
			++dCnt[s][ds];
		}

		fillPrefixFromDSCount(dCnt[s], prefixD[s]);
	}

	int maxAdd[2] = {
		maxPrefixCount(prefixD[0], kRemainLL),
		maxPrefixCount(prefixD[1], kRemainLL)};

	int need[2] = {
		std::max(0, lb[0] - (int)S[0].size()),
		std::max(0, lb[1] - (int)S[1].size())};

	if (maxAdd[0] < need[0] || maxAdd[1] < need[1])
		return false;

	const long long pairCount =
		(long long)(maxAdd[0] - need[0] + 1) *
		(long long)(maxAdd[1] - need[1] + 1);

	for (int s = 0; s <= 1; ++s)
	{
		int maxCap = maxAdd[s ^ 1];
		degHist[s].assign(maxCap + 1, 0);
		for (const auto &x : cand[s])
		{
			int dc = std::min(x.second, maxCap);
			++degHist[s][dc];
		}
	}

	if (pairCount <= SHORTAGE_EXACT_PAIR_CAP)
	{

		std::vector<std::vector<long long>> eccURow(maxAdd[1] + 1);
		std::vector<std::vector<long long>> eccVRow(maxAdd[0] + 1);
		std::vector<char> eccUBuilt(maxAdd[1] + 1, 0);
		std::vector<char> eccVBuilt(maxAdd[0] + 1, 0);

		auto getEccU = [&](int capJ, int cntI) -> long long
		{
			if (cntI <= 0)
				return 0;
			if (!eccUBuilt[capJ])
			{
				eccURow[capJ] =
					buildTopMinRowFromDegreeHist(degHist[0], maxAdd[0], capJ);
				eccUBuilt[capJ] = 1;
			}
			return eccURow[capJ][std::min(cntI, maxAdd[0])];
		};

		auto getEccV = [&](int capI, int cntJ) -> long long
		{
			if (cntJ <= 0)
				return 0;
			if (!eccVBuilt[capI])
			{
				eccVRow[capI] =
					buildTopMinRowFromDegreeHist(degHist[1], maxAdd[1], capI);
				eccVBuilt[capI] = 1;
			}
			return eccVRow[capI][std::min(cntJ, maxAdd[1])];
		};

		auto eccUB = [&](int i, int j) -> long long
		{
			long long fromU = getEccU(j, i);
			long long fromV = getEccV(i, j);
			long long ub = std::min(fromU, fromV);
#ifdef DEBUG_SHORTAGE_UB
			long long totalPairs = (long long)i * j;
			assert(ub <= totalPairs);
			assert(totalPairs - ub >= 0);
#endif
			return ub;
		};

		bool jointReady = false;
		std::vector<std::vector<int>> jointHist[2];
		std::vector<std::vector<long long>> suffixU, suffixV;
		std::vector<std::vector<long long>> jointURow, jointVRow;
		std::vector<char> builtU, builtV;

		auto ensureJoint = [&]()
		{
			if (jointReady)
				return;

			for (int s = 0; s <= 1; ++s)
			{
				int maxCap = maxAdd[s ^ 1];
				jointHist[s].assign(
					kRemain + 1,
					std::vector<int>(maxCap + 1, 0));

				for (const auto &x : cand[s])
				{
					int ds = x.first;
					int dc = std::min(x.second, maxCap);
					++jointHist[s][ds][dc];
				}
			}

			suffixU = buildJointSuffix(jointHist[0], kRemain, maxAdd[1]);
			suffixV = buildJointSuffix(jointHist[1], kRemain, maxAdd[0]);

			jointURow.resize(maxAdd[1] + 1);
			jointVRow.resize(maxAdd[0] + 1);
			builtU.assign(maxAdd[1] + 1, 0);
			builtV.assign(maxAdd[0] + 1, 0);
			jointReady = true;
		};

		auto getJointU = [&](int capJ, int cntI) -> long long
		{
			if (cntI <= 0)
				return 0;
			ensureJoint();
			if (!builtU[capJ])
			{
				const int JSB_BATCH_ATTACK_CAP_LOCAL = 8;
				const int rowCount =
					std::min((int)cand[0].size(),
							 maxAdd[0] + JSB_BATCH_ATTACK_CAP_LOCAL);
				jointURow[capJ] =
					buildJointRowFromBucket(
						jointHist[0], suffixU,
						kRemain, rowCount, maxAdd[1], capJ);
				builtU[capJ] = 1;
			}
			return jointURow[capJ][std::min(cntI, maxAdd[0])];
		};

		auto getJointV = [&](int capI, int cntJ) -> long long
		{
			if (cntJ <= 0)
				return 0;
			ensureJoint();
			if (!builtV[capI])
			{
				const int JSB_BATCH_ATTACK_CAP_LOCAL = 8;
				const int rowCount =
					std::min((int)cand[1].size(),
							 maxAdd[1] + JSB_BATCH_ATTACK_CAP_LOCAL);
				jointVRow[capI] =
					buildJointRowFromBucket(
						jointHist[1], suffixV,
						kRemain, rowCount, maxAdd[0], capI);
				builtV[capI] = 1;
			}
			return jointVRow[capI][std::min(cntJ, maxAdd[1])];
		};

		for (int i = need[0]; i <= maxAdd[0]; ++i)
		{
			if (prefixD[0][i] > kRemainLL)
				break;

			for (int j = need[1]; j <= maxAdd[1]; ++j)
			{
				long long baseDef = prefixD[0][i] + prefixD[1][j];
				if (baseDef > kRemainLL)
					break;

				long long ij = (long long)i * j;

				long long oldPairUB =
					edgesS +
					(long long)i * sV - prefixD[0][i] +
					(long long)j * sU - prefixD[1][j] +
					ij;
				if (oldPairUB <= numEdgesSs)
					continue;

				long long candidateEccUB = eccUB(i, j);
				long long candidateShortage = ij - candidateEccUB;

				if (baseDef + candidateShortage > kRemainLL)
					break;

				long long cheapUB = oldPairUB - candidateShortage;
				if (cheapUB <= numEdgesSs)
					continue;

				long long compBase =
					candidateEccUB - prefixD[0][i] - prefixD[1][j];
				long long compJointU =
					getJointU(j, i) - prefixD[1][j];
				long long compJointV =
					getJointV(i, j) - prefixD[0][i];

				const long long rawComp =
					std::min(compBase, std::min(compJointU, compJointV));

				unsigned char activeMask = 0;
				if (compBase == rawComp)
					activeMask |= 1u;
				if (compJointU == rawComp)
					activeMask |= 2u;
				if (compJointV == rawComp)
					activeMask |= 4u;

				long long comp = rawComp;
				if (ij < comp)
				{
					comp = ij;

					activeMask = 8u;
				}

				long long defectLB = ij - comp;
				if (defectLB > kRemainLL)
					continue;

				long long ub =
					edgesS +
					(long long)i * sV +
					(long long)j * sU +
					comp;

				if (ub > numEdgesSs)
				{
					jsbcrit::witnessValid = true;
					jsbcrit::witnessI = i;
					jsbcrit::witnessJ = j;
					jsbcrit::activeMask = activeMask;

					const int JSB_BATCH_ATTACK_CAP_LOCAL = 8;
					const long long witnessGap =
						std::max(1LL, ub - (long long)numEdgesSs);

					auto computeAttackH =
						[&](int side,
							const std::vector<long long> &row,
							int selected,
							int candCount)
					{
						if (selected <= 0)
							return;

						const unsigned char jointBit =
							(side == 0 ? 2u : 4u);
						if ((activeMask & jointBit) == 0)
							return;

						const int hMax =
							std::min(JSB_BATCH_ATTACK_CAP_LOCAL, candCount);

						for (int h = 1; h <= hMax; ++h)
						{

							if (candCount - h < selected)
							{
								jsbcrit::attackH[side] = h;
								jsbcrit::attackDrop[side] = witnessGap;
								return;
							}

							if (selected + h >= (int)row.size())
								break;

							const long long newTop =
								row[selected + h] - row[h];
							const long long drop =
								row[selected] - newTop;

							if (drop >= witnessGap)
							{
								jsbcrit::attackH[side] = h;
								jsbcrit::attackDrop[side] = drop;
								return;
							}
						}
					};

					if (i > 0 &&
						jointReady &&
						j >= 0 && j < (int)jointURow.size() &&
						j < (int)builtU.size() && builtU[j])
					{
						computeAttackH(
							0, jointURow[j], i, (int)cand[0].size());
					}

					if (j > 0 &&
						jointReady &&
						i >= 0 && i < (int)jointVRow.size() &&
						i < (int)builtV.size() && builtV[i])
					{
						computeAttackH(
							1, jointVRow[i], j, (int)cand[1].size());
					}

					return true;
				}
			}
		}

		return false;
	}

	auto eccUBCheap = [&](int i, int j) -> long long
	{
		long long fromU = sumTopMinFromDegreeHist(degHist[0], i, j);
		long long fromV = sumTopMinFromDegreeHist(degHist[1], j, i);
		long long ub = std::min(fromU, fromV);
#ifdef DEBUG_SHORTAGE_UB
		long long totalPairs = (long long)i * j;
		assert(ub <= totalPairs);
		assert(totalPairs - ub >= 0);
#endif
		return ub;
	};

	long long bestCrossToVS = 0;
	for (int i = need[0]; i <= maxAdd[0]; ++i)
		bestCrossToVS =
			std::max(bestCrossToVS,
					 (long long)i * sV - prefixD[0][i]);

	long long bestCrossToUS = 0;
	for (int j = need[1]; j <= maxAdd[1]; ++j)
		bestCrossToUS =
			std::max(bestCrossToUS,
					 (long long)j * sU - prefixD[1][j]);

	long long newUB =
		edgesS +
		bestCrossToVS +
		bestCrossToUS +
		eccUBCheap(maxAdd[0], maxAdd[1]);

#ifdef DEBUG_SHORTAGE_UB
	assert(newUB >= 0);
#endif
	return newUB > numEdgesSs;
}

MDB::BakPos MDB::backup(const VertexSet V[2])
{
	MDB::BakPos pos;
	pos.backup(V);
	return pos;
}

MDB::BakPos MDB::update(int uSide, int u)
{
	using namespace traversal;

	auto pos = backup(C);

	moveC2S(uSide, u);

	for (int s = 0; s <= 1; ++s)
	{
		for (int v : C[s])
		{

			if (numNnbS   + nnbS(s, v)   > k)
			{
				subC(s, v);
			}
		}
	}

	for (int s = 0; s <= 1; ++s)
	{
		for (int v : C[s])
			if (degSub(s, v) < lb[s ^ 1] - k + numNnbS)
			{
				subC(s, v);

			}
	}

	for (int s = 0; s <= 1; ++s)
	{
		for (int v : reverse(C[s]))
			if (nnbSub(s, v) == 0)
			{
				moveC2S(s, v);
			}
	}

	return pos;
}

MDB::BakPos MDB::minus(int uSide, int u)
{
	using namespace traversal;

	auto pos = backup(C);

	subC(uSide, u);

	if ((flags & FLAG_1NN) && nnbS(uSide, u) + nnbC(uSide, u) == 1)
	{

		int uNnbC = -1;
		if (nnbC(uSide, u) == 1)
		{
			for (int v : C[uSide ^ 1])
				if (!G.connect(uSide, u, v))
				{
					uNnbC = v;
					break;
				}
		}

		for (int v : C[uSide])
			if (nnbS(uSide, v) > 0 || (uNnbC != -1 && !G.connect(uSide, v, uNnbC)))
			{
				subC(uSide, v);
			}
	}

	for (int s = 0; s <= 1; ++s)
	{
		for (int v : C[s])
			if (degSub(s, v) < lb[s ^ 1] - k + numNnbS)
			{
				subC(s, v);

			}
	}

	for (int s = 0; s <= 1; ++s)
	{
		for (int v : reverse(C[s]))
			if (nnbSub(s, v) == 0)
			{
				moveC2S(s, v);
			}
	}

	return pos;
}

void MDB::restore(BakPos &pos)
{
	for (int s = 0; s <= 1; ++s)
	{
		for (int i = C[s].backPos(); i < pos.p[1][s]; ++i)
		{
			moveS2C(s, C[s][i]);
		}
	}

	for (int s = 0; s <= 1; ++s)
	{
		for (int i = C[s].frontPos() - 1; i >= pos.p[0][s]; --i)
		{
			addC(s, C[s][i]);
		}
	}
}
