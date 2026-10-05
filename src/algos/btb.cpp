#include "btb.h"

#include "../utils/log.hpp"
#include <algorithm>
#include <vector>
#include <queue>
#include <cstdint>
#include <functional>
#include <utility>

using namespace logging;

namespace jsbcrit
{
	extern thread_local bool witnessValid;
	extern thread_local int witnessI;
	extern thread_local int witnessJ;
	extern thread_local unsigned char activeMask;
	extern thread_local int attackH[2];
	extern thread_local long long attackDrop[2];
}

namespace
{
	const int JSB_BSET_ATTACK_CAP = 8;
	const int PROBE_BUDGET_BETA = 16;
	const double CD2_RISK_THRESHOLD = 1.0;
	const int SELECTOR_EXTRA_ORIGINAL_LEVELS = 1;
	const int PARTIAL_RESIDUAL_ORIGINAL_LEVELS = 1;

	thread_local int selectorOriginalLevels = 0;
	thread_local int residualOriginalLevels = 0;

	struct ScopedIntRestore
	{
		int &value;
		int oldValue;

		explicit ScopedIntRestore(int &x)
			: value(x), oldValue(x)
		{
		}

		~ScopedIntRestore()
		{
			value = oldValue;
		}
	};

	struct WitnessSnapshot
	{
		bool valid = false;
		int i = 0;
		int j = 0;
		unsigned char activeMask = 0;
		int attackH[2] = {-1, -1};
		long long attackDrop[2] = {0, 0};
	};

	WitnessSnapshot saveWitness()
	{
		WitnessSnapshot x;
		x.valid = jsbcrit::witnessValid;
		x.i = jsbcrit::witnessI;
		x.j = jsbcrit::witnessJ;
		x.activeMask = jsbcrit::activeMask;
		x.attackH[0] = jsbcrit::attackH[0];
		x.attackH[1] = jsbcrit::attackH[1];
		x.attackDrop[0] = jsbcrit::attackDrop[0];
		x.attackDrop[1] = jsbcrit::attackDrop[1];
		return x;
	}

	void restoreWitness(const WitnessSnapshot &x)
	{
		jsbcrit::witnessValid = x.valid;
		jsbcrit::witnessI = x.i;
		jsbcrit::witnessJ = x.j;
		jsbcrit::activeMask = x.activeMask;
		jsbcrit::attackH[0] = x.attackH[0];
		jsbcrit::attackH[1] = x.attackH[1];
		jsbcrit::attackDrop[0] = x.attackDrop[0];
		jsbcrit::attackDrop[1] = x.attackDrop[1];
	}
}

SolverResult BTB::run(const std::string &dataPath, int q[2], int k, int flags, int numThreads)
{
	selectorOriginalLevels = 0;
	residualOriginalLevels = 0;

	BTB solver;
	solver.findMDB(dataPath, q, k, flags, numThreads);
	return solver.getResult();
}

void BTB::branch(int dep)
{
	if (flags & FLAG_DEBUG)
	{
		log("\n---------- dep=%d ----------", dep);
		logSet(S[0]);
		logSet(S[1]);
		logSet(C[0]);
		logSet(C[1]);
		for (int s = 0; s <= 1; ++s)
		{
			for (int u : G.V[s])
			{
				int degSu = 0, degCu = 0;
				for (int v : G.nbr[s][u])
				{
					if (S[s ^ 1].inside(v))
						++degSu;
					if (C[s ^ 1].inside(v))
						++degCu;
				}
				if (degSu != degS[s][u] || degCu != degC[s][u])
					log("Wrong degree of %d-%d", s, u);
			}
		}
	}

	if (C[0].size() == 0 && C[1].size() == 0)
	{
		if (S[0].size() >= lb[0] && S[1].size() >= lb[1])
		{
#pragma omp critical(update)
			{
				if (numEdgesS > numEdgesSs)
				{
					for (int s = 0; s <= 1; ++s)
					{
						Ss[s].clear();
						for (int v : S[s])
							Ss[s].push(v);
					}
					numNnbSs = numNnbS;
					log("New S*! |E|=%d", numEdgesSs);
				}
			}
		}
		return;
	}

	if (S[0].size() + C[0].size() < lb[0] ||
		S[1].size() + C[1].size() < lb[1])
		return;

	if (flags & FLAG_UB)
	{
		int numEdgesSub = 0;
		int s = (int)(S[0].size() + C[0].size() >
					  S[1].size() + C[1].size());

		for (int v : S[s])
			numEdgesSub += degSub(s, v);
		for (int v : C[s])
			numEdgesSub += degSub(s, v);

		if (numEdgesSub <= numEdgesSs)
			return;

		if (!upperbound())
		{
			++numUbPruned;
			return;
		}
	}
	++numBranches;

	struct PivotChoice
	{
		int vertex = -1;
		int side = -1;
		int nnbC = 0;
	};

	auto findTruePivot = [&]() -> PivotChoice
	{
		PivotChoice choice;
		choice.nnbC = G.n[0] + G.n[1];

		for (int s = 0; s <= 1; ++s)
		{
			for (int v : C[s])
			{
				if (nnbS(s, v) == 0 && nnbC(s, v) < choice.nnbC)
				{
					choice.vertex = v;
					choice.side = s;
					choice.nnbC = nnbC(s, v);
				}
			}
		}

		return choice;
	};

	auto branchOriginalOneLayer = [&]() -> void
	{
		const PivotChoice pivot = findTruePivot();

		if (pivot.vertex == -1)
		{
			++numBipartite;

			int u = -1;
			int uSide = -1;
			int bestDS = -1;

			for (int s = 0; s <= 1; ++s)
			{
				for (int v : C[s])
				{
					const int ds = nnbS(s, v);
					if (ds > bestDS)
					{
						bestDS = ds;
						u = v;
						uSide = s;
					}
				}
			}

			if (u == -1)
				return;

			auto includePos = update(uSide, u);
			branch(dep + 1);
			restore(includePos);

			auto excludePos = minus(uSide, u);
			branch(dep + 1);
			restore(excludePos);
			return;
		}

		++numPivoting;

		BakPos parentPos;
		parentPos.backup(C);

		std::vector<int> nonNeighbors;
		for (int v : C[pivot.side ^ 1])
		{
			if (!G.connect(pivot.side, pivot.vertex, v))
				nonNeighbors.push_back(v);
		}

		auto pivotPos = update(pivot.side, pivot.vertex);
		branch(dep + 1);
		restore(pivotPos);

		minus(pivot.side, pivot.vertex);

		bool alreadyForced = false;
		for (int v : nonNeighbors)
		{
			if (S[pivot.side ^ 1].inside(v))
			{
				alreadyForced = true;
				branch(dep + 1);
				break;
			}
		}

		if (!alreadyForced)
		{
			for (int v : nonNeighbors)
			{
				if (!C[pivot.side ^ 1].inside(v))
					continue;

				auto posV = update(pivot.side ^ 1, v);
				branch(dep + 1);
				restore(posV);
				minus(pivot.side ^ 1, v);
			}
		}

		restore(parentPos);
	};

	if (residualOriginalLevels > 0)
	{
		ScopedIntRestore restoreResidualState(residualOriginalLevels);
		--residualOriginalLevels;
		branchOriginalOneLayer();
		return;
	}

	bool selectorCooldownCancelled = false;
	ScopedIntRestore restoreSelectorState(selectorOriginalLevels);

	if (selectorOriginalLevels > 0)
	{
		const PivotChoice forcedPivot = findTruePivot();

		if (forcedPivot.vertex != -1)
		{
			--selectorOriginalLevels;
			branchOriginalOneLayer();
			return;
		}

		selectorOriginalLevels = 0;
		selectorCooldownCancelled = true;
	}

	if (!selectorCooldownCancelled)
	{
		const PivotChoice pivot = findTruePivot();

		if (pivot.vertex != -1)
		{
			const int originalFanout = std::max(1, 1 + pivot.nnbC);

			int attackH = -1;
			long long bestDrop = -1;

			if ((flags & FLAG_UB) && jsbcrit::witnessValid)
			{
				for (int s = 0; s <= 1; ++s)
				{
					const int h = jsbcrit::attackH[s];
					if (h <= 0 || h > JSB_BSET_ATTACK_CAP)
						continue;

					if (attackH == -1 ||
						h < attackH ||
						(h == attackH && jsbcrit::attackDrop[s] > bestDrop))
					{
						attackH = h;
						bestDrop = jsbcrit::attackDrop[s];
					}
				}
			}

			if (attackH > 0)
			{
				const double risk =
					(double)attackH / (double)originalFanout;

				if (risk >= CD2_RISK_THRESHOLD)
				{
					selectorOriginalLevels = SELECTOR_EXTRA_ORIGINAL_LEVELS;
					branchOriginalOneLayer();
					return;
				}
			}
		}
	}

	const std::uint64_t probeParentC =
		(std::uint64_t)C[0].size() + (std::uint64_t)C[1].size();
	const std::uint64_t probeBudget =
		(std::uint64_t)PROBE_BUDGET_BETA * probeParentC;
	std::uint64_t probeWork = 0;
	bool probeBudgetExceeded = false;

	const WitnessSnapshot parentWitness = saveWitness();

	BakPos shadow0;
	shadow0.backup(C);

	std::vector<std::pair<int, int>> branchSet;
	bool certified = false;

	auto updateIncumbentFromCurrentS = [&]() -> void
	{
		if (S[0].size() < lb[0] ||
			S[1].size() < lb[1] ||
			numEdgesS <= numEdgesSs)
			return;

#pragma omp critical(update)
		{
			if (numEdgesS > numEdgesSs)
			{
				for (int s = 0; s <= 1; ++s)
				{
					Ss[s].clear();
					for (int v : S[s])
						Ss[s].push(v);
				}
				numNnbSs = numNnbS;
			}
		}
	};

	auto appendMaxDSVertex = [&]() -> int
	{
		int bestSide = -1;
		int bestV = -1;
		int bestDS = -1;

		for (int s = 0; s <= 1; ++s)
		{
			for (int v : C[s])
			{
				const int ds = nnbS(s, v);
				if (ds > bestDS)
				{
					bestDS = ds;
					bestSide = s;
					bestV = v;
				}
			}
		}

		if (bestV == -1)
			return 0;

		branchSet.push_back({bestSide, bestV});
		subC(bestSide, bestV);
		return 1;
	};

	auto appendTopPhiBatch =
		[&](int side, int h, int cap) -> int
	{
		if (h <= 0)
			return 0;

		using HeapItem = std::pair<int, int>;
		std::priority_queue<
			HeapItem,
			std::vector<HeapItem>,
			std::greater<HeapItem>>
			top;

		const int kRemain = k - numNnbS;

		for (int v : C[side])
		{
			const int ds = nnbS(side, v);
			if (ds < 0 || ds > kRemain)
				continue;

			const int phi =
				std::min(degC[side][v], cap) - ds;

			const HeapItem item(phi, v);
			if ((int)top.size() < h)
				top.push(item);
			else if (item > top.top())
			{
				top.pop();
				top.push(item);
			}
		}

		if (top.empty())
			return 0;

		std::vector<HeapItem> selected;
		selected.reserve(top.size());
		while (!top.empty())
		{
			selected.push_back(top.top());
			top.pop();
		}

		std::sort(
			selected.begin(),
			selected.end(),
			[](const HeapItem &a, const HeapItem &b)
			{
				if (a.first != b.first)
					return a.first > b.first;
				return a.second > b.second;
			});

		int added = 0;
		for (const auto &x : selected)
		{
			const int v = x.second;
			if (!C[side].inside(v))
				continue;

			branchSet.push_back({side, v});
			subC(side, v);
			++added;
		}

		return added;
	};

	auto appendFallbackBatch = [&]() -> int
	{
		int side = -1;
		int cap = 0;

		if (jsbcrit::witnessValid)
		{
			const bool activeU = (jsbcrit::activeMask & 2u) != 0;
			const bool activeV = (jsbcrit::activeMask & 4u) != 0;

			if (activeU && C[0].size() != 0)
			{
				side = 0;
				cap = jsbcrit::witnessJ;
			}

			if (activeV && C[1].size() != 0)
			{
				if (side == -1 || C[1].size() < C[0].size())
				{
					side = 1;
					cap = jsbcrit::witnessI;
				}
			}
		}

		if (side != -1)
		{
			const int h =
				std::min(JSB_BSET_ATTACK_CAP, (int)C[side].size());
			const int added = appendTopPhiBatch(side, h, cap);
			if (added > 0)
				return added;
		}

		return appendMaxDSVertex();
	};

	while (!certified)
	{
		const std::uint64_t currentProbeC =
			(std::uint64_t)C[0].size() + (std::uint64_t)C[1].size();

		if (probeWork + currentProbeC > probeBudget)
		{
			probeBudgetExceeded = true;
			break;
		}

		probeWork += currentProbeC;

		int attackSide = -1;
		int attackH = -1;
		long long bestDrop = -1;

		if ((flags & FLAG_UB) && jsbcrit::witnessValid)
		{
			for (int s = 0; s <= 1; ++s)
			{
				const int h = jsbcrit::attackH[s];
				if (h <= 0 || h > JSB_BSET_ATTACK_CAP)
					continue;

				if (attackSide == -1 ||
					h < attackH ||
					(h == attackH && jsbcrit::attackDrop[s] > bestDrop))
				{
					attackSide = s;
					attackH = h;
					bestDrop = jsbcrit::attackDrop[s];
				}
			}
		}

		int added = 0;
		if (attackSide != -1)
		{
			const int cap =
				(attackSide == 0 ? jsbcrit::witnessJ : jsbcrit::witnessI);
			added = appendTopPhiBatch(attackSide, attackH, cap);
		}
		else
		{
			added = appendFallbackBatch();
		}

		if (added <= 0)
			added = appendMaxDSVertex();

		if (C[0].size() == 0 && C[1].size() == 0)
		{
			updateIncumbentFromCurrentS();
			certified = true;
			continue;
		}

		if (S[0].size() + C[0].size() < lb[0] ||
			S[1].size() + C[1].size() < lb[1])
		{
			certified = true;
			continue;
		}

		int numEdgesSub = 0;
		const int side =
			(int)(S[0].size() + C[0].size() >
				  S[1].size() + C[1].size());

		for (int v : S[side])
			numEdgesSub += degSub(side, v);
		for (int v : C[side])
			numEdgesSub += degSub(side, v);

		if (numEdgesSub <= numEdgesSs)
		{
			certified = true;
			continue;
		}

		if (flags & FLAG_UB)
			certified = !upperbound();
	}

	restore(shadow0);
	restoreWitness(parentWitness);

	std::sort(
		branchSet.begin(), branchSet.end(),
		[&](const std::pair<int, int> &a, const std::pair<int, int> &b)
		{
			const int sa = a.first, va = a.second, sb = b.first, vb = b.second;
			const int an = nnbC(sa, va), bn = nnbC(sb, vb);
			if (an != bn)
				return an > bn;
			const int as = nnbS(sa, va), bs = nnbS(sb, vb);
			if (as != bs)
				return as > bs;
			const int ad = degC[sa][va], bd = degC[sb][vb];
			if (ad != bd)
				return ad < bd;
			if (sa != sb)
				return sa < sb;
			return va < vb;
		});

	BakPos pos0;
	pos0.backup(C);

	bool residualPossible = true;
	for (int idx = 0; idx < (int)branchSet.size(); ++idx)
	{
		const int s = branchSet[idx].first;
		const int v = branchSet[idx].second;

		if (S[s].inside(v))
		{

			branch(dep + 1);
			residualPossible = false;
			break;
		}

		if (!C[s].inside(v))
			continue;

		auto posV = update(s, v);
		branch(dep + 1);
		restore(posV);

		if (!certified || idx + 1 < (int)branchSet.size())
			minus(s, v);
	}

	if (!certified && probeBudgetExceeded && residualPossible)
	{
		ScopedIntRestore restoreResidualState(residualOriginalLevels);
		residualOriginalLevels = PARTIAL_RESIDUAL_ORIGINAL_LEVELS;
		branch(dep + 1);
	}

	restore(pos0);
	return;
}
