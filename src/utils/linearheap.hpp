#ifndef LINEAR_HEAP
#define LINEAR_HEAP

#include "vertexset.hpp"
#include <exception>
#include <vector>
#include <algorithm>
#include <cassert>
#include <functional>

class LinearHeap {
	std::vector<int> rk, pos, bin;
	std::vector<int> keys, vals;
	std::vector<bool> vis;
	int pl, pr;
	int capKey, capVal;

	void resizeKey(int capacity) {
		capKey = capacity;
		pos.resize(capacity + 1);
		rk.resize(capacity);

		vals.resize(capacity);
		vis.resize(capacity);
	}
	void resizeVal(int capacity) {
		capVal = capacity;
		bin.resize(capacity + 1);
	}
public:
	LinearHeap(): pl(0), pr(0), capKey(0), capVal(0) {}
	LinearHeap(const std::vector<int> &values): capKey(0), capVal(0) {
		build(values);
	}
	LinearHeap(const std::vector<int> &keys, const std::vector<int> &values): capKey(0), capVal(0) {
		build(keys, values);
	}
	LinearHeap(const VertexSet &V, const std::vector<int> &values): capKey(0), capVal(0) {
		build(V, values);
	}
	LinearHeap(const VertexSet &V, const std::function<int(int)> &valueFunc): capKey(0), capVal(0) {
		build(V, valueFunc);
	}

	void build(const std::vector<int> &values) {
		if (values.size() > capKey) resizeKey(values.size());
		pl = pr = 0;
		int maxValue = 0;
		for (int i = 0; i < values.size(); ++i) {
			keys[i] = i;
			vis[i] = true;
			vals[i] = values[i];
			maxValue = std::max(maxValue, vals[i]);
		}
		if (maxValue >= capVal) resizeVal(maxValue+1);
		bin.assign(maxValue + 2, 0);
		for (int v : values) ++bin[v+1];
		for (int v = 1; v <= maxValue; ++v) bin[v] += bin[v-1];
		for (int i = 0; i < vals.size(); ++i) {
			rk[i] = bin[values[i]]++;
			pos[rk[i]] = i;
		}
	}

	void build(const std::vector<int> &keys, const std::vector<int> &values) {
		this->keys.assign(keys.begin(), keys.end());
		pl = pr = 0;
		int maxValue = 0, maxKey = 0;
		for (int k : keys) {
			maxValue = std::max(maxValue, vals[k]);
			maxKey = std::max(maxKey, k);
		}
		if (maxKey >= capKey) resizeKey(maxKey+1);
		if (maxValue >= capVal) resizeVal(maxValue+1);
		bin.assign(maxValue + 2, 0);
		for (int k : keys) ++bin[values[k]+1];
		for (int v = 1; v <= maxValue; ++v) bin[v] += bin[v-1];
		for (int k : keys) {
			vis[k] = true;
			vals[k] = values[k];
			rk[k] = bin[values[k]]++;
			pos[rk[k]] = k;
		}
		pr = bin[maxValue] - 1;
	}

	void build(const VertexSet &V, const std::vector<int> &values) {
		keys.assign(V.begin(), V.end());
		pl = pr = 0;
		int maxValue = 0, maxKey = 0;
		for (int k : keys) {
			maxValue = std::max(maxValue, values[k]);
			maxKey = std::max(maxKey, k);
		}
		if (maxKey >= capKey) resizeKey(maxKey+1);
		if (maxValue >= capVal) resizeVal(maxValue+1);
		bin.assign(maxValue + 2, 0);
		for (int k : V) ++bin[values[k]+1];
		for (int v = 1; v <= maxValue; ++v) bin[v] += bin[v-1];
		for (int k : V) {
			vis[k] = true;
			vals[k] = values[k];
			rk[k] = bin[values[k]]++;
			pos[rk[k]] = k;
		}
		pr = bin[maxValue] - 1;
	}

	void build(const VertexSet &V, const std::function<int(int)> &valueFunc) {
		keys.assign(V.begin(), V.end());
		pl = pr = 0;
		int maxValue = 0, maxKey = 0;
		for (int k : V) { maxKey = std::max(maxKey, k); }
		if (maxKey >= capKey) resizeKey(maxKey+1);
		for (int k : keys) {
			vals[k] = valueFunc(k);
			maxValue = std::max(maxValue, vals[k]);
		}
		if (maxValue >= capVal) resizeVal(maxValue+1);
		bin.assign(maxValue + 2, 0);
		for (int k : V) ++bin[vals[k]+1];
		for (int v = 1; v <= maxValue; ++v) bin[v] += bin[v-1];
		for (int k : V) {
			vis[k] = true;
			rk[k] = bin[vals[k]]++;
			pos[rk[k]] = k;
		}
		pr = bin[maxValue] - 1;
	}

	void inc(int k) {
		if (!inside(k)) return;
		if (vals[k]+1 >= capVal) { resizeVal(vals[k]+2); bin[vals[k]+1] = bin[vals[k]]; }
		if (bin[vals[k]] > pr+1) bin[vals[k]] = pr+1;
		int k2 = pos[--bin[vals[k]]];
		rk[k2] = rk[k];
		rk[k] = bin[vals[k]++];
		pos[rk[k]] = k;
		pos[rk[k2]] = k2;

	}

	void dec(int k) {
		if (!inside(k)) return;
		assert(vals[k] > 0);
		if (bin[vals[k]-1] < pl) bin[vals[k]-1] = pl;
		int k2 = pos[bin[--vals[k]]];
		rk[k2] = rk[k];
		rk[k] = bin[vals[k]]++;
		pos[rk[k]] = k;
		pos[rk[k2]] = k2;
	}

	bool empty() {
		return pl > pr;
	}

	bool inside(int k) {
		return k < capKey && vis[k] && rk[k] >= pl && rk[k] <= pr;
	}

	int size() {
		return pr - pl + 1;
	}

	int value(int k) {
		assert(k < capKey);
		return vals[k];
	}

	int rank(int k) {
		return rk[k];
	}

	int max(int i=0) {
		return pos[pr-i];
	}

	int min(int i=0) {
		return pos[pl+i];
	}

	void popMin() {
		if (pl <= pr) ++pl;
	}

	void popMax() {
		if (pl <= pr) --pr;
	}

};

#endif
