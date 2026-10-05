#ifndef VERTEX_SET
#define VERTEX_SET

#include <cstdint>
#include <vector>
#include <algorithm>
#include <cassert>
#include <string>
#include <sstream>

class VertexSet {
	int lp, rp;
	std::vector<int> s, pos;

  void swapByPos(int i, int j) {
		std::swap(s[i], s[j]);
		pos[s[i]] = i;
		pos[s[j]] = j;
	}

	int capacity() const {
		return s.size();
	}

public:
	VertexSet(int n) {
		reserve(n);
		clear();
	}

	VertexSet() {
		clear();
	}

  void swapByVal(int u, int v) {
    swapByPos(pos[u], pos[v]);
  }

	void reserve(int new_capacity) {
		int old_capacity = capacity();
		if (old_capacity >= new_capacity) return;
		s.resize(new_capacity);
		pos.resize(new_capacity);
		for (int i = old_capacity; i < new_capacity; ++i)
			s[i] = pos[i] = i;

	}

	void push(int v) {
		if (v >= capacity()) reserve((v+1) << 1);
		if (pos[v] < lp) pushFront(v);
		else if (pos[v] >= rp) pushBack(v);
	}

	void pushFront(int v) {
		if (v >= capacity()) reserve((v+1) << 1);
		if (pos[v] < lp) swapByPos(pos[v], --lp);
	}

	void pushBack(int v) {
		if (v >= capacity()) reserve((v+1) << 1);
		if (pos[v] >= rp) swapByPos(pos[v], rp++);
	}

	void pop(int v) {
		popFront(v);
	}

	void popFront(int v) {
		if (inside(v)) swapByPos(pos[v], lp++);
	}

	void popBack(int v) {
		if (inside(v)) swapByPos(pos[v], --rp);
	}

	bool inside(int v) const {
		return v < capacity() && pos[v] >= lp && pos[v] < rp;
	}

  int size() const {
		return rp - lp;
	}

	int frontPos() const {
		return lp;
	}

	int backPos() const {
		return rp;
	}

	void restore(int pos) {
		restoreFront(pos);
	}

	void restoreFront(int pos) {
		lp = pos;
	}

	void restoreBack(int pos) {
		rp = pos;
	}

	void clear() {
		lp = rp = 0;
	}

	int operator [] (int index) const {
		return s[index];
	}

	int vertexPos(int v) const {
		return pos[v];
	}

	std::vector<int>::iterator begin() {
		return s.begin() + lp;
	}

	std::vector<int>::const_iterator begin() const {
		return s.begin() + lp;
	}

	std::vector<int>::iterator end() {
		return s.begin() + rp;
	}

	std::vector<int>::const_iterator end() const {
		return s.begin() + rp;
	}

	std::vector<int>::reverse_iterator rbegin() {
		return std::vector<int>::reverse_iterator(s.begin() + rp);
	}

	std::vector<int>::const_reverse_iterator rbegin() const {
		return std::vector<int>::const_reverse_iterator(s.begin() + rp);
	}

	std::vector<int>::reverse_iterator rend() {
		return std::vector<int>::reverse_iterator(s.begin() + lp);
	}

	std::vector<int>::const_reverse_iterator rend() const {
		return std::vector<int>::const_reverse_iterator(s.begin() + lp);
	}

	VertexSet& operator += (VertexSet &V) {
		for (int v : V)
			push(v);
		return *this;
	}

	friend VertexSet operator + (VertexSet a, VertexSet b) {
		int cap_a = a.capacity(), cap_b = b.capacity();
		VertexSet c = cap_a > cap_b ? std::move(a) : std::move(b);
		if (cap_a > cap_b) {
			for (int v : b)
				c.push(v);
		}
		else {
			for (int v : a)
				c.push(v);
		}
		return c;
	}

	std::string toString(const std::string &sep = ",") {
		std::stringstream ss;
		ss << "{";
		for (int i = lp; i < rp; ++i) {
			if (i > lp) ss << sep;
			ss << s[i];
		}
		ss << "}";
		return ss.str();
	}

};

#endif
