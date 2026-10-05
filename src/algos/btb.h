#ifndef BTB_H
#define BTB_H

#include "mdb.h"

class BTB : public MDB {
public:
	static SolverResult run(const std::string &dataPath, int q[2], int k, int flags, int numThreads);
private:
	void branch(int dep) override;
};

#endif
