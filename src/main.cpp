#include <chrono>
#include <algorithm>
#include <csignal>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include "utils/cmdline.hpp"
#include "utils/log.hpp"
#include "algos/btb.h"

using namespace logging;

namespace
{
	constexpr const char *RESULT_FILE = "result.txt";
	char timeoutLine[1024];
	size_t timeoutLineLen = 0;

	bool fileExists(const std::string &path)
	{
		struct stat buffer;
		return stat(path.c_str(), &buffer) == 0 && S_ISREG(buffer.st_mode);
	}

	std::string baseName(const std::string &path)
	{
		size_t end = path.find_last_not_of('/');
		if (end == std::string::npos)
			return path;
		size_t start = path.find_last_of('/', end);
		return path.substr(start == std::string::npos ? 0 : start + 1, end - (start == std::string::npos ? 0 : start + 1) + 1);
	}

	std::string dirName(const std::string &path)
	{
		size_t end = path.find_last_not_of('/');
		if (end == std::string::npos)
			return "";
		size_t slash = path.find_last_of('/', end);
		if (slash == std::string::npos)
			return "";
		return path.substr(0, slash);
	}

	std::string instanceNameFromInput(const std::string &inputPath)
	{
		if (baseName(inputPath) == "edges.txt")
			return baseName(dirName(inputPath));
		return baseName(inputPath);
	}

	std::string resolveDataPath(const std::string &inputPath)
	{
		if (fileExists(inputPath))
			return inputPath;
		return inputPath + "/edges.txt";
	}

	void readGraphHeader(const std::string &dataPath, int &numEdges, int &numLeft, int &numRight)
	{
		std::ifstream in(dataPath);
		if (!(in >> numEdges >> numLeft >> numRight))
		{
			log("Failed to read graph header from %s", dataPath.c_str());
			exit(1);
		}
	}

	void writeLine(const std::string &line)
	{
		std::string text = line + '\n';
		int fd = open(RESULT_FILE, O_WRONLY | O_CREAT | O_APPEND, 0644);
		if (fd >= 0)
		{
			(void)write(fd, text.c_str(), text.size());
			close(fd);
		}
	}

	void timeoutHandler(int)
	{
		int fd = open(RESULT_FILE, O_WRONLY | O_CREAT | O_APPEND, 0644);
		if (fd >= 0)
		{
			(void)write(fd, timeoutLine, timeoutLineLen);
			close(fd);
		}
		_exit(124);
	}

	void installTimeoutHandlers(const std::string &instanceName, int numEdges, int numLeft, int numRight, int k, int q)
	{
		std::ostringstream line;
		line << instanceName << '\t'
			 << numEdges << '\t'
			 << numLeft << '\t'
			 << numRight << '\t'
			 << k << '\t'
			 << q << '\t'
			 << 0;
		std::string text = line.str() + '\n';
		timeoutLineLen = std::min(text.size(), sizeof(timeoutLine) - 1);
		std::memcpy(timeoutLine, text.c_str(), timeoutLineLen);
		timeoutLine[timeoutLineLen] = '\0';

		std::signal(SIGTERM, timeoutHandler);
		std::signal(SIGINT, timeoutHandler);
		std::signal(SIGALRM, timeoutHandler);
	}
}

int main(int argc, char *argv[])
{
	cmdline::parser args;
	args.add<std::string>("data", 'd', "dataset path", true, "");
	args.add<int>("key", 'k', "value of k", true, 0);
	args.add<int>("lb", 'q', "lower bound size", false, 0);
	args.add<int>("jobs", 'j', "number of parallel jobs", false, 0);
	args.add("no-ub", '\0', "disable upper bound techniques");
	args.add("no-core", '\0', "disable core reduction");
	args.add("no-cn", '\0', "disable common neighbor reduction");
	args.add("no-1nn", '\0', "disable one non-neighbor reduction");
	args.add("no-pb", '\0', "disable progressive bounding");
	args.add("no-order", '\0', "disable ordering reduction");
	args.add("no-heu", '\0', "disable heuristic algorithm");
	args.add("no-br", '\0', "disable branching rules");
	args.add("debug", '\0', "output intermediate data");
	args.parse_check(argc, argv);

	auto inputDataPath = args.get<std::string>("data");
	auto dataPath = resolveDataPath(inputDataPath);
	auto instanceName = instanceNameFromInput(dataPath);
	int k = args.get<int>("key");
	int qValue = args.get<int>("lb");
	int lb[2] = {qValue, qValue};
	int numEdges = 0, numLeft = 0, numRight = 0;
	readGraphHeader(dataPath, numEdges, numLeft, numRight);
	installTimeoutHandlers(instanceName, numEdges, numLeft, numRight, k, qValue);

	int flags = !args.exist("no-order") * FLAG_ORDER | !args.exist("no-core") * FLAG_CORE | !args.exist("no-cn") * FLAG_CN | !args.exist("no-1nn") * FLAG_1NN | !args.exist("no-pb") * FLAG_PB | !args.exist("no-heu") * FLAG_HEU | !args.exist("no-br") * FLAG_BR | !args.exist("no-ub") * FLAG_UB | args.exist("debug") * FLAG_DEBUG;

	int numThreads = std::max(args.get<int>("jobs"), 1);

	auto startTimePoint = std::chrono::steady_clock::now();
	SolverResult result = BTB::run(dataPath, lb, k, flags, numThreads);

	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
		std::chrono::steady_clock::now() - startTimePoint);
	double durationSeconds = duration.count() / 1000.0;

	log("Total time spent: %lld ms", (long long)duration.count());

	std::ostringstream line;
	line << std::fixed << std::setprecision(3);
	line << instanceName << '\t'
		 << numEdges << '\t'
		 << numLeft << '\t'
		 << numRight << '\t'
		 << k << '\t'
		 << qValue << '\t'
		 << 1 << '\t'
		 << result.optimalEdges << '\t'
		 << result.numBranches << '\t'
		 << result.numPruned << '\t'
		 << durationSeconds;
	writeLine(line.str());
}
