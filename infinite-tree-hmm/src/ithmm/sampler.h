#pragma once
#include <random>

namespace ithmm {
	namespace sampler {
		extern std::mt19937 mt;
		double gamma(double a, double b);
		double beta(double a, double b);
		double bernoulli(double p);
		double uniform(double min, double max);
		double uniform_int(int min, int max);
		double normal(double mean, double stddev);
	}
}