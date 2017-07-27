#pragma once
#include <random>
#include <iostream>
#include <chrono>

namespace {
	class Sampler{
	public:
		static std::mt19937 mt;
		static std::default_random_engine rand_gen;

		static double gamma(double a, double b){
			std::gamma_distribution<double> distribution(a, 1.0 / b);
			return distribution(rand_gen);
		}

		static double beta(double a, double b){
			double ga = gamma(a, 1.0);
			double gb = gamma(b, 1.0);
			return ga / (ga + gb);
		}

		static double bernoulli(double p){
			std::uniform_real_distribution<double> rand(0, 1);
			double r = rand(mt);
			if(r > p){
				return 0;
			}
			return 1;
		}
		static double uniform(double min = 0, double max = 0){
			std::uniform_real_distribution<double> rand(min, max);
			return rand(mt);
		}
		static double uniform_int(int min = 0, int max = 0){
			std::uniform_int_distribution<> rand(min, max);
			return rand(mt);
		}
		static double normal(double mean = 0, double stddev = 1){
			std::normal_distribution<double> rand(mean, stddev);
			return rand(mt);
		}
	};

	int seed = std::chrono::system_clock::now().time_since_epoch().count();
	// int seed = 0;
	std::mt19937 Sampler::mt(seed);
	std::default_random_engine Sampler::rand_gen(seed);
}
