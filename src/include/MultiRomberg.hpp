#pragma once

#include <thread>
#include <map>

#include <Communication.hpp>

class MultiRomberg {
public:
	MultiRomberg(
		double (*f)(double), 
		const double a, 
		const double b, 
		const uint64_t _sweeps
	);

	void wait();
	double get_ans();
	double* get_table();

private:
	Communication comm;
	uint64_t sweeps;

	double *table;
	std::vector<std::thread> threads;
};