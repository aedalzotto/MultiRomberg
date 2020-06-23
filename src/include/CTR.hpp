#pragma once

#include <map>

#include <Communication.hpp>

class CTR {
public:
	CTR(Communication &_comm, double (*_f)(double), const double _a, const double _b, const uint64_t _sweeps, double* _table);

	void run(uint64_t i);

private:
	Communication &comm;
	double (*f)(double);
	const double a;
	const double b;
	const uint64_t sweeps;
	double* table;
};