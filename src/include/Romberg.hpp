#pragma once

#include <map>

#include <Communication.hpp>

class Romberg {
public:
	Romberg(Communication &_comm, const uint64_t sweeps, double* _table);

	void run(uint64_t i);

private:
	Communication &comm;
	const uint64_t sweeps;
	double* table;
	
};