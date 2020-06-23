#include <cmath>

#include <MultiRomberg.hpp>
#include <Romberg.hpp>
#include <CTR.hpp>

#include <iostream>

MultiRomberg::MultiRomberg(
	double (*f)(double), 
	double a, 
	double b, 
	uint64_t _sweeps
) : 
	comm(_sweeps),
	sweeps(_sweeps),
	table(nullptr)
{
	/* Create the table */
	table = (double*)malloc(sizeof(double) * sweeps * sweeps);

	/* Create the threads */
	threads.reserve(sweeps * 2 - 1);

	for(uint64_t i = 0; i < sweeps; i++)
		threads.push_back(std::thread(&CTR::run, CTR(comm, f, a, b, sweeps, table), i));

	for(uint64_t i = sweeps; i < sweeps * 2 - 1; i++)
		threads.push_back(std::thread(&Romberg::run, Romberg(comm, sweeps, table), i));
}

void MultiRomberg::wait()
{
	for(auto& thread : threads)
		thread.join();
}

double MultiRomberg::get_ans()
{
	return table[sweeps*sweeps - (sweeps > 1 ? 2 : 1)];
}

double* MultiRomberg::get_table()
{
	return table;
}