#include <Communication.hpp>

#include <iostream>

Communication::Communication(uint64_t sweeps)
{
	/* Reserve and create communication queue, synchronization mechanism and message buffers for all processes */
	comm.reserve(sweeps * 2 - 1);
	sync.reserve(sweeps * 2 - 1);

	for(uint64_t i = 0; i < sweeps * 2 - 1; i++){
		comm.push_back(std::queue<double>());

		sem_t s;
		sem_init(&s, 0, 0);
		sync.push_back(s);
	}
}

void Communication::send(double value, uint64_t channel)
{
	comm[channel].push(value);

	/* @todo change to 'sync[channel].acquire()' when C++20 implements it */
	sem_post(&(sync[channel]));
}

double Communication::receive(uint64_t channel)
{
	/* @todo change to 'sync[channel].acquire()' when C++20 implements it */
	sem_wait(&sync[channel]);

	double value = comm[channel].front();
	comm[channel].pop();
	return value;
}