#pragma once

#include <vector>
#include <queue>
#include <sstream>

#include <semaphore.h>

class Communication {
public:
	Communication(uint64_t sweeps);
	void send(double value, uint64_t channel);
	double receive(uint64_t channel);

private:
	std::vector<std::queue<double>> comm;
	std::vector<sem_t> sync;

};