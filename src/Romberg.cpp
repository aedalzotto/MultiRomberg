#include <cmath>

#include <Romberg.hpp>

Romberg::Romberg(
	Communication &_comm, 
	const uint64_t _sweeps,
	double* _table
) :
	comm(_comm),
	sweeps(_sweeps),
	table(_table)
{

}

void Romberg::run(uint64_t i)
{
	/* Calculate the weight of this depth */
	double weight = pow(4, i - sweeps + 1);

	/* Receive the first CTR/Romberg sent to this depth */
	double r0 = comm.receive(i);

	for(uint64_t r = i - sweeps + 1; r < sweeps; r++){
		/* Receive another CTR/Romberg sent to this depth */
		double r1 = comm.receive(i);

		/* Apply the Romberg */
		double romberg = (weight * r1 - r0) / (weight - 1);

		/* Send to the next depth */
		if(i < sweeps * 2 - 2)
			comm.send(romberg, i + 1);

		/* Prepare for next Sweep */
		r0 = r1;

		/* Save result */
		table[r*sweeps + (i - sweeps)] = romberg;
	}
}