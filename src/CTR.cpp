#include <cmath>

#include <CTR.hpp>

CTR::CTR(
	Communication &_comm, 
	double (*_f)(double),
	const double _a,
	const double _b,
	const uint64_t _sweeps,
	double* _table
) :
	comm(_comm),
	f(_f),
	a(_a),
	b(_b),
	sweeps(_sweeps),
	table(_table)
{

}

void CTR::run(uint64_t i)
{
	/* Recalculate this each thread, is faster than waiting for a 'master' */
	double interval = b - a;		/* Integration interval */
	double limit_sum = f(a) + f(b);	/* Sum of the integration of the limits */

	uint64_t n = 1 << i;							/* Number of trapezoids = pow(2, i-1) */
	double h = interval / static_cast<double>(n);	/* Size of trapezoids */

	/* Apply the CTR (only the odd ones) */
	double trapezoids_odd = 0;
	for(uint64_t k = 1; k < n; k += 2)
		trapezoids_odd += f(a + k*h);


	/* Receive the even trapezoids */
	double trapezoids_even = i ? comm.receive(i) : 0;

	double trapezoids = trapezoids_even + trapezoids_odd;

	/* Send this thread's trapezoids to the next */
	if(i < sweeps - 1)
		comm.send(trapezoids, i + 1);

	/* Calculate the CTR */
	double ctr = h/2.0 * (limit_sum + 2.0 * trapezoids);

	/* Send this thread's CTR to the first depth of Romberg */
	if(sweeps > 1)
		comm.send(ctr, sweeps);

	/* Save result */
	table[i] = ctr;

}