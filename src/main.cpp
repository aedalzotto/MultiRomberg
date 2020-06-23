#include <iostream>

#include <cmath>
#include <argp.h>

#include <MultiRomberg.hpp>

double f1(double x)
{
	/* Avoid calculating ln(0) */
	return x ? log(exp(x * log(x)))/log(4) * sin(x) : 0;
}

double f2(double x)
{
	return x * x;
}

const char *argp_program_bug_address = "<angelo.dalzotto@edu.pucrs.br>";

static char args_doc[] = "A B SWEEPS";

static struct argp_option options[] = {
	{"verbose", 'v', 0, 0, "Produce verbose output" },
	{ 0 }
};

struct arguments
{
	bool verbose;
	double a, b;
	uint64_t sweeps;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
	struct arguments *args = static_cast<struct arguments *>((state->input));

	switch (key){
		case 'v':
			args->verbose = true;
			break;
		case ARGP_KEY_ARG:
			switch(state->arg_num){
				case 0:
					args->a = std::stod(arg);
					break;
				case 1:
					args->b = std::stod(arg);
					break;
				case 2:
					args->sweeps = std::stoull(arg);
					if(args->sweeps < 1){
						std::cerr << "É necessário ao menos uma iteração." <<std::endl;
						argp_usage(state);
					}
					break;
				default:
					/* Too many arguments. */
					argp_usage(state);
					break;
			}
			break;
		case ARGP_KEY_END:
			if (state->arg_num < 3)
				/* Not enough arguments. */
				argp_usage(state);

			break;
		default:
			return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

static struct argp argp = {options, parse_opt, args_doc, 0};

int main(int argc, char **argv)
{
	struct arguments args;
	args.verbose = 0;
	argp_parse (&argp, argc, argv, 0, 0, &args);

	std::cout << "Integrating log4(x * ln(x)) * sin(x) with [" << args.a << "," << args.b << "] limits with " << args.sweeps << " sweeps\n\n";
	
	MultiRomberg app(f1, args.a, args.b, args.sweeps);

	app.wait();

	/* Print table */
	if(args.verbose){
		double* table = app.get_table();
		for(uint64_t i = 0; i < args.sweeps; i++){
			std::cout << "R" << i << ",0" << " = " << table[i];
			for(uint64_t j = 1; j <= i; j++)
				std::cout << "\tR" << i << "," << j << "=" << table[i*args.sweeps + j - 1];
			std::cout << std::endl;
		}
	}

	std::cout << "\nAnswer = " << app.get_ans() << std::endl;
	return 0;
}
