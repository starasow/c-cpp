all:
	g++ -O2 -openmp main.c set_init_cond.c calc.c gauss.c print.c -o calc