#include "main_lib.h"

double v0 = 0.002;
int M  = 100;              		/// -number of nodes in grid
int N = 1000000;		   	  	/// -number of MAX iteration over time (although payload should reach the other end earlier)
double e0 = 0.0024;        		/// -initial deformation from thw small satellite
double eps = 0.000000001;  		/// -criterion for newton's method |s_n+1 - s_n| < eps
double dimless_omega = 0.0035;	/// -dimension less angular velocity of larger sat (check out the paper)
double rho_S_L_over_M = 0.015;  ///  0.015 = 100kg; mass of small satellite in dimensionless units
double rho_S_L_over_m = 1.5;    ///  0.15 = 10kg; mass of payload  in dimensionless units

//g++ -O2 -openmp main.c set_init_cond.c calc.c gauss.c print.c -o calc