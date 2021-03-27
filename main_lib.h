#pragma once

///GLOBAL VARIABLES
extern double v0;
extern int M;              		/// -number of nodes in grid
extern int N;		   	  		/// -number of MAX iterations over time (although payload should reach the other end earlier)
extern double e0;        		/// -initial deformation from thw small satellite
extern double eps;  		    /// -criterion for newton's method |s_n+1 - s_n| < eps
extern double dimless_omega;	/// -dimension less angular velocity of larger sat (check out the paper)
extern double rho_S_L_over_M;  	/// -mass of small satellite in dimensionless units
extern double rho_S_L_over_m;	/// -mass of payload  in dimensionless units

typedef struct point {
	double coord;
	double x;
	double y;
	double x_s;
	double y_s;
	double x_t;
	double y_t;
	double e;
	double lambda;
	double COS_Phi;
	double SIN_Phi;
}*p[];
 
typedef struct num_triangle {
	point* A1;
	point* P;
	point* A2;
	point* next_P;
}*tr;

typedef struct PAYLOAD_POINT {
	double s;
	double s_t;
	double x_t;
	double y_t;
	double x;
	double y;
	double next_s_t;
	double next_x_t;
	double next_y_t;
	double next_x;
	double next_y;
	double next_s;
}*m_pt;

//initiators, used only once at the start of the program;
double set_init_y_cond (register int m);
double set_init_x_cond (register int m);
double set_init_y_s_cond (register int m, double h);
double set_init_x_s_cond (register int m, double h);  
double set_init_y_t_cond (register int m);
double set_init_x_t_cond (register int m);

//used multiple times to set all the parameters, calculated from y_s and x_s
double set_COS_Phi(point* p);
double set_SIN_Phi(point* p);
double set_e(point* p);
double set_lambda(point* p, double E_over_rho);

void init_payload (PAYLOAD_POINT* m_pt, point** p, double h, int lft_pt, int rt_pt);
void payload_calc (point** p, point** next_p, PAYLOAD_POINT* m_pt, double h, double tau, int lft_pt, int rt_pt, double E_over_rho,  double rho_S_L_over_m);
void right_boundary_calc(double** MatrBnd, point** p, point** next_p, double tau, double h, double a, double* solution_boundary, int M, double E_over_rho, double rho_S_L_over_M);
void left_boundary_calc(double** MatrBnd, point** p, point** next_p, double tau, double h, double a, double* solution_boundary, double E_over_rho);
void calculation(double** Matr, point** p, num_triangle* tr, double h, double tau, register int m, double a, double* solution, double E_over_rho);
void init_triangle(register int m, num_triangle* tr, point** p, point** next_p, point* temp_next);
void swap (register int m, point** p, point** next_p, point* temp_for_swap);

double linearB2(double s_B2, double lftPoint, double func_lftpnt, double func_rtpnt, double h);
double linearB1(double s_B1, double rtPoint, double func_lftpnt, double func_rtpnt, double h);
double COS_Phi (double x_s, double y_s);
double SIN_Phi (double x_s, double y_s);
double get_module(double* vec, int length);
double e (double x_s, double y_s);
double set_e1(double nxt_s, point** p, point** next_p, PAYLOAD_POINT* m_pt, double h, int lft_pt, int rt_pt, double c);
double set_e2(double nxt_s, point** p, point** next_p, PAYLOAD_POINT* m_pt, double h, int lft_pt, int rt_pt, double c);
double calculate_s(point** p, point** next_p, PAYLOAD_POINT* m_pt, double h, int lft_pt, int rt_pt, double c, double tau);
double newtons_method(point** p, point** next_p, PAYLOAD_POINT* m_pt, double h, int lft_pt, int rt_pt, double c, double s);
double energy_calc (point** p, PAYLOAD_POINT* m_pt,  double rho_S_L_over_M, double rho_S_L_over_m, int M);

void print_in_file (point** p, PAYLOAD_POINT* m_pt, int M, int n, int lft_pt);
void printVector(double* vec, int N);
void printMatrix(double** matrix, int Height, int Width);
void print_struct(point* p);
