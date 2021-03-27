#include "main_lib.h"
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include "gauss.h" 

double get_module(double* vec, int length){
    int i; double res = 0.0; 
    for(i = 0; i < length; ++i){res += vec[i]*vec[i];}
    return sqrt(res);
}

double linearB2(double s_B2, double lftPoint, double func_lftpnt, double func_rtpnt, double h){
  ///linear interpolation between two points from left to right
  return func_lftpnt + (s_B2 - lftPoint)*(func_rtpnt - func_lftpnt)/h;
}  

double linearB1(double s_B1, double rtPoint, double func_lftpnt, double func_rtpnt, double h){
  ///linear interpolation between two points from right to left
  return func_rtpnt - (rtPoint - s_B1)*(func_rtpnt - func_lftpnt)/h;
}

//this functions are duplicating funcs from set_init_cond.c, carrying the same sense, the get different arguments
double e(double x_s, double y_s){
	return sqrt((x_s + 1.0 + e0)*(x_s + 1.0 + e0) + y_s*y_s) - 1.0;
}

double set_e1(double nxt_s, point** p, point** next_p, PAYLOAD_POINT* m_pt, double h, int lft_pt, int rt_pt, double c)
{ 
	double next_lft_x_s, next_lft_y_s, e1;
	
	next_lft_x_s = (c - nxt_s - next_p[lft_pt-1]->x)/(nxt_s - next_p[lft_pt-1]->coord); //(c - nxt_s - next_p[lft_pt]->x)/h + (h - nxt_s + next_p[lft_pt]->coord)*(next_p[lft_pt]->x - next_p[lft_pt-1]->x)/h/h;
	next_lft_y_s = (m_pt->next_y - next_p[lft_pt-1]->y)/(nxt_s - next_p[lft_pt-1]->coord); // + (h - nxt_s + next_p[lft_pt]->coord)*(next_p[lft_pt]->y 
	e1 = sqrt((1.0 + e0 + next_lft_x_s)*(1 + e0 + next_lft_x_s) + (next_lft_y_s)*(next_lft_y_s)) - 1.0; 
	//printf("e1 = %lf, nxt_s = %lf, next_lft_x_s = %lf, next_lft_y_s = %lf\n", e1, nxt_s, next_lft_x_s, next_lft_y_s);
	return (e1) > (0.0) ? (e1) : (0.0);
}

double set_e2(double nxt_s, point** p, point** next_p, PAYLOAD_POINT* m_pt, double h, int lft_pt, int rt_pt, double c){
	double e2, next_rt_x_s, next_rt_y_s;
	
	next_rt_x_s = (next_p[rt_pt + 1]->x - c + nxt_s)/(next_p[rt_pt + 1]->coord - nxt_s); //(next_p[rt_pt]->x - c + nxt_s)/h + (h + nxt_s - next_p[rt_pt]->coord)*(next_p[rt_pt+1]->x - next_p[rt_pt]->x)/h/h;
	next_rt_y_s = (next_p[rt_pt + 1]->y - m_pt->next_y)/(next_p[rt_pt + 1]->coord - nxt_s); //(next_p[rt_pt]->y - m_pt->next_y)/h + (h + nxt_s - next_p[rt_pt]->coord)*(next_p[rt_pt+1]->y - next_p[rt_pt]->y)/h/h;
	e2 = sqrt((1.0 + e0 + next_rt_x_s)*(1.0 + e0 + next_rt_x_s) + (next_rt_y_s)*(next_rt_y_s)) - 1.0;
	return ((e2) > (0.0) ? (e2) : (0.0));
}

double COS_Phi (double x_s, double y_s){
	return (x_s + 1.0 + e0)/sqrt((x_s + 1.0 + e0)*(x_s + 1.0 + e0) + y_s*y_s); 
}

double  SIN_Phi (double x_s, double y_s){
	return y_s/sqrt((x_s + 1.0 + e0)*(x_s + 1.0 + e0) + y_s*y_s);
}

void calculation(double** Matr, point** p, num_triangle* tr, double h, double tau, register int m, double a, double* solution, double E_over_rho)
{	
	double s_B1, s_B2, lambda_B1, lambda_B2, x_s_B1, x_s_B2, x_t_B1, x_t_B2, y_s_B1, y_s_B2, y_t_B1, y_t_B2, x_B1, x_B2, y_B1, y_B2;
	
	//printf("started calculations...mesh_size = %lf, tau = %lf m = %d\n", h, tau, m);	
	///A1 B1 P B2 A2
	//linearisation;
	s_B1 = tr->P->coord - (tr->P->lambda)*tau*h/(h + tau*(tr->P->lambda - tr->A1->lambda)); 
	s_B2 = tr->P->coord + (tr->P->lambda)*tau*h/(h - tau*(tr->A2->lambda - tr->P->lambda)); 
	lambda_B1 = linearB1(s_B1, tr->P->coord, tr->A1->lambda, tr->P->lambda, h);		//s_A[m], lambda_cur[m-1], lambda_cur[m], h);
	lambda_B2 = linearB2(s_B2, tr->P->coord, tr->P->lambda, tr->A2->lambda, h);		//s_A[m], lambda_cur[m], lambda_cur[m+1], h);  
    
	x_s_B1 = linearB1(s_B1, tr->P->coord, tr->A1->x_s, tr->P->x_s, h);      //s_A[m], x_s_cur_A[m-1], x_s_cur_A[m], h);
	x_s_B2 = linearB2(s_B2, tr->P->coord, tr->P->x_s,  tr->A2->x_s, h);     //s_A[m], x_s_cur_A[m], x_s_cur_A[m+1], h);
	y_s_B1 = linearB1(s_B1, tr->P->coord, tr->A1->y_s, tr->P->y_s, h); 		//s_A[m], y_s_cur_A[m-1], y_s_cur_A[m], h);
	y_s_B2 = linearB2(s_B2, tr->P->coord, tr->P->y_s,  tr->A2->y_s, h);     //s_A[m], y_s_cur_A[m], y_s_cur_A[m+1], h);

	x_t_B1 = linearB1(s_B1, tr->P->coord, tr->A1->x_t, tr->P->x_t, h); 		//s_A[m], x_t_cur_A[m-1], x_t_cur_A[m], h);
	x_t_B2 = linearB2(s_B2, tr->P->coord, tr->P->x_t,  tr->A2->x_t, h); 	//s_A[m], x_t_cur_A[m], x_t_cur_A[m+1], h);
	y_t_B1 = linearB1(s_B1, tr->P->coord, tr->A1->y_t, tr->P->y_t, h); 	    //s_A[m], y_t_cur_A[m-1], y_t_cur_A[m], h);
	y_t_B2 = linearB2(s_B2, tr->P->coord, tr->P->y_t,  tr->A2->y_t, h);		//s_A[m], y_t_cur_A[m], y_t_cur_A[m+1], h);
	
	x_B1 = linearB1(s_B1, tr->P->coord, tr->A1->x, tr->P->x, h);
	x_B2 = linearB2(s_B2, tr->P->coord, tr->P->x,  tr->A2->x, h); 
	y_B1 = linearB1(s_B1, tr->P->coord, tr->A1->y, tr->P->y, h);
	y_B2 = linearB2(s_B2, tr->P->coord, tr->P->y,  tr->A2->y, h);	
		
	///printf("Start filling in the matrix...\n");
	Matr[0][0] = tr->A1->COS_Phi;						///A1
	Matr[0][1] = tr->A1->SIN_Phi;          				///A1		
	Matr[0][2] = -a*(tr->A1->COS_Phi);  				///A1
	Matr[0][3] = -a*(tr->A1->SIN_Phi);  				///A1
	Matr[0][4] = (tr->A1->COS_Phi)*(tr->A1->x_t) + (tr->A1->SIN_Phi)*(tr->A1->y_t) - a*(tr->A1->COS_Phi)*(tr->A1->x_s) - a*(tr->A1->SIN_Phi)*tr->A1->y_s  + tau*((tr->A1->COS_Phi)*(-2.0*dimless_omega*(tr->A1->y_t) + 3.0*dimless_omega*dimless_omega*(tr->A1->x + tr->A1->coord)) + (tr->A1->SIN_Phi)*(2.0*dimless_omega*(tr->A1->x_t) - dimless_omega*dimless_omega*(tr->A1->y)));
	
	Matr[1][0] = tr->A2->COS_Phi;         				///A2
	Matr[1][1] = tr->A2->SIN_Phi;         				///A2
	Matr[1][2] = a*tr->A2->COS_Phi;						///A2
	Matr[1][3] = a*tr->A2->SIN_Phi; 					///A2
	Matr[1][4] = tr->A2->COS_Phi*tr->A2->x_t + tr->A2->SIN_Phi*tr->A2->y_t + tr->A2->COS_Phi*tr->A2->x_s + tr->A2->SIN_Phi*tr->A2->y_s + tau*((tr->A2->COS_Phi)*(-2.0*dimless_omega*tr->A2->y_t + 3.0*dimless_omega*dimless_omega*(tr->A2->x + tr->A2->coord)) + (tr->A2->SIN_Phi)*(2.0*dimless_omega*(tr->A2->x_t) - dimless_omega*dimless_omega*(tr->A2->y)));

	Matr[2][0] = -SIN_Phi(x_s_B1, y_s_B1);				///B1 
	Matr[2][1] = COS_Phi(x_s_B1, y_s_B1); 				///B1		
	Matr[2][2] = lambda_B1*SIN_Phi(x_s_B1, y_s_B1);   	///B1
	Matr[2][3] = -lambda_B1*COS_Phi(x_s_B1, y_s_B1); 	///B1		
	Matr[2][4] = -SIN_Phi(x_s_B1, y_s_B1)*x_t_B1 + COS_Phi(x_s_B1, y_s_B1)*y_t_B1 + lambda_B1*SIN_Phi(x_s_B1, y_s_B1)*x_s_B1 - lambda_B1*COS_Phi(x_s_B1, y_s_B1)*y_s_B1 + tau*(COS_Phi(x_s_B1, y_s_B1)*(2.0*dimless_omega*x_t_B1  - dimless_omega*dimless_omega*y_B1)) - SIN_Phi(x_s_B1, y_s_B1)*(-2.0*dimless_omega*y_t_B1 + 3.0*dimless_omega*dimless_omega*(s_B1 + x_B1)); 

	Matr[3][0] = -SIN_Phi(x_s_B2, y_s_B2); 				///B2		
	Matr[3][1] = COS_Phi(x_s_B2, y_s_B2);				///B2
	Matr[3][2] = -lambda_B2*SIN_Phi(x_s_B2, y_s_B2);	///B2	
	Matr[3][3] = lambda_B2*COS_Phi(x_s_B2, y_s_B2);		///B2		
	Matr[3][4] = -SIN_Phi(x_s_B2, y_s_B2)*x_t_B2 + COS_Phi(x_s_B2, y_s_B2)*y_t_B2 - lambda_B2*SIN_Phi(x_s_B2, y_s_B2)*x_s_B2 + lambda_B2*COS_Phi(x_s_B2, y_s_B2)*y_s_B2 + tau*(COS_Phi(x_s_B2, y_s_B2)*(2.0*dimless_omega*x_t_B2 - dimless_omega*dimless_omega*y_B2) - SIN_Phi(x_s_B2, y_s_B2)*(-2.0*dimless_omega*y_t_B2 + 3.0*dimless_omega*dimless_omega*(s_B2 + x_B2))); 

	GaussElimination(Matr, 4);
	SolveforDiagonal(solution, Matr, 4);
	
	tr->next_P->x_t = solution[0]; //x_t_nxt[m] = solution[0]; u
	tr->next_P->y_t = solution[1]; //y_t_nxt[m] = solution[1]; nu 
	tr->next_P->x_s = solution[2]; //x_s_nxt[m] = solution[2]; mu
	tr->next_P->y_s = solution[3]; //y_s_nxt[m] = solution[3]; chi	
	tr->next_P->e = set_e(tr->next_P);
	tr->next_P->lambda = set_lambda(tr->next_P, E_over_rho);
	tr->next_P->x = tr->P->x + (tr->P->x_t)*tau;
	tr->next_P->y = tr->P->y + (tr->P->y_t)*tau;
	tr->next_P->coord = tr->P->coord;
	tr->next_P->COS_Phi = set_COS_Phi(tr->next_P);
	tr->next_P->SIN_Phi = set_SIN_Phi(tr->next_P);
	//print_struct(tr->next_P);
	return;	
}

void left_boundary_calc(double** MatrBnd, point** p, point** next_p, double tau, double h, double a, double* solution_boundary, double E_over_rho)
{
	   ///left end
	double s_B2, lambda_B2, x_s_B2, y_s_B2, x_t_B2, y_t_B2;
	
	///check for positive deformation, if not positive no forces applied, x_t, y_t - same, x_s, y_s recalculate through x,  y
	if(p[0]->e > eps && p[1]->e > eps)
	{		
		s_B2 = p[0]->coord + (p[0]->lambda)*tau*h/(h - tau*(p[1]->lambda - p[0]->lambda));
		lambda_B2 = linearB2(s_B2, p[0]->coord, p[0]->lambda, p[1]->lambda, h);
		x_s_B2 = linearB2(s_B2, p[0]->coord, p[0]->x_s, p[1]->x_s, h);
		y_s_B2 = linearB2(s_B2, p[0]->coord, p[0]->y_s, p[1]->y_s, h);
		x_t_B2 = linearB2(s_B2, p[0]->coord, p[0]->x_t, p[1]->x_t, h);
		y_t_B2 = linearB2(s_B2, p[0]->coord, p[0]->y_t, p[1]->y_t, h);
	
		MatrBnd[0][0] = a*p[1]->COS_Phi; 	//COS_Phi(p[1]->x_s, p[1]->y_s);    ///A2
		MatrBnd[0][1] = a*p[1]->SIN_Phi;	//SIN_Phi(p[1]->x_s, p[1]->y_s);    ///A2
		MatrBnd[0][2] = (p[1]->COS_Phi)*(p[1]->x_t + a*(p[1]->x_s)) + (p[1]->SIN_Phi)*(p[1]->y_t + a*p[1]->y_s) + tau*((p[1]->COS_Phi)*(-2.0*dimless_omega*(p[1]->y_t) + 3.0*dimless_omega*dimless_omega*(p[1]->x + p[1]->coord)) + (p[1]->SIN_Phi)*(2.0*dimless_omega*(p[1]->x_t))); 

		MatrBnd[1][0] = -lambda_B2*SIN_Phi(x_s_B2, y_s_B2);  ///B2
		MatrBnd[1][1] = lambda_B2*COS_Phi(x_s_B2, y_s_B2);   ///B2
		MatrBnd[1][2] = -SIN_Phi(x_s_B2, y_s_B2)*(x_t_B2 + lambda_B2*x_s_B2) + COS_Phi(x_s_B2, y_s_B2)*(y_t_B2 + lambda_B2*y_s_B2) + tau*(COS_Phi(x_s_B2, y_s_B2)*(2.0*dimless_omega*x_t_B2) - SIN_Phi(x_s_B2, y_s_B2)*(-2.0*dimless_omega*y_t_B2 + 3.0*dimless_omega*dimless_omega*s_B2)); 

		GaussElimination(MatrBnd, 2);
		SolveforDiagonal(solution_boundary, MatrBnd, 2);

		next_p[0]->x_s = solution_boundary[0];
		next_p[0]->y_s = solution_boundary[1];  
		next_p[0]->e = set_e(next_p[0]);
		next_p[0]->lambda = set_lambda(next_p[0], E_over_rho);
		next_p[0]->coord = 0.0;
		next_p[0]->SIN_Phi = set_SIN_Phi(next_p[0]);
		next_p[0]->COS_Phi = set_COS_Phi(next_p[0]);
	} else {
		printf("e dropped to 0 at left end\n"); 
		next_p[0]->x_s = (-(next_p[2]->x) + 4*(next_p[1]->x) - 3*(next_p[0]->x))/2.0/h;
		next_p[0]->y_s = (-(next_p[2]->y) + 4*(next_p[1]->y) - 3*(next_p[0]->y))/2.0/h;
		next_p[0]->e = set_e(next_p[0]);
		next_p[0]->lambda = set_lambda(next_p[0], E_over_rho);
		next_p[0]->coord = 0.0;
		next_p[0]->SIN_Phi = set_SIN_Phi(next_p[0]);
		next_p[0]->COS_Phi = set_COS_Phi(next_p[0]);
	}
	
	next_p[0]->x_t = 0.0;
	next_p[0]->y_t = 0.0;
	next_p[0]->x = p[0]->x + tau*(p[0]->x_t);
	next_p[0]->y = p[0]->y + tau*(p[0]->y_t);
	//free(solution_boundary); */
	return;
}

void right_boundary_calc(double** MatrBnd, point** p, point** next_p, double tau, double h, double a, double* solution_boundary, int M, double E_over_rho, double rho_S_L_over_M)
{
	///right end
	double s_B1, lambda_B1, x_s_B1, y_s_B1, x_t_B1, y_t_B1, f_x, f_y;
	//double rho_S_L_over_M = 0.015; //0.015 = 100kg;
	//g++ -O2 -openmp main.c set_init_cond.c calc.c gauss.c print.c -o calc
	
	f_x = -2.0*dimless_omega*(p[M-1]->y_t) + 3.0*dimless_omega*dimless_omega*(1.0 + p[M-1]->x);
	f_y = 2.0*dimless_omega*(p[M-1]->x_t) - dimless_omega*dimless_omega*(p[M-1]->y);
	printf("coriolis = %lf, 3omega^2(1 + x_L) = %lf, f_y = %lf\n", -2.0*dimless_omega*(p[M-1]->y_t), 3.0*dimless_omega*dimless_omega*(1.0 + p[M-1]->x), f_y);
	
	next_p[M-1]->x_t = p[M-1]->x_t + tau*(-rho_S_L_over_M*(p[M-1]->e)*(p[M-1]->COS_Phi) + f_x);
	printf("x_t_M = %lf, -rhoSLecos(Phi)/M = %lf, f_x = %lf\n", next_p[M-1]->x_t, -rho_S_L_over_M*(p[M-1]->e)*(p[M-1]->COS_Phi), f_x);
	
	next_p[M-1]->y_t = p[M-1]->y_t + tau*(-rho_S_L_over_M*(p[M-1]->e)*(p[M-1]->SIN_Phi) + f_y);
	next_p[M-1]->x = p[M-1]->x + tau*(next_p[M-1]->x_t);
	next_p[M-1]->y = p[M-1]->y + tau*(next_p[M-1]->y_t);  
	
	s_B1 = p[M-1]->coord - (p[M-1]->lambda)*tau*h/(h + tau*(p[M-1]->lambda - p[M-2]->lambda));
	lambda_B1 = linearB1(s_B1, p[M-1]->coord, p[M-2]->lambda, p[M-1]->lambda, h);		
	x_s_B1 = linearB1(s_B1, p[M-1]->coord, p[M-2]->x_s, p[M-1]->x_s, h);
	y_s_B1 = linearB1(s_B1, p[M-1]->coord, p[M-2]->y_s, p[M-1]->y_s, h);
	x_t_B1 = linearB1(s_B1, p[M-1]->coord, p[M-2]->x_t, p[M-1]->x_t, h);
	y_t_B1 = linearB1(s_B1, p[M-1]->coord, p[M-2]->y_t, p[M-1]->y_t, h);
	
	MatrBnd[0][0] = -a*p[M-2]->COS_Phi;
	MatrBnd[0][1] = -a*p[M-2]->SIN_Phi;
	MatrBnd[0][2] = (p[M-2]->COS_Phi)*(p[M-2]->x_t - next_p[M-1]->x_t - a*p[M-2]->x_s + f_x*tau) + (p[M-2]->SIN_Phi)*(-next_p[M-1]->y_t - a*p[M-2]->y_s + f_y*tau + p[M-2]->y_t);
	
	MatrBnd[1][0] = lambda_B1*SIN_Phi(x_s_B1, y_s_B1);
	MatrBnd[1][1] = -lambda_B1*COS_Phi(x_s_B1, y_s_B1);
	MatrBnd[1][2] = SIN_Phi(x_s_B1, y_s_B1)*(-x_t_B1 + lambda_B1*x_s_B1 - f_x*tau + next_p[M-1]->x_t) + COS_Phi(x_s_B1, y_s_B1)*(y_t_B1 - lambda_B1*y_s_B1 + f_y*tau - next_p[M-1]->y_t);

	GaussElimination(MatrBnd, 2);
	SolveforDiagonal(solution_boundary, MatrBnd, 2);
	printf("\n"); 
	 
	next_p[M-1]->x_s = solution_boundary[0]; //+ 0.5*(p[M-1]->x - p[M-2]->x)/(p[M-1]->coord - p[M-2]->coord); // + p[M-1]->x - p[M-1]->x); 	//solution_boundary[0]; 
	next_p[M-1]->y_s = solution_boundary[1]; //+ 0.5*(p[M-1]->y - p[M-2]->y)/(p[M-1]->coord - p[M-2]->coord); // + p[M-1]->x - p[M-1]->x);//-solution_boundary[1];
	next_p[M-1]->e = set_e(next_p[M-1]); 
	next_p[M-1]->lambda = set_lambda(next_p[M-1], E_over_rho);
	next_p[M-1]->coord = 1.0;
	next_p[M-1]->SIN_Phi = set_SIN_Phi(next_p[M-1]);
	next_p[M-1]->COS_Phi = set_COS_Phi(next_p[M-1]);

	return;
}
 
 void payload_calc (point** p, point** next_p, PAYLOAD_POINT* m_pt, double h, double tau, int lft_pt, int rt_pt, double E_over_rho, double rho_S_L_over_m) 
{
	//double rho_S_L_over_m = 1.5; //0.15 = 10kg; comes from dimensionless calculations
	double c, k, e; 
	double next_lft_x_s, next_lft_y_s, next_rt_x_s, next_rt_y_s;
	//g++ -O2 -openmp main.c set_init_cond.c calc.c gauss.c print.c -o calc
	
	//x and y on t+1 level
	next_p[lft_pt]->x = p[lft_pt]->x + (p[lft_pt]->x_t)*tau;
	next_p[rt_pt]->x = p[rt_pt]->x + (p[rt_pt]->x_t)*tau;
	next_p[lft_pt]->y = p[lft_pt]->y + (p[lft_pt]->y_t)*tau;
	next_p[rt_pt]->y = p[rt_pt]->y + (p[rt_pt]->y_t)*tau;		
	next_p[lft_pt]->coord = p[lft_pt]->coord;
	next_p[rt_pt]->coord = p[rt_pt]->coord; 
	
	//printf("e1 - e2 = %lf\n", (p[lft_pt]->e) - (p[rt_pt]->e));
	e = 0.5*(p[lft_pt]->e + p[rt_pt]->e);
	k = rho_S_L_over_m*(e/(e + 1.0) + (m_pt->s_t)*(m_pt->s_t));
	
	c = m_pt->s + m_pt->x + (m_pt->s_t + m_pt->x_t)*tau + (k*(p[rt_pt]->x_s - p[lft_pt]->x_s) - 2.0*dimless_omega*(m_pt->y_t) + 3.0*dimless_omega*dimless_omega*(m_pt->s + m_pt->x))*tau*tau;
	
	m_pt->next_y = m_pt->y + (m_pt->y_t)*tau + tau*tau*(k*(p[rt_pt]->y_s - p[lft_pt]->y_s) + 2.0*dimless_omega*(m_pt->x_t) - dimless_omega*dimless_omega*(m_pt->y));
	
	///calculate next position of the mass - s coordinate
	//m_pt->next_s = calculate_s(p, next_p, m_pt, h, lft_pt, rt_pt, c, tau);
	m_pt->next_s = newtons_method(p, next_p, m_pt, h, lft_pt, rt_pt, c, m_pt->s);
	///IMPORTANT///
	
	///now that we have s, we calculate rest of the stuff for t+1
	m_pt->next_s_t = (m_pt->next_s - m_pt->s)/tau;
	m_pt->next_x = c - m_pt->next_s;
	m_pt->next_y_t = (m_pt->next_y - m_pt->y)/tau; 
	m_pt->next_x_t = (m_pt->next_x - m_pt->x)/tau;
	
 	 
	//these equations are on the same level tau + 1!!!
	next_lft_x_s = (m_pt->next_x - next_p[lft_pt-1]->x)/(m_pt->next_s - next_p[lft_pt-1]->coord); 
	next_lft_y_s = (m_pt->next_y - next_p[lft_pt-1]->y)/(m_pt->next_s - next_p[lft_pt-1]->coord); 
	next_rt_x_s = (next_p[rt_pt+1]->x - m_pt->next_x)/(next_p[rt_pt + 1]->coord - m_pt->next_s); 
	next_rt_y_s = (next_p[rt_pt+1]->y - m_pt->next_y)/(next_p[rt_pt + 1]->coord - m_pt->next_s); 
	
	//find x_t, y_t on both sides and rest
	next_p[lft_pt]->x_t = m_pt->next_x_t - (m_pt->next_s_t)*next_lft_x_s;
	next_p[lft_pt]->y_t = m_pt->next_y_t - (m_pt->next_s_t)*next_lft_y_s;
	next_p[lft_pt]->x_s = next_lft_x_s;
	next_p[lft_pt]->y_s = next_lft_y_s;
	next_p[lft_pt]->e = set_e(next_p[lft_pt]);
	next_p[lft_pt]->lambda = set_lambda(next_p[lft_pt], E_over_rho);
	next_p[lft_pt]->COS_Phi = set_COS_Phi(next_p[lft_pt]);
	next_p[lft_pt]->SIN_Phi = set_SIN_Phi(next_p[lft_pt]);

	next_p[rt_pt]->x_t = m_pt->next_x_t - (m_pt->next_s_t)*next_rt_x_s;
	next_p[rt_pt]->y_t = m_pt->next_y_t - (m_pt->next_s_t)*next_rt_y_s;
	next_p[rt_pt]->x_s = next_rt_x_s;
	next_p[rt_pt]->y_s = next_rt_y_s;
	next_p[rt_pt]->e = set_e(next_p[rt_pt]);
	next_p[rt_pt]->lambda = set_lambda(next_p[rt_pt], E_over_rho);
	next_p[rt_pt]->COS_Phi = set_COS_Phi(next_p[rt_pt]);
	next_p[rt_pt]->SIN_Phi = set_SIN_Phi(next_p[rt_pt]); 	
	return;
} 
 
double newtons_method(point** p, point** next_p, PAYLOAD_POINT* m_pt, double h, int lft_pt, int rt_pt, double c, double s)
{
	///s - is the pre corrected value, for first iteration (m_pt->s) is taken;
	double nxt_s, best_nxt_s;	
	double f, df_ds;
	double c1, c2, e1, e2;
	// 
	
	e1 = set_e1(s, p, next_p, m_pt, h, lft_pt, rt_pt, c);
	e2 = set_e2(s, p, next_p, m_pt, h, lft_pt, rt_pt, c);
	c1 = next_p[lft_pt-1]->coord - c + next_p[lft_pt-1]->x - m_pt->next_y + next_p[lft_pt-1]->y;
	c2 = next_p[rt_pt+1]->coord - c + next_p[rt_pt+1]->x - m_pt->next_y + next_p[rt_pt+1]->y; 
	//printf("difference = %e, error = %e\n", e1 - e2, abs(e1 - e2)/e1);
	
	///Newton's method
	///f(x_n):
	f = e1 - e2;
	
	///f'(x_n):
	df_ds = ((s - next_p[lft_pt-1]->coord)*(s - next_p[lft_pt-1]->coord)*(1.0 + e0) + c1)/(s - next_p[lft_pt-1]->coord)/(s - next_p[lft_pt-1]->coord)/(e1 + 1.0) - ((next_p[rt_pt+1]->coord - s)*(next_p[rt_pt+1]->coord - s)*(1.0 + e0) + c2)/(next_p[rt_pt+1]->coord - s)/(next_p[rt_pt+1]->coord - s)/(e2 + 1.0);
	
	///x_n+1 = x_n - f(x_n)/f'(x_n);
	nxt_s = s - f/df_ds;
	
	if(abs(nxt_s - s) > eps){
		nxt_s = newtons_method(p, next_p, m_pt, h, lft_pt, rt_pt, c, nxt_s);
	}
	//printf("difference = %e, error = %e\n", e1 - e2, abs(e1 - e2)/e1);
	return nxt_s;
}

double energy_calc (point** p, PAYLOAD_POINT* m_pt, double rho_S_L_over_M, double rho_S_L_over_m, int M){
	double h = 1.0 / M;
	double tether_energy = 0.0;
	double M_end_over_rho_S_L = 1.0 / rho_S_L_over_M;
	double M_payload_over_rho_S_L = 1.0 / rho_S_L_over_m;
	double payload_energy, endmass_energy, payload_energy_kinetic;
	double payload_energy_potential; // = M_payload_over_rho_S_L*3.0*dimless_omega*dimless_omega;
	register int m;
	#pragma omp parallel for
	for (m = 0; m < M; m++){
		tether_energy += 0.5*(p[m]->x_t)*(p[m]->x_t)*h + 0.5*(p[m]->y_t)*(p[m]->y_t)*h + 0.5*(p[m]->e)*(p[m]->e)*h;
	}
	
	payload_energy_kinetic = 0.5*M_payload_over_rho_S_L*((m_pt->x_t + m_pt->s_t)*(m_pt->x_t + m_pt->s_t) + (m_pt->y_t)*(m_pt->y_t));
	payload_energy_potential = 0.5*M_payload_over_rho_S_L*3.0*dimless_omega*dimless_omega*(1.0 - (m_pt->s + m_pt->x)*(m_pt->s + m_pt->x));
	
	payload_energy = payload_energy_kinetic + payload_energy_potential;
	endmass_energy = 0.5*M_end_over_rho_S_L*((p[M-1]->x_t)*(p[M-1]->x_t) + (p[M-1]->y_t)*(p[M-1]->y_t) - 3.0*dimless_omega*dimless_omega*(2.0*(p[M-1]->x) + (p[M-1]->x)*(p[M-1]->x)));
	
	return payload_energy + tether_energy + endmass_energy;
}