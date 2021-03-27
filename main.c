#include "main_lib.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>


int main(void)
{
	double time_spent;
	
	//each point in the grid is a structure, m is its' position in the array from 0 to M-1
	//which represents its' Lagrangian coordinate on the dimensionless tether: m/M-1 -- between 0 and 1;
	clock_t begin = clock();	
	FILE* command;
	FILE* payload; 
	FILE* end_mass;
	
	point** p;
	point** next_p;
	PAYLOAD_POINT* m_pt; 
	point* temp_next;
	point* temp_for_swap;
	num_triangle* tr;
	
	double** Matr;
	double** MatrBnd; 
	double* solution;
	double* solution_boundary;
	register int m, n;
	int lft_pt = 3; //0.35
	int rt_pt = lft_pt + 1;
	//double a = 1.0;
	double E_over_rho = 1.0;
	double a;
	double tau;
	double h = 1.0 / M;
	
	a = sqrt(E_over_rho);
	tau = h / a;
	
	p = (point**)malloc(M*sizeof(point*));
		for(m = 0; m < M; m++)
			p[m] = (point*)malloc(M*sizeof(point)); 
	
	next_p = (point**)malloc(M*sizeof(point*));
		for(m = 0; m < M; m++)
			next_p[m] = (point*)malloc(M*sizeof(point));	
	
	temp_for_swap = (point*)malloc(sizeof(point));
	m_pt = (PAYLOAD_POINT*)malloc(sizeof(PAYLOAD_POINT));
	temp_next = (point*)malloc(sizeof(point));
	 
	tr = (num_triangle*)malloc(sizeof(num_triangle));
	solution = (double*)malloc(4*sizeof(double));
	solution_boundary = (double*)malloc(2*sizeof(double));
	
	//setting initial deformations of the tether x and y
	for(m = 0; m < M; m++){
		p[m]->coord = (double)m/(M-1.0);
		p[m]->x = set_init_x_cond(m);
		p[m]->y = set_init_y_cond(m);
		p[m]->x_t = set_init_x_t_cond (m);
		p[m]->y_t = set_init_y_t_cond (m);
	}	
	
	//setting initial derivatives of deformations x_t, y_t, x_s, y_s; 
	for(m = 1; m < M-1; m++){
		p[m]->x_s = set_init_x_s_cond (m, h);
		p[m]->y_s = set_init_y_s_cond (m, h);
		//init_derivatives(next_p[m], m);
	}
	//boundary derivatives
	p[0]->x_s = (-(p[2]->x) + 4*(p[1]->x) - 3*(p[0]->x))/2.0/h;
	p[0]->y_s = (-(p[2]->y) + 4*(p[1]->y) - 3*(p[0]->y))/2.0/h;
	p[M-1]->x_s = (3*(p[M-1]->x) - 4*(p[M-2]->x) + p[M-3]->x)/2.0/h;
	p[M-1]->y_s = (3*(p[M-1]->y) - 4*(p[M-2]->y) + p[M-3]->y)/2.0/h;
	
	//setting  e, lambda, SIN, COS
	for(m = 0; m < M; m++){
		p[m]->e = set_e(p[m]);
		p[m]->lambda = set_lambda(p[m], E_over_rho);
		p[m]->SIN_Phi = set_SIN_Phi(p[m]);
		p[m]->COS_Phi = set_COS_Phi(p[m]);
	}
	
	Matr = (double**)malloc(4*sizeof(double*));
		for(m = 0; m < 4; m++) ///lines
			Matr[m] = (double*)malloc(5*sizeof(double));

	MatrBnd = (double**)malloc(2*sizeof(double*));
		for(m = 0; m < 2; m++) ///lines
			MatrBnd[m] = (double*)malloc(3*sizeof(double));
	
	for(m = 0; m < M; m++)
		print_struct(p[m]);
	
	command = fopen("command.txt", "w");
	fprintf(command, "set xrange[-0.05:1.05];\n");
    fprintf(command, "set yrange[-0.5:0.5];\n");
	fclose(command);
	
	payload = fopen("payload.txt", "w");
	end_mass = fopen("end_mass.txt", "w");
	
	init_payload (m_pt, p, h, lft_pt, rt_pt);
	///main cycle over time
	for(n = 0; n < N; n++)
	{	
			printf("payload is between %d and %d points\n", lft_pt, rt_pt);			
			#pragma omp parallel for
			for (m = 1; m < lft_pt; m++){
					init_triangle(m, tr, p, next_p, temp_next);
					calculation(Matr, p, tr, h, tau, m, a, solution, E_over_rho);
				}
				
			#pragma omp parallel for
			for (m = rt_pt + 1; m < M - 1; m++){
					init_triangle(m, tr, p, next_p, temp_next);
					calculation(Matr, p, tr, h, tau, m, a, solution, E_over_rho);
			}
		 
			left_boundary_calc(MatrBnd, p, next_p, tau, h, a, solution_boundary, E_over_rho);
			right_boundary_calc(MatrBnd, p, next_p, tau, h, a, solution_boundary, M, E_over_rho, rho_S_L_over_M);
			payload_calc (p, next_p, m_pt, h, tau, lft_pt, rt_pt, E_over_rho, rho_S_L_over_m);
			
			m_pt->s_t = m_pt->next_s_t;
			m_pt->x_t = m_pt->next_x_t;
			m_pt->y_t = m_pt->next_y_t;
			m_pt->x = m_pt->next_x;
			m_pt->y = m_pt->next_y;
			m_pt->s = m_pt->next_s;
			
			//swapping levels  
			for(m = 0; m < M; m++){
				swap (m, p, next_p, temp_for_swap);
			}
			
			print_in_file (p, m_pt, M, n, lft_pt);

			if(m_pt->next_s > p[rt_pt]->coord){
				//m_pt->s = p[rt_pt]->coord;
				lft_pt++;
				rt_pt++;
			}
			
			if(m_pt->next_s < p[lft_pt]->coord){
				//m_pt->s = p[lft_pt]->coord;
				lft_pt--;
				rt_pt--;
			}
			
			if(lft_pt == 2 || rt_pt == M - 3){
				printf("reached one end\n");
				break;
			} 
			//print_in_file (p, M, n);
	} 
	
	fclose(payload);
	fclose(end_mass);
	
	for (m = 0; m < M-1; m++){
		free(p[m]);
		free(next_p[m]);
	}
	
	free(temp_for_swap);
	free(temp_next);
	free(m_pt);
	free(tr); 
	free(solution_boundary);
	free(solution);
	
	clock_t end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("time spent = %lf sec\n", time_spent);
	return 0;
}