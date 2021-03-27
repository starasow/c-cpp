#include "main_lib.h"
#include "math.h"


//g++ -O2 -openmp main.c set_init_cond.c calc.c gauss.c print.c -o calc

void init_payload (PAYLOAD_POINT* m_pt, point** p, double h, int lft_pt, int rt_pt)
{
	m_pt->x_t = 0.0;
	m_pt->y_t = 0.0;
	m_pt->s_t = v0; //0.2; //initial velocity along the tether;	//0.045 is the speed of transversal wave
	m_pt->s = (p[lft_pt]->coord + p[rt_pt]->coord)/2.0;		
	m_pt->x = (p[lft_pt]->x + p[rt_pt]->x)/2.0;
	m_pt->y = (p[lft_pt]->y + p[rt_pt]->y)/2.0;
	p[M-1]->y_t = 0.0; 
	p[M-1]->x_t = 0.0;
	return;		
}

 void init_triangle(register int m, num_triangle* tr, point** p, point** next_p, point* temp_next)
{
	//printf("initializing triangle scheme m = %d\n", m);
	tr->A1 = p[m-1];
	tr->P = p[m];
	tr->A2 = p[m+1];
	tr->next_P = next_p[m];
	//printf("initializing triangle scheme m = %d, success\n", m);
	return;
}

void swap (register int m, point** p, point** next_p, point* temp_for_swap)
{
	temp_for_swap = p[m];
	p[m] = next_p[m];
	next_p[m] = temp_for_swap;
}

double set_init_x_cond (register int m){
	double h = 0.00000001;
	return h*sin(M_PI*m/(M - 1.0)); //1.0*h*m/(M - 1.0);
}

double set_init_y_cond (register int m){
	double h = 0.00000001;//-0.25;
	return h*sin(2.0*M_PI*m/(M - 1.0));
	/*double s_0 = 0.35; //0.35
    double sum = 0;
    register int k;
    for(k = 1; k < M; k++) {
        sum = sum + h*sin(M_PI*k*s_0)*sin(M_PI*k*(1.0*m/(M - 1.0)))/M_PI/M_PI/k/k/s_0/(1.0 - s_0);
    }
	return sum;*/ 
}

double set_init_x_t_cond (register int m){
	return 0.0;
}

double set_init_y_t_cond (register int m){
	return 0.0;
}

//such scheme is employed: y'(m) = (y[m+1] - y[m-1])/2/h; 
double set_init_x_s_cond (register int m, double h){
	return (set_init_x_cond(m + 1) - set_init_x_cond(m - 1))/h/2.0;
}

double set_init_y_s_cond (register int m, double h){
	return (set_init_y_cond(m + 1) - set_init_y_cond(m - 1))/h/2.0;;
}

double set_SIN_Phi(point* p){
	return (p->y_s)/sqrt((p->x_s + 1.0 + e0)*(p->x_s + 1.0 + e0) + (p->y_s)*(p->y_s));
}
 
double set_COS_Phi(point* p){  
	return (p->x_s + 1.0 + e0)/sqrt((p->x_s + 1.0 + e0)*(p->x_s + 1.0 + e0) + (p->y_s)*(p->y_s));
}

double set_e(point* p)
{
	double e = sqrt((p->x_s + 1.0 + e0)*(p->x_s + 1.0 + e0) + (p->y_s)*(p->y_s)) - 1.0;
	return ((e) > (0.0) ? (e) : (0.0)); 
}
double set_lambda(point* p, double E_over_rho)
{
	return sqrt(E_over_rho*(p->e)/(1.0 + p->e));
}
