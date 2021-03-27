#include "main_lib.h"
#include <stdio.h>
#include <stdlib.h>

void print_in_file (point** p, PAYLOAD_POINT* m_pt, int M, int n, int lft_pt)
{
	if(n % 100 == 0){
	int m;
	FILE* command;
	FILE* payload;
	FILE* end_mass;
	FILE* plot;
	char name[100];
	char fileNameLetter[] = "T_";
		
	sprintf(name, "%s%d.txt", fileNameLetter, n);
    plot = fopen(name, "w");
	
	payload = fopen("payload.txt", "a");
	end_mass = fopen("end_mass.txt", "a");
    printf("level: %d, energy =  %e\n", n, energy_calc(p, m_pt,rho_S_L_over_M, rho_S_L_over_m, M));
	command = fopen("command.txt", "a");
	
	//fprintf(test, "set xrange[-0.05:1.05];\n");
    //fprintf(test, "set yrange[-0.05:0.05];\n");
	fprintf(payload, "%d %lf %lf %lf %lf %lf %e\n", n, m_pt->s, m_pt->s_t, m_pt->x, p[lft_pt]->e, p[lft_pt]->y, energy_calc(p, m_pt,rho_S_L_over_M, rho_S_L_over_m, M));
	fprintf(end_mass, "%d %lf %lf %lf %lf %lf %lf %lf\n", n, p[M-1]->x, p[M-1]->y, p[M-1]->x_s, p[M-1]->y_s, p[M-1]->x_t, p[M-1]->y_t, p[M-1]->e);
	
	fprintf(command, "plot \"%s\" u 1:2 with lines, \"%s\" u 1:3 with lines, \"%s\" u 1:4 with lines, \"%s\" u 1:5 with lines, \"%s\" u 1:6 with lines, \"%s\" u 1:7 with lines, \"%s\" u 1:8 with lines\n", name, name, name, name, name, name, name);
	fprintf(command, "pause 0.001;\n");
	
	for(m = 0;  m < M; m++){
		fprintf(plot, "%lf %lf %lf %lf %lf %lf %lf %lf\n", p[m]->coord, p[m]->x, p[m]->y, p[m]->x_s, p[m]->y_s, p[m]->x_t, p[m]->y_t, p[m]->e);
	}
		
	fclose(plot);
	fclose(payload);
	fclose(end_mass);
	fclose(command);
	}
	return;
}

void printVector(double* vec, int N){
   int i;
   for (i = 0; i < N; ++i) {
       printf("\n");
       printf ("%lf ", vec[i]);
   }
printf("\n");
return;
}

void printMatrix(double** matrix, int Height, int Width){
	int i, j;
	for (i = 0; i < Height; i++){
		for(j = 0; j < Width; j++){
			printf("%lf ", matrix[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void print_struct(point* p){
	printf("printing struct: s = %lf, x = %lf, y = %lf, x_s = %lf, y_s = %lf, x_t = %lf, y_t = %lf, lambda = %lf, e = %lf\n",p->coord, p->x, p->y, p->x_s, p->y_s, p->x_t, p->y_t, p->lambda, p->e);
	return;
}