#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gauss.h"
#include "main_lib.h"

#define eps 0.00000000001

double BackStep(int Row, int N, double** A, double* x)
{
  register int j;
  double Sum = 0.0;
  for(j = 0; j < Row; j++)
  {
    if(j > 0)
	{
		Sum = Sum + A[N - Row][N-j]*x[N-j];
    }
  }
  x[N - Row] = (A[N - Row][N] - Sum)/A[N-Row][N-Row];
  Sum = 0.0;
  return x[N - Row];
}

void SetZeroRow(int Row, int N, double** A)
{
register int i, j;
double buf = fabs(A[0][0]);
  for(i = Row +1; i < N; i++)
   {
    buf = A[i][Row];
    for(j = Row ; j < N + 1; j++)
        {
             if(fabs(A[Row][Row]) > pow(10, -20))
            {
                A[i][j] = A[i][j] - A[Row][j]*(buf)/(A[Row][Row]);
            }
            else
			 printf("Error. Can't divide by zero\n");
			
        }
   }
}

double** GaussElimination (double** A, int N)
{
    register int j = 0;
     for(j = 0; j < N-1; j++)
     {
        SetZeroRow(j, N, A);
     }
     return A;
}

double* SolveforDiagonal(double* x, double** A, int N)
{
    register int i;
    if(get_module(A[N-1], N+1) > eps) 
        {
        for(i = 1; i < N+1; i++)
            {
                BackStep(i, N, A, x);
            }
        }
    else
    {
        printf("system is linearly dependent\n");
        printMatrix(A, N,  N+1);
        printVector(x, N);
        printVector(A[N-1], N+1);
        printf("\n");
	}
    return x;
}