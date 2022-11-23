#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define EPS 1e-4
#define tests_num 5
#define y0 0
#define y1 0
#define dy0 30
#define ddy0 0.5

double f(double x);
double solution(double x);
double Abs(double x);
double f_vect(double x, double y_1, double y_2, double y_3, int num);
double find_ddy0(double* x, double n);
double* RungeKutta(double* x, double* res_y, double n);

double f(double x)
{
    double res;
    res = - 27000 * cos(30 * x);
    res = res + 30 * cos(x) * cos(30 * x);
    res = res + sin(30 * x);
    return res;
}

double solution(double x)
{
    return sin(30 * x);
}

double Abs(double x)
{
    return (x > 0) ? x : -x;
}

double f_vect(double x, double y_1, double y_2, double y_3, int num)
{
    switch(num){
        case 0:
            return y_2; 
            break;
        case 1:
            return y_3;
            break;
        case 2:
            return f(x) - y1 - cos(x) * y_2;
            break;
    }
}

double find_ddy0(double *x, double n)
{
    double up_border = -100;
    double down_border = 100;
    double res = 0;
    double *u;

    u = (double *)malloc(n * sizeof(double));

    while(Abs(u[1] - y1) > EPS)
    {
        RungeKutta(x, u, n); //убрать define'ы
    }


    free(u);
    return res;
}

double* RungeKutta(double* x, double* res_y, double n)
{
    double h = 1 / (n - 1);
    double k1[3];
    double k2[3];
    double k3[3];
    double k4[3];
    double** y;

    y = (double **)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++)
        y[i] = (double *)malloc(3 * sizeof(double));

    y[0][0] = y0;
    y[0][1] = dy0;
    y[0][2] = ddy0;

    for (int i = 0; i < n - 1; i++)
    {
        for (int j = 0; j < 3; j++)
            k1[j] = f_vect(x[i], y[i][0], y[i][1], y[i][2], j);
        for (int j = 0; j < 3; j++)  
            k2[j] = f_vect(x[i] + (h / 2), y[i][0] + (h / 2) * k1[0], y[i][1] + (h / 2) * k1[1], y[i][2] + (h / 2) * k1[2], j);
        for (int j = 0; j < 3; j++)
            k3[j] = f_vect(x[i] + (h / 2), y[i][0] + (h / 2) * k2[0], y[i][1] + (h / 2) * k2[1], y[i][2] + (h / 2) * k2[2], j);
        for (int j = 0; j < 3; j++)
            k4[j] = f_vect(x[i] + h, y[i][0] + h * k3[0], y[i][1] + h * k3[1], y[i][2] + h * k3[2], j);
        for (int j = 0; j < 3; j++)    
            y[i + 1][j] = y[i][j] + (h / 6) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
    }
    for (int i = 0; i < n; i++)
        res_y[i] = y[i][0];
    for(int i = 0; i < n; i++)
        free(y[i]);
    free(y);

}

int main()
{
    int n;
    int k;
    double h;
    double* x;
    double* y;
    double error[tests_num];

    printf("Enter initial number of points: ");
    scanf("%d", &n);
    printf("Initial h = %lf\n", 1. / (n - 1));
    printf("Enter integer coefficient k (h -> h / k): ");
    scanf("%d", &k);

    for (int j = 0; j < tests_num; j++)
    {
        h = 1. / (n - 1);
        x = (double *)malloc(n * sizeof(double));
        y = (double *)malloc(n * sizeof(double));
        for (int i = 0; i < n; i++)
            x[i] = i * h;
        RungeKutta(x, y, n);

        FILE *fp;
        switch(j){
            case 0:
                fp = fopen("test1.txt", "w");
                break;
            case 1:
                fp = fopen("test2.txt", "w");
                break;
            case 2:
                fp = fopen("test3.txt", "w");
                break;
            case 3:
                fp = fopen("test4.txt", "w");
                break;
            case 4:
                fp = fopen("test5.txt", "w");
                break;
        }
        error[j] = 0;
        for(int i = 0; i < n; i++)
        {
            if (error[j] < Abs(y[i] - solution(x[i])))
                error[j] = Abs(y[i] - solution(x[i]));
            fprintf(fp, "%lf %lf\n", x[i], y[i]);
        }
        fclose(fp);
        printf("h = %.8f: error = %le\n", h, error[j]);
        free(x);
        free(y);
        n = (n - 1) * k + 1;
    }

    printf("Ratio error[i] / error[i + 1]:\n");
    for (int i = 0; i < tests_num - 1; i++)
        printf("%.3f\n", error[i] / error[i + 1]);
    
    return 0;
}