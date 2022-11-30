#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define EPS 1e-10
#define tests_num 5

double f(double x);
double solution(double x);
double Abs(double x);
double f_vect(double x, double y_1, double y_2, double y_3, int num);
double find_ddu0(double* x, int n, double y0, double y1, double dy0);
void RungeKutta(double* x, double* res_y, int n, double y0, double dy0, double ddy0);

double f(double x)
{
    double res;
    res = 100 * (10 * x * sin(10 * x) - 3 * cos(10 * x));
    res = res + cos(x) * (cos(10 * x) - 10 * x * sin(10 * x));
    res = res + x * cos(10 * x);
    return (1 / cos(10)) * res;
}

double solution(double x)
{
    return (1 / cos(10)) * x * cos(10 * x);
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
            return f(x) - y_1 - cos(x) * y_2;
            break;
    }
}

double find_ddu0(double *x, int n, double y0, double y1, double dy0)
{
    double f_cent, f_right;
    double right_diff;
    double result = 0;
    double *u;

    u = (double *)malloc(n * sizeof(double));
    
    RungeKutta(x, u, n, y0, dy0, result);
    f_cent = u[n - 1];
    RungeKutta(x, u, n, y0, dy0, result + 1);
    f_right = u[n - 1];
    right_diff = (f_right - f_cent);
    result += (y1 - f_cent) / right_diff;  

    free(u);
    return result;
}

void RungeKutta(double* x, double* res_y, int n, double y0, double dy0, double ddy0)
{
    double h = 1. / (n - 1);
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
    {
        res_y[i] = y[i][0];
        free(y[i]);
    }
    free(y);

}

int main()
{
    int n, k;
    double h;
    double u0, u1, du0, ddu0;
    double *x, *y;
    double error[tests_num];

    printf("Enter u(0): ");
    scanf("%lf", &u0);
    printf("Enter u(1): ");
    scanf("%lf", &u1);
    printf("Enter u'(0): ");
    scanf("%lf", &du0);
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
        x[0] = 0;
        for (int i = 1; i < n; i++)
            x[i] = x[i - 1] + h;
        ddu0 = find_ddu0(x, n, u0, u1, du0);
        RungeKutta(x, y, n, u0, du0, ddu0);

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