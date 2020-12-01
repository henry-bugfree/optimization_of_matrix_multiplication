/*
 * author: Henry Liang
 * date: 11/30/2020
 */

#include <iostream>
#include <cstring>
#include <ctime>
using namespace std;

#define MATRIX_SIZE 1024
#define BLOCK_SIZE 32

#define NAIVE_MUL
#define ADVANCED_MUL
#define ANOTHER_ADVANCED_MUL
#define BLOCKED_MUL
#define ADVANCED_BLOCKED_MUL

double matrix_a[MATRIX_SIZE*MATRIX_SIZE];
double matrix_b[MATRIX_SIZE*MATRIX_SIZE];
double matrix_c[MATRIX_SIZE*MATRIX_SIZE];

int initial_matrix(double* matrix);
int naive_mul(double* matrix_des, const double* matrix_1, const double* matrix_2);
int advanced_mul(double* matrix_des, const double* matrix_1, const double* matrix_2);
int another_advanced_mul(double* matrix_des, const double* matrix_1, const double* matrix_2);
int cal_one_block(double* matrix_des,const double* matrix_1,const double* matrix_2);
int advanced_cal_one_block(double* matrix_des,const double* matrix_1,const double* matrix_2);
int blocked_mul(double* matrix_des, double* matrix_1, double* matrix_2);
int advanced_blocked_mul(double* matrix_des, double* matrix_1, double* matrix_2);

int main() {
    clock_t start_time, end_time;

    initial_matrix(matrix_a);
    initial_matrix(matrix_b);
#ifdef NAIVE_MUL
    //naive_mul
    memset(matrix_c,0,MATRIX_SIZE*MATRIX_SIZE*sizeof(double));

    start_time = clock();
    naive_mul(matrix_c, matrix_a, matrix_b);
    end_time = clock();
    cout<<"naive_mul runtime:"<<(double)(end_time-start_time)/CLOCKS_PER_SEC<<"s"<<endl;

    cout<<"random check:"<<matrix_c[5*MATRIX_SIZE+24]<<endl;
#endif
#ifdef ADVANCED_MUL
    //advanced_mul
    memset(matrix_c,0,MATRIX_SIZE*MATRIX_SIZE*sizeof(double));

    start_time = clock();
    advanced_mul(matrix_c, matrix_a, matrix_b);
    end_time = clock();
    cout<<"advanced_mul runtime:"<<(double)(end_time-start_time)/CLOCKS_PER_SEC<<"s"<<endl;

    cout<<"random check:"<<matrix_c[5*MATRIX_SIZE+24]<<endl;
#endif
#ifdef ANOTHER_ADVANCED_MUL
    //another_advanced_mul
    memset(matrix_c,0,MATRIX_SIZE*MATRIX_SIZE*sizeof(double));

    start_time = clock();
    another_advanced_mul(matrix_c, matrix_a, matrix_b);
    end_time = clock();
    cout<<"another_advanced_mul runtime:"<<(double)(end_time-start_time)/CLOCKS_PER_SEC<<"s"<<endl;

    cout<<"random check:"<<matrix_c[5*MATRIX_SIZE+24]<<endl;
#endif
#ifdef BLOCKED_MUL
    //blocked_mul
    memset(matrix_c,0,MATRIX_SIZE*MATRIX_SIZE*sizeof(double));

    start_time = clock();
    blocked_mul(matrix_c, matrix_a, matrix_b);
    end_time = clock();
    cout<<"blocked_mul runtime:"<<(double)(end_time-start_time)/CLOCKS_PER_SEC<<"s"<<endl;

    cout<<"random check:"<<matrix_c[5*MATRIX_SIZE+24]<<endl;
#endif
#ifdef ADVANCED_BLOCKED_MUL
    //advanced_blocked_mul
    memset(matrix_c,0,MATRIX_SIZE*MATRIX_SIZE*sizeof(double));

    start_time = clock();
    advanced_blocked_mul(matrix_c, matrix_a, matrix_b);
    end_time = clock();
    cout<<"advanced_blocked_mul runtime:"<<(double)(end_time-start_time)/CLOCKS_PER_SEC<<"s"<<endl;

    cout<<"random check:"<<matrix_c[5*MATRIX_SIZE+24]<<endl;
#endif
    return 0;
}

int initial_matrix(double* matrix)
{
    for(int i=0;i<MATRIX_SIZE;i++)
        for(int j=0;j<MATRIX_SIZE;j++)
            matrix[i*MATRIX_SIZE+j] = 2.0;

    return 0;
}

int naive_mul(double* matrix_des, const double* matrix_1, const double* matrix_2)
{
    for(int i=0;i<MATRIX_SIZE;i++)
        for(int j=0;j<MATRIX_SIZE;j++)
            for(int k=0;k<MATRIX_SIZE;k++)
                matrix_des[i*MATRIX_SIZE+j] += matrix_1[i*MATRIX_SIZE+k]*matrix_2[k*MATRIX_SIZE+j];

    return 0;
}

int advanced_mul(double* matrix_des, const double* matrix_1, const double* matrix_2)
{
    for(int k=0;k<MATRIX_SIZE;k++)
        for(int i=0;i<MATRIX_SIZE;i++)
            for(int j=0;j<MATRIX_SIZE;j++)
                matrix_des[i*MATRIX_SIZE+j] += matrix_1[i*MATRIX_SIZE+k]*matrix_2[k*MATRIX_SIZE+j];

     return 0;
}

int another_advanced_mul(double* matrix_des, const double* matrix_1, const double* matrix_2)
{
    for(int k=0;k<MATRIX_SIZE;k++)
        for(int i=0;i<MATRIX_SIZE;i++)
        {
            double temp = matrix_1[i*MATRIX_SIZE+k];
            for(int j=0;j<MATRIX_SIZE;j++)
                matrix_des[i*MATRIX_SIZE+j] += temp*matrix_2[k*MATRIX_SIZE+j];
        }

    return 0;
}

int cal_one_block(double* matrix_des,const double* matrix_1,const double* matrix_2)
{
    for(int i=0;i<BLOCK_SIZE;i++)
        for(int j=0;j<BLOCK_SIZE;j++)
            for(int k=0;k<BLOCK_SIZE;k++)
                matrix_des[i*MATRIX_SIZE+j] += matrix_1[i*MATRIX_SIZE+k]*matrix_2[k*MATRIX_SIZE+j];

    return 0;
}

int advanced_cal_one_block(double* matrix_des,const double* matrix_1,const double* matrix_2)
{
    for(int k=0;k<BLOCK_SIZE;k++)
        for(int i=0;i<BLOCK_SIZE;i++)
            for(int j=0;j<BLOCK_SIZE;j++)
                matrix_des[i*MATRIX_SIZE+j] += matrix_1[i*MATRIX_SIZE+k]*matrix_2[k*MATRIX_SIZE+j];

    return 0;
}

int blocked_mul(double* matrix_des, double* matrix_1, double* matrix_2)
{
    for(int i=0;i<MATRIX_SIZE;i+=BLOCK_SIZE)
        for(int j=0;j<MATRIX_SIZE;j+=BLOCK_SIZE)
            for(int k=0;k<MATRIX_SIZE;k+=BLOCK_SIZE)
                cal_one_block(matrix_des+i*MATRIX_SIZE+j,
                              matrix_1+i*MATRIX_SIZE+k,matrix_2+k*MATRIX_SIZE+j);

    return 0;
}

int advanced_blocked_mul(double* matrix_des, double* matrix_1, double* matrix_2)
{
    for(int i=0;i<MATRIX_SIZE;i+=BLOCK_SIZE)
        for(int j=0;j<MATRIX_SIZE;j+=BLOCK_SIZE)
            for(int k=0;k<MATRIX_SIZE;k+=BLOCK_SIZE)
                advanced_cal_one_block(matrix_des+i*MATRIX_SIZE+j,
                              matrix_1+i*MATRIX_SIZE+k,matrix_2+k*MATRIX_SIZE+j);

    return 0;
}

