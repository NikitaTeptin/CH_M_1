#ifndef GAUSS
#define GAUSS

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>

#define REAL // смена точности
#define REAL_SUM
#define EPS 1e-13

#ifdef REAL
typedef float real;
#define REALIN "%f"
#define REALOUT "%f\n"
#else
typedef double real;
#define REALIN "%lf"
#define REALOUT "%.15f\n"
#endif

#ifdef REAL_SUM
typedef float realsum;
#else
typedef double realsum;
#endif

using namespace std;

void LLT(real* DI, real* AL, int* IA, int& n, int& n_profile);
void CalcX(real* AL, real* DI, int* IA, real* Y, real* X, int& n);
void CalcY(real* AL, real* DI, int* IA, real* F, real* Y, int& n);
void Input_size(const char* s, int& n);
real* GAUSS_METHOD(vector<vector<real>> A, real* y, int& n);
void Input_arr(real*& DI, int*& IA, real*& AL, int& n, int& n_profile);
void Input_vector(real*& F, int& n);
vector<vector<real>> Profile_To_Dense(real* DI, real* AL, int* IA, real* F, int& n);
void Output_Vector(real* x, int& n);
void Hilbert_Matrix(real*& DI, real*& AL, int*& IA, real*& F, int n, int& n_profile);
void different(real* x, int n);
#endif