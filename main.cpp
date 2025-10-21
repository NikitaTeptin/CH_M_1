#include "GAUSS.h"

int main()
{
    int n = 0;
    int n_profile = 0;
    int* IA = NULL;
    real* DI = NULL, * X = NULL, * Y = NULL, * F = NULL, * AL = NULL;
    setlocale(LC_ALL, "ru");
    int choice = 0;
    cout << "Введите команду: 1 - метод Гаусса\n2 - LLT-разложение\n3 - Генарация матрицы Гильберта" << endl;
    cin >> choice;
    switch ( choice )
    {
    case 1:
    {
        Input_size("SIZE.txt", n);
        Input_arr(DI, IA, AL, n, n_profile);
        Input_vector(F, n);
        vector<vector<real>>A = Profile_To_Dense(DI, AL, IA, F, n);
        for ( int i = 0; i < n; i++ )
        {
            for ( int j = 0; j < n; j++ )
                cout << A[i][j] << ' ';
            cout << endl;
        }
        X = GAUSS_METHOD(A, F, n);
        //Output_Vector(GAUSS_METHOD(A, F, n), n);
        different(X, n);
        cout << "Решение записано в файл X.txt";
        break;
    }
    case 2:
    {
        Input_size("SIZE.txt", n);
        Input_arr(DI, IA, AL, n, n_profile);
        Input_vector(F, n);
        Y = new real[n](), X = new real[n]();
        vector<vector<real>>A = Profile_To_Dense(DI, AL, IA, F, n);
        for ( int i = 0; i < n; i++ )
        {
            for ( int j = 0; j < n; j++ )
                cout << A[i][j] << ' ';
            cout << endl;
        }
        LLT(DI, AL, IA, n, n_profile);
        CalcY(AL, DI, IA, F, Y, n);
        CalcX(AL, DI, IA, Y, X, n);
        Output_Vector(X, n);
        different(X, n);
        cout << "Решение записано в файл X.txt";
        break;
    }
    case 3:
    {
        cout << "Введите размер матрицы Гильберта: ";
        cin >> n;
        Y = new real[n], X = new real[n];
        Hilbert_Matrix(DI, AL, IA, F, n, n_profile);
        vector<vector<real>>A = Profile_To_Dense(DI, AL, IA, F, n);
        for ( int i = 0; i < n; i++ )
        {
            for ( int j = 0; j < n; j++ )
                cout << A[i][j] << ' ';
            cout << endl;
        }
        LLT(DI, AL, IA, n, n_profile);
        CalcY(AL, DI, IA, F, Y, n);
        CalcX(AL, DI, IA, Y, X, n);
        Output_Vector(X, n);
        different(X, n);
        cout << "Решение записано в файл X.txt";
        break;
    }
    default:
    {
        cout << "Некорректный ввод команды.";
        break;
    }
    }

    delete[] DI;
    delete[] IA;
    delete[] AL;
    delete[] F;
    delete[] X;
    delete[] Y;

    return 0;
}