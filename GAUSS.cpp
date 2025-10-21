#include "GAUSS.h"

void Input_size(const char* s, int& n)
{
    FILE* SIZE;
    fopen_s(&SIZE, s, "r");
    fscanf_s(SIZE, "%d", &n);
    fclose(SIZE);

}

void Input_arr(real*& DI, int*& IA, real*& AL, int& n, int& n_profile)
{
    FILE* DI_F, * IA_F, * AL_F;

    // Чтение IA
    fopen_s(&IA_F, "IA.txt", "r");
    IA = new int[n + 1];
    for ( int i = 0; i <= n; i++ )
        fscanf_s(IA_F, "%d", &IA[i]);
    fclose(IA_F);

    if ( IA[0] )
        for ( int i = 0; i <= n; i++ ) IA[i]--;

    n_profile = IA[n];

  // Чтение DI
    fopen_s(&DI_F, "DI.txt", "r");
    DI = new real[n];
    for ( int i = 0; i < n; i++ )
        fscanf_s(DI_F, REALIN, &DI[i]);
    fclose(DI_F);

    // Чтение AL
    fopen_s(&AL_F, "AL.txt", "r");
    AL = new real[n_profile];
    for ( int i = 0; i < n_profile; i++ )
        fscanf_s(AL_F, REALIN, &AL[i]);
    fclose(AL_F);
}

void Input_vector(real*& F, int& n)
{
    // Чтение F
    FILE* F_f;
    fopen_s(&F_f, "F.txt", "r");
    F = new real[n];
    for ( int i = 0; i < n; i++ )
        fscanf_s(F_f, REALIN, &F[i]);
    fclose(F_f);
}

vector<vector<real>> Profile_To_Dense(real* DI, real* AL, int* IA, real* F, int& n)
{
    vector<vector<real>> A(n, vector<real>(n, 0.0));

    // Заполняем диагональ
    for ( int i = 0; i < n; i++ )
        A[i][i] = DI[i];

    // Заполняем внедиагональные элементы
    for ( int i = 0; i < n; i++ )
    {
        int start = IA[i];
        int num_elem = IA[i + 1] - IA[i];

        for ( int k = 0; k < num_elem; k++ )
        {
            int j = i - num_elem + k;
            A[i][j] = AL[start + k];
            A[j][i] = AL[start + k];
        }
    }
    return A;
}


void Output_Vector(real* x, int& n)
{
    FILE* X_f;
    fopen_s(&X_f, "X.txt", "w");
    for ( int i = 0; i < n; i++ )
        fprintf_s(X_f, REALOUT, x[i]);
    fclose(X_f);
}

real* GAUSS_METHOD(vector<vector<real>> A, real* y, int& n)
{
    vector<vector<real>> new_A(n, vector<real>(n + 1));

    // Инициализизация расширенной матрицы
    for ( int i = 0; i < n; i++ )
    {
        for ( int j = 0; j < n; j++ )
            new_A[i][j] = A[i][j];
        new_A[i][n] = y[i];
    }

    // Прямой ход метода Гаусса
    for ( int k = 0; k < n; k++ )
    {
        // Поиск главного элемента
        int max_row = k;
        for ( int i = k + 1; i < n; i++ )
            if ( fabs(new_A[i][k]) > fabs(new_A[max_row][k]) )
                max_row = i;

        // Обмен строк
        if ( max_row != k )
            swap(new_A[k], new_A[max_row]);

        // Нормализация текущей строки
        double divisor = new_A[k][k];
        if ( fabs(divisor) < EPS )
        {
            cout << "Matrix is singular or nearly singular" << endl;
            return NULL;
        }

        for ( int j = k; j <= n; j++ )
            new_A[k][j] /= divisor;

        // Исключение переменной из последующих строк
        for ( int i = k + 1; i < n; i++ )
        {
            double factor = new_A[i][k];
            for ( int j = k; j <= n; j++ )
                new_A[i][j] -= new_A[k][j] * factor;
        }
    }

    // Обратный ход метода Гаусса
    real* x = new real[n];
    for ( int k = n - 1; k >= 0; k-- )
    {
        x[k] = new_A[k][n];
        for ( int i = k + 1; i < n; i++ )
            x[k] -= new_A[k][i] * x[i];
        x[k] /= new_A[k][k];
    }

    return x;
}

void Hilbert_Matrix(real*& DI, real*& AL, int*& IA, real*& F, int n, int& n_profile)
{
    IA = new int[n + 1];
    DI = new real[n];
    F = new real[n];

    IA[0] = 0;
    for ( int i = 0; i < n; i++ )
        IA[i + 1] = IA[i] + i;

    n_profile = IA[n];
    AL = new real[n_profile];

    real* x = new real[n];
    for ( int i = 0; i < n; i++ )
        x[i] = i + 1.0;

    for ( int i = 0; i < n; i++ )
    {
        DI[i] = 1.0 / (2.0 * i + 1.0);

        F[i] = 0.0;
        for ( int j = 0; j < n; j++ )
            F[i] += (1.0 / (i + j + 1.0)) * x[j];

        for ( int j = 0; j < i; j++ )
        {
            int idx = IA[i] + j;
            AL[idx] = 1.0 / (i + j + 1.0);
        }
    }
    delete[] x;
}

void LLT(real* DI, real* AL, int* IA, int& n, int& n_profile)
{
    for ( int i = 0; i < n; i++ )
    {
        int i0 = IA[i];
        int i1 = IA[i + 1];

        // Количество элементов в строке i
        int profile = i1 - i0;

        int first_col_i = i - profile;

        // Обработка внедиагональных элементов
        for ( int j = first_col_i; j < i; j++ )
        {
            realsum sum = 0.0;

            int j0 = IA[j];
            int j1 = IA[j + 1];
            int profile_j = j1 - j0;
            int first_col_j = j - profile_j;

            int common_start = max(first_col_i, first_col_j);

            for ( int k = common_start; k < j; k++ )
            {
                int ik = i0 + (k - first_col_i);
                int jk = j0 + (k - first_col_j);
                sum += AL[ik] * AL[jk];
            }

            int ij = i0 + (j - first_col_i);
            if ( ij < i1 )
                AL[ij] = (AL[ij] - sum) / DI[j];
        }
        realsum sum_diag = 0.0;
        for ( int j = first_col_i; j < i; j++ )
        {
            int ij = i0 + (j - first_col_i);
            if ( ij < i1 )
                sum_diag += AL[ij] * AL[ij];
        }
        DI[i] = sqrt(DI[i] - sum_diag);
    }
}

void CalcY(real* AL, real* DI, int* IA, real* F, real* Y, int& n)
{
    for ( int i = 0; i < n; i++ )
    {
        realsum sum = 0.0;
        int i0 = IA[i];
        int i1 = IA[i + 1];
        int first_col_i = i - (i1 - i0);

        for ( int j = first_col_i; j < i; j++ )
        {
            int ij = i0 + (j - first_col_i);
            if ( ij < i1 )
                sum += AL[ij] * Y[j];
        }
        Y[i] = (F[i] - sum) / DI[i];
    }
}

void CalcX(real* AL, real* DI, int* IA, real* Y, real* X, int& n)
{
    for ( int i = n - 1; i >= 0; i-- )
    {
        realsum sum = 0.0;
        for ( int j = i + 1; j < n; j++ )
        {
            int j0 = IA[j];
            int j1 = IA[j + 1];
            int first_col_j = j - (j1 - j0);

            if ( i >= first_col_j && i < j )
            {
                int ji = j0 + (i - first_col_j);
                if ( ji < j1 )
                    sum += AL[ji] * X[j];
            }
        }
        X[i] = (Y[i] - sum) / DI[i];

    }
}

void different(real* x, int n)
{
    FILE* DIFF;
    fopen_s(&DIFF, "DIFF.txt", "w");
    for ( int i = 0; i < n; i++ )
        fprintf_s(DIFF, "%E\n", fabs((real)(i + 1) - x[i]));
    fclose(DIFF);
}