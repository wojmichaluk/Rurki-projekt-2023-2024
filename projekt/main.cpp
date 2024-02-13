#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

typedef double (*Func1vFp)(double);
typedef double (*Func3vFp)(double, double, double);

// funkcje pojawiajace sie w operatorach - liniowym i biliniowym
double k(double);
double ei(double, double, double);
double eprim(double, double, double);

// funkcje pomocnicze do rozwiazywania ukladu macierzowego
void init_matrix(double [], int);
void init_result_column(double [], int);
void build_matrix(double [], int, double);
double coefficient(int, int, double, int);
double quad_simpson(Func1vFp, Func3vFp, int, int, double, int);
void echelon_form(double [], double [], int);
void solve(double [], const double [], double [], int);

// minimalna liczba punktow w calkowaniu po pojedynczym przedziale
const int MIN_INTEGRATION_POINTS = 100;

// dziedzina funkcji
const double LOWER_BOUND = 0;
const double UPPER_BOUND = 2;

// glowna funkcja programu
int main(int argc, char *argv[]) {
    int n;

    // zakladam, ze dane sa poprawne - poprawnego typu i sensowne
    if (argc >= 2) {
        n = (int) *argv[1];
    } else {
        cout << "Podaj n - ilosc elementow w MES: ";
        cin >> n;
    }

    double* matrix = new double [n * n];
    double* f_column = new double [n];
    double* u_column = new double [n];
    ofstream file_writer;
    string filename = "results.csv";

    // inicjalizacja macierzy
    init_matrix(matrix, n);
    init_result_column(f_column, n);

    // rozwiazywanie macierzowego ukladu rownan
    build_matrix(matrix, n, (UPPER_BOUND - LOWER_BOUND) / n);
    echelon_form(matrix, f_column, n);
    solve(matrix, f_column, u_column, n);

    // zapisywanie wartosci do pliku
    file_writer.open(filename);
    file_writer << "i,ui" << endl;

    for (int i = 0; i < n; i++) {
        file_writer << i << "," << u_column[i] << endl;
    }

    file_writer.close();
    delete [] matrix;
    delete [] f_column;
    delete [] u_column;
    return 0;
}

// funkcja k wystepujaca w rownaniu transferu ciepla
double k(double x) {
    return x > 1.0 ? 2.0 : 1.0;
}

// funkcja bazowa ei
double ei(double xi, double h, double x) {
    if (x < xi - h || x > xi + h) {
        return 0.0;
    } else if (x < xi) {
        return (x - xi + h) / h;
    } else {
        return (xi + h - x) / h;
    }
}

// pochodna funkcji bazowej ei
double eprim(double xi, double h, double x) {
    if (x < xi - h || x > xi + h) {
        return 0.0;
    } else if (x < xi) {
        return 1.0 / h;
    } else {
        return -1.0 / h;
    }
}

// funkcja inicjalizuje macierz - wypelnia ja zerami
void init_matrix(double A[], int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[n * i + j] = 0.0;
        }
    }
}

// funkcja uzupelnia macierz (kolumne) wystepujaca po prawej
// stronie rownania - wartosci podane dla konkretnego zadania
void init_result_column(double f[], int n) {
    f[0] = -20.0;

    for (int i = 1; i < n; i++) {
        f[i] = 0.0;
    }
}

// funkcja buduje macierz z odpowiednimi wspolczynnikami, korzystajac
// z wlasnosci funkcji bazowych ei - macierz wstegowa (trojdiagonalna)
void build_matrix(double A[], int n, double h) {
    A[0] = coefficient(0, 0, h, n);
    A[1] = coefficient(0, 1, h, n);

    for (int i = 1; i < n - 1; i++) {
        for (int j = i - 1; j < i + 2; j++){
            A[n * i + j] = coefficient(i, j, h, n);
        }
    }

    A[n * n - 2] = coefficient(n - 1, n - 2, h, n);
    A[n * n - 1] = coefficient(n - 1, n - 1, h, n);
}

// funkcja oblicza wspolczynnik w podanej komorce macierzy
double coefficient(int row, int column, double h, int n) {
    return -ei(column * h, h, 0.0) * ei(row * h, h, 0.0) + quad_simpson(k, eprim, row, column, h, n);
}

// funkcja oblicza calke z przekazanej funkcji f na zbiorze,
// gdzie jest niezerowa, korzystajac z kwadratury Simpsona
double quad_simpson(Func1vFp f1, Func3vFp f3, int row, int column, double h, int n) {
    double sum = 0.0;
    int lower = max(row, column) == 0 ? 0 : max(row, column) - 1;
    int upper = min(row, column) + 1;
    double a = lower * h;
    double b = upper * h;
    int integration_points = max(MIN_INTEGRATION_POINTS, n);
    double ih = (b - a) / (integration_points + 2);
    double first = f1(a + ih) * f3(column * h, h, a + ih) * f3(row * h, h, a + ih);
    double mid;
    double last;

    for (int i = 2; i < integration_points + 2; i++){
        last = f1(a + i * ih) * f3(column * h, h, a + i * ih) * f3(row * h, h, a + i * ih);
        mid = f1(a + (i - 0.5) * ih) * f3(column * h, h, a + (i - 0.5) * ih) * f3(row * h, h, a + (i - 0.5) * ih);
        sum += 1.0 / 6.0 * ih * (first + 4.0 * mid + last);
        first = last;
    }

    return sum;
}

// funkcja przeksztalca macierz w postac schodkowa,
// przygotowujac do rozwiazania ukladu metoda Gaussa
void echelon_form(double A[], double f[], int n) {
    double element;
    double proportion;
    int j;

    for (int i = 0; i < n - 1; i++) {
        element = A[i * n + i];
        j = i + 1;
        proportion = A[j * n + i] / element;
        A[j * n + i] = 0.0;

        for (int k = j; k < n; k++)
            A[j * n + k] -= proportion * A[i * n + k];
        f[j] -= proportion * f[i];
    }
}

// funkcja oblicza rozwiazanie ukladu rownan metoda Gaussa
// zapisujac wartosci kolejnych wspolczynnikow w tablicy u
void solve(double A[], const double f[], double u[], int n) {
    double sum = 0.0;

    for (int i = n - 1; i >= 0; i--){
        for (int j = 0; j < n - 1 - i; j++)
            sum -= A[i * n + n - j - 1] * u[n - j - 1];

        sum += f[i];
        u[i] = sum / A[i * n + i];
        sum = 0.0;
    }
}
