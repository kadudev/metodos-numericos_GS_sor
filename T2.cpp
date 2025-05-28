#include <iostream>
#include <vector>
#include <cmath>
#include <locale.h>

using namespace std;

//função para o método de Gauss-Seidel
vector<double> gaussSeidel(const vector<vector<double>>& A, const vector<double>& b, vector<double> x0, int max_inter, double eps) {
    int n = b.size(); //tamanho do sistema

    for (int k = 0; k < max_inter; k++) {
        vector<double> x_old = x0;

        for (int i = 0; i < n; i++) {
            double sum1 = 0.0;
            double sum2 = 0.0;

            for (int j = 0; j < i; j++) {
                sum1 += A[i][j] * x0[j];  //usa valor atualizado
            }

            for (int j = i + 1; j < n; j++) {
                sum2 += A[i][j] * x0[j];  //usa valor anterior
            }

            x0[i] = (b[i] - sum1 - sum2) / A[i][i];
        }

        //cálculo da variação de x
        double dx = 0.0; //dx indica quanto a solução está mudando a cada iteração.
        for (int i = 0; i < n; i++) {
            dx += pow(x0[i] - x_old[i], 2);
        }
        dx = sqrt(dx);

        //cálculo do erro
        double erro = 0.0;
        for (int i = 0; i < n; i++) {
            double Ax = 0.0;
            for (int j = 0; j < n; j++) {
                Ax += A[i][j] * x0[j];
            }
            erro += pow(Ax - b[i], 2);
        }
        erro = sqrt(erro);

        cout << "Iteração " << k << ": dx = " << dx << ", erro = " << erro << ", x = [";
for (int i = 0; i < n; i++) {
    cout << x0[i]; //impressão de x1,x2,x3
    if (i < n - 1) cout << ", ";
}
cout << "]" << endl;

        if (dx < eps) { //criterio de parada
            break;
        }
    }

    return x0;
}

//função para o método SOR-Gauss-Seidel
vector<double> SORGaussSeidel(const vector<vector<double>>& A, const vector<double>& b, vector<double> x0, double lambda, int max_inter, double eps) {
    int n = b.size();

    for (int k = 0; k < max_inter; k++) { //loop pricipal
        vector<double> x_old = x0;

        for (int i = 0; i < n; i++) {
            double sum1 = 0.0;
            double sum2 = 0.0;

            for (int j = 0; j < i; j++) {
                sum1 += A[i][j] * x0[j];
            }

            for (int j = i + 1; j < n; j++) {
                sum2 += A[i][j] * x0[j];
            }

            double x_new = (b[i] - sum1 - sum2) / A[i][i];
            x0[i] = lambda * x_new + (1 - lambda) * x0[i];
        }

        double dx = 0.0;
        for (int i = 0; i < n; i++) {
            dx += pow(x0[i] - x_old[i], 2);
        }
        dx = sqrt(dx);

        double erro = 0.0;
        for (int i = 0; i < n; i++) {
            double Ax = 0.0;
            for (int j = 0; j < n; j++) {
                Ax += A[i][j] * x0[j];
            }
            erro += pow(Ax - b[i], 2);
        }
        erro = sqrt(erro);

        cout << "Iteração " << k << ": dx = " << dx << ", erro = " << erro << ", x = [";
for (int i = 0; i < n; i++) {
    cout << x0[i];
    if (i < n - 1) cout << ", ";
}
cout << "]" << endl;

        if (dx < eps) {
            break;
        }
    }

    return x0;
}

int main() {
    setlocale(LC_ALL, "Portuguese");

    //definindo o sistema linear
    vector<vector<double>> A = {
        {4, -1, 0},
        {-1, 4, -1},
        {0, -1, 3}
    };

    //definição das contantes b
    vector<double> b = {15, 10, 10};
    vector<double> x0(3, 0.0); //todos os x's como 0
    int max_inter = 20; //maximo de 20 iterações
    double eps = 1e-6; //definição do eps
    double lambda = 1.05; //fator lambda para sor

    cout << "Solução Gauss-Seidel:" << endl;
    vector<double> solucao_gs = gaussSeidel(A, b, x0, max_inter, eps);
    cout << "Resultado: ";
    for (double x : solucao_gs) {
        cout << x << " ";
    }
    cout << endl;

    cout << "\nSolução SOR-Gauss-Seidel:" << endl;
    vector<double> solucao_sor = SORGaussSeidel(A, b, x0, lambda, max_inter, eps);
    cout << "Resultado: ";
    for (double x : solucao_sor) {
        cout << x << " ";
    }
    cout << endl;

    return 0;
}
