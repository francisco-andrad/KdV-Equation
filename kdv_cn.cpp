#include <cmath>
#include <fstream>
#include <iostream>
const int space_steps = 20001;                            // número de passos no espaço
const int time_steps = 1000000;                           // número de passos no tempo
const double x_init = -1000.0;                            // início do intervalo no espaço
const double x_final = 0.0;                               // fim do intervalo no espaço
const double dx = (x_final - x_init) / (space_steps - 1); // incremento espaço
const double dt = 0.00001;                                // incremento tempo

void general_initial_conditions(double *ic, double c, double t, double x0) // ta ok
{
    double aux;
    int i;
    for (i = 0; i < space_steps; i++)
    {
        aux = cosh((0.5 * sqrt(c) * (ic[i] - (c * t) - x0)));
        ic[i] = c / (2.0 * aux * aux);
    }
}

void soliton_initial_conditions(double *ic, int n)
{
    double aux;
    int i;
    for (i = 0; i < space_steps; i++)
    {
        aux = cosh(ic[i]);
        ic[i] = (double)(n) * (n + 1.0) / (aux * aux);
    }
}

void discretize_axis(double *x) // ta ok
{
    int i = 0;
    x[0] = x_init;
    for (i = 1; i < space_steps; i++)
    {
        x[i] = x[i - 1] + dx;
        // std::cout << x[i] << "\n";
    }
}

// Procedimento que calcula a combinação linear de dois vetores.
// No momento só funciona para esse caso específico, mas pode ser
// facilmente portada.
void linear_combination(double alpha, double *v1, double beta, double *v2, double *combination)
{
    for (int i = 0; i < space_steps; i++)
    {
        combination[i] = (alpha * v1[i]) + (beta * v2[i]);
    }
}

// Norma L2 da função
double mass_conservation(double *x)
{
    double sum = 0.0;
    sum += pow(x[0], 2.0);
    for (int i = 1; i < space_steps - 1; i++)
    {
        if (i % 2 != 0)
        {
            sum += 4 * pow(x[i], 2);
        }
        else
        {
            sum += 2 * pow(x[i], 2);
        }
        // sum += 2 * pow(x[i], 2.0);
    }
    sum += pow(x[space_steps - 1], 2);
    sum *= dx / 3;
    return sqrt(sum);
}

// Procedimento que discretiza e aproxima as derivadas espaciais, transformando
// o problema de EDP inicial em um sistema de EDO's com relação ao tempo.
void space_finite_diff(double *xt, double *xtplus1, double *aux, double dx)
{

    int i = 0;
    int ultimo = space_steps - 1;
    // passo zero: condição de fronteira periódica. Tudo indica que esse
    // valor possa ser substituido por 0.
    aux[i] = -6.0 *
             ((xtplus1[i + 1] + xt[i + 1] + xtplus1[i] + xt[i] + 0.0 + 0.0) *
              (xtplus1[i + 1] + xt[i + 1] - xtplus1[ultimo] - xt[ultimo])) /
             (24.0 * dx); // ta ok
    // derivada terceira:
    // aux[0] = 0;
    aux[i] += -((xtplus1[i + 2] + xt[i + 2]) - 2.0 * (xtplus1[i + 1] + xt[i + 1]) + 2.0 * (0.0) - (0.0)) /
              (4.0 * dx * dx * dx);
    ; // ta ok

    i = 1;

    // termo não linear:
    aux[i] = -6.0 *
             ((xtplus1[i + 1] + xt[i + 1] + xtplus1[i] + xt[i] + xtplus1[i - 1] + xt[i - 1]) *
              (xtplus1[i + 1] + xt[i + 1] - xtplus1[i - 1] - xt[i - 1])) /
             (24.0 * dx);

    // derivada terceira:
    aux[i] += -((xtplus1[i + 2] + xt[i + 2]) - 2.0 * (xtplus1[i + 1] + xt[i + 1]) + 2.0 * (xtplus1[i - 1] + xt[i - 1]) -
                (0.0)) /
              (4.0 * dx * dx * dx);

    // dentro do initervalo:
    for (i = 2; i < space_steps - 2; i++)
    {
        // termo não-linear:
        aux[i] = -6.0 *
                 ((xtplus1[i + 1] + xt[i + 1] + xtplus1[i] + xt[i] + xtplus1[i - 1] + xt[i - 1]) *
                  (xtplus1[i + 1] + xt[i + 1] - xtplus1[i - 1] - xt[i - 1])) /
                 (24.0 * dx);
        // derivada terceira:
        aux[i] += -((xtplus1[i + 2] + xt[i + 2]) - 2.0 * (xtplus1[i + 1] + xt[i + 1]) +
                    2.0 * (xtplus1[i - 1] + xt[i - 1]) - (xtplus1[i - 2] + xt[i - 2])) /
                  (4.0 * dx * dx * dx);
    }
    // para os últimos dois pontos, para os quais é necessário pontos
    // fora do intervalo:
    i = space_steps - 2;
    // termo não-linear:
    // aux[i] = 0;
    aux[i] = -6.0 *
             ((xtplus1[i + 1] + xt[i + 1] + xtplus1[i] + xt[i] + xtplus1[i - 1] + xt[i - 1]) *
              (xtplus1[i + 1] + xt[i + 1] - xtplus1[i - 1] - xt[i - 1])) /
             (24.0 * dx); // ta ok
    // derivada terceira:
    // aux[i] = 0;
    aux[i] += -((0.0) - 2.0 * (xtplus1[i + 1] + xt[i + 1]) + 2.0 * (xtplus1[i - 1] + xt[i - 1]) -
                (xtplus1[i - 2] + xt[i - 2])) /
              (4.0 * dx * dx * dx); // ta ok
    i = space_steps - 1;
    // termo não-linear:
    // aux[i] = 0;
    aux[i] = -6.0 *
             ((0.0 + 0.0 + xtplus1[i] + xt[i] + xtplus1[i - 1] + xt[i - 1]) *
              (xtplus1[0] + xt[0] - xtplus1[i - 1] - xt[i - 1])) /
             (24.0 * dx); // ta ok
    // derivada terceira:
    // aux[i] = 0;
    aux[i] += -((0.0) - 2.0 * (0.0) + 2.0 * (xtplus1[i - 1] + xt[i - 1]) - (xtplus1[i - 2] + xt[i - 2])) /
              (4.0 * dx * dx * dx);
}

// Resolve o vetor de EDO's com relação ao tempo
void time_rkf()
{

    std::fstream mass_file, ic_file, kdv_file;

    double *f1 = new double[space_steps];
    double *f2 = new double[space_steps];
    double *f3 = new double[space_steps];
    double *f4 = new double[space_steps];
    double *ic = new double[space_steps];
    double *aux = new double[space_steps];
    double *aux_ic = new double[space_steps];
    double mass;

    discretize_axis(ic);
    discretize_axis(aux_ic);

    // soliton_initial_conditions(ic, 2);
    general_initial_conditions(ic, 16.0, 0.0, -90.0); // condição inicial, pode-se alterar
    general_initial_conditions(aux_ic, 16.0, 0.0, -90.0);
    // linear_combination(1.0, ic, 1.0, aux_ic, ic);
    // general_initial_conditions(aux_ic, 8., 0, -11);
    // linear_combination(1.0, ic, 1.0, aux_ic, ic);
    ic_file.open("ic_data.txt", std::ios::out);
    for (int i = 0; i < space_steps; i++)
    {
        ic_file << ic[i] << std::endl;
    }
    ic_file.close();

    int cont_conv;
    kdv_file.open("kdv_data.txt", std::ios::out);
    mass_file.open("mass_data.txt", std::ios::out);
    for (int i = 0; i < time_steps; i++)
    {
        // 15 iterações do método de ponto fixo costuma funcionar para esse caso
        for (int m = 0; m < 15; m++)
        {
            space_finite_diff(ic, aux_ic, f1, dx);
            linear_combination(1.0, ic, dt, f1, aux_ic);
        }
        for (int j = 0; j < space_steps; j++)
        {
            ic[j] = aux_ic[j];
        }
        if (i == 0)
            mass = mass_conservation(ic);
        if (i % 1000 == 0)
        {
            double aux_mass = mass;
            mass = mass_conservation(ic);
            mass_file << std::scientific << mass << std::endl;
            // mass_file << std::endl;

            for (int k = 0; k < space_steps; k++)
            {
                kdv_file << ic[k] << std::endl;
            }
            kdv_file << std::endl;
        }
    }
    mass_file.close();
    kdv_file.close();

    delete[] f1;
    delete[] f2;
    delete[] f3;
    delete[] f4;
    delete[] aux;
    delete[] ic;
    delete[] aux_ic;
}

int main(void)
{
    time_rkf();
    return 0;
}
