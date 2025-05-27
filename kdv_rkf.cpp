#include <cmath>
#include <fstream>
#include <iostream>
const int space_steps = 2001;                             // número de passos no espaço
const int time_steps = 100000;                            // número de passos no tempo
const double x_init = -100.0;                             // início do intervalo no espaço
const double x_final = 0.0;                               // fim do intervalo no espaço
const double dx = (x_final - x_init) / (space_steps - 1); // incremento espaço
const double dt = 0.0001;                                 // incremento tempo

void general_initial_conditions(double *ic, double c, double t, double x0)
{
    double aux;
    int i;
    for (i = 0; i < space_steps; i++)
    {
        aux = cosh(0.5 * sqrt(c) * (ic[i] - (c * t) - x0));
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

void discretize_axis(double *x)
{
    int i = 0;
    x[0] = x_init;
    for (i = 1; i < space_steps; i++)
    {
        x[i] = x[i - 1] + dx;
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
// Para aproximar a derivada primeira, foi usado o esquema de diferenças centrais
// u_x(x_i,t_i) = [u_(x+1) - u_(x-1)]/2h; para aproximar a derivada terceira também foi
// usado o esquema de diferenças centrais, dado por u_xxx(x_i, t_i) = [u_xx_(i + 1) - u_xx_(i-1)]/2h.
// A derivada segunda foi aproximada usando o esquema u_xx(x_i, t_i) = [u_(i+1) - 2*u_i + u_(i-1)]/h².
// Note que não é necessário calcular a derivada primeira para calcular a derivada segunda.
void space_finite_diff(double *x, double *aux, double dx)
{
    int i = 0;
    // passo zero: condição de fronteira periódica. Tudo indica que esse
    // valor possa ser substituido por 0.
    // talvez o sinal do primeiro termo esteja errado
    // termo não-linear:
    aux[0] = 0;
    // aux[0] = -6.0 * x[0] * ((x[1] - x[space_steps - 1]) / (2.0 * dx));
    // derivada terceira:
    aux[0] = 0;
    // aux[0] += -((x[2] - 2.0 * x[1] + 2.0 * x[space_steps - 1] - x[space_steps - 2])) / (2.0 * dx * dx * dx);

    // passo 1:
    // termo não-linear:
    aux[0] = 0;
    // aux[1] = -6.0 * x[1] * ((x[2] - x[0]) / (2.0 * dx));
    // derivada terceira:
    aux[0] = 0;
    // aux[1] += -((x[3] - 2.0 * x[2] + 2.0 * x[0] - x[space_steps - 1])) / (2.0 * dx * dx * dx);

    // dentro do initervalo:
    for (i = 2; i < space_steps - 2; i++)
    {
        // termo não-linear:
        aux[i] = -6.0 * x[i] * ((x[i + 1] - x[i - 1]) / (2.0 * dx));
        // derivada terceira:
        aux[i] += -((x[i + 2] - (2.0 * x[i + 1]) + (2.0 * x[i - 1]) - x[i - 2])) / (2.0 * dx * dx * dx);
    }

    // para os últimos dois pontos, para os quais é necessário pontos
    // fora do intervalo:
    i = space_steps - 2;
    // termo não-linear:
    aux[0] = 0;
    // aux[i] = -6.0 * x[i] * ((x[i + 1] - x[i - 1]) / (2.0 * dx));
    // derivada terceira:
    aux[0] = 0;
    // aux[i] += -((x[i + 2 - space_steps] - (2.0 * x[i + 1]) + (2.0 * x[i - 1]) - x[i - 2])) / (2.0 * dx * dx * dx);
    i = space_steps - 1;
    // termo não-linear:
    aux[0] = 0;
    // aux[i] = -6.0 * x[i] * ((x[i + 1 - space_steps] - x[i - 1]) / (2.0 * dx));
    // derivada terceira:
    aux[0] = 0;
    // aux[i] += -((x[i + 2 - space_steps] - (2.0 * x[i + 1 - space_steps]) + (2.0 * x[i - 1]) - x[i - 2])) /
    //           (2.0 * dx * dx * dx);
}

// Resolve o vetor de EDO's com relação ao tempo usando um método de
// Runge-Kutta com convergência O(h⁴).
// O método é dado por y_(k+1) = y_(k) + h/6[f1 + 2*f2 + 2*f3 + f4], onde
// f1 = f(x_k,y_k) * h; f2 = (x_k + h/2, y_k + h/2*f1);
// f3 = f( x_k + h/2, y_k + h/2*f2) e f4 = (x_k + h, y_k + h*f3).

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
    general_initial_conditions(ic, 16., 0, -90.0); // condição inicial, pode-se alterar
    general_initial_conditions(aux_ic, 4., 0, -85);
    linear_combination(1.0, ic, 1.0, aux_ic, ic);
    //
    // discretize_axis(aux_ic);
    // general_initial_conditions(aux_ic, 8., 0, -11);
    // linear_combination(1.0, ic, 1.0, aux_ic, ic);
    ic_file.open("ic_data.txt", std::ios::out);
    for (int i = 0; i < space_steps; i++)
    {
        ic_file << ic[i] << std::endl;
    }
    ic_file.close();

    kdv_file.open("kdv_data.txt", std::ios::out);
    mass_file.open("mass_data.txt", std::ios::out);
    for (int i = 0; i < time_steps; i++)
    {
        space_finite_diff(ic, f1, dx);
        linear_combination(1.0, ic, dt / 2, f1, aux);
        space_finite_diff(aux, f2, dx);
        linear_combination(1.0, ic, dt / 2, f2, aux);
        space_finite_diff(aux, f3, dx);
        linear_combination(1.0, ic, dt, f3, aux);
        space_finite_diff(aux, f4, dx);

        for (int j = 0; j < space_steps; j++)
        {
            ic[j] = ic[j] + ((dt / 6.0) * (f1[j] + 2.0 * f2[j] + 2.0 * f3[j] + f4[j]));
        }

        if (i == 0)
            mass = mass_conservation(ic);
        if (i % 100 == 0)
        {
            double aux_mass = mass;
            mass = mass_conservation(ic);
            mass_file << std::scientific << mass << std::endl;
            // mass_file << std::endl;

            if (std::abs(mass - aux_mass) > 0.1)
            {
                // std::cout << "Conservação de massa violada.\n";
                // return;
            }

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
