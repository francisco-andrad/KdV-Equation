
#include <cmath>
#include <fstream>
#include <iostream>
const double mu = 1;
const int space_steps = 20001;                            // Mantido
const int time_steps = 10000000;                          // AVISO: Aumento significativo!
const double x_init = -200.0;                             // Novo domínio
const double x_final = 200.0;                             // Novo domínio
const double dx = (x_final - x_init) / (space_steps - 1); // Será ~0.02
const double dt = 0.000001; // dt drasticamente menor                            // incremento tempo

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

void testing_initial_conditions(double *ic, double x, double a, double b)
{
    for (int i = 0; i < space_steps; i++)
    {
        if ((ic[i] < b) && (ic[i] > a))
        {
            ic[i] = x;
        }
        else
        {
            ic[i] = 0;
        }
    }
}
// u(x,0) = A \cdot e^{-(x-x_0)^2 / w^2}
void gaussian_pulse_initial_conditions(double *ic, double A, double w, double x0)
{
    double x_val = x_init;
    for (int i = 0; i < space_steps; i++)
    {
        double arg = (x_val - x0) / w;
        ic[i] = A * exp(-arg * arg);
        x_val += dx;
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
    return sum;
}

double energy_conservation(double *x, double *x_prime)
{
    double sum = 0.0;
    // E(u) = 0.5 * (u_x)^2 + (mu/3.0) * |u|^3
    sum = (0.5 * pow(x_prime[0], 2)) + ((mu / 3.0) * pow(std::abs(x[0]), 3.0));
    for (int i = 1; i < space_steps - 1; i++)
    {
        double term = (0.5 * pow(x_prime[i], 2)) + ((mu / 3.0) * pow(std::abs(x[i]), 3.0));
        if (i % 2 != 0)
        {
            sum += 4 * term;
        }
        else
        {
            sum += 2 * term;
        }
    }
    sum += (0.5 * pow(x_prime[space_steps - 1], 2)) + ((mu / 3.0) * pow(std::abs(x[space_steps - 1]), 3.0));
    sum *= dx / 3.0;
    return sum;
}

void calculate_first_x_derivative(double *x, double *x_prime)
{
    x_prime[0] = (x[1] - x[space_steps - 1]) / (2.0 * dx);
    for (int i = 1; i < space_steps - 1; i++)
    {
        x_prime[i] = (x[i + 1] - x[i - 1]) / (2.0 * dx);
    }
    x_prime[space_steps - 1] = (x[0] - x[space_steps - 2]) / (2.0 * dx);
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
    // Ponto i = 0
    // termo não-linear: esquema de Zabusky e Kruskal
    aux[i] = (2.0 * mu * (std::abs(x[i + 1]) + std::abs(x[i]) + std::abs(x[space_steps - 1])) / 3.0) *
             ((x[i + 1] - x[space_steps - 1]) / (2.0 * dx));
    // derivada terceira:
    aux[i] += -((x[2] - 2.0 * x[1] + 2.0 * x[space_steps - 1] - x[space_steps - 2])) / (2.0 * dx * dx * dx);

    i = 1;
    // termo não-linear:
    aux[i] = (2.0 * mu * (std::abs(x[i + 1]) + std::abs(x[i]) + std::abs(x[i - 1])) / 3.0) *
             ((x[i + 1] - x[i - 1]) / (2.0 * dx));
    // derivada terceira:
    aux[1] += -((x[3] - 2.0 * x[2] + 2.0 * x[0] - x[space_steps - 1])) / (2.0 * dx * dx * dx);

    // Dentro do intervalo:
    for (i = 2; i < space_steps - 2; i++)
    {
        // termo não-linear:
        aux[i] = (2.0 * mu * (std::abs(x[i + 1]) + std::abs(x[i]) + std::abs(x[i - 1])) / 3.0) *
                 ((x[i + 1] - x[i - 1]) / (2.0 * dx));

        // derivada terceira:
        aux[i] += -((x[i + 2] - (2.0 * x[i + 1]) + (2.0 * x[i - 1]) - x[i - 2])) / (2.0 * dx * dx * dx);
    }

    i = space_steps - 2;
    // termo não-linear:
    aux[i] = (2.0 * mu * (std::abs(x[i + 1]) + std::abs(x[i]) + std::abs(x[i - 1])) / 3.0) *
             ((x[i + 1] - x[i - 1]) / (2.0 * dx));
    // derivada terceira:
    aux[i] +=
        -((x[0] - (2.0 * x[space_steps - 1]) + (2.0 * x[space_steps - 3]) - x[space_steps - 4])) / (2.0 * dx * dx * dx);

    i = space_steps - 1;
    // termo não-linear:
    aux[i] =
        (2.0 * mu * (std::abs(x[0]) + std::abs(x[i]) + std::abs(x[i - 1])) / 3.0) * ((x[0] - x[i - 1]) / (2.0 * dx));
    // derivada terceira:
    aux[i] += -((x[1] - (2.0 * x[0]) + (2.0 * x[space_steps - 2]) - x[space_steps - 3])) / (2.0 * dx * dx * dx);
}

// Resolve o vetor de EDO's com relação ao tempo usando um método de
// Runge-Kutta com convergência O(h⁴).
// O método é dado por y_(k+1) = y_(k) + h/6[f1 + 2*f2 + 2*f3 + f4], onde
// f1 = f(x_k,y_k) * h; f2 = (x_k + h/2, y_k + h/2*f1);
// f3 = f( x_k + h/2, y_k + h/2*f2) e f4 = (x_k + h, y_k + h*f3).

void time_rkf()
{

    std::fstream mass_file, ic_file, energy_file, kdv_file;

    double *f1 = new double[space_steps];
    double *f2 = new double[space_steps];
    double *f3 = new double[space_steps];
    double *f4 = new double[space_steps];
    double *ic = new double[space_steps];
    double *aux = new double[space_steps];
    double *aux_ic = new double[space_steps];
    double *ic_prox = new double[space_steps];
    double *ic_prime = new double[space_steps];
    double mass, energy;

    discretize_axis(ic);
    discretize_axis(aux_ic);

    // soliton_initial_conditions(ic, 2);
    // testing_initial_conditions(ic, 1.0, -5.0, 5.0); // condição inicial, pode-se alterar
    general_initial_conditions(ic, 13., 0, 0.);
    // linear_combination(1.0, ic, 1.0, aux_ic, ic);
    // gaussian_pulse_initial_conditions(ic, 6.0, 2.0, 0.0);
    //
    // discretize_axis(aux_ic);
    general_initial_conditions(aux_ic, 13., 0, 0.);
    // linear_combination(1.0, ic, 1.0, aux_ic, ic);
    ic_file.open("ic_data.txt", std::ios::out);
    for (int i = 0; i < space_steps; i++)
    {
        ic_file << ic[i] << std::endl;
    }
    ic_file.close();

    kdv_file.open("kdv_data.txt", std::ios::out);
    mass_file.open("mass_data.txt", std::ios::out);
    energy_file.open("energy_data.txt", std::ios::out);
    for (int i = 0; i < time_steps; i++)
    {
        // space_finite_diff(ic, f1, dx);
        // linear_combination(1.0, ic, dt / 2, f1, aux);
        // space_finite_diff(aux, f2, dx);
        // linear_combination(1.0, ic, dt / 2, f2, aux);
        // space_finite_diff(aux, f3, dx);
        // linear_combination(1.0, ic, dt, f3, aux);
        // space_finite_diff(aux, f4, dx);

        // for (int j = 0; j < space_steps; j++)
        // {
        //     ic[j] = ic[j] + ((dt / 6.0) * (f1[j] + 2.0 * f2[j] + 2.0 * f3[j] + f4[j]));
        // }

        if (i == 0)
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
            mass = mass_conservation(ic);
            calculate_first_x_derivative(ic, ic_prime);
            energy = energy_conservation(ic, ic_prime);
        }
        space_finite_diff(ic, f1, dx);
        linear_combination(1.0, aux_ic, 2.0 * dt, f1, ic_prox);
        linear_combination(1.0, ic, 0.0, f1, aux_ic);
        linear_combination(1.0, ic_prox, 0.0, f1, ic);

        if (i % 1000 == 0)
        {
            double aux_mass = mass;
            mass = mass_conservation(ic);
            mass_file << std::scientific << mass << std::endl;
            calculate_first_x_derivative(ic, ic_prime);
            energy = energy_conservation(ic, ic_prime);
            energy_file << std::scientific << energy << std::endl;

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
    energy_file.close();

    delete[] f1;
    delete[] f2;
    delete[] f3;
    delete[] f4;
    delete[] aux;
    delete[] ic;
    delete[] aux_ic;
    delete[] ic_prox;
    delete[] ic_prime;
}

int main(void)
{
    time_rkf();
    return 0;
}
