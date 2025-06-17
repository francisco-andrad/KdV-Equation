#include <cmath>
#include <fstream>
#include <iostream>
const int max_iterations = 100;
const double tolerance = 1e-6;
const double mu = 1;
const int space_steps = 20001;
const int time_steps = 20000;
const double x_init = -200.0;
const double x_final = 200.0;
const double dx = (x_final - x_init) / (space_steps - 1); // Será ~0.02
const double dt = 0.000001;                               // dt ultra-seguro

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

// Recebe u^n (xt) e a estimativa de u^{n+1} (xtplus1)
// Retorna em aux o valor de 0.5 * (F(xt) + F(xtplus1))
void space_finite_diff_CN(double *xt, double *xtplus1, double *aux, double dx)
{
    // Vetores temporários para F(u^n) e F(u^{n+1})
    double *F_n = new double[space_steps];
    double *F_n_plus_1 = new double[space_steps];

    // Calcula F(u^n) usando a sua função Leapfrog/RK4 já corrigida
    // O primeiro argumento é a onda, o segundo é onde o resultado é salvo
    space_finite_diff(xt, F_n, dx);

    // Calcula F(u^{n+1}) usando a mesma função
    space_finite_diff(xtplus1, F_n_plus_1, dx);

    // Calcula a média e salva em aux
    for (int i = 0; i < space_steps; i++)
    {
        aux[i] = 0.5 * (F_n[i] + F_n_plus_1[i]);
    }

    delete[] F_n;
    delete[] F_n_plus_1;
}

void time_crank_nicolson()
{
    std::fstream mass_file, energy_file, kdv_file;

    // --- Vetores necessários ---
    double *ic = new double[space_steps];           // Guarda u^n (u_atual)
    double *u_next_guess = new double[space_steps]; // Guarda a estimativa para u^{n+1}
    double *u_prev_guess = new double[space_steps]; // Guarda a estimativa anterior para a verificação
    double *F_average = new double[space_steps];    // Guarda a média de F
    double *ic_prime = new double[space_steps];
    double mass, energy;

    // --- Condição Inicial ---
    discretize_axis(ic);
    general_initial_conditions(ic, 13., 0, 0.);

    // --- Arquivos de Saída ---
    kdv_file.open("kdv_data.txt", std::ios::out);
    mass_file.open("mass_data.txt", std::ios::out);
    energy_file.open("energy_data.txt", std::ios::out);

    // --- Loop de Tempo Principal ---
    for (int i = 0; i < time_steps; i++)
    {
        // Chute inicial para a iteração: u^{n+1, k=0} = u^n
        for (int k = 0; k < space_steps; k++)
        {
            u_next_guess[k] = ic[k];
        }

        // --- Laço da Iteração de Ponto Fixo com Verificação ---
        for (int j = 0; j < max_iterations; j++)
        {
            // PASSO 1: Salvar a estimativa atual em 'u_prev_guess' ANTES de modificá-la.
            // ESTA É A CORREÇÃO CRUCIAL.
            for (int k = 0; k < space_steps; k++)
            {
                u_prev_guess[k] = u_next_guess[k];
            }

            // PASSO 2: Calcular a próxima estimativa e salvá-la em 'u_next_guess'
            space_finite_diff_CN(ic, u_next_guess, F_average, dx);
            linear_combination(1.0, ic, dt, F_average, u_next_guess);

            // PASSO 3: Verificar a convergência
            double error_norm = 0.0;
            for (int k = 0; k < space_steps; k++)
            {
                error_norm += pow(u_next_guess[k] - u_prev_guess[k], 2);
            }
            error_norm = sqrt(error_norm / space_steps); // Norma RMS - mais robusta

            if (error_norm < tolerance)
            {
                break; // Convergiu, podemos sair do laço de iteração
            }

            if (j == max_iterations - 1)
            {
                std::cout << "AVISO: Solver não convergiu no passo de tempo " << i << std::endl;
            }
        }
        // --- Fim do Laço de Iteração ---

        // A iteração terminou. Atualizamos a solução principal 'ic' para o próximo passo.
        for (int j = 0; j < space_steps; j++)
        {
            ic[j] = u_next_guess[j];
        }

        // --- Bloco para salvar os dados ---
        if (i % 1000 == 0)
        {
            mass = mass_conservation(ic);
            mass_file << std::scientific << mass << std::endl;
            calculate_first_x_derivative(ic, ic_prime);
            energy = energy_conservation(ic, ic_prime);
            energy_file << std::scientific << energy << std::endl;
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

    // --- Limpeza de memória ---
    delete[] ic;
    delete[] u_next_guess;
    delete[] u_prev_guess;
    delete[] F_average;
    delete[] ic_prime;
}

int main(void)
{
    time_crank_nicolson();
    return 0;
}
