#include <cmath>
#include <fstream>
#include <iostream>

// --- PARÂMETROS DA SIMULAÇÃO (ALTA RESOLUÇÃO COM CRANK-NICOLSON) ---
const double mu = 1.0;
const int space_steps = 20001;
const int time_steps = 100000; // Para T=10 com dt=0.0001. Pode ajustar conforme necessário.
const double x_init = -200.0;
const double x_final = 200.0;
const double dx = (x_final - x_init) / (space_steps - 1); // Aprox. 0.02
const double dt = 0.0001;                                 // dt relativamente grande, graças ao método implícito

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

// Esta função calcula apenas o termo não-linear F_conv(u) = -2*mu*|u|*u_x
// Ela usa o esquema conservativo de Zabusky & Kruskal.
// u_in: o estado da onda de entrada
// rhs_out: o vetor de saída onde o resultado será salvo
void calculate_nonlinear_rhs(double *u_in, double *rhs_out)
{
    int i = 0;
    // Ponto i = 0 (Contorno Periódico)
    rhs_out[i] = -(2.0 * mu * (std::abs(u_in[i + 1]) + std::abs(u_in[i]) + std::abs(u_in[space_steps - 1])) / 3.0) *
                 ((u_in[i + 1] - u_in[space_steps - 1]) / (2.0 * dx));

    // Ponto i = 1
    i = 1;
    rhs_out[i] = -(2.0 * mu * (std::abs(u_in[i + 1]) + std::abs(u_in[i]) + std::abs(u_in[i - 1])) / 3.0) *
                 ((u_in[i + 1] - u_in[i - 1]) / (2.0 * dx));

    // Dentro do intervalo
    for (i = 2; i < space_steps - 1; i++) // Otimizado para ir até o penúltimo ponto
    {
        rhs_out[i] = -(2.0 * mu * (std::abs(u_in[i + 1]) + std::abs(u_in[i]) + std::abs(u_in[i - 1])) / 3.0) *
                     ((u_in[i + 1] - u_in[i - 1]) / (2.0 * dx));
    }

    // Ponto i = space_steps - 1 (Contorno Periódico)
    i = space_steps - 1;
    rhs_out[i] = -(2.0 * mu * (std::abs(u_in[0]) + std::abs(u_in[i]) + std::abs(u_in[i - 1])) / 3.0) *
                 ((u_in[0] - u_in[i - 1]) / (2.0 * dx));
}

// Recebe um estado u_in e avança no tempo por 'dt_step' usando Euler Explícito
// u_in: o estado inicial do passo
// u_out: o vetor de saída para salvar o resultado
// dt_step: a duração do passo (pode ser dt ou dt/2)
void solve_nonlinear_step(double *u_in, double *u_out, double dt_step)
{
    // Aloca memória para o resultado do lado direito (RHS)
    double *rhs = new double[space_steps];

    // 1. Calcula o lado direito F_conv(u_in)
    calculate_nonlinear_rhs(u_in, rhs);

    // 2. Aplica a fórmula de Euler: u_out = u_in + dt_step * rhs
    for (int i = 0; i < space_steps; i++)
    {
        u_out[i] = u_in[i] + dt_step * rhs[i];
    }

    // Libera a memória temporária
    delete[] rhs;
}

// Resolve um sistema linear pentadiagonal M*x = r, onde M tem 5 diagonais.
// a: diagonal inferior (a_i = M_{i, i-2})
// b: diagonal sub-principal (b_i = M_{i, i-1})
// c: diagonal principal (c_i = M_{i, i})
// d: diagonal super-principal (d_i = M_{i, i+1})
// e: diagonal superior (e_i = M_{i, i+2})
// r: o vetor do lado direito (RHS)
// x: o vetor da solução (saída)
// n: o tamanho dos vetores (space_steps)
void solve_pentadiagonal_system(double *a, double *b, double *c, double *d, double *e, double *r, double *x, int n)
{
    // Aloca vetores temporários para a eliminação
    double *c_prime = new double[n];
    double *d_prime = new double[n];
    double *r_prime = new double[n];

    // --- PASSO 1: Eliminação Progressiva ---

    // Primeira linha
    double m = b[0] / c[0];
    c_prime[0] = c[0];
    d_prime[0] = d[0];
    r_prime[0] = r[0];

    // Segunda linha
    c_prime[1] = c[1] - m * d_prime[0];
    m = b[1] / c_prime[0];
    d_prime[1] = d[1] - m * e[0];
    r_prime[1] = r[1] - m * r_prime[0];

    // Linhas do meio (de i=2 até n-1)
    for (int i = 2; i < n; i++)
    {
        double m1 = a[i] / c_prime[i - 2];
        double m2 = b[i] - m1 * d_prime[i - 2];
        c_prime[i] = c[i] - m1 * e[i - 2] - m2 / c_prime[i - 1] * d_prime[i - 1];
        if (i < n - 1)
        {
            d_prime[i] = d[i] - m2 / c_prime[i - 1] * e[i - 1];
        }
        r_prime[i] = r[i] - m1 * r_prime[i - 2] - m2 / c_prime[i - 1] * r_prime[i - 1];
    }

    // --- PASSO 2: Substituição Regressiva ---

    // Última linha
    x[n - 1] = r_prime[n - 1] / c_prime[n - 1];

    // Penúltima linha
    x[n - 2] = (r_prime[n - 2] - d_prime[n - 2] * x[n - 1]) / c_prime[n - 2];

    // Linhas restantes (de n-3 até 0)
    for (int i = n - 3; i >= 0; i--)
    {
        x[i] = (r_prime[i] - d_prime[i] * x[i + 1] - e[i] * x[i + 2]) / c_prime[i];
    }

    // Libera memória
    delete[] c_prime;
    delete[] d_prime;
    delete[] r_prime;
}

// Resolve o passo de dispersão u_t = -u_xxx usando Crank-Nicolson.
// u_in: o estado inicial do passo
// u_out: o vetor de saída para salvar o resultado
// dt_step: a duração do passo (será dt)
void solve_dispersion_step(double *u_in, double *u_out, double dt_step)
{
    // --- 1. Definir constantes e alocar memória ---
    const int n = space_steps;
    const double k = dt_step / (4.0 * dx * dx * dx);

    double *a = new double[n]; // diagonal M(i, i-2)
    double *b = new double[n]; // diagonal M(i, i-1)
    double *c = new double[n]; // diagonal M(i, i)
    double *d = new double[n]; // diagonal M(i, i+1)
    double *e = new double[n]; // diagonal M(i, i+2)
    double *r = new double[n]; // vetor do lado direito (RHS)

    // --- 2. Montar as 5 diagonais da matriz M ---
    // M = I + (dt/2)*D_xxx, onde D_xxx é o operador da derivada terceira
    for (int i = 0; i < n; i++)
    {
        // Coeficientes para a j-ésima linha da matriz M
        // M(i, i-2) = -k
        // M(i, i-1) = +2k
        // M(i, i)   = 1
        // M(i, i+1) = -2k
        // M(i, i+2) = +k
        a[i] = -k;
        b[i] = 2.0 * k;
        c[i] = 1.0;
        d[i] = -2.0 * k;
        e[i] = k;
    }
    // Tratamento de fronteiras (simplificação para matriz não-circular)
    // Para as primeiras duas linhas, não há termos u_{i-1} ou u_{i-2}
    a[0] = 0;
    b[0] = 0;
    a[1] = 0;
    // Para as duas últimas, não há termos u_{i+1} ou u_{i+2}
    e[n - 1] = 0;
    d[n - 1] = 0;
    e[n - 2] = 0;

    // --- 3. Calcular o vetor do lado direito (RHS), r ---
    // r = (I - (dt/2)*D_xxx) * u_in
    for (int i = 0; i < n; i++)
    {
        // u_{xxx} em u_in[i], com cuidado nas fronteiras (u=0 fora do domínio)
        double u_im2 = (i > 1) ? u_in[i - 2] : 0.0;
        double u_im1 = (i > 0) ? u_in[i - 1] : 0.0;
        double u_ip1 = (i < n - 1) ? u_in[i + 1] : 0.0;
        double u_ip2 = (i < n - 2) ? u_in[i + 2] : 0.0;

        double u_xxx_n = (u_ip2 - 2.0 * u_ip1 + 2.0 * u_im1 - u_im2) / (2.0 * dx * dx * dx);

        r[i] = u_in[i] - (dt_step / 2.0) * u_xxx_n;
    }

    // --- 4. Chamar o solver pentadiagonal ---
    // A solução será salva em u_out
    solve_pentadiagonal_system(a, b, c, d, e, r, u_out, n);

    // --- 5. Liberar memória ---
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
    delete[] e;
    delete[] r;
}

void time_operator_splitting()
{
    std::fstream mass_file, energy_file, kdv_file;

    double *u_current = new double[space_steps];
    double *u_temp1 = new double[space_steps];
    double *u_temp2 = new double[space_steps];
    double *u_prime = new double[space_steps];
    double mass, energy;

    // Condição Inicial
    discretize_axis(u_current);
    gaussian_pulse_initial_conditions(u_current, 15.0, 5.0, 0.0);

    // Abrir arquivos
    kdv_file.open("kdv_data.txt", std::ios::out);
    mass_file.open("mass_data.txt", std::ios::out);
    energy_file.open("energy_data.txt", std::ios::out);

    // Loop de tempo principal usando Strang Splitting
    for (int i = 0; i < time_steps; i++)
    {
        // A(dt/2) - Primeiro passo de convecção
        solve_nonlinear_step(u_current, u_temp1, dt / 2.0);

        // B(dt) - Passo inteiro de dispersão
        solve_dispersion_step(u_temp1, u_temp2, dt);

        // A(dt/2) - Segundo passo de convecção
        solve_nonlinear_step(u_temp2, u_current, dt / 2.0);
        // u_current agora contém a solução no tempo n+1

        // Salva os dados a cada N passos
        if (i % 1000 == 0)
        {
            mass = mass_conservation(u_current);
            mass_file << std::scientific << mass << std::endl;
            calculate_first_x_derivative(u_current, u_prime);
            energy = energy_conservation(u_current, u_prime);
            energy_file << std::scientific << energy << std::endl;

            for (int k = 0; k < space_steps; k++)
            {
                kdv_file << u_current[k] << std::endl;
            }
            kdv_file << std::endl;
        }
    }

    mass_file.close();
    kdv_file.close();
    energy_file.close();

    delete[] u_current;
    delete[] u_temp1;
    delete[] u_temp2;
    delete[] u_prime;
}

int main(void)
{
    time_operator_splitting();
    return 0;
}
