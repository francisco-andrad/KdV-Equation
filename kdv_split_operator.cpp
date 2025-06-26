#include <cmath>
#include <fstream>
#include <iostream>

const double mu = 0.0;
const int space_steps = 20001;
const int time_steps = 100000; // Para T=10 com dt=0.0001
const double x_init = -200.0;
const double x_final = 200.0;
const double dx = (x_final - x_init) / (space_steps - 1); //  0.02
const double dt = 0.0001;

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

// Aplica um filtro de baixa-passagem explícito de 3 pontos na solução u
// alpha: força do filtro
void apply_filter(double *u, double alpha)
{
    double *u_original = new double[space_steps];
    for (int i = 0; i < space_steps; i++)
    {
        u_original[i] = u[i];
    }

    for (int i = 0; i < space_steps; i++)
    {
        // Contorno periódico
        int im1 = (i == 0) ? space_steps - 1 : i - 1;
        int ip1 = (i == space_steps - 1) ? 0 : i + 1;

        u[i] = u_original[i] + alpha * (u_original[ip1] - 2.0 * u_original[i] + u_original[im1]);
    }

    delete[] u_original;
}

// Esta função calcula apenas o termo não-linear F_conv(u) = -2*mu*|u|*u_x
void calculate_nonlinear_rhs(double *u_in, double *rhs_out)
{
    for (int i = 0; i < space_steps; i++)
    {
        int im1 = (i == 0) ? space_steps - 1 : i - 1;
        int ip1 = (i == space_steps - 1) ? 0 : i + 1;

        rhs_out[i] = -(2.0 * mu * (std::abs(u_in[ip1]) + std::abs(u_in[i]) + std::abs(u_in[im1])) / 3.0) *
                     ((u_in[ip1] - u_in[im1]) / (2.0 * dx));
    }
}

// Resolve u_t = F_conv(u) por um passo dt_step usando RK4
void solve_nonlinear_step_rk4(double *u_in, double *u_out, double dt_step)
{
    double *k1 = new double[space_steps];
    double *k2 = new double[space_steps];
    double *k3 = new double[space_steps];
    double *k4 = new double[space_steps];
    double *u_temp = new double[space_steps];

    // k1 = F(u_in)
    calculate_nonlinear_rhs(u_in, k1);

    // k2 = F(u_in + (dt/2)*k1)
    for (int i = 0; i < space_steps; i++)
        u_temp[i] = u_in[i] + (dt_step / 2.0) * k1[i];
    calculate_nonlinear_rhs(u_temp, k2);

    // k3 = F(u_in + (dt/2)*k2)
    for (int i = 0; i < space_steps; i++)
        u_temp[i] = u_in[i] + (dt_step / 2.0) * k2[i];
    calculate_nonlinear_rhs(u_temp, k3);

    // k4 = F(u_in + dt*k3)
    for (int i = 0; i < space_steps; i++)
        u_temp[i] = u_in[i] + dt_step * k3[i];
    calculate_nonlinear_rhs(u_temp, k4);

    for (int i = 0; i < space_steps; i++)
    {
        u_out[i] = u_in[i] + (dt_step / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }

    delete[] k1;
    delete[] k2;
    delete[] k3;
    delete[] k4;
    delete[] u_temp;
}

// Resolve um sistema linear pentadiagonal M*x = r, onde M tem 5 diagonais.
void solve_pentadiagonal_system(double *a, double *b, double *c, double *d, double *e, double *r, double *x, int n)
{
    double *c_prime = new double[n];
    double *d_prime = new double[n];
    double *r_prime = new double[n];

    double m = b[0] / c[0];
    c_prime[0] = c[0];
    d_prime[0] = d[0];
    r_prime[0] = r[0];

    c_prime[1] = c[1] - m * d_prime[0];
    m = b[1] / c_prime[0];
    d_prime[1] = d[1] - m * e[0];
    r_prime[1] = r[1] - m * r_prime[0];

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

    // Última linha
    x[n - 1] = r_prime[n - 1] / c_prime[n - 1];

    // Penúltima linha
    x[n - 2] = (r_prime[n - 2] - d_prime[n - 2] * x[n - 1]) / c_prime[n - 2];

    // Linhas restantes (de n-3 até 0)
    for (int i = n - 3; i >= 0; i--)
    {
        x[i] = (r_prime[i] - d_prime[i] * x[i + 1] - e[i] * x[i + 2]) / c_prime[i];
    }

    delete[] c_prime;
    delete[] d_prime;
    delete[] r_prime;
}

// Resolve o passo de dispersão u_t = -u_xxx usando Crank-Nicolson.
void solve_dispersion_step(double *u_in, double *u_out, double dt_step)
{
    const int n = space_steps;
    const double k = dt_step / (4.0 * dx * dx * dx);

    double *a = new double[n]; // diagonal M(i, i-2)
    double *b = new double[n]; // diagonal M(i, i-1)
    double *c = new double[n]; // diagonal M(i, i)
    double *d = new double[n]; // diagonal M(i, i+1)
    double *e = new double[n]; // diagonal M(i, i+2)
    double *r = new double[n]; // vetor do lado direito

    for (int i = 0; i < n; i++)
    {
        a[i] = -k;
        b[i] = 2.0 * k;
        c[i] = 1.0;
        d[i] = -2.0 * k;
        e[i] = k;
    }

    a[0] = 0;
    b[0] = 0;
    a[1] = 0;
    e[n - 1] = 0;
    d[n - 1] = 0;
    e[n - 2] = 0;

    for (int i = 0; i < n; i++)
    {
        double u_im2 = (i > 1) ? u_in[i - 2] : 0.0;
        double u_im1 = (i > 0) ? u_in[i - 1] : 0.0;
        double u_ip1 = (i < n - 1) ? u_in[i + 1] : 0.0;
        double u_ip2 = (i < n - 2) ? u_in[i + 2] : 0.0;

        double u_xxx_n_stencil = (u_ip2 - 2.0 * u_ip1 + 2.0 * u_im1 - u_im2);

        r[i] = u_in[i] - k * u_xxx_n_stencil;
    }

    solve_pentadiagonal_system(a, b, c, d, e, r, u_out, n);

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

    discretize_axis(u_current);
    // gaussian_pulse_initial_conditions(u_current, 15.0, 5.0, 0.0);

    general_initial_conditions(u_current, 13., 0, 0.);

    kdv_file.open("kdv_data.txt", std::ios::out);
    mass_file.open("mass_data.txt", std::ios::out);
    energy_file.open("energy_data.txt", std::ios::out);

    for (int i = 0; i < time_steps; i++)
    {
        // A(dt/2) - Primeiro passo não-linear
        solve_nonlinear_step_rk4(u_current, u_temp1, dt / 2.0);

        // B(dt) - Passo inteiro de dispersão
        solve_dispersion_step(u_temp1, u_temp2, dt);

        // A(dt/2) - Segundo passo não-linear
        solve_nonlinear_step_rk4(u_temp2, u_current, dt / 2.0);

        // Aplica o filtro de frquências
        if (i > 0 && i % 100 == 0)
        {
            apply_filter(u_current, 0.5);
        }

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
