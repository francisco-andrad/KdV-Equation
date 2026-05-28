#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring> // Incluído para memcpy
#include <fftw3.h>
#include <fstream>
#include <iostream>
#include <vector>

using Complex = std::complex<double>;

// --- Parâmetros da Simulação 2D ---
const int Nx = std::pow(2, 14);
// const int Ny = 1024;
const double Lx = 100.0;
// const double Ly = 100.0;
const double T_FINAL = 0.01;
const double DT = 0.00001;

// --- Parâmetros da Equação ZK (Generalizada) ---
const double P = 2.0;
const double MU = 1.0;

// --- Parâmetros da Condição Inicial ---
const double SOLITON_AMP = 12.0;
const double SOLITON_X0 = 0.0;
const double GAUSS_AMP = 10.0;    // Amplitude do pulso
const double GAUSS_X0 = 0.0;      // Posição central em x
const double GAUSS_Y0 = 0.0;      // Posição central em y
const double GAUSS_WIDTH_X = 5.0; // Largura do pulso em x
const double GAUSS_WIDTH_Y = 5.0; // Largura do pulso em y
// --- Parâmetros da Condição Inicial ("Lump" Soliton) ---
const double LUMP_AMP = 24.0; // Amplitude
const double LUMP_X0 = 0.0;   // Posição central em x
const double LUMP_Y0 = 0.0;   // Posição central em y
// Parâmetros que definem a forma do "lump"
const double LUMP_A = 3.0;
const double LUMP_B = 12.0;
const double LUMP_C = 1.0;

// --- Constantes Derivadas ---
const double DX = Lx / Nx;

// a biblioteca lida com a simetria
const int NUM_COEFFS = (Nx / 2) + 1;

const int max_iter = 10000;
const double tol = 1e-12;

void calculate_nonlinear_power_fourier(fftw_complex *u_hat_in, fftw_complex *n_hat_out, fftw_plan plan_bwd,
                                       fftw_plan plan_fwd)
{
    double *u_real = (double *)fftw_malloc(sizeof(double) * Nx);

    fftw_execute_dft_c2r(plan_bwd, u_hat_in, u_real);

    double norm_factor = 1.0 / (double)(Nx);
    for (int i = 0; i < Nx; ++i)
    {
        // é necessário uma constante de normalização, particularidade da biblioteca
        double u_val = u_real[i] * norm_factor;
        u_real[i] = pow(u_val, P);
    }

    fftw_execute_dft_r2c(plan_fwd, u_real, n_hat_out);
    fftw_free(u_real);
}

void perform_irk2_startup_step(fftw_complex *u_hat_0, fftw_complex *u_hat_1, const std::vector<double> &k,
                               fftw_plan plan_bwd, fftw_plan plan_fwd)
{
    fftw_complex *u_hat_half_guess = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NUM_COEFFS);
    fftw_complex *u_hat_half_new = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NUM_COEFFS);
    fftw_complex *n_hat_half = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NUM_COEFFS);
    memcpy(u_hat_half_guess, u_hat_0, sizeof(fftw_complex) * NUM_COEFFS);
    int iter = 0;

    std::cout << "  Iniciando iteracão de ponto fixo..." << std::endl;
    while (iter < max_iter)
    {
        calculate_nonlinear_power_fourier(u_hat_half_guess, n_hat_half, plan_bwd, plan_fwd);
        for (int i = 0; i < NUM_COEFFS; ++i)
        {
            Complex non_lin_hat_half = Complex(n_hat_half[i][0], n_hat_half[i][1]);
            Complex ik = Complex(0.0, k[i]);
            Complex u_hat_0_i = Complex(u_hat_0[i][0], u_hat_0[i][1]);

            Complex RHSh = u_hat_0_i + (MU * ((DT / 2.0) * (1.0) * ik) * non_lin_hat_half);
            Complex one = Complex(1, 0);
            Complex two = Complex(2, 0);
            Complex dt = Complex(DT, 0);

            RHSh = RHSh / (one + (ik * ik * ik) * (dt / two));

            u_hat_half_new[i][0] = RHSh.real();
            u_hat_half_new[i][1] = RHSh.imag();
        }

        double error = 0.0;

        for (int i = 0; i < NUM_COEFFS; ++i)
        {
            double diff_real = u_hat_half_new[i][0] - u_hat_half_guess[i][0];
            double diff_imag = u_hat_half_new[i][1] - u_hat_half_guess[i][1];

            double modulo = (diff_real * diff_real) + (diff_imag * diff_imag);

            if (i == 0 || i == NUM_COEFFS - 1)
            {
                error += modulo;
            }
            else
            {
                error += 2 * modulo;
            }
        }

        // Teorema de Plancherel discreto
        error = sqrt(error * (1.0 / (double)Nx));

        memcpy(u_hat_half_guess, u_hat_half_new, sizeof(fftw_complex) * NUM_COEFFS);

        if (error < tol)
        {
            std::cout << "  Ponto fixo convergiu em " << iter + 1 << " iteracoes." << std::endl;
            break;
        }
        iter++;
        if (iter >= max_iter)
        {
            std::cout << "  AVISO: Iteracao de ponto fixo nao convergiu! Erro = " << error << std::endl;
            break;
        }
    }

    for (int i = 0; i < NUM_COEFFS; ++i)
    {
        u_hat_1[i][0] = 2.0 * u_hat_half_new[i][0] - u_hat_0[i][0];
        u_hat_1[i][1] = 2.0 * u_hat_half_new[i][1] - u_hat_0[i][1];
    }

    fftw_free(u_hat_half_guess);
    fftw_free(u_hat_half_new);
    fftw_free(n_hat_half);
}

void perform_cn_leapfrog_step(fftw_complex *u_hat_antigo, fftw_complex *u_hat_atual, fftw_complex *u_hat_novo,
                              const std::vector<double> &k, fftw_plan plan_bwd, fftw_plan plan_fwd)
{
    fftw_complex *n_hat = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NUM_COEFFS);
    calculate_nonlinear_power_fourier(u_hat_atual, n_hat, plan_bwd, plan_fwd);
    for (int i = 0; i < NUM_COEFFS; ++i)
    {

        Complex non_lin_hat_half = Complex(n_hat[i][0], n_hat[i][1]);
        Complex ik = Complex(0.0, k[i]);
        Complex u_hat_ntime = Complex(u_hat_antigo[i][0], u_hat_antigo[i][1]);

        Complex RHSh = u_hat_ntime + (MU * (DT * (1.0) * ik) * non_lin_hat_half);
        Complex one = Complex(1, 0);
        Complex dt = Complex(DT, 0);

        RHSh = RHSh / (one + (ik * ik * ik) * dt);

        u_hat_atual[i][0] = RHSh.real();
        u_hat_atual[i][1] = RHSh.imag();
    }

    for (int i = 0; i < NUM_COEFFS; ++i)
    {
        u_hat_novo[i][0] = 2.0 * u_hat_atual[i][0] - u_hat_antigo[i][0];
        u_hat_novo[i][0] = 2.0 * u_hat_atual[i][1] - u_hat_antigo[i][1];
    }

    fftw_free(n_hat);
}

double mass_integral(const double *u_real)
{
    double sum = 0.0;
    // A normalização é necessária porque a IFFT da FFTW não é normalizada. ????????
    // double norm_factor = 1.0 / (double)(Nx);

    for (int i = 0; i < Nx; ++i)
    {
        double u_val = u_real[i];
        sum += u_val * u_val;
    }

    return sum * DX;
}

// Calcula a derivada parcial em X usando o método espectral
void calculate_derivative_fourier(const fftw_complex *u_hat, double *u_x_real, const std::vector<double> &k,
                                  fftw_plan plan_bwd)
{

    fftw_complex *ux_hat = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NUM_COEFFS);

    for (int i = 0; i < NUM_COEFFS; ++i)
    {
        Complex u_hat_comp(u_hat[i][0], u_hat[i][1]);
        Complex ikx(0.0, k[i]);
        Complex ux_hat_comp = ikx * u_hat_comp;
        ux_hat[i][0] = ux_hat_comp.real();
        ux_hat[i][1] = ux_hat_comp.imag();
    }

    fftw_execute_dft_c2r(plan_bwd, ux_hat, u_x_real);

    double norm_factor = 1.0 / (double)(Nx);
    for (int i = 0; i < Nx; ++i)
    {
        u_x_real[i] = u_x_real[i] * norm_factor;
    }

    fftw_free(ux_hat);
}

double energy_integral(const double *u_real, const double *u_x_real)
{
    double sum = 0.0;

    for (int i = 0; i < Nx; ++i)
    {
        double u_val = u_real[i];
        double ux_val = u_x_real[i];

        sum += 0.5 * (ux_val)*ux_val + (MU / (P + 1.0)) * pow(std::abs(u_val), P + 1.0);
    }
    return sum * DX;
}
double mass_center(const double *u_real)
{
    double sum = 0.0;

    for (int i = 0; i < Nx; ++i)
    {
        double x = -(Lx / 2.0) + (i * DX);
        int index = i * Nx + i;

        double u_val;
        double mass_density = u_val * u_val;

        sum += x * mass_density;
    }

    return (sum * DX) / mass_integral(u_real);
}

double energy_center(const double *u_real, const double *u_x_real)
{
    double sum = 0.0;

    for (int i = 0; i < Nx; ++i)
    {

        double x = -Lx / 2.0 + i * DX;

        double u_val = u_real[i];
        double ux_val = u_x_real[i];

        double energy_density = 0.5 * (ux_val * ux_val) + (MU / (P + 1.0)) * pow(std::abs(u_val), P + 1.0);
        sum += x * energy_density;
    }
    return sum * DX / energy_integral(u_real, u_x_real);
}

void save_solution(std::ofstream &solution_stream, std::ofstream &mass_stream, std::ofstream &energy_stream,
                   std::ofstream &mass_energy_error, std::ofstream &mon_stream, std::vector<double> &k, double *u_real,
                   fftw_plan fft_fwd, fftw_plan fft_bwd, double initial_mass, double initial_energy,
                   double current_time)
{
    int i = 0;
    for (i = 0; i < Nx; i++)
    {
        double x = -(Lx / 2.0) + (i * DX);
        solution_stream << current_time << "," << x << "," << u_real[i] << std::endl;
    }

    fftw_complex *u_hat = (fftw_complex *)malloc(NUM_COEFFS * sizeof(fftw_complex));
    fftw_complex *ux_hat = (fftw_complex *)malloc(NUM_COEFFS * sizeof(fftw_complex));
    double *ux_real = (double *)malloc(Nx * sizeof(double));
    fftw_execute_dft_r2c(fft_fwd, u_real, u_hat);

    calculate_derivative_fourier(u_hat, ux_real, k, fft_bwd);

    double current_mass = mass_integral(u_real);
    double current_energy = energy_integral(u_real, ux_real);

    mass_stream << std::scientific << current_mass << std::endl;
    energy_stream << std::scientific << current_energy << std::endl;

    std::cout << "Massa e energia atuais: " << current_mass << " " << current_energy << std::endl;
    double mass_error = std::abs(initial_mass - current_mass);
    double energy_error = std::abs(initial_energy - current_energy);

    mass_energy_error << current_time << "," << mass_error << "," << energy_error << std::endl;
    std::cout << mass_error << " " << energy_error << std::endl;
    double mass_c = mass_center(u_real);
    double energy_c = energy_center(u_real, ux_real);
    mon_stream << current_time << energy_c - mass_c << std::endl;

    // std::cout << "Tempo: " << current_time << "\n";

    free(u_hat);
    free(ux_hat);
    free(ux_real);
}

double *gaussian_initial_conditions(double *axis, double x0, double var, double c)
{
    for (int j = 0; j < Nx; ++j)
    {
        double x = -Lx / 2.0 + j * DX;
        double arg_x = (x - x0) / var;
        axis[j] = c * exp(-(arg_x * arg_x));
    }
    return axis;
}

int main(void)
{
    std::cout << "Iniciando simulação: KdV fft-leap-frog \n";

    double *u_real = (double *)fftw_malloc(sizeof(double) * Nx);
    fftw_complex *u_hat_antigo = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NUM_COEFFS);
    fftw_complex *u_hat_atual = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NUM_COEFFS);
    fftw_complex *u_hat_novo = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NUM_COEFFS);

    // Não preserva o input por definição
    fftw_plan plan_fwd = fftw_plan_dft_r2c_1d(Nx, u_real, u_hat_atual, FFTW_MEASURE | FFTW_PRESERVE_INPUT);
    fftw_plan plan_bwd = fftw_plan_dft_c2r_1d(Nx, u_hat_atual, u_real, FFTW_MEASURE | FFTW_PRESERVE_INPUT);

    std::cout << "Planos de transformada escolhidos! \n ";

    std::ofstream mass_file("mass_data.txt");
    std::ofstream energy_file("energy_data.txt");
    std::ofstream kdv_file("kdv_data.csv");
    kdv_file << "time,x,u\n";

    std::ofstream centers_file("centers_data.csv");
    centers_file << "time,mon,\n";

    std::ofstream error_file("error_kdv_data.csv");
    error_file << "time,mass_error,energy_error\n";

    std::vector<double> k(NUM_COEFFS);
    for (int i = 0; i < NUM_COEFFS; ++i)
    {
        k[i] = 2.0 * M_PI * ((double)i / (double)Nx);
    }

    u_real = gaussian_initial_conditions(u_real, 0.0, 1.0, 2.0);

    //

    fftw_execute_dft_r2c(plan_fwd, u_real, u_hat_atual);
    double initial_mass, initial_energy;
    double *ux_real_init = (double *)fftw_malloc(sizeof(double) * Nx);

    // Precisamos do u_real inicial normalizado para os cálculos
    // isso tudo é surto
    // fftw_complex *u_hat_temp = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NUM_COEFFS);
    // memcpy(u_hat_temp, u_hat_atual, sizeof(fftw_complex) * NUM_COEFFS);
    // fftw_execute_dft_c2r(plan_bwd, u_hat_temp, u_real);

    calculate_derivative_fourier(u_hat_atual, ux_real_init, k, plan_bwd);

    initial_mass = mass_integral(u_real);
    initial_energy = energy_integral(u_real, ux_real_init);

    std::cout << "Massa Inicial: " << std::scientific << initial_mass << std::endl;
    std::cout << "Energia Inicial: " << std::scientific << initial_energy << std::endl;

    fftw_free(ux_real_init);

    std::cout << "Executando passo inicial com IRK2 para encontrar u^1..." << std::endl;
    perform_irk2_startup_step(u_hat_atual, u_hat_novo, k, plan_bwd, plan_fwd);
    std::cout << "Passo inicial concluido." << std::endl;

    memcpy(u_hat_antigo, u_hat_atual, sizeof(fftw_complex) * NUM_COEFFS);
    memcpy(u_hat_atual, u_hat_novo, sizeof(fftw_complex) * NUM_COEFFS);

    //

    int num_steps = T_FINAL / DT;

    for (int i = 1; i < num_steps; ++i)
    {
        if (i % 100 == 0)
        {
            fftw_complex *u_hat_temp = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NUM_COEFFS);
            memcpy(u_hat_temp, u_hat_atual, sizeof(fftw_complex) * NUM_COEFFS);
            fftw_execute_dft_c2r(plan_bwd, u_hat_temp, u_real);
            double norm_factor = 1.0 / (double)(Nx);
            for (int j = 0; j < Nx; ++j)
            {
                u_real[j] = u_real[j] * norm_factor;
            }

            double current_time = i * DT;
            save_solution(kdv_file, mass_file, energy_file, error_file, centers_file, k, u_real, plan_fwd, plan_bwd,
                          initial_mass, initial_energy, current_time);

            free(u_hat_temp);

            std::cout << "Passo " << i << " / " << num_steps << std::endl;
        }

        perform_cn_leapfrog_step(u_hat_antigo, u_hat_atual, u_hat_novo, k, plan_bwd, plan_fwd);

        memcpy(u_hat_antigo, u_hat_atual, sizeof(fftw_complex) * NUM_COEFFS);
        memcpy(u_hat_atual, u_hat_novo, sizeof(fftw_complex) * NUM_COEFFS);
    }

    free(u_real);
    free(u_hat_antigo);
    free(u_hat_atual);
    free(u_hat_novo);
    return 0;
}
