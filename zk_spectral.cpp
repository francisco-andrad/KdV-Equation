// #include <algorithm>
#include <cmath>
#include <complex>
#include <cstring> // Incluído para memcpy
#include <fftw3.h>
#include <fstream>
#include <iostream>
#include <vector>

using Complex = std::complex<double>;

// --- Parâmetros da Simulação 2D ---
const int Nx = 1024;
const int Ny = 1024;
const double Lx = 100.0;
const double Ly = 100.0;
const double T_FINAL = 0.4;
const double DT = 0.0001;

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
const double DY = Ly / Ny;
const int NUM_COEFFS = Ny * (Nx / 2 + 1);

// Protótipos de Funções
void save_solution(std::ofstream &file, const double *u_real, double current_time);
void calculate_nonlinear_rhs_fourier(fftw_complex *u_hat_in, fftw_complex *n_hat_out,
                                     const std::vector<std::vector<double>> &kx,
                                     const std::vector<std::vector<double>> &ky, fftw_plan plan_bwd,
                                     fftw_plan plan_fwd);
void perform_irk2_startup_step(fftw_complex *u_hat_0, fftw_complex *u_hat_1, const std::vector<std::vector<double>> &kx,
                               const std::vector<std::vector<double>> &ky, fftw_plan plan_bwd, fftw_plan plan_fwd);
void perform_cn_leapfrog_step(fftw_complex *u_hat_antigo, fftw_complex *u_hat_atual, fftw_complex *u_hat_novo,
                              const std::vector<std::vector<double>> &kx, const std::vector<std::vector<double>> &ky,
                              fftw_plan plan_bwd, fftw_plan plan_fwd);
void calculate_x_derivative_fourier(const fftw_complex *u_hat, double *u_x_real,
                                    const std::vector<std::vector<double>> &kx, fftw_plan plan_bwd);
void calculate_y_derivative_fourier(const fftw_complex *u_hat, double *u_y_real,
                                    const std::vector<std::vector<double>> &kx, fftw_plan plan_bwd);
double mass_conservation_2d(const double *u_real);
double energy_conservation_2d(const double *u_real, const double *u_x_real, const double *u_y_real);
double energy_center_y_2d(const double *u_real, const double *u_x_real, const double *u_y_real);
double mass_center_x_2d(const double *u_real);
double mass_center_y_2d(const double *u_real);
double energy_center_x_2d(const double *u_real, const double *u_x_real, const double *u_y_real);

int main()
{
    std::cout << "Iniciando simulacao ZK com Crank-Nicolson-Leapfrog..." << std::endl;

    double *u_real = (double *)fftw_malloc(sizeof(double) * Nx * Ny);
    fftw_complex *u_hat_antigo = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NUM_COEFFS);
    fftw_complex *u_hat_atual = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NUM_COEFFS);
    fftw_complex *u_hat_novo = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NUM_COEFFS);

    fftw_plan plan_fwd = fftw_plan_dft_r2c_2d(Ny, Nx, u_real, u_hat_atual, FFTW_ESTIMATE);
    fftw_plan plan_bwd = fftw_plan_dft_c2r_2d(Ny, Nx, u_hat_atual, u_real, FFTW_ESTIMATE);

    std::ofstream mass_file("mass_data_2d.txt");
    std::ofstream energy_file("energy_data_2d.txt");
    std::ofstream kdv_file("zk_data.csv");
    kdv_file << "time,x,y,u\n";

    std::ofstream centers_file("centers_data.csv");
    centers_file << "time,monoton_x,monoton_y\n";

    std::ofstream error_file("error_zk_data.txt");
    error_file << "time,mass_error,energy_error\n";

    std::vector<std::vector<double>> kx(Ny, std::vector<double>(Nx / 2 + 1));
    std::vector<std::vector<double>> ky(Ny, std::vector<double>(Nx / 2 + 1));
    for (int i = 0; i < Ny; ++i)
    {
        for (int j = 0; j < Nx / 2 + 1; ++j)
        {
            kx[i][j] = 2.0 * M_PI * j / Lx;
            if (i < Ny / 2)
            {
                ky[i][j] = 2.0 * M_PI * i / Ly;
            }
            else
            {
                ky[i][j] = 2.0 * M_PI * (i - Ny) / Ly;
            }
        }
    }

    // Usando a condição inicial Gaussiana
    for (int i = 0; i < Ny; ++i)
    {
        double y = -Ly / 2.0 + i * DY;
        for (int j = 0; j < Nx; ++j)
        {
            double x = -Lx / 2.0 + j * DX;
            double arg_x = (x - GAUSS_X0) / 5.0; // Usando largura 5 como exemplo
            double arg_y = (y - 0.0) / 5.0;
            u_real[i * Nx + j] = SOLITON_AMP * exp(-(arg_x * arg_x + arg_y * arg_y));
        }
    }

    // for (int i = 0; i < Ny; ++i)
    // {
    //     double y = -Ly / 2.0 + i * DY;
    //     for (int j = 0; j < Nx; ++j)
    //     {
    //         double x = -Lx / 2.0 + j * DX;
    //         double sech_arg = 0.5 * sqrt(SOLITON_AMP) * (x - SOLITON_X0);
    //         u_real[i * Nx + j] = 0.5 * SOLITON_AMP / (cosh(sech_arg) * cosh(sech_arg));
    //     }
    // }

    // for (int i = 0; i < Ny; ++i)
    // {
    //     double y = -Ly / 2.0 + i * DY;
    //     for (int j = 0; j < Nx; ++j)
    //     {
    //         double x = -Lx / 2.0 + j * DX;
    //
    //         double term_x = x - LUMP_X0;
    //         double term_y = y - LUMP_Y0;
    //
    //         double denominator_part1 = LUMP_C + (term_x * term_x) + LUMP_A * (term_y * term_y);
    //         double denominator_part2 = LUMP_B * (term_y * term_y);
    //
    //         double numerator = LUMP_C - (term_x * term_x) + LUMP_A * (term_y * term_y);
    //         double denominator = (denominator_part1 * denominator_part1) + denominator_part2;
    //
    //         int index = i * Nx + j;
    //         if (denominator > 1e-9)
    //         { // Evitar divisão por zero no centro exato
    //             u_real[index] = LUMP_AMP * (numerator / denominator);
    //         }
    //         else
    //         {
    //             u_real[index] = 0.0;
    //         }
    //     }
    // }if (
    // ) {

    fftw_execute_dft_r2c(plan_fwd, u_real, u_hat_atual);

    double initial_mass, initial_energy;
    {
        double *u_x_real_init = (double *)fftw_malloc(sizeof(double) * Nx * Ny);
        double *u_y_real_init = (double *)fftw_malloc(sizeof(double) * Nx * Ny);

        // Precisamos do u_real inicial normalizado para os cálculos
        fftw_complex *u_hat_temp = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NUM_COEFFS);
        memcpy(u_hat_temp, u_hat_atual, sizeof(fftw_complex) * NUM_COEFFS);
        fftw_execute_dft_c2r(plan_bwd, u_hat_temp, u_real);

        calculate_x_derivative_fourier(u_hat_atual, u_x_real_init, kx, plan_bwd);
        calculate_y_derivative_fourier(u_hat_atual, u_y_real_init, ky, plan_bwd);

        initial_mass = mass_conservation_2d(u_real);
        initial_energy = energy_conservation_2d(u_real, u_x_real_init, u_y_real_init);

        std::cout << "Massa Inicial: " << std::scientific << initial_mass << std::endl;
        std::cout << "Energia Inicial: " << std::scientific << initial_energy << std::endl;

        fftw_free(u_x_real_init);
        fftw_free(u_y_real_init);
        fftw_free(u_hat_temp);
    }

    std::cout << "Executando passo inicial com IRK2 para encontrar u^1..." << std::endl;
    perform_irk2_startup_step(u_hat_atual, u_hat_novo, kx, ky, plan_bwd, plan_fwd);
    std::cout << "Passo inicial concluido." << std::endl;

    memcpy(u_hat_antigo, u_hat_atual, sizeof(fftw_complex) * NUM_COEFFS);
    memcpy(u_hat_atual, u_hat_novo, sizeof(fftw_complex) * NUM_COEFFS);

    int num_steps = T_FINAL / DT;
    for (int i = 1; i < num_steps; ++i)
    {
        if (i % 40 == 0)
        {
            // --- BLOCO DE SALVAMENTO SEGURO ---
            fftw_complex *u_hat_temp = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NUM_COEFFS);
            memcpy(u_hat_temp, u_hat_atual, sizeof(fftw_complex) * NUM_COEFFS);
            fftw_execute_dft_c2r(plan_bwd, u_hat_temp, u_real);

            double current_time = i * DT;
            save_solution(kdv_file, u_real, current_time);

            fftw_free(u_hat_temp);
            // --- FIM DO BLOCO DE SALVAMENTO SEGURO ---

            // Alocar memória para as derivadas
            double *u_x_real = (double *)fftw_malloc(sizeof(double) * Nx * Ny);
            double *u_y_real = (double *)fftw_malloc(sizeof(double) * Nx * Ny);

            // Calcular as derivadas a partir do estado atual no espaço de Fourier
            calculate_x_derivative_fourier(u_hat_atual, u_x_real, kx, plan_bwd);
            calculate_y_derivative_fourier(u_hat_atual, u_y_real, ky, plan_bwd);

            double current_mass = mass_conservation_2d(u_real);
            double current_energy = energy_conservation_2d(u_real, u_x_real, u_y_real);

            mass_file << std::scientific << current_mass << std::endl;
            energy_file << std::scientific << current_energy << std::endl;

            double mass_error = std::abs(current_mass - initial_mass);
            double energy_error = std::abs(current_energy - initial_energy);
            error_file << current_time << "," << mass_error << "," << energy_error << "\n";

            std::cout << current_time << "," << mass_error << "," << energy_error << "\n";

            double mass_cx = mass_center_x_2d(u_real);
            double mass_cy = mass_center_y_2d(u_real);
            double energy_cx = energy_center_x_2d(u_real, u_x_real, u_y_real);
            double energy_cy = energy_center_y_2d(u_real, u_x_real, u_y_real);

            centers_file << current_time << "," << energy_cx - mass_cx << "," << energy_cy - mass_cy << "\n";

            // Liberar memória das derivadas
            fftw_free(u_x_real);
            fftw_free(u_y_real);

            std::cout << "Passo " << i << "/" << num_steps << ", Tempo = " << current_time << std::endl;
        }

        perform_cn_leapfrog_step(u_hat_antigo, u_hat_atual, u_hat_novo, kx, ky, plan_bwd, plan_fwd);

        memcpy(u_hat_antigo, u_hat_atual, sizeof(fftw_complex) * NUM_COEFFS);
        memcpy(u_hat_atual, u_hat_novo, sizeof(fftw_complex) * NUM_COEFFS);
    }

    fftw_destroy_plan(plan_fwd);
    fftw_destroy_plan(plan_bwd);
    fftw_free(u_real);
    fftw_free(u_hat_antigo);
    fftw_free(u_hat_atual);
    fftw_free(u_hat_novo);
    kdv_file.close();
    mass_file.close();
    energy_file.close();
    centers_file.close();
    error_file.close();

    std::cout << "Simulacao concluida." << std::endl;
    return 0;
}

// =============================================================================
//               *** FUNÇÕES IMPLEMENTADAS / CORRIGIDAS ***
// =============================================================================

void save_solution(std::ofstream &file, const double *u_real, double current_time)
{
    for (int i = 0; i < Ny; ++i)
    {
        double y = -Ly / 2.0 + i * DY;
        for (int j = 0; j < Nx; ++j)
        {
            double x = -Lx / 2.0 + j * DX;
            int index = i * Nx + j;
            double u_val = u_real[index] / (double)(Nx * Ny); // Normaliza IFFT da FFTW
            file << current_time << "," << x << "," << y << "," << u_val << "\n";
        }
    }
}

void calculate_nonlinear_rhs_fourier(fftw_complex *u_hat_in, fftw_complex *n_hat_out,
                                     const std::vector<std::vector<double>> &kx,
                                     const std::vector<std::vector<double>> &ky, fftw_plan plan_bwd, fftw_plan plan_fwd)
{
    double *u_real = (double *)fftw_malloc(sizeof(double) * Ny * Nx);
    fftw_complex *u_hat_dealiased = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NUM_COEFFS);

    memcpy(u_hat_dealiased, u_hat_in, sizeof(fftw_complex) * NUM_COEFFS);

    int k_cutoff_x = Nx / 3;
    int k_cutoff_y = Ny / 3;
    for (int i = 0; i < Ny; ++i)
    {
        for (int j = 0; j < Nx / 2 + 1; ++j)
        {
            if (j > k_cutoff_x || std::abs(ky[i][j]) > (2.0 * M_PI * k_cutoff_y / Ly))
            {
                int index = i * (Nx / 2 + 1) + j;
                u_hat_dealiased[index][0] = 0.0;
                u_hat_dealiased[index][1] = 0.0;
            }
        }
    }

    fftw_execute_dft_c2r(plan_bwd, u_hat_dealiased, u_real);

    // MELHORIA DE CLAREZA: Separar normalização e cálculo da potência
    double norm_factor = 1.0 / (double)(Nx * Ny);
    for (int i = 0; i < Ny * Nx; ++i)
    {
        double u_val = u_real[i] * norm_factor;
        u_real[i] = pow(u_val, P);
    }

    fftw_execute_dft_r2c(plan_fwd, u_real, n_hat_out);

    for (int i = 0; i < Ny; ++i)
    {
        for (int j = 0; j < Nx / 2 + 1; ++j)
        {
            int index = i * (Nx / 2 + 1) + j;
            Complex nonlin_term = Complex(n_hat_out[index][0], n_hat_out[index][1]);
            Complex ikx = Complex(0.0, kx[i][j]);
            Complex result = -MU * ikx * nonlin_term;
            n_hat_out[index][0] = result.real();
            n_hat_out[index][1] = result.imag();
        }
    }

    fftw_free(u_real);
    fftw_free(u_hat_dealiased);
}

void perform_irk2_startup_step(fftw_complex *u_hat_0, fftw_complex *u_hat_1, const std::vector<std::vector<double>> &kx,
                               const std::vector<std::vector<double>> &ky, fftw_plan plan_bwd, fftw_plan plan_fwd)
{
    fftw_complex *u_hat_half_guess = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NUM_COEFFS);
    fftw_complex *u_hat_half_new = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NUM_COEFFS);
    fftw_complex *n_hat_half = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NUM_COEFFS);
    memcpy(u_hat_half_guess, u_hat_0, sizeof(fftw_complex) * NUM_COEFFS);
    int iter = 0;
    const int max_iter = 100;
    const double tol = 1e-10;

    std::cout << "  Iniciando iteracao de ponto fixo..." << std::endl;
    while (iter < max_iter)
    {
        calculate_nonlinear_rhs_fourier(u_hat_half_guess, n_hat_half, kx, ky, plan_bwd, plan_fwd);
        for (int i = 0; i < Ny; ++i)
        {
            for (int j = 0; j < Nx / 2 + 1; ++j)
            {
                int index = i * (Nx / 2 + 1) + j;
                double kx_val = kx[i][j];
                double ky_val = ky[i][j];
                Complex L_op = Complex(0.0, kx_val * (kx_val * kx_val + ky_val * ky_val));
                Complex u0 = Complex(u_hat_0[index][0], u_hat_0[index][1]);
                Complex n_half = Complex(n_hat_half[index][0], n_hat_half[index][1]);
                Complex numerator = u0 + (DT / 2.0) * n_half;
                Complex denominator = 1.0 + (DT / 2.0) * L_op;
                Complex u_new = numerator / denominator;
                u_hat_half_new[index][0] = u_new.real();
                u_hat_half_new[index][1] = u_new.imag();
            }
        }

        // *** BUG CORRIGIDO AQUI: Calcular o erro no espaço de Fourier ***
        double error = 0.0;
        for (int i = 0; i < NUM_COEFFS; ++i)
        {
            double diff_real = u_hat_half_new[i][0] - u_hat_half_guess[i][0];
            double diff_imag = u_hat_half_new[i][1] - u_hat_half_guess[i][1];
            error += diff_real * diff_real + diff_imag * diff_imag;
        }
        error = sqrt(error / (double)NUM_COEFFS); // Erro médio

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
                              const std::vector<std::vector<double>> &kx, const std::vector<std::vector<double>> &ky,
                              fftw_plan plan_bwd, fftw_plan plan_fwd)
{
    fftw_complex *n_hat = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NUM_COEFFS);
    calculate_nonlinear_rhs_fourier(u_hat_atual, n_hat, kx, ky, plan_bwd, plan_fwd);
    for (int i = 0; i < Ny; ++i)
    {
        for (int j = 0; j < Nx / 2 + 1; ++j)
        {
            int index = i * (Nx / 2 + 1) + j;
            double kx_val = kx[i][j];
            double ky_val = ky[i][j];
            Complex L_op = Complex(0.0, kx_val * (kx_val * kx_val + ky_val * ky_val));
            Complex u_antigo = Complex(u_hat_antigo[index][0], u_hat_antigo[index][1]);
            Complex n_atual = Complex(n_hat[index][0], n_hat[index][1]);
            Complex numerator = u_antigo * (1.0 / (2.0 * DT) - 0.5 * L_op) + n_atual;
            Complex denominator = (1.0 / (2.0 * DT) + 0.5 * L_op);
            Complex u_novo = numerator / denominator;
            u_hat_novo[index][0] = u_novo.real();
            u_hat_novo[index][1] = u_novo.imag();
        }
    }
    fftw_free(n_hat);
}

// Calcula a massa (norma L2 ao quadrado) para a solução 2D, com a normalização correta
double mass_conservation_2d(const double *u_real)
{
    double sum = 0.0;
    // A normalização é necessária porque a IFFT da FFTW não é normalizada.
    double norm_factor = 1.0 / (double)(Nx * Ny);

    // O loop pode ser simplificado para percorrer o array 1D diretamente
    for (int i = 0; i < Ny * Nx; ++i)
    {
        // Normaliza o valor ANTES de elevá-lo ao quadrado
        double u_val = u_real[i] * norm_factor;
        sum += u_val * u_val;
    }

    // Multiplica pelo elemento de área para completar a integral
    return sum * DX * DY;
}

// Calcula a derivada parcial em X usando o método espectral
void calculate_x_derivative_fourier(const fftw_complex *u_hat, double *u_x_real,
                                    const std::vector<std::vector<double>> &kx, fftw_plan plan_bwd)
{

    fftw_complex *ux_hat = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NUM_COEFFS);

    // Multiplicar por ikx no espaço de Fourier
    for (int i = 0; i < Ny; ++i)
    {
        for (int j = 0; j < Nx / 2 + 1; ++j)
        {
            int index = i * (Nx / 2 + 1) + j;
            Complex u_hat_comp(u_hat[index][0], u_hat[index][1]);
            Complex ikx(0.0, kx[i][j]);
            Complex ux_hat_comp = ikx * u_hat_comp;
            ux_hat[index][0] = ux_hat_comp.real();
            ux_hat[index][1] = ux_hat_comp.imag();
        }
    }

    // Transformar de volta para o espaço real
    fftw_execute_dft_c2r(plan_bwd, ux_hat, u_x_real);
    fftw_free(ux_hat);
}

// Calcula a derivada parcial em Y usando o método espectral
void calculate_y_derivative_fourier(const fftw_complex *u_hat, double *u_y_real,
                                    const std::vector<std::vector<double>> &ky, fftw_plan plan_bwd)
{

    fftw_complex *uy_hat = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NUM_COEFFS);

    // Multiplicar por iky no espaço de Fourier
    for (int i = 0; i < Ny; ++i)
    {
        for (int j = 0; j < Nx / 2 + 1; ++j)
        {
            int index = i * (Nx / 2 + 1) + j;
            Complex u_hat_comp(u_hat[index][0], u_hat[index][1]);
            Complex iky(0.0, ky[i][j]);
            Complex uy_hat_comp = iky * u_hat_comp;
            uy_hat[index][0] = uy_hat_comp.real();
            uy_hat[index][1] = uy_hat_comp.imag();
        }
    }

    // Transformar de volta para o espaço real
    fftw_execute_dft_c2r(plan_bwd, uy_hat, u_y_real);
    fftw_free(uy_hat);
}

double energy_conservation_2d(const double *u_real, const double *u_x_real, const double *u_y_real)
{
    double sum = 0.0;
    double norm_factor = 1.0 / (double)(Nx * Ny);

    for (int i = 0; i < Ny * Nx; ++i)
    {
        // Normaliza os valores (pois eles vêm da IFFT da FFTW)
        double u_val = u_real[i] * norm_factor;
        double ux_val = u_x_real[i] * norm_factor;
        double uy_val = u_y_real[i] * norm_factor;

        sum += 0.5 * (ux_val * ux_val + uy_val * uy_val) + (MU / (P + 1.0)) * pow(std::abs(u_val), P + 1.0);
    }
    return sum * DX * DY;
}
// double mass_motion_integral(const std::vector<double> &u, const std::vector<double> &u_x)
// {
//     double sum = 0.0;
//     for (int i = 0; i < N; i++)
//     {
//         // Correção: Adicionar parênteses em (P + 1.0)
//         sum += 3.0 * pow(u_x[i], 2.0) + (2.0 * MU * P / (P + 1.0)) * pow(u[i], P + 1.0);
//     }
//     return sum * DX * (-1.0) * (1.0 / mass_conservation(u));
// }
// double energy_motion_integral(const std::vector<double> &u, const std::vector<double> &u_x,
//                               const std::vector<double> &u_xx)
// {
//     double sum = 0.0;
//     for (int i = 0; i < N; i++)
//     {
//         // Correção: O termo deve ser 2*mu*u^(p-1)*ux^2
//         sum += 1.5 * pow(u_xx[i], 2.0) + 2.0 * MU * pow(u[i], P - 1.0) * pow(u_x[i], 2.0) +
//                (MU * MU * 0.5) * pow(u[i], 2.0 * P);
//     }
//     return sum * DX * (-1.0) * (1.0 / energy_conservation(u, u_x));
// }
//
// double mass_center(const std::vector<double> &u, const std::vector<double> &x_values)
// {
//     double sum = 0.0;
//     for (int i = 0; i < N; i++)
//     {
//         sum += x_values[i] * (pow(u[i], 2.0));
//     }
//     return sum * DX * (1.0 / mass_conservation(u));
// }
// double energy_center(const std::vector<double> &u, const std::vector<double> &u_x, const std::vector<double>
// &x_values)
// {
//     double sum = 0.0;
//     for (int i = 0; i < N; i++)
//     {
//         sum += x_values[i] * (0.5 * u_x[i] * u_x[i] + (MU / (P + 1.0)) * pow(u[i], P + 1.0));
//     }
//     return sum * DX * (1.0 / mass_conservation(u));
// }
//
// =============================================================================
//               *** NOVAS FUNÇÕES DE DIAGNÓSTICO 2D ***
// =============================================================================

// Calcula a componente X do centro de massa
double mass_center_x_2d(const double *u_real)
{
    double sum = 0.0;
    double norm_factor = 1.0 / (double)(Nx * Ny);

    for (int i = 0; i < Ny; ++i)
    {
        for (int j = 0; j < Nx; ++j)
        {
            double x = -Lx / 2.0 + j * DX;
            int index = i * Nx + j;

            double u_val = u_real[index] * norm_factor;
            double mass_density = u_val * u_val;

            sum += x * mass_density;
        }
    }
    // Normaliza pelo elemento de área e pela massa total
    return sum * DX * DY / mass_conservation_2d(u_real);
}

// Calcula a componente Y do centro de massa
double mass_center_y_2d(const double *u_real)
{
    double sum = 0.0;
    double norm_factor = 1.0 / (double)(Nx * Ny);

    for (int i = 0; i < Ny; ++i)
    {
        double y = -Ly / 2.0 + i * DY;
        for (int j = 0; j < Nx; ++j)
        {
            int index = i * Nx + j;

            double u_val = u_real[index] * norm_factor;
            double mass_density = u_val * u_val;

            sum += y * mass_density;
        }
    }
    return sum * DX * DY / mass_conservation_2d(u_real);
}

// Calcula a componente X do centro de energia
double energy_center_x_2d(const double *u_real, const double *u_x_real, const double *u_y_real)
{
    double sum = 0.0;
    double norm_factor = 1.0 / (double)(Nx * Ny);

    for (int i = 0; i < Ny; ++i)
    {
        for (int j = 0; j < Nx; ++j)
        {
            double x = -Lx / 2.0 + j * DX;
            int index = i * Nx + j;

            double u_val = u_real[i * Nx + j] * norm_factor;
            double ux_val = u_x_real[i * Nx + j] * norm_factor;
            double uy_val = u_y_real[i * Nx + j] * norm_factor;

            double energy_density =
                0.5 * (ux_val * ux_val + uy_val * uy_val) + (MU / (P + 1.0)) * pow(std::abs(u_val), P + 1.0);

            sum += x * energy_density;
        }
    }
    return sum * DX * DY / energy_conservation_2d(u_real, u_x_real, u_y_real);
}

// Calcula a componente Y do centro de energia
double energy_center_y_2d(const double *u_real, const double *u_x_real, const double *u_y_real)
{
    double sum = 0.0;
    double norm_factor = 1.0 / (double)(Nx * Ny);

    for (int i = 0; i < Ny; ++i)
    {
        double y = -Ly / 2.0 + i * DY;
        for (int j = 0; j < Nx; ++j)
        {
            int index = i * Nx + j;

            double u_val = u_real[i * Nx + j] * norm_factor;
            double ux_val = u_x_real[i * Nx + j] * norm_factor;
            double uy_val = u_y_real[i * Nx + j] * norm_factor;

            double energy_density =
                0.5 * (ux_val * ux_val + uy_val * uy_val) + (MU / (P + 1.0)) * pow(std::abs(u_val), P + 1.0);

            sum += y * energy_density;
        }
    }
    return sum * DX * DY / energy_conservation_2d(u_real, u_x_real, u_y_real);
}
