#include <cmath>
#include <complex>
#include <fftw3.h>
#include <fstream>
#include <iostream>
#include <vector>

using Complex = std::complex<double>;

// --- Parâmetros da Simulação Refinada ---
const int N = 2048;     // Podemos usar mais pontos agora!
const double L = 200.0; // Domínio de -100 a 100
const double T_FINAL = 25.0;
const double DT = 0.005; // DT muito maior graças ao Fator de Integração!

// --- Parâmetros da KdV (de Tao) ---
const double P = 3.0;
const double MU = -1.0;

// --- Constantes Derivadas ---
const double DX = L / (N - 1);

// Função para salvar a solução (x, u)
void save_solution(const std::vector<double> &x, const std::vector<double> &u, int step)
{
    char filename[100];
    sprintf(filename, "solution_%04d.txt", step);
    std::ofstream file(filename);
    for (int i = 0; i < N; ++i)
    {
        file << x[i] << " " << u[i] << std::endl;
    }
    file.close();
}

// Calcula o termo não-linear N(u) = -mu * d/dx(u^p) no espaço de Fourier
void calculate_nonlinear_term_fourier(fftw_complex *u_hat, fftw_complex *n_hat, const std::vector<Complex> &k_vals,
                                      fftw_plan plan_bwd, fftw_plan plan_fwd)
{

    double *u_real = (double *)fftw_malloc(sizeof(double) * N);

    // De-aliasing (Regra dos 2/3)
    fftw_complex *u_hat_dealiased = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1));
    for (int i = 0; i < N / 2 + 1; ++i)
        u_hat_dealiased[i][0] = u_hat[i][0], u_hat_dealiased[i][1] = u_hat[i][1];
    int k_cutoff = N / 3;
    for (int i = k_cutoff; i < N / 2 + 1; ++i)
        u_hat_dealiased[i][0] = 0.0, u_hat_dealiased[i][1] = 0.0;

    // IFFT para obter u no espaço real
    fftw_execute_dft_c2r(plan_bwd, u_hat_dealiased, u_real);

    // Calcular u^p no espaço real
    for (int i = 0; i < N; ++i)
    {
        u_real[i] = pow(u_real[i] / N, P); // Normaliza a IFFT e eleva à potência
    }

    // FFT(u^p)
    fftw_execute_dft_r2c(plan_fwd, u_real, n_hat);

    // Multiplicar por -mu*ik para obter F(-mu * d/dx(u^p))
    for (int i = 0; i < N / 2 + 1; ++i)
    {
        Complex nonlin_term = Complex(n_hat[i][0], n_hat[i][1]);
        Complex result = -MU * k_vals[i] * nonlin_term;
        n_hat[i][0] = result.real();
        n_hat[i][1] = result.imag();
    }

    fftw_free(u_real);
    fftw_free(u_hat_dealiased);
}

int main()
{
    std::cout << "Iniciando simulacao KdV Refinada (FFTW + Fator de Integracao)..." << std::endl;

    // --- Alocação de Memória Otimizada para FFTW ---
    double *u_real = (double *)fftw_malloc(sizeof(double) * N);
    fftw_complex *u_hat = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1));
    fftw_complex *n_hat = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1));

    // --- Planos da FFTW ---
    fftw_plan plan_fwd = fftw_plan_dft_r2c_1d(N, u_real, u_hat, FFTW_ESTIMATE);
    fftw_plan plan_bwd = fftw_plan_dft_c2r_1d(N, u_hat, u_real, FFTW_ESTIMATE);

    // --- Domínio e Números de Onda ---
    std::vector<double> x(N);
    for (int i = 0; i < N; i++)
        x[i] = -L / 2.0 + i * DX;

    std::vector<Complex> k_vals(N / 2 + 1);
    for (int i = 0; i < N / 2 + 1; i++)
        k_vals[i] = Complex(0.0, 2.0 * M_PI * i / L);

    // --- NOVO: Fator de Integração Pré-calculado ---
    std::vector<Complex> integrating_factor(N / 2 + 1);
    for (int i = 0; i < N / 2 + 1; ++i)
    {
        Complex ik3 = k_vals[i] * k_vals[i] * k_vals[i];
        integrating_factor[i] = exp(-ik3 * DT);
    }

    // --- NOVO: Condição Inicial com Dois Sólitons ---
    double c1 = 16.0, x1 = -40.0;
    double c2 = 4.0, x2 = -10.0;
    for (int i = 0; i < N; i++)
    {
        double sech_arg1 = 0.5 * sqrt(c1) * (x[i] - x1);
        double soliton1 = 0.5 * c1 / (cosh(sech_arg1) * cosh(sech_arg1));
        double sech_arg2 = 0.5 * sqrt(c2) * (x[i] - x2);
        double soliton2 = 0.5 * c2 / (cosh(sech_arg2) * cosh(sech_arg2));
        u_real[i] = soliton1 + soliton2;
    }

    // Transformada inicial
    fftw_execute_dft_r2c(plan_fwd, u_real, u_hat);

    // --- Loop Temporal com Fator de Integração ---
    int num_steps = T_FINAL / DT;
    for (int i = 0; i <= num_steps; ++i)
    {
        if (i % 10 == 0)
        { // Salvar com mais frequência, pois DT é maior
            fftw_execute_dft_c2r(plan_bwd, u_hat, u_real);
            std::vector<double> u_current(N);
            for (int j = 0; j < N; ++j)
                u_current[j] = u_real[j] / N; // Normalizar
            save_solution(x, u_current, i / 10);
            std::cout << "Passo " << i << "/" << num_steps << ", Tempo = " << i * DT << std::endl;
        }

        // Calcular o termo não-linear N(u_n)
        calculate_nonlinear_term_fourier(u_hat, n_hat, k_vals, plan_bwd, plan_fwd);

        // Aplicar a fórmula do Fator de Integração de 1ª ordem
        for (int j = 0; j < N / 2 + 1; ++j)
        {
            Complex u_hat_comp(u_hat[j][0], u_hat[j][1]);
            Complex n_hat_comp(n_hat[j][0], n_hat[j][1]);

            Complex updated_u_hat = integrating_factor[j] * (u_hat_comp + DT * n_hat_comp);

            u_hat[j][0] = updated_u_hat.real();
            u_hat[j][1] = updated_u_hat.imag();
        }
    }

    // --- Limpeza ---
    fftw_destroy_plan(plan_fwd);
    fftw_destroy_plan(plan_bwd);
    fftw_free(u_real);
    fftw_free(u_hat);
    fftw_free(n_hat);

    std::cout << "Simulacao refinada concluida." << std::endl;
    return 0;
}
