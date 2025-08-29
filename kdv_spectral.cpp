#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <vector>

using Complex = std::complex<double>;

// --- Parâmetros da Simulação ---
const int N = 1024;
const double L = 100.0;
const double T_FINAL = 1.0;
const double DT = 0.0001;

// --- Parâmetros da KdV (de Tao) ---
const double P = 3.0;
const double MU = 1.0;

// --- Parâmetros da Condição Inicial ---
const double SOLITON_AMP = 6.0;
const double SOLITON_X0 = 40.0;

// --- Constantes Derivadas ---
const double DX = L / N;

// Funções de DFT/IDFT
void dft(const std::vector<double> &in, std::vector<Complex> &out);
void idft(const std::vector<Complex> &in, std::vector<double> &out);

double mass_conservation(const std::vector<double> &u);
double energy_conservation(const std::vector<double> &u, const std::vector<double> &u_x);
void calculate_derivative_fourier(const std::vector<double> &u, std::vector<double> &u_x,
                                  const std::vector<Complex> &k_vals)
{
    std::vector<Complex> u_hat(N);
    std::vector<Complex> ux_hat(N);

    // Passo 1: Transformar u para o espaço de Fourier
    dft(u, u_hat);

    // Passo 2: Multiplicar por ik no espaço de Fourier
    for (int i = 0; i < N; ++i)
    {
        ux_hat[i] = k_vals[i] * u_hat[i];
    }

    // Passo 3: Transformar de volta para o espaço real para obter u_x
    idft(ux_hat, u_x);
}

// Calcula o RHS
void kdv_rhs_fourier(const std::vector<Complex> &u_hat_in, std::vector<Complex> &rhs_hat_out,
                     const std::vector<Complex> &k_vals)
{
    std::vector<Complex> u_hat_truncated = u_hat_in;

    // --- PASSO 1: truncar
    int k_cutoff = N / 3;
    for (int i = k_cutoff; i < N - k_cutoff; ++i)
    {
        u_hat_truncated[i] = Complex(0.0, 0.0);
    }

    // --- PASSO 2: Calcular o termo não-linear a partir do espectro truncado
    std::vector<double> u_real;
    idft(u_hat_truncated, u_real);

    std::vector<double> nonlin_term_real(N);
    for (int i = 0; i < N; ++i)
    {
        nonlin_term_real[i] = pow(u_real[i], P);
    }

    std::vector<Complex> nonlin_term_hat;
    dft(nonlin_term_real, nonlin_term_hat);

    // --- PASSO 3:
    rhs_hat_out.resize(N);
    for (int i = 0; i < N; ++i)
    {
        Complex ik = k_vals[i];
        Complex ik3 = ik * ik * ik;

        Complex linear_part = -ik3 * u_hat_truncated[i];

        Complex nonlin_part = -MU * ik * nonlin_term_hat[i];

        rhs_hat_out[i] = linear_part + nonlin_part;
    }
}

int main()
{
    std::cout << "Iniciando simulacao KdV" << std::endl;

    // Setup dos vetores de números de onda (ik) para N pontos
    std::vector<Complex> k_vals(N);
    for (int i = 0; i < N / 2; ++i)
    {
        k_vals[i] = Complex(0.0, 2.0 * M_PI * i / L);
    }
    for (int i = N / 2; i < N; ++i)
    {
        k_vals[i] = Complex(0.0, 2.0 * M_PI * (i - N) / L);
    }

    // Condição Inicial
    std::vector<double> u_current(N);
    for (int i = 0; i < N; i++)
    {
        double x = i * DX;
        double sech_arg = 0.5 * sqrt(SOLITON_AMP) * (x - SOLITON_X0);
        u_current[i] = 0.5 * SOLITON_AMP / (cosh(sech_arg) * cosh(sech_arg));
    }
    u_current[0] = 0.0;
    u_current[N - 1] = 0.0;

    std::vector<Complex> u_hat;
    dft(u_current, u_hat);

    std::ofstream kdv_file("kdv_data.txt");
    std::ofstream mass_file("mass_data.txt");
    std::ofstream energy_file("energy_data.txt");

    // Loop Temporal (RK4)
    int num_steps = T_FINAL / DT;
    std::vector<Complex> k1(N), k2(N), k3(N), k4(N), u_temp_hat(N);

    for (int i = 0; i <= num_steps; ++i)
    {
        if (i % 100 == 0)
        {
            idft(u_hat, u_current);
            u_current[0] = 0.0;
            u_current[N - 1] = 0.0;

            for (int k = 0; k < N; ++k)
                kdv_file << u_current[k] << "\n";
            kdv_file << "\n";

            mass_file << std::scientific << mass_conservation(u_current) << "\n";
            std::cout << "massa: " << std::scientific << mass_conservation(u_current) << "\n";
            std::vector<double> u_x(N);
            calculate_derivative_fourier(u_current, u_x, k_vals);
            energy_file << std::scientific << energy_conservation(u_current, u_x) << "\n";
            std::cout << "energia: " << std::scientific << energy_conservation(u_current, u_x) << "\n";

            std::cout << "Passo " << i << "/" << num_steps << ", Tempo = " << i * DT << std::endl;
        }

        dft(u_current, u_hat); // Esta linha pode ser removida se a de baixo for descomentada

        // Com esta filosofia, o truncamento deve acontecer no final de cada passo completo
        // Descomente a linha abaixo e remova a DFT acima para um truncamento mais 'puro'
        // for(int j=N/3; j < N - N/3; ++j) u_hat[j] = Complex(0.0, 0.0);

        kdv_rhs_fourier(u_hat, k1, k_vals);
        for (int j = 0; j < N; ++j)
            u_temp_hat[j] = u_hat[j] + 0.5 * DT * k1[j];
        kdv_rhs_fourier(u_temp_hat, k2, k_vals);
        for (int j = 0; j < N; ++j)
            u_temp_hat[j] = u_hat[j] + 0.5 * DT * k2[j];
        kdv_rhs_fourier(u_temp_hat, k3, k_vals);
        for (int j = 0; j < N; ++j)
            u_temp_hat[j] = u_hat[j] + DT * k3[j];
        kdv_rhs_fourier(u_temp_hat, k4, k_vals);
        for (int j = 0; j < N; ++j)
        {
            u_hat[j] += (DT / 6.0) * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
        }
    }

    kdv_file.close();
    mass_file.close();
    std::cout << "Simulacao concluida." << std::endl;
    return 0;
}

// Implementações das funções restantes
void dft(const std::vector<double> &in, std::vector<Complex> &out)
{
    int size = in.size();
    out.resize(size);
    for (int k = 0; k < size; ++k)
    {
        out[k] = Complex(0.0, 0.0);
        for (int n = 0; n < size; ++n)
        {
            double angle = -2.0 * M_PI * k * n / size;
            out[k] += in[n] * Complex(cos(angle), sin(angle));
        }
    }
}
void idft(const std::vector<Complex> &in, std::vector<double> &out)
{
    int size = in.size();
    out.resize(size);
    std::vector<Complex> temp_out(size);
    for (int n = 0; n < size; ++n)
    {
        temp_out[n] = Complex(0.0, 0.0);
        for (int k = 0; k < size; ++k)
        {
            double angle = 2.0 * M_PI * k * n / size;
            temp_out[n] += in[k] * Complex(cos(angle), sin(angle));
        }
        out[n] = temp_out[n].real() / size;
    }
}
double mass_conservation(const std::vector<double> &u)
{
    double sum = 0.0;
    for (int i = 0; i < N; ++i)
        sum += u[i] * u[i];
    return sum * DX;
}

// Definição de energy_conservation que estava faltando no seu arquivo original
double energy_conservation(const std::vector<double> &u, const std::vector<double> &u_x)
{
    double sum = 0.0;
    for (int i = 0; i < N; i++)
    {
        // E(u) = integral( 1/2 * (u_x)^2 - mu/(p+1) * u^(p+1) )
        // Note que o sinal de MU é importante aqui.
        sum += 0.5 * u_x[i] * u_x[i] - (MU / (P + 1.0)) * pow(u[i], P + 1.0);
    }
    return sum * DX;
}
