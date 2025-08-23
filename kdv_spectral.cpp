#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <vector>

using Complex = std::complex<double>;

// --- Parâmetros da Simulação ---
const int N = 512;
const double L = 100.0;
const double T_FINAL = 2.0;
const double DT = 0.0001;

// --- Parâmetros da KdV (de Tao) ---
const double P = 3.0;
const double MU = -1.0;

// --- Parâmetros da Condição Inicial ---
const double SOLITON_AMP = 12.0;
const double SOLITON_X0 = 25.0;

// --- Constantes Derivadas ---
const double DX = L / N;
const int N_EXT = 2 * N;

void discretize_axis(std::vector<double> &x_axis)
{
    x_axis.resize(N);
    x_axis[0] = 0.0; // Domínio começa em 0 para a Transformada Seno
    for (int i = 1; i < N; i++)
    {
        x_axis[i] = x_axis[i - 1] + DX;
    }
}

void general_initial_conditions(std::vector<double> &ic, double c, double x0)
{
    std::vector<double> x_coords = ic; // Salva as coordenadas
    for (int i = 0; i < N; i++)
    {
        double sech_arg = 0.5 * sqrt(c) * (x_coords[i] - x0);
        ic[i] = 0.5 * c / (cosh(sech_arg) * cosh(sech_arg));
    }
}

// DFT: Transformada Discreta de Fourier
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

// IDFT: Transformada Discreta de Fourier Inversa
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

void calculate_derivative_fourier(const std::vector<double> &u, std::vector<double> &u_x,
                                  const std::vector<Complex> &k_vals);
double mass_conservation(const std::vector<double> &u);
double energy_conservation(const std::vector<double> &u, const std::vector<double> &u_x);

// Calcula o lado direito (RHS) da KdV no espaço de Fourier
void kdv_rhs_fourier(const std::vector<Complex> &u_hat_in, std::vector<Complex> &rhs_hat_out,
                     const std::vector<Complex> &k_vals)
{

    std::vector<Complex> u_hat_dealiased = u_hat_in;

    // Filtro de frequências:
    //  Zera o terço superior das frequências para evitar poluição espectral.
    int k_cutoff = N_EXT / 3;
    for (int i = k_cutoff; i < N_EXT - k_cutoff; ++i)
    {
        u_hat_dealiased[i] = Complex(0.0, 0.0);
    }

    std::vector<double> u_real_ext;
    idft(u_hat_dealiased, u_real_ext); // Usa a versão sem aliasing

    // Calcular o termo não linear u^p no espaço real
    std::vector<double> nonlin_term_real(N_EXT);
    for (int i = 0; i < N_EXT; ++i)
    {
        nonlin_term_real[i] = pow(u_real_ext[i], P);
    }

    std::vector<Complex> nonlin_term_hat;
    dft(nonlin_term_real, nonlin_term_hat);

    // Combinar termos
    rhs_hat_out.resize(N_EXT);
    for (int i = 0; i < N_EXT; ++i)
    {
        Complex ik = k_vals[i];
        Complex ik3 = ik * ik * ik;

        Complex linear_part = -ik3 * u_hat_in[i]; // O termo linear usa o u_hat original
        Complex nonlin_part = -MU * ik * nonlin_term_hat[i];

        rhs_hat_out[i] = linear_part + nonlin_part;
    }
}

int main()
{
    std::cout << "Iniciando simulacao KdV com DFT" << std::endl;

    // --- Setup dos vetores de números de onda (ik)
    std::vector<Complex> k_vals(N_EXT);
    for (int i = 0; i < N_EXT / 2; ++i)
    {
        k_vals[i] = Complex(0.0, 2.0 * M_PI * i / (2.0 * L));
    }
    for (int i = N_EXT / 2; i < N_EXT; ++i)
    {
        k_vals[i] = Complex(0.0, 2.0 * M_PI * (i - N_EXT) / (2.0 * L));
    }

    // --- CONDIÇÃO INICIAL
    std::vector<double> u_initial;
    discretize_axis(u_initial);
    general_initial_conditions(u_initial, SOLITON_AMP, SOLITON_X0);

    // --- Extensão Ímpar
    std::vector<double> u_ext_real(N_EXT, 0.0);
    for (int i = 0; i < N; ++i)
        u_ext_real[i] = u_initial[i];
    for (int i = 1; i < N; ++i)
        u_ext_real[N_EXT - i] = -u_initial[i];

    std::vector<Complex> u_hat;
    dft(u_ext_real, u_hat);

    // Arquivos de Saída
    std::ofstream kdv_file("kdv_data.txt");
    std::ofstream mass_file("mass_data.txt");
    std::ofstream energy_file("energy_data.txt");

    // Loop Temporal (RK4)
    int num_steps = T_FINAL / DT;
    std::vector<Complex> k1(N_EXT), k2(N_EXT), k3(N_EXT), k4(N_EXT), u_temp_hat(N_EXT);

    for (int i = 0; i <= num_steps; ++i)
    {
        if (i % 100 == 0)
        {
            idft(u_hat, u_ext_real);
            std::vector<double> u_current(N);
            for (int j = 0; j < N; ++j)
                u_current[j] = u_ext_real[j];

            for (int k = 0; k < N; ++k)
                kdv_file << u_current[k] << "\n";
            kdv_file << "\n";

            std::vector<double> u_x;
            calculate_derivative_fourier(u_current, u_x, k_vals);
            mass_file << std::scientific << mass_conservation(u_current) << "\n";
            energy_file << std::scientific << energy_conservation(u_current, u_x) << "\n";

            std::cout << "Passo " << i << "/" << num_steps << ", Tempo = " << i * DT << std::endl;
        }

        // RK4
        kdv_rhs_fourier(u_hat, k1, k_vals);
        for (int j = 0; j < N_EXT; ++j)
            u_temp_hat[j] = u_hat[j] + 0.5 * DT * k1[j];
        kdv_rhs_fourier(u_temp_hat, k2, k_vals);
        for (int j = 0; j < N_EXT; ++j)
            u_temp_hat[j] = u_hat[j] + 0.5 * DT * k2[j];
        kdv_rhs_fourier(u_temp_hat, k3, k_vals);
        for (int j = 0; j < N_EXT; ++j)
            u_temp_hat[j] = u_hat[j] + DT * k3[j];
        kdv_rhs_fourier(u_temp_hat, k4, k_vals);
        for (int j = 0; j < N_EXT; ++j)
        {
            u_hat[j] += (DT / 6.0) * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
        }
    }

    kdv_file.close();
    mass_file.close();
    energy_file.close();
    std::cout << "Simulacao concluida." << std::endl;

    return 0;
}

void calculate_derivative_fourier(const std::vector<double> &u, std::vector<double> &u_x,
                                  const std::vector<Complex> &k_vals)
{
    std::vector<Complex> u_hat(N_EXT), ux_hat(N_EXT);
    std::vector<double> u_ext(N_EXT, 0.0);
    for (int i = 0; i < N; ++i)
        u_ext[i] = u[i];
    for (int i = 1; i < N; ++i)
        u_ext[N_EXT - i] = -u[i];
    dft(u_ext, u_hat);
    for (int i = 0; i < N_EXT; ++i)
        ux_hat[i] = k_vals[i] * u_hat[i];
    std::vector<double> ux_ext;
    idft(ux_hat, ux_ext);
    u_x.resize(N);
    for (int i = 0; i < N; ++i)
        u_x[i] = ux_ext[i];
}

double mass_conservation(const std::vector<double> &u)
{
    double sum = 0.0;
    for (int i = 0; i < N; ++i)
        sum += u[i] * u[i];
    return sum * DX;
}

double energy_conservation(const std::vector<double> &u, const std::vector<double> &u_x)
{
    double sum = 0.0;
    for (int i = 0; i < N; ++i)
    {
        sum += 0.5 * u_x[i] * u_x[i] - (MU / (P + 1.0)) * pow(u[i], P + 1.0);
    }
    return sum * DX;
}
