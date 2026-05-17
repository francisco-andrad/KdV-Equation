#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

using Complex = std::complex<double>;

// --- Parâmetros da Simulação ---
const int Nx = 512;
const int Ny = 512;
const double Lx = 100.0;
const double Ly = 100.0;
const double T_FINAL = 10.0;
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
void dft_2d(const std::vector<std::vector<double>> &in, std::vector<std::vector<Complex>> &out);
void idft_1d_complex(const std::vector<Complex> &in, std::vector<Complex> &out);
void idft_2d(const std::vector<std::vector<Complex>> &in, std::vector<std::vector<double>> &out);

double mass_conservation(const std::vector<double> &u);
double energy_conservation(const std::vector<double> &u, const std::vector<double> &u_x);
double mass_motion_integral(const std::vector<double> &u, const std::vector<double> &u_x);
double energy_motion_integral(const std::vector<double> &u, const std::vector<double> &u_x,
                              const std::vector<double> &u_xx);
double mass_center(const std::vector<double> &u, const std::vector<double> &x_values);
double energy_center(const std::vector<double> &u, const std::vector<double> &u_x, const std::vector<double> &x_values);
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

void calculate_second_derivative_fourier(const std::vector<double> &u, std::vector<double> &u_xx,
                                         const std::vector<Complex> &k_vals)
{
    std::vector<Complex> u_hat(N);
    std::vector<Complex> uxx_hat(N);

    // Passo 1: Transformar u para o espaço de Fourier
    dft(u, u_hat);

    // Passo 2: Multiplicar por ik no espaço de Fourier
    for (int i = 0; i < N; ++i)
    {
        uxx_hat[i] = k_vals[i] * k_vals[i] * u_hat[i];
    }

    // Passo 3: Transformar de volta para o espaço real para obter u_x
    idft(uxx_hat, u_xx);
}

// Calcula o RHS
void kdv_rhs_fourier_2d(const std::vector<std::vector<Complex>> &u_hat_in,
                        std::vector<std::vector<Complex>> &rhs_hat_out, const std::vector<std::vector<double>> &kx,
                        const std::vector<std::vector<double>> &ky)
{
    int Ny = u_hat_in.size();
    int Nx = u_hat_in[0].size();

    std::vector<std::vector<Complex>> u_hat_truncated = u_hat_in;

    // --- PASSO 1: Truncamento (De-aliasing) 2D ---
    // A regra dos 2/3 é aplicada em cada direção independentemente.
    int k_cutoff_x = Nx / 3;
    int k_cutoff_y = Ny / 3;
    for (int i = 0; i < Ny; ++i)
    {
        for (int j = 0; j < Nx; ++j)
        {
            // Zera se a frequência em X OU em Y for muito alta
            if (std::abs(kx[i][j]) > (2.0 * M_PI * k_cutoff_x / Lx) ||
                std::abs(ky[i][j]) > (2.0 * M_PI * k_cutoff_y / Ly))
            {
                u_hat_truncated[i][j] = Complex(0.0, 0.0);
            }
        }
    }

    // --- PASSO 2: Calcular o termo não-linear a partir do espectro truncado ---
    // Usar vetores e funções 2D
    std::vector<std::vector<double>> u_real;
    idft_2d(u_hat_truncated, u_real);

    std::vector<std::vector<double>> nonlin_term_real(Ny, std::vector<double>(Nx));
    for (int i = 0; i < Ny; ++i)
    {
        for (int j = 0; j < Nx; ++j)
        {
            nonlin_term_real[i][j] = pow(u_real[i][j], P);
        }
    }

    std::vector<std::vector<Complex>> nonlin_term_hat;
    dft_2d(nonlin_term_real, nonlin_term_hat);

    // --- PASSO 3: Combinar os termos Linear e Não-Linear ---
    rhs_hat_out.assign(Ny, std::vector<Complex>(Nx));
    for (int i = 0; i < Ny; ++i)
    {
        for (int j = 0; j < Nx; ++j)
        {
            double kx_val = kx[i][j];
            double ky_val = ky[i][j];

            // Termo Linear: ikx * (kx^2 + ky^2) * u_hat
            Complex linear_factor = Complex(0.0, kx_val * (kx_val * kx_val + ky_val * ky_val));
            Complex linear_part = linear_factor * u_hat_truncated[i][j];

            // Termo Não-Linear: -MU * ikx * F(u^p)
            Complex nonlin_factor = Complex(0.0, -MU * kx_val);
            Complex nonlin_part = nonlin_factor * nonlin_term_hat[i][j];

            // A equação é u_t = - (parte_linear) - (parte_nao_linear)
            rhs_hat_out[i][j] = -linear_part - nonlin_part;
        }
    }
}

int main()
{
    std::cout << "Iniciando simulacao ZK" << std::endl;

    std::vector<std::vector<double>> kx(Ny, std::vector<double>(Nx));
    std::vector<std::vector<double>> ky(Ny, std::vector<double>(Nx));

    // Loop sobre cada ponto (i, j) da grade de Fourier
    for (int i = 0; i < Ny; ++i)
    {
        for (int j = 0; j < Nx; ++j)
        {
            // Calcular kx baseado no índice da coluna (j)
            if (j < Nx / 2)
            {
                kx[i][j] = 2.0 * M_PI * j / Lx;
            }
            else
            {
                kx[i][j] = 2.0 * M_PI * (j - Nx) / Lx;
            }

            // Calcular ky baseado no índice da linha (i)
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

    std::vector<std::vector<double>> u_current(Ny, std::vector<double>(Nx));

    // Percorrer a grade 2D com loops aninhados
    for (int i = 0; i < Ny; ++i)
    { // Loop sobre as linhas (eixo y)
        // double y = -Ly / 2.0 + i * DY;

        for (int j = 0; j < Nx; ++j)
        { // Loop sobre as colunas (eixo x)
            double x = -Lx / 2.0 + j * DX;

            // A fórmula do sóliton depende apenas de x para criar a "linha"
            double sech_arg = 0.5 * sqrt(SOLITON_AMP) * (x - SOLITON_X0);
            u_current[i][j] = 0.5 * SOLITON_AMP / (cosh(sech_arg) * cosh(sech_arg));
        }
    }

    std::vector<std::vector<Complex>> u_hat;
    dft_2d(u_current, u_hat);

    std::ofstream kdv_file("kdv_data.txt");
    std::ofstream mass_file("mass_data.txt");
    std::ofstream energy_file("energy_data.txt");
    std::ofstream monotonicity_file("monotonicidade.txt");
    std::ofstream motion_file("movimento.txt");
    std::ofstream mass_centre_file("mass_centre.txt");
    std::ofstream energy_centre_file("energy_centre.txt");
    std::ofstream mass_error("mass_error.txt");
    std::ofstream energy_error("energy_error.txt");

    // Loop Temporal (RK4)
    int num_steps = T_FINAL / DT;
    std::vector<Complex> k1(N), k2(N), k3(N), k4(N), u_temp_hat(N);
    double e_i, m_i;
    m_i = mass_conservation(u_current);
    calculate_derivative_fourier(u_current, u_x, k_vals);
    e_i = energy_conservation(u_current, u_x);

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
            double monotone, e_error, m_error, motion, curr_mass, curr_en;
            curr_mass = mass_conservation(u_current);
            mass_file << std::scientific << curr_mass << "\n";
            std::cout << "massa: " << std::scientific << curr_mass << "\n";
            mass_error << std::scientific << std::abs(m_i - curr_mass) << "\n";
            std::vector<double> u_x(N), u_xx(N);
            calculate_derivative_fourier(u_current, u_x, k_vals);
            calculate_second_derivative_fourier(u_current, u_xx, k_vals);
            curr_en = energy_conservation(u_current, u_x);
            motion = mass_motion_integral(u_current, u_x) - energy_motion_integral(u_current, u_x, u_xx);
            monotone = energy_center(u_current, u_x, x_values) - mass_center(u_current, x_values);
            // std::cout << "movimento massa: " << std::scientific << mass_motion_integral(u_current, u_x) << "\n";
            // std::cout << "movimento energia: " << std::scientific << energy_motion_integral(u_current, u_x, u_xx);
            energy_file << std::scientific << curr_en << "\n";
            energy_error << std::scientific << std::abs(e_i - curr_en) << "\n";
            std::cout << "energia: " << std::scientific << curr_en << "\n";
            std::cout << "monotonicidade: " << std::scientific << monotone << "\n";
            motion_file << std::scientific << motion << "\n";

            monotonicity_file << std::scientific << monotone << "\n";
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
    energy_file.close();
    monotonicity_file.close();
    mass_centre_file.close();
    energy_centre_file.close();
    mass_error.close();
    energy_error.close();
    motion_file.close();
    std::cout << "Simulacao concluida." << std::endl;
    return 0;
}

// Função auxiliar para calcular a DFT 1D em um vetor
void dft_1d(const std::vector<Complex> &in, std::vector<Complex> &out)
{
    int N = in.size();
    out.resize(N);
    for (int k = 0; k < N; ++k)
    {
        out[k] = Complex(0.0, 0.0);
        for (int n = 0; n < N; ++n)
        {
            double angle = -2.0 * M_PI * k * n / N;
            out[k] += in[n] * Complex(cos(angle), sin(angle));
        }
    }
}

// Função principal para a DFT 2D
void dft_2d(const std::vector<std::vector<double>> &in, std::vector<std::vector<Complex>> &out)
{
    int N_rows = in.size();
    if (N_rows == 0)
        return;
    int N_cols = in[0].size();

    // Matriz temporária para armazenar o resultado da primeira etapa
    std::vector<std::vector<Complex>> temp_matrix(N_rows, std::vector<Complex>(N_cols));

    // --- ETAPA 1: Aplicar DFT 1D em cada LINHA ---
    for (int i = 0; i < N_rows; ++i)
    {
        // Converte a linha de double para complex para a função dft_1d
        std::vector<Complex> row_in(N_cols);
        for (int j = 0; j < N_cols; ++j)
        {
            row_in[j] = Complex(in[i][j], 0.0);
        }
        dft_1d(row_in, temp_matrix[i]);
    }

    // --- ETAPA 2: Aplicar DFT 1D em cada COLUNA da matriz temporária ---
    out.assign(N_rows, std::vector<Complex>(N_cols));
    for (int j = 0; j < N_cols; ++j)
    {
        // Extrai a coluna da matriz temporária
        std::vector<Complex> col_in(N_rows);
        for (int i = 0; i < N_rows; ++i)
        {
            col_in[i] = temp_matrix[i][j];
        }

        // Calcula a DFT da coluna
        std::vector<Complex> col_out;
        dft_1d(col_in, col_out);

        // Armazena o resultado na coluna correspondente da matriz de saída
        for (int i = 0; i < N_rows; ++i)
        {
            out[i][j] = col_out[i];
        }
    }
}

void idft_1d_complex(const std::vector<Complex> &in, std::vector<Complex> &out)
{
    int N = in.size();
    out.resize(N);
    for (int n = 0; n < N; ++n)
    {
        out[n] = Complex(0.0, 0.0);
        for (int k = 0; k < N; ++k)
        {
            double angle = 2.0 * M_PI * k * n / N;
            out[n] += in[k] * Complex(cos(angle), sin(angle));
        }
    }
}

// Função principal e CORRIGIDA para a IDFT 2D
// A saída agora é uma matriz de 'double'
void idft_2d(const std::vector<std::vector<Complex>> &in, std::vector<std::vector<double>> &out)
{
    int N_rows = in.size();
    if (N_rows == 0)
        return;
    int N_cols = in[0].size();

    // Matriz temporária para armazenar o resultado da primeira etapa (ainda complexa)
    std::vector<std::vector<Complex>> temp_matrix(N_rows, std::vector<Complex>(N_cols));

    // --- ETAPA 1: Aplicar IDFT 1D (Complex -> Complex) em cada LINHA ---
    for (int i = 0; i < N_rows; ++i)
    {
        // A entrada da linha 'i' já é complexa, não precisa converter
        idft_1d_complex(in[i], temp_matrix[i]);
    }

    // Matriz temporária final para o resultado completo (ainda complexa)
    std::vector<std::vector<Complex>> final_complex_matrix(N_rows, std::vector<Complex>(N_cols));

    // --- ETAPA 2: Aplicar IDFT 1D (Complex -> Complex) em cada COLUNA da matriz temporária ---
    for (int j = 0; j < N_cols; ++j)
    {
        // Extrai a coluna da matriz temporária
        std::vector<Complex> col_in(N_rows);
        for (int i = 0; i < N_rows; ++i)
        {
            col_in[i] = temp_matrix[i][j];
        }

        // Calcula a IDFT da coluna
        std::vector<Complex> col_out;
        idft_1d_complex(col_in, col_out);

        // Armazena o resultado na coluna correspondente da matriz complexa final
        for (int i = 0; i < N_rows; ++i)
        {
            final_complex_matrix[i][j] = col_out[i];
        }
    }

    // --- ETAPA 3: Copiar a parte real para a saída e aplicar a NORMALIZAÇÃO 2D ---
    out.assign(N_rows, std::vector<double>(N_cols));
    double normalization_factor = 1.0 / (N_rows * N_cols);
    for (int i = 0; i < N_rows; ++i)
    {
        for (int j = 0; j < N_cols; ++j)
        {
            out[i][j] = final_complex_matrix[i][j].real() * normalization_factor;
        }
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
        sum += 0.5 * u_x[i] * u_x[i] + (MU / (P + 1.0)) * pow(u[i], P + 1.0);
    }
    return sum * DX;
}

double mass_motion_integral(const std::vector<double> &u, const std::vector<double> &u_x)
{
    double sum = 0.0;
    for (int i = 0; i < N; i++)
    {
        // Correção: Adicionar parênteses em (P + 1.0)
        sum += 3.0 * pow(u_x[i], 2.0) + (2.0 * MU * P / (P + 1.0)) * pow(u[i], P + 1.0);
    }
    return sum * DX * (-1.0) * (1.0 / mass_conservation(u));
}
double energy_motion_integral(const std::vector<double> &u, const std::vector<double> &u_x,
                              const std::vector<double> &u_xx)
{
    double sum = 0.0;
    for (int i = 0; i < N; i++)
    {
        // Correção: O termo deve ser 2*mu*u^(p-1)*ux^2
        sum += 1.5 * pow(u_xx[i], 2.0) + 2.0 * MU * pow(u[i], P - 1.0) * pow(u_x[i], 2.0) +
               (MU * MU * 0.5) * pow(u[i], 2.0 * P);
    }
    return sum * DX * (-1.0) * (1.0 / energy_conservation(u, u_x));
}

double mass_center(const std::vector<double> &u, const std::vector<double> &x_values)
{
    double sum = 0.0;
    for (int i = 0; i < N; i++)
    {
        sum += x_values[i] * (pow(u[i], 2.0));
    }
    return sum * DX * (1.0 / mass_conservation(u));
}
double energy_center(const std::vector<double> &u, const std::vector<double> &u_x, const std::vector<double> &x_values)
{
    double sum = 0.0;
    for (int i = 0; i < N; i++)
    {
        sum += x_values[i] * (0.5 * u_x[i] * u_x[i] + (MU / (P + 1.0)) * pow(u[i], P + 1.0));
    }
    return sum * DX * (1.0 / mass_conservation(u));
}
