#include <cstdlib>
#include <fftw3.h>
#include <iostream>

int main(void)
{

    double *u_real = (double *)fftw_malloc(sizeof(double) * 100);
    fftw_complex *u_hat_antigo = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * 100);
    fftw_complex *u_hat_atual = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * 100);
    fftw_complex *u_hat_novo = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * 100);
    fftw_plan plan_fwd = fftw_plan_dft_r2c_1d(100, u_real, u_hat_atual, FFTW_MEASURE | FFTW_PRESERVE_INPUT);
    fftw_plan plan_bwd = fftw_plan_dft_c2r_1d(100, u_hat_atual, u_real, FFTW_MEASURE | FFTW_PRESERVE_INPUT);

    for (int i = 0; i < 100; i++)
    {
        u_real[i] = (double)i;
        std::cout << u_real[i] << " ";
    }
    std::cout << "\n";

    fftw_execute_dft_r2c(plan_fwd, u_real, u_hat_novo);
    fftw_execute_dft_c2r(plan_bwd, u_hat_novo, u_real);
    for (int i = 0; i < 100; i++)
    {
        // u_real[i] = (double)i;
        std::cout << u_real[i] << " ";
    }
    std::cout << "\n";
    free(u_real);
    free(u_hat_novo);
    free(u_hat_atual);
    free(u_hat_antigo);
    fftw_destroy_plan(plan_bwd);
    fftw_destroy_plan(plan_fwd);

    return 0;
}
