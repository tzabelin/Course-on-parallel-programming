/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns: 0 <= x < nx
- element at row y and column x is stored in data[x + y*nx]
- correlation between rows i and row j has to be stored in result[i + j*ny]
- only parts with 0 <= j <= i < ny need to be filled
*/

#include <cmath>
#include <vector>
#include <x86intrin.h>

typedef float double4_t __attribute__((vector_size(8 * sizeof(float))));

static inline double4_t swap4(double4_t x) { return _mm256_permute2f128_ps(x, x, 0b00000001); }
static inline double4_t swap2(double4_t x) { return _mm256_permute_ps(x, 0b01001110); }
static inline double4_t swap1(double4_t x) { return _mm256_permute_ps(x, 0b10110001); }


constexpr double4_t d4zeros
{
    0.0, 0.0, 0.0, 0.0
};


void norm_rows(const int& ny, const int& nx, const float* data, std::vector<double4_t>& vd, const int& na, const int& nb)
{
#pragma omp parallel for
    for (int block = 0; block < na; block++)
    {
        for (int i = 0; i < nb; i++)
        {
            if (block * nb + i < ny)
            {
                std::vector<double> norm_data(nx, 0);
                double sum = 0.0;
                for (int x = 0; x < nx; x++)
                {
                    sum += data[x + (nb * block + i) * nx];
                }
                double mean = sum / nx;

                double sqr_sum = 0.0;
                for (int x = 0; x < nx; x++)
                {
                    norm_data[x] = data[x + (nb * block + i) * nx] - mean;
                    sqr_sum += norm_data[x] * norm_data[x];
                }

                double norm = std::sqrt(sqr_sum);

                for (int x = 0; x < nx; x++)
                {
                    vd[nx * block + x][i] = norm_data[x] / norm;
                }
            }
        }
    }
}

void correlate(int ny, int nx, const float* data, float* result) {

    int na = (ny + 4 - 1) / 4;

    std::vector<double4_t> vd(na * nx, d4zeros);


    norm_rows(ny, nx, data, vd, na, 4);


}