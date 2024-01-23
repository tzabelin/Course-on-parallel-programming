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

void norm_rows(const int& ny, const int& nx, const float* data, std::vector<double>& norm_data)
{
    for (int y = 0; y < ny; y++)
    {
        double sum = 0.0;
        for (int x = 0; x < nx; x++)
        {
            sum += data[x + y * nx];
        }
        double mean = sum / nx;

        double sqr_sum = 0.0;
        for (int x = 0; x < nx; x++)
        {
            norm_data[x + y * nx] = data[x + y * nx] - mean;
            sqr_sum += norm_data[x + y * nx] * norm_data[x + y * nx];
        }

        double norm = std::sqrt(sqr_sum);
        for (int x = 0; x < nx; x++)
        {
            norm_data[x + y * nx] /= norm;
        }
    }
}

void correlate(int ny, int nx, const float* data, float* result)
{
    std::vector<double> norm_data(ny * nx);
    norm_rows(ny, nx, data, norm_data);
    double corr = 0.0;

    for (int i = 0; i < ny; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            corr = 0.0;
            for (int x = 0; x < nx; x++)
            {
                corr += norm_data[x + i * nx] * norm_data[x + j * nx];
            }
            result[i + j * ny] = corr;
        }
    }
}

