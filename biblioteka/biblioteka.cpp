#include "biblioteka.hpp"
#include <matplot/matplot.h>
#include <ctime>
#include <cstdlib>
#include <cstdio>

std::vector<point> sinus(int sampleNumber, double timeStart, double timeEnd, double frequency)
{
    std::vector<point> p(sampleNumber);
    double dt = (timeEnd - timeStart) / (sampleNumber - 1);
    for (int i = 0; i < sampleNumber; ++i)
    {
        p[i].x = timeStart + i * dt;
        p[i].y = sin(2 * matplot::pi * frequency * p[i].x);
    }
    return p;
}

std::vector<point> cosinus(int sampleNumber, double timeStart, double timeEnd, double frequency)
{
    std::vector<point> p(sampleNumber);
    double dt = (timeEnd - timeStart) / (sampleNumber - 1);
    for (int i = 0; i < sampleNumber; ++i)
    {
        p[i].x = timeStart + i * dt;
        p[i].y = cos(2 * matplot::pi * frequency * p[i].x);
    }
    return p;
}

std::vector<point> rectangle(int sampleNumber, double timeStart, double timeEnd, double frequency)
{
    std::vector<point> p(sampleNumber);
    double dt = (timeEnd - timeStart) / (sampleNumber - 1);
    for (int i = 0; i < sampleNumber; ++i)
    {
        p[i].x = timeStart + i * dt;
        double s = sin(2 * matplot::pi * frequency * p[i].x);
        p[i].y = s >= 0 ? 1.0 : -1.0;
    }
    return p;
}

std::vector<point> saw(int sampleNumber, double timeStart, double timeEnd, double frequency)
{
    std::vector<point> p(sampleNumber);
    double dt = (timeEnd - timeStart) / (sampleNumber - 1);
    for (int i = 0; i < sampleNumber; ++i)
    {
        p[i].x = timeStart + i * dt;
        double t = p[i].x;
        p[i].y = 2.0 * (frequency * t - floor(frequency * t + 0.5));
    }
    return p;
}

void draw(std::vector<point> data, std::string descritpion)
{
    std::vector<double> x(data.size()), y(data.size());
    for (int i = 0; i < data.size(); i++)
    {
        x[i] = data[i].x;
        y[i] = data[i].y;
    }

    double min = 0;
    double max = 0;
    for (int i = 0; i < y.size(); i++)
    {
        if (y[i] > max)
            max = y[i];
        if (y[i] < min)
            min = y[i];
    }
    matplot::ylim({min - 0.5, max + 0.5});
    min = 0;
    max = 0;
    for (int i = 0; i < x.size(); i++)
    {
        if (x[i] > max)
            max = x[i];
        if (x[i] < min)
            min = x[i];
    }
    matplot::plot(x, y);
    matplot::xlim({min - (0.05 * max), max + (0.05 * max)});
    matplot::title(descritpion);
    matplot::xlabel("Czas [s]");
    matplot::ylabel("Amplituda");
    matplot::grid(true);
    matplot::show();
}

std::vector<std::complex<double>> DFT(std::vector<point> &inputSignal)
{
    int N = inputSignal.size();
    std::vector<std::complex<double>> dft(N);

    for (int k = 0; k < N; ++k)
    {
        std::complex<double> sum(0.0, 0.0);
        for (int n = 0; n < N; ++n)
        {
            double angle = -2.0 * matplot::pi * k * n / N;
            std::complex<double> w(cos(angle), sin(angle));
            sum += inputSignal[n].y * w;
        }
        dft[k] = sum;
    }

    return dft;
}

std::vector<point> spectrumDFT(std::vector<std::complex<double>> &dft)
{
    int N = dft.size();
    std::vector<point> spectrum(N);

    for (int k = 0; k < N; ++k)
    {
        spectrum[k].x = k;
        spectrum[k].y = std::abs(dft[k]) / N;
    }

    return spectrum;
}

std::vector<point> IDFT(std::vector<std::complex<double>> &dft, double timeStart, double timeEnd)
{
    int N = dft.size();

    std::vector<point> signal(N);
    double dt = (timeEnd - timeStart) / (N - 1);

    for (int n = 0; n < N; ++n)
    {
        std::complex<double> sum(0.0, 0.0);
        for (int k = 0; k < N; ++k)
        {
            double angle = 2.0 * matplot::pi * k * n / N;
            std::complex<double> w(cos(angle), sin(angle));
            sum += dft[k] * w;
        }

        signal[n].x = timeStart + n * dt;
        signal[n].y = sum.real() / N;
    }

    return signal;
}

std::vector<std::complex<double>> filter1D(std::vector<std::complex<double>> &dft, double cutoffRatio)
{
    int N = dft.size();
    int cutoffIndex = static_cast<int>(cutoffRatio * N / 2.0);

    std::vector<std::complex<double>> filteredDFT(N, {0.0, 0.0});

    for (int k = 0; k <= cutoffIndex; ++k)
        filteredDFT[k] = dft[k];

    for (int k = N - cutoffIndex; k < N; ++k)
        filteredDFT[k] = dft[k];

    return filteredDFT;
}

std::vector<point> generateMatrix(int rows, int cols)
{
    std::vector<point> matrix(cols);
    std::srand(std::time(0));
    for (int x = 0; x < cols; x++)
    {

        matrix[x].y = std::rand() % rows;

        matrix[x].x = x;
    }
    return matrix;
}

std::vector<std::vector<double>> convert1Dto2D(std::vector<point> &inputSignal, int scaleToRows)
{
    int cols = inputSignal.size();
    int rows = cols;
    double sr = (double)scaleToRows * 1.0;

    std::vector<std::vector<double>> result(rows, std::vector<double>(cols));
    for (int x = 0; x < cols; x++)
    {
        double oy = inputSignal[x].y;
        if (scaleToRows > 0)
            oy = oy * (sr / 2.0) + (sr / 2.0);
        for (int y = 0; y < rows; y++)
        {
            if (y == static_cast<int>(oy))
                result[y][x] = 1;
            else
                result[y][x] = 0;
        }
    }
    return result;
}

std::vector<std::vector<double>> gaussianKernel2D(int size, double sigma)
{
    std::vector<std::vector<double>> kernel(size, std::vector<double>(size));
    double sum = 0.0;
    int center = size / 2;

    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            int x = i - center;
            int y = j - center;
            kernel[i][j] = std::exp(-(x * x + y * y) / (2 * sigma * sigma));
            sum += kernel[i][j];
        }
    }
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            kernel[i][j] /= sum;

    return kernel;
}

std::vector<std::vector<double>> applyGaussianBlur2D(std::vector<point> &inputSignal, int size, double sigma, int scaleToRows)
{

    std::vector<std::vector<double>> &kernel = gaussianKernel2D(size, sigma);
    std::vector<std::vector<double>> &input = convert1Dto2D(inputSignal, scaleToRows);
    int rows = input.size();
    int cols = input[0].size();
    int kSize = kernel.size();
    int kCenter = kSize / 2;

    std::vector<std::vector<double>> output(rows, std::vector<double>(cols, 0.0));

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            double sum = 0.0;
            for (int ki = 0; ki < kSize; ++ki)
            {
                for (int kj = 0; kj < kSize; ++kj)
                {
                    int ii = i + ki - kCenter;
                    int jj = j + kj - kCenter;
                    if (ii >= 0 && ii < rows && jj >= 0 && jj < cols)
                    {
                        sum += input[ii][jj] * kernel[ki][kj];
                    }
                }
            }
            output[i][j] = sum;
        }
    }

    return output;
}

void plotMatrix(std::vector<std::vector<double>> &mat, std::string &titleStr)
{
    int size = mat.size();
    std::vector<double> ticks(11);

    for (int i = 0; i < 11; i++)
    {
        ticks[i] = 0.1 * i * size;
    }

    matplot::heatmap(mat);
    matplot::xlim({0.0, size * 1.0});
    matplot::ylim({0.0, size * 1.0});
    matplot::title(titleStr);
    matplot::xlabel("Czas [s]");
    matplot::ylabel("Amplituda");
    matplot::xticks(ticks);
    matplot::yticks(ticks);
    matplot::grid(true);

    matplot::show();
}
