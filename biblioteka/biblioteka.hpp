#pragma once
#include "point.hpp"
#include <vector>
#include <complex>
#include <string>

std::vector<point> sinus(int sampleNumber, double timeStart, double timeEnd, double frequency);
std::vector<point> cosinus(int sampleNumber, double timeStart, double timeEnd, double frequency);
std::vector<point> rectangle(int sampleNumber, double timeStart, double timeEnd, double frequency);
std::vector<point> saw(int sampleNumber, double timeStart, double timeEnd, double frequency);
void draw(std::vector<point> data, std::string descritpion);
std::vector<std::complex<double>> DFT(std::vector<point> &inputSignal);
std::vector<point> spectrumDFT(std::vector<std::complex<double>> &dft);
std::vector<point> IDFT(std::vector<std::complex<double>> &dft, double timeStart, double timeEnd);
std::vector<std::complex<double>> filter1D(std::vector<std::complex<double>> &dft, double cutoffRatio);
std::vector<point> generateMatrix(int rows, int cols);
std::vector<std::vector<double>> applyGaussianBlur2D(std::vector<point> &inputSignal, int size, double sigma, int scaleToRows);
void plotMatrix(std::vector<std::vector<double>> &mat, std::string &titleStr);
std::vector<std::vector<double>> convert1Dto2D(std::vector<point> &inputSignal, int scaleToRows);