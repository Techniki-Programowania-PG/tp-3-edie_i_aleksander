#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <complex>
#include <filesystem>

constexpr double PI = 3.14159265358;
constexpr int N = 100;
constexpr double T = 5.0;

class Fourier;

class Signal {
public:
    double f;
    std::string name;
    std::vector<double> samples;

    Signal(double frequency, double faza, std::string signal_name);
    Signal(std::string signal_name, const Fourier& transformata);
    Signal(const std::vector<double>& samples);
    ~Signal();
};

class Fourier {
public:
    std::vector<std::complex<double>> X;

    Fourier(const Signal& signal);
    ~Fourier();
};

std::vector<double> generate(double f, double faza, std::string name);
std::vector<std::complex<double>> dft(Signal signal);
std::vector<double> idft(Fourier fourier);
Signal filter(Signal signal, double cutoff_hz);
std::vector<double> correlation(const Signal& signal1, const Signal& signal2);
void plot_signal(Signal signal);
void plot_fourier(Fourier fourier);
