#include "naglowek.h"
#include <matplot/matplot.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

// Implementacja klasy Signal
Signal::Signal(double frequency, double faza, std::string signal_name) : f(frequency), name(signal_name) {
    samples = generate(f, faza, name);
}
Signal::Signal(std::string signal_name, const Fourier& transformata) : name(signal_name) {
    samples = idft(transformata);
}
Signal::~Signal() {
    std::cout << "Signal destructor called for " << name << std::endl;
}

// Implementacja klasy Fourier
Fourier::Fourier(const Signal& signal) {
    X = dft(signal);
}
Fourier::~Fourier() {}

// Funkcje przetwarzania
std::vector<double> generate(double f, double faza, std::string name) {
    int sf = N;
    int n = N * T;
    std::vector<double> signal(n);

    for (int i = 0; i < n; i++) {
        double phase = (2 * PI * f * i) / sf + faza * PI / 180;
        if (name == "sin") signal[i] = sin(phase);
        else if (name == "cos") signal[i] = cos(phase);
        else if (name == "kwadrat") signal[i] = sin(phase) >= 0 ? 1 : -1;
        else if (name == "pila") signal[i] = (1 / PI) * (phase - PI * floor(phase / PI));
        else return {};
    }
    return signal;
}

std::vector<std::complex<double>> dft(Signal signal) {
    std::vector<std::complex<double>> X(signal.samples.size());
    const double scale = 1.0 / std::sqrt(signal.samples.size());
    for (int k = 0; k < signal.samples.size(); k++) {
        X[k] = std::complex<double>(0.0, 0.0);
        for (int n = 0; n < signal.samples.size(); n++) {
            double angle = 2.0 * PI * k * n / signal.samples.size();
            X[k] += signal.samples[n] * std::polar(1.0, -angle);
        }
        X[k] *= scale;
    }
    return X;
}

std::vector<double> idft(Fourier fourier) {
    if (fourier.X.empty()) return {};
    const std::vector<std::complex<double>> X = fourier.X;
    std::vector<double> samples(X.size());
    for (int n = 0; n < X.size(); n++) {
        std::complex<double> sum(0.0, 0.0);
        for (int k = 0; k < X.size(); k++) {
            sum += X[k] * std::polar(1.0, 2.0 * PI * k * n / X.size());
        }
        samples[n] = std::real(sum);
    }
    double max_val = *std::max_element(samples.begin(), samples.end());
    for (double& val : samples) val /= max_val;
    return samples;
}

Signal filter(Signal signal, double cutoff_hz) {
    Fourier transformata(signal);
    int M = signal.samples.size();
    double fs = N;
    double df = fs / M;
    for (int k = 0; k < M; ++k) {
        double freq = (k >= M / 2) ? (k - M) * df : k * df;
        if (std::abs(freq) > cutoff_hz) transformata.X[k] = { 0.0, 0.0 };
    }
    return Signal("Filtered " + signal.name, transformata);
}
// Rysowanie wykresów
void plot_signal(Signal signal) {
    using namespace matplot;
    std::filesystem::create_directory("../raport");
    std::vector<double> t(signal.samples.size());
    double dt = T / signal.samples.size();
    for (size_t i = 0; i < t.size(); i++) t[i] = i * dt;
    xlim({ 0, T - 1 });
    ylim({ -1.2, 1.2 });
    plot(t, signal.samples);
    xlabel("Time [s]");
    ylabel("Amplitude");
    title("Signal: " + signal.name);
    save("../raport/" + signal.name + "_signal.png");
    show();
}

void plot_fourier(Fourier fourier) {
    using namespace matplot;
    std::filesystem::create_directory("../raport");
    size_t M = fourier.X.size();
    double fs = N, df = fs / M;
    size_t H = M / 2 + 1;
    std::vector<double> f(H), mag(H);
    for (size_t k = 0; k < H; ++k) {
        f[k] = k * df;
        mag[k] = 2.0 * std::abs(fourier.X[k]) / M;
    }
    xlim({ 0, *std::max_element(f.begin(), f.end()) * 1.2 });
    ylim({ 0, *std::max_element(mag.begin(), mag.end()) * 1.2 });
    plot(f, mag);
    xlabel("Frequency [Hz]");
    ylabel("Magnitude");
    title("Fourier Transform");
    save("../raport/Fourier_transform.png");
    show();
}

// PYBIND11
namespace py = pybind11;
PYBIND11_MODULE(example, m) {
    m.doc() = "Signal processing module";

    py::class_<Signal>(m, "Signal")
        .def(py::init<double, double, std::string>())
        .def(py::init<std::string, Fourier>())
        .def(py::init<const std::vector<double>&>())
        .def_readwrite("frequency", &Signal::f)
        .def_readwrite("name", &Signal::name)
        .def_readwrite("samples", &Signal::samples);

    py::class_<Fourier>(m, "Fourier")
        .def(py::init<Signal>())
        .def_readwrite("X", &Fourier::X);

    m.attr("N") = N;
    m.attr("PI") = PI;

    m.def("generate", &generate);
    m.def("plot_signal", &plot_signal);
    m.def("plot_fourier", &plot_fourier);
    m.def("dft", &dft);
    m.def("idft", &idft);
    m.def("filter", &filter);
}
