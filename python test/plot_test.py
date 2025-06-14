import example
import sys

# dodanie katalogu z .pyd do ścieżki
sys.path.append(".")

def main():
    freq = float(input("Podaj częstotliwość sygnału: "))
    phase = float(input("Podaj fazę sygnału (w stopniach): "))
    name = input("Podaj typ sygnału (sin, cos, kwadrat, pila): ")

    sig = example.Signal(freq, phase, name)
    example.plot_signal(sig)

    filtruj = input("Czy chcesz wykonać filtrację sygnału? (t/n): ").strip().lower()
    if filtruj == 't':
        cutoff = float(input("Podaj częstotliwość odcięcia [Hz]: "))
        filtered = example.filter(sig, cutoff)
        example.plot_signal(filtered)

    transformuj = input("Czy chcesz wykonać transformatę Fouriera? (t/n): ").strip().lower()
    if transformuj == 't':
        trans = example.Fourier(sig)
        example.plot_fourier(trans)

if __name__ == '__main__':
    main()
