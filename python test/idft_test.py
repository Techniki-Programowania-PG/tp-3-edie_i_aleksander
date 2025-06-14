import example
import matplotlib.pyplot as plt

# Parametry
freq = float(input("Podaj częstotliwość (Hz): "))
signal_type = input("Podaj typ sygnału (sin/cos/kwadrat/pila): ")

# Generuj sygnał
sig = example.Signal(freq, 0.0, signal_type)

# Transformata Fouriera
fourier = example.Fourier(sig)

# Odwrotna transformata
samples = example.idft(fourier)
reconstructed = example.Signal(samples)

# Oś czasu
T = 5.0  # Zdefiniowane lokalnie, zamiast korzystać z example.T
t = [i * (T / len(sig.samples)) for i in range(len(sig.samples))]

# Wykresy
plt.figure(figsize=(10, 5))
plt.plot(t, sig.samples, label="Oryginalny")
plt.plot(t, reconstructed.samples, label="Po IDFT", linestyle="--")
plt.xlabel("Czas [s]")
plt.ylabel("Amplituda")
plt.title("Porównanie sygnału oryginalnego i po IDFT")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
