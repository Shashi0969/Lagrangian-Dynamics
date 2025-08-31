# Apply principle of Least Action and obtain the path that minimizes the action relative to all other possible paths between the same points for a Harmonic Oscillator System with potential defined by V(x)= 1/2 kx^2

import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import simpson
# Constants
m = 1.0      # mass
L = 1.0      # final position
T = 1.0      # total time
hbar = 1.0      # reduced Planck's constant
k = 5.0      # spring constant (non-zero potential)
num_paths = 100  # number of paths to simulate

# Time array
t = np.linspace(0, T, 1000)
dt = t[1] - t[0]

# Initialize lists
amplitudes = []
paths = []

def Potenial_Energy(x):
    return 0.5 * k * x**2
    
# Generate multiple paths
for i in range(num_paths):
    epsilon = np.random.uniform(-0.2, 0.2)  # small wiggle amplitude
    wiggle = epsilon * np.sin(2 * np.pi * t / T)  # sinusoidal perturbation
    x_path = (L / T) * t + wiggle  # perturbed path
    v_path = np.gradient(x_path, dt) # velocity
    T_path = 0.5 * m * v_path**2    # Kinetic Energy
    V_path = Potenial_Energy(x_path) # Potential Energy
    L_path = T_path - V_path    # Lagrangian (kinetic energy only)
    S_path = simpson(L_path, t)  # action
    phase = np.exp(1j * S_path / hbar)  # complex amplitude
    amplitudes.append(phase)
    paths.append(x_path)

# Sum of amplitudes
total_amplitude = sum(amplitudes)
probability = abs(total_amplitude)**2

# Plotting a subset of paths
plt.figure(figsize=(12, 6))
for i in range(0, num_paths, 10):  # Plot every 10th path
    plt.plot(t, paths[i], alpha=0.6, linewidth=1.2)

# Classical path (straight line)
x_classical = (L / T) * t
plt.plot(t, x_classical, 'k--', linewidth=2, label="Classical Path (Least Action)")

# Labels and display
plt.title("Many Quantum Paths from x=0 to x=L (Feynman Path Integral Picture)")
plt.xlabel("Time")
plt.ylabel("Position")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# Output the total quantum probability
print(f"Total probability amplitude (|Σ $e^(iS/ħ)$|²): {probability:.4f}")
