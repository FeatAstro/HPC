"""Figure for Step 3c: whistler frequency against the wavenumber and resistive decay of a wave,
computed by the code (coloured) against the theory (black, dashed).

Run after building:  python3 plots/plot_step3c_faraday.py
"""
import math

import numpy as np

from style import BLUE, INK, ORANGE, plt, save
from common import follow_amplitude, pic, plasma_parameters, transverse_wave

GRID = pic.Grid(128, 2.0 * math.pi / 128)
WAVE_AMPLITUDE = 1e-3
EXPECTED = dict(color=INK, lw=1, ls="--")


def whistler_frequency(mode, Bx=1.0, dt=5e-4):
    fields, k = transverse_wave(GRID, mode, WAVE_AMPLITUDE, Bx)
    steps = round(2 * 2.0 * math.pi / (k**2 * Bx) / dt)
    times, amplitudes = follow_amplitude(fields, GRID, plasma_parameters(), mode, dt, max(steps // 40, 1), 40)
    return k, abs(np.polyfit(times, np.unwrap(np.angle(amplitudes)), 1)[0])


def plot_whistler(ax):
    modes = np.arange(1, 9)
    wavenumbers, measured = np.array([whistler_frequency(mode) for mode in modes]).T
    expected = wavenumbers**2
    fine_k = np.linspace(0.0, wavenumbers[-1], 200)
    ax.plot(wavenumbers, measured, "o", color=BLUE, ms=9, label="measured: rotation frequency of the wave pattern")
    ax.plot(fine_k, fine_k**2, label="expected: ω = k² Bx / n", **EXPECTED)
    ax.set_xlabel("wavenumber k")
    ax.set_ylabel("frequency ω")
    ax.set_title("Whistler waves (Hall term only, ions at rest, Bx = n = 1)\n"
                 f"largest relative difference = {np.max(np.abs(measured - expected) / expected):.1%} "
                 "(at the shortest wavelength)")
    ax.legend(loc="upper left", fontsize=8)


def plot_decay(ax):
    mode, dt = 2, 2e-4
    parameters = plasma_parameters(resistivity=0.05, hyper_resistivity=0.001)
    fields, k = transverse_wave(GRID, mode, WAVE_AMPLITUDE, Bx=0.0)
    rate = parameters.resistivity * k**2 + parameters.hyper_resistivity * k**4
    times, amplitudes = follow_amplitude(fields, GRID, parameters, mode, dt, 700, 50)
    measured = np.abs(amplitudes) / WAVE_AMPLITUDE
    expected = np.exp(-rate * times)
    ax.semilogy(times, measured, color=ORANGE, lw=3, label="measured: amplitude of the wave")
    ax.semilogy(times, expected, label="expected: exp(−(η k² + ν k⁴) t)", **EXPECTED)
    ax.set_xlabel("time t")
    ax.set_ylabel("amplitude / initial amplitude")
    ax.set_title(f"Resistive decay (no Bx, η = {parameters.resistivity}, ν = {parameters.hyper_resistivity}, k = {k:g})\n"
                 f"largest relative difference = {np.max(np.abs(measured - expected) / expected):.1e}")
    ax.legend(loc="upper right", fontsize=8)


def plot_faraday():
    fig, (ax_whistler, ax_decay) = plt.subplots(1, 2, figsize=(13, 5))
    plot_whistler(ax_whistler)
    plot_decay(ax_decay)
    save(fig, "step3c_faraday.png")


if __name__ == "__main__":
    plot_faraday()
