"""Figures for Step 1 (Boris pusher), one per test: what the code gives (thick, coloured)
against what is expected (thin, black, dashed).

Run after building:  python3 plots/plot_step1_boris.py
"""
import numpy as np

from style import AQUA, BLUE, INK, INK_SECONDARY, ORANGE, plt, save
from common import NO_FIELD, cyclotron_period, field_along_z, larmor_radius, make_particle, pic

B = 1.0
STEPS_PER_PERIOD = 200

# The electron gets m = 1 (with its real mass its orbit would be invisible next to the ions),
# and every particle has its own speed so that no curve hides another.
PARTICLES = [
    # label,                    q,    m,   v,   colour
    ("ion, m = 1, v = 1", 1.0, 1.0, 1.0, BLUE),
    ("ion, m = 2, v = 0.8", 1.0, 2.0, 0.8, AQUA),
    ("electron (m = 1), v = 0.6", -1.0, 1.0, 0.6, ORANGE),
]

MEASURED = dict(lw=3)
EXPECTED = dict(color=INK, lw=1, ls="--")


def arrow_between(ax, xs, ys, index, color):
    ax.annotate("", xy=(xs[index + 1], ys[index + 1]), xytext=(xs[index], ys[index]),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=2, mutation_scale=18))


def plot_larmor_orbits():
    fig, (ax_x, ax_v, ax_radius) = plt.subplots(1, 3, figsize=(15, 5.2), gridspec_kw={"width_ratios": [1.5, 1, 1]})

    dt = cyclotron_period(1.0, 1.0, B) / STEPS_PER_PERIOD
    largest_difference = 0.0
    for index, (label, q, m, v, color) in enumerate(PARTICLES):
        particle = make_particle(v=(v, 0.0, 0.0), q=q, m=m)
        trajectory = pic.run_particle(particle, NO_FIELD, field_along_z(B), dt, 2 * STEPS_PER_PERIOD)
        t = trajectory.t
        cyclotron_frequency = abs(q) * B / m
        expected_x = larmor_radius(v, q, m, B) * np.sin(cyclotron_frequency * t)
        largest_difference = max(largest_difference, np.max(np.abs(trajectory.x - expected_x)))

        ax_x.plot(t, trajectory.x, color=color, label=f"measured: {label}", **MEASURED)
        ax_x.plot(t, expected_x, label="expected: r_L sin(ω_c t),  r_L = m v / |q| B" if index == len(PARTICLES) - 1 else None,
                  **EXPECTED)

        one_period = round(cyclotron_period(q, m, B) / dt)
        angle = np.linspace(0.0, 2.0 * np.pi, 200)
        ax_v.plot(trajectory.vx[:one_period + 1], trajectory.vy[:one_period + 1], color=color, **MEASURED)
        ax_v.plot(v * np.cos(angle), v * np.sin(angle),
                  label="expected: circle of radius v" if index == 0 else None, **EXPECTED)
        ax_v.plot(v, 0, "o", color=INK, ms=4)
        arrow_between(ax_v, trajectory.vx, trajectory.vy, one_period // 8, color)
        arrow_between(ax_v, trajectory.vx, trajectory.vy, 5 * one_period // 8, color)

    ax_x.set_xlabel("time t")
    ax_x.set_ylabel("position x")
    ax_x.set_title("Position: x oscillates with amplitude r_L\n"
                   f"largest |measured − expected| = {largest_difference:.0e} (Boris phase lag)")
    ax_x.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), fontsize=8, ncol=2)

    ax_v.set_aspect("equal")
    ax_v.set_xlim(-1.3, 1.3)
    ax_v.set_ylim(-1.3, 1.6)
    ax_v.set_xlabel("v_x")
    ax_v.set_ylabel("v_y")
    ax_v.set_title("Velocity (B out of the page), dots = start\n"
                   "expected: ions clockwise, electron anticlockwise")
    ax_v.legend(loc="upper center", fontsize=8)

    masses = np.array([1, 2, 4, 8, 16, 32, 64], dtype=float)
    measured_radius = []
    for m in masses:
        particle = make_particle(v=(1.0, 0.0, 0.0), m=m)
        trajectory = pic.run_particle(particle, NO_FIELD, field_along_z(B),
                                      cyclotron_period(1.0, m, B) / STEPS_PER_PERIOD, STEPS_PER_PERIOD)
        measured_radius.append((trajectory.x.max() - trajectory.x.min()) / 2)
    expected_radius = larmor_radius(1.0, 1.0, masses, B)
    relative_difference = np.max(np.abs(measured_radius - expected_radius) / expected_radius)
    ax_radius.plot(masses, measured_radius, "o", color=BLUE, ms=9, label="measured: (x_max − x_min) / 2")
    ax_radius.plot(masses, expected_radius, label="expected: r_L = m v / |q| B", **EXPECTED)
    ax_radius.set_xlabel("mass m  (q = 1, v = 1, B = 1)")
    ax_radius.set_ylabel("Larmor radius r_L")
    ax_radius.set_title("Larmor radius against mass\n"
                        f"largest relative difference = {relative_difference:.0e}")
    ax_radius.legend(loc="upper left", fontsize=8)

    save(fig, "step1_larmor_orbits.png")


def plot_error_vs_dt():
    period = cyclotron_period(1.0, 1.0, B)
    steps_per_period = np.array([16, 32, 64, 128, 256, 512, 1024])
    dts = period / steps_per_period
    measured_error = []
    for steps in steps_per_period:
        particle = make_particle(v=(1.0, 0.0, 0.0))
        pic.run_particle(particle, NO_FIELD, field_along_z(B), period / steps, int(steps))
        measured_error.append(abs(particle.x))
    measured_slope = np.polyfit(np.log(dts), np.log(measured_error), 1)[0]
    # Boris turns v by 2 atan(dt/2) per step instead of dt: after one period it is late by this angle
    phase_lag = 2.0 * np.pi - steps_per_period * 2.0 * np.arctan(np.pi / steps_per_period)
    expected_error = np.sin(phase_lag)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(dts, measured_error, "o-", color=BLUE, ms=9, label="measured: |x| after one period", **MEASURED)
    ax.loglog(dts, expected_error, label="expected: r_L sin(phase lag of Boris) ≈ π dt² / 6", **EXPECTED)
    ax.set_xlabel("time step dt  (one period = 2π)")
    ax.set_ylabel("|x| after one period  (exact motion: 0)")
    ax.set_title("Error of the Boris pusher against the time step\n"
                 f"measured slope = {measured_slope:.3f}  (expected: 2)")
    ax.legend(loc="upper left", fontsize=8)
    save(fig, "step1_error_vs_dt.png")


def plot_ExB_drift():
    Ey = 0.2
    drift_speed = Ey / B  # v_E = E x B / B^2 = (Ey / B, 0, 0)
    periods = 8
    dt = cyclotron_period(1.0, 1.0, B) / STEPS_PER_PERIOD

    started_at_drift = ("particle started at v_E", 1.0, 1.0, drift_speed, AQUA)
    gyrating = [particle for particle in PARTICLES if particle[2] == 1.0]

    fig, ax = plt.subplots(figsize=(11, 5.2))
    largest_difference = 0.0
    for index, (label, q, m, v, color) in enumerate(gyrating + [started_at_drift]):
        particle = make_particle(v=(v, 0.0, 0.0), q=q, m=m)
        trajectory = pic.run_particle(particle, pic.Vec3(0.0, Ey, 0.0), field_along_z(B), dt,
                                      periods * STEPS_PER_PERIOD)
        t = trajectory.t
        expected_x = drift_speed * t + (v - drift_speed) * np.sin(t)
        largest_difference = max(largest_difference, np.max(np.abs(trajectory.x - expected_x)))
        ax.plot(t, trajectory.x, color=color, label=f"measured: {label}", **MEASURED)
        ax.plot(t, expected_x, label="expected: v_E t + (v − v_E) sin(ω_c t)" if index == len(gyrating) else None, **EXPECTED)

    ax.set_xlabel("time t  (one period = 2π)")
    ax.set_ylabel("position x")
    ax.set_title(f"E×B drift, E = (0, {Ey}, 0), B = (0, 0, {B:g}): every particle advances at v_E = {drift_speed}, "
                 "whatever its charge\n"
                 f"largest |measured − expected| = {largest_difference:.0e} (Boris phase lag)")
    ax.legend(loc="upper left", fontsize=8)
    save(fig, "step1_ExB_drift.png")


if __name__ == "__main__":
    plot_larmor_orbits()
    plot_error_vs_dt()
    plot_ExB_drift()
