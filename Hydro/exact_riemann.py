"""Exact solution of the 1-D Riemann problem for an ideal gas (Toro, chapter 4)."""
import numpy as np

NEWTON_TOLERANCE = 1e-10
NEWTON_MAX_ITERATIONS = 100


class RiemannState:
    """Uniform gas state on one side of the initial discontinuity."""

    def __init__(self, density, velocity, pressure, gamma):
        self.density = density
        self.velocity = velocity
        self.pressure = pressure
        self.sound_speed = np.sqrt(gamma * pressure / density)


def velocity_jump(p, state, gamma):
    """f_K(p) and df_K/dp: the velocity change across the wave between state K and pressure p.

    p > P_K gives a shock, p <= P_K a rarefaction.
    """
    if p > state.pressure:
        scale = 2 / ((gamma + 1) * state.density)
        pressure_shift = (gamma - 1) / (gamma + 1) * state.pressure
        root = np.sqrt(scale / (p + pressure_shift))
        jump = (p - state.pressure) * root
        derivative = root * (1 - (p - state.pressure) / (2 * (p + pressure_shift)))
    else:
        ratio = p / state.pressure
        jump = 2 * state.sound_speed / (gamma - 1) * (ratio ** ((gamma - 1) / (2 * gamma)) - 1)
        derivative = ratio ** (-(gamma + 1) / (2 * gamma)) / (state.density * state.sound_speed)
    return jump, derivative


def star_region(left, right, gamma):
    """P_star and u_star: f_L(P_star) + f_R(P_star) + u_R - u_L = 0, solved with Newton."""
    p = 0.5 * (left.pressure + right.pressure)
    for iteration in range(NEWTON_MAX_ITERATIONS):
        jump_left, derivative_left = velocity_jump(p, left, gamma)
        jump_right, derivative_right = velocity_jump(p, right, gamma)
        residual = jump_left + jump_right + right.velocity - left.velocity
        p_new = max(p - residual / (derivative_left + derivative_right), 1e-12)
        converged = abs(p_new - p) < NEWTON_TOLERANCE * p
        p = p_new
        if converged:
            break
    u_star = 0.5 * (left.velocity + right.velocity) + 0.5 * (jump_right - jump_left)
    return p, u_star


def inside_rarefaction_fan(state, speed, gamma, side):
    """Gas state inside a rarefaction fan at speed = x/t; side is -1 for the left wave, +1 for the right."""
    c = 2 / (gamma + 1) * (state.sound_speed - side * (gamma - 1) / 2 * (state.velocity - speed))
    u = 2 / (gamma + 1) * (-side * state.sound_speed + (gamma - 1) / 2 * state.velocity + speed)
    density = state.density * (c / state.sound_speed) ** (2 / (gamma - 1))
    pressure = state.pressure * (c / state.sound_speed) ** (2 * gamma / (gamma - 1))
    return density, u, pressure


def sample_side(state, p_star, u_star, speed, gamma, side):
    """Gas state at speed = x/t on one side of the contact; side is -1 (left) or +1 (right)."""
    ratio = p_star / state.pressure
    if p_star > state.pressure:
        density_star = state.density * (ratio + (gamma - 1) / (gamma + 1)) / ((gamma - 1) / (gamma + 1) * ratio + 1)
        shock_speed = state.velocity + side * state.sound_speed * np.sqrt(
            (gamma + 1) / (2 * gamma) * ratio + (gamma - 1) / (2 * gamma))
        if side * (speed - shock_speed) > 0:
            return state.density, state.velocity, state.pressure
        return density_star, u_star, p_star

    density_star = state.density * ratio ** (1 / gamma)
    sound_speed_star = state.sound_speed * ratio ** ((gamma - 1) / (2 * gamma))
    head_speed = state.velocity + side * state.sound_speed
    tail_speed = u_star + side * sound_speed_star
    if side * (speed - head_speed) > 0:
        return state.density, state.velocity, state.pressure
    if side * (speed - tail_speed) < 0:
        return density_star, u_star, p_star
    return inside_rarefaction_fan(state, speed, gamma, side)


def exact_solution(x, x_interface, t, left, right, gamma):
    """Density, velocity and pressure at the positions x, at time t > 0."""
    p_star, u_star = star_region(left, right, gamma)
    density = np.zeros(len(x))
    velocity = np.zeros(len(x))
    pressure = np.zeros(len(x))
    for i in range(len(x)):
        speed = (x[i] - x_interface) / t
        if speed < u_star:
            values = sample_side(left, p_star, u_star, speed, gamma, side=-1)
        else:
            values = sample_side(right, p_star, u_star, speed, gamma, side=+1)
        density[i], velocity[i], pressure[i] = values
    return density, velocity, pressure
