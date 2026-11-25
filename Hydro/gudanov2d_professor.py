"""A documented 2-D Euler finite-volume solver (Godunov / MUSCL-Hancock style).

Variables
---------
Primitive : W = (rho, u, v, p)
Conserved : U = (rho, rho*u, rho*v, E)
E = p/(gamma-1) + 0.5*rho*(u**2+v**2)

Array convention
----------------
All state arrays have shape (4, nx+2, ny+2).  The first and last
indices in each spatial direction are one layer of ghost cells.  Thus
(nx, ny) passed to Godunov2D are the numbers of PHYSICAL cells.
The physical domain is represented by indices 1:-1 in both directions.

Numerical sequence for one full step
------------------------------------
1. compute_dt: determine a CFL-limited dt from the current primitive state.
2. advance: evaluate the conservative spatial operator at U^n, build U^(1/2),
   evaluate the operator again, and form U^(n+1) (SSP-RK2 / midpoint form).
3. step: impose boundary conditions, recover primitive variables, check floors.

The implementation is deliberately explicit and educational.  It uses a
single ghost-cell layer and first-order boundary-face reconstruction; interior
slopes are limited with minmod.  For demanding production calculations, use
at least two ghost layers, a more systematic positivity-preserving method,
and convergence tests.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

import numpy as np

try:
    import matplotlib.pyplot as plt
except ImportError:  # plotting is optional for running the solver core
    plt = None

Array = np.ndarray
TINY = 1.0e-12


def prim_to_cons(W: Array, gamma: float) -> Array:
    """Convert primitive variables ``(rho,u,v,p)`` to conserved variables.

    Parameters
    ----------
    W : ndarray, shape (4, ...)
        Primitive state.
    gamma : float
        Ideal-gas adiabatic index.

    Returns
    -------
    U : ndarray, shape (4, ...)
        Conserved state ``(rho,rho*u,rho*v,E)``.
    """
    rho, u, v, p = W
    energy = p / (gamma - 1.0) + 0.5 * rho * (u*u + v*v)
    return np.stack((rho, rho*u, rho*v, energy), axis=0)


def cons_to_prim(U: Array, gamma: float, floor: float = TINY) -> Array:
    """Convert conserved to primitive variables, applying density/pressure floors.

    The pressure floor is a numerical safeguard, not a physical heating term.
    A production code should also report where and how often it is activated.
    """
    rho = np.maximum(U[0], floor)
    u = U[1] / rho
    v = U[2] / rho
    kinetic = 0.5 * rho * (u*u + v*v)
    p = (gamma - 1.0) * (U[3] - kinetic)
    p = np.maximum(p, floor)
    return np.stack((rho, u, v, p), axis=0)


def sound_speed(p: Array, rho: Array, gamma: float) -> Array:
    """Return the ideal-gas adiabatic sound speed ``sqrt(gamma*p/rho)``."""
    return np.sqrt(gamma * np.maximum(p, TINY) / np.maximum(rho, TINY))


def minmod(a: Array, b: Array) -> Array:
    """Minmod limiter: retain the smaller same-sign one-sided gradient."""
    same_sign = (a * b) > 0.0
    return np.where(same_sign, np.sign(a) * np.minimum(np.abs(a), np.abs(b)), 0.0)


def reconstruct(W: Array, idir: int) -> Tuple[Array, Array]:
    """Reconstruct left/right primitive states at faces in one direction.

    ``idir=0`` reconstructs x-faces; ``idir=1`` reconstructs y-faces.
    For face index i, WL[...,i] is the state extrapolated from the cell
    on the left/below, and WR[...,i] is the state extrapolated from the
    cell on the right/above.  Boundary ghost-cell slopes are set to zero.
    This avoids using np.roll's wrap-around values as physical neighbours.
    """
    if idir not in (0, 1):
        raise ValueError("idir must be 0 (x) or 1 (y)")
    axis = idir + 1
    d_minus = W - np.roll(W, 1, axis=axis)
    d_plus = np.roll(W, -1, axis=axis) - W
    slope = minmod(d_minus, d_plus)
    # Ghost-zone slopes are not reconstructed: boundary faces are first order.
    sl = [slice(None)] * W.ndim
    sl[axis] = 0
    slope[tuple(sl)] = 0.0
    sl[axis] = -1
    slope[tuple(sl)] = 0.0
    W_left_from_cell = W + 0.5 * slope
    W_right_from_cell = W - 0.5 * slope
    W_L = W_left_from_cell
    W_R = np.roll(W_right_from_cell, -1, axis=axis)
    return W_L, W_R


def physical_flux(U: Array, gamma: float, idir: int) -> Tuple[Array, Array]:
    """Return primitive state and Euler flux in x (idir=0) or y (idir=1)."""
    W = cons_to_prim(U, gamma)
    rho, u, v, p = W
    energy = U[3]
    F = np.zeros_like(U)
    if idir == 0:
        F[0] = rho*u
        F[1] = rho*u*u + p
        F[2] = rho*u*v
        F[3] = u*(energy+p)
    elif idir == 1:
        F[0] = rho*v
        F[1] = rho*u*v
        F[2] = rho*v*v + p
        F[3] = v*(energy+p)
    else:
        raise ValueError("idir must be 0 (x) or 1 (y)")
    return W, F


def hll_flux(U_L: Array, U_R: Array, idir: int, gamma: float) -> Array:
    """HLL approximate Riemann flux using bounding acoustic wave speeds."""
    W_L, F_L = physical_flux(U_L, gamma, idir)
    W_R, F_R = physical_flux(U_R, gamma, idir)
    rhoL, uL, vL, pL = W_L
    rhoR, uR, vR, pR = W_R
    vnL, vnR = (uL, uR) if idir == 0 else (vL, vR)
    cL = sound_speed(pL, rhoL, gamma)
    cR = sound_speed(pR, rhoR, gamma)
    sL = np.minimum(vnL-cL, vnR-cR)
    sR = np.maximum(vnL+cL, vnR+cR)
    sLb, sRb = sL[None, ...], sR[None, ...]
    denom = np.where(np.abs(sRb-sLb) > TINY, sRb-sLb, TINY)
    f_hll = (sRb*F_L - sLb*F_R + sLb*sRb*(U_R-U_L)) / denom
    return np.where(sLb >= 0.0, F_L,
                    np.where(sRb <= 0.0, F_R, f_hll))


def hllc_flux(U_L: Array, U_R: Array, idir: int, gamma: float) -> Array:
    """HLLC flux; restores the contact wave omitted by the HLL model."""
    W_L, F_L = physical_flux(U_L, gamma, idir)
    W_R, F_R = physical_flux(U_R, gamma, idir)
    rhoL, uL, vL, pL = W_L
    rhoR, uR, vR, pR = W_R
    if idir == 0:
        wL, wtL, wR, wtR = uL, vL, uR, vR
    elif idir == 1:
        wL, wtL, wR, wtR = vL, uL, vR, uR
    else:
        raise ValueError("idir must be 0 (x) or 1 (y)")
    cL = sound_speed(pL, rhoL, gamma)
    cR = sound_speed(pR, rhoR, gamma)
    sL = np.minimum(wL-cL, wR-cR)
    sR = np.maximum(wL+cL, wR+cR)
    den = rhoL*(sL-wL) - rhoR*(sR-wR)
    den = np.where(np.abs(den) > TINY, den, np.sign(den)*TINY + (den == 0)*TINY)
    sM = (pR-pL + rhoL*wL*(sL-wL) - rhoR*wR*(sR-wR)) / den
    pM = np.maximum(pL + rhoL*(sL-wL)*(sM-wL), TINY)
    rM_L = rhoL*(sL-wL) / np.where(np.abs(sL-sM)>TINY, sL-sM, TINY)
    rM_R = rhoR*(sR-wR) / np.where(np.abs(sR-sM)>TINY, sR-sM, TINY)

    U_star_L = np.zeros_like(U_L)
    U_star_R = np.zeros_like(U_R)
    U_star_L[0], U_star_R[0] = rM_L, rM_R
    if idir == 0:
        U_star_L[1], U_star_L[2] = rM_L*sM, rM_L*wtL
        U_star_R[1], U_star_R[2] = rM_R*sM, rM_R*wtR
    else:
        U_star_L[1], U_star_L[2] = rM_L*wtL, rM_L*sM
        U_star_R[1], U_star_R[2] = rM_R*wtR, rM_R*sM
    U_star_L[3] = pM/(gamma-1.0) + 0.5*rM_L*(sM*sM+wtL*wtL)
    U_star_R[3] = pM/(gamma-1.0) + 0.5*rM_R*(sM*sM+wtR*wtR)

    sLb, sRb, sMb = sL[None,...], sR[None,...], sM[None,...]
    F = np.zeros_like(F_L)
    left = sL >= 0.0
    left_star = (sL < 0.0) & (sM >= 0.0)
    right_star = (sM < 0.0) & (sR > 0.0)
    right = sR <= 0.0
    F = np.where(left[None,...], F_L, F)
    F = np.where(left_star[None,...], F_L+sLb*(U_star_L-U_L), F)
    F = np.where(right_star[None,...], F_R+sRb*(U_star_R-U_R), F)
    F = np.where(right[None,...], F_R, F)
    return F


@dataclass
class RunConfig:
    nx: int = 200
    ny: int = 100
    xmin: float = 0.0
    xmax: float = 10.0
    ymin: float = 0.0
    ymax: float = 5.0
    gamma: float = 1.4
    final_time: float = 5.0
    cfl: float = 0.45
    initial_condition: str = "riemann"
    scheme: str = "hllc"
    boundary: str = "outflow"


class Godunov2D:
    """Finite-volume solver for the two-dimensional compressible Euler equations."""

    def __init__(self, nx: int, ny: int, xmin: float, xmax: float,
                 ymin: float, ymax: float, gamma: float = 1.4,
                 initial_condition: str = "two_clouds", scheme: str = "hll",
                 bc: str = "outflow") -> None:
        if nx < 3 or ny < 3:
            raise ValueError("nx and ny must be at least 3 physical cells")
        if gamma <= 1.0 or xmax <= xmin or ymax <= ymin:
            raise ValueError("Require gamma>1 and increasing domain bounds")
        if scheme.lower() not in ("hll", "hllc"):
            raise ValueError("scheme must be 'hll' or 'hllc'")
        if bc not in ("outflow", "periodic"):
            raise ValueError("bc must be 'outflow' or 'periodic'")

        self.nx, self.ny = nx, ny  # physical cell counts
        self.nxt, self.nyt = nx+2, ny+2  # include one ghost layer
        self.xmin, self.xmax = xmin, xmax
        self.ymin, self.ymax = ymin, ymax
        self.gamma = gamma
        self.initial_condition = initial_condition
        self.scheme = scheme.lower()
        self.bc = bc
        self.dx = (xmax-xmin)/nx
        self.dy = (ymax-ymin)/ny
        # Cell-centred coordinates; indices 1:-1 are the physical cells.
        self.x = xmin + (np.arange(self.nxt)-0.5)*self.dx
        self.y = ymin + (np.arange(self.nyt)-0.5)*self.dy
        self.xx, self.yy = np.meshgrid(self.x, self.y, indexing="ij")
        self.U = np.zeros((4, self.nxt, self.nyt), dtype=float)
        self.W = np.zeros_like(self.U)
        self.t = 0.0
        self._initial_conditions()
        self.U = self.apply_bc(self.U)
        self.W = cons_to_prim(self.U, self.gamma)
        self.fignames: List[str] = []

    def _set_region(self, mask: Array, W_state: Array) -> None:
        """Set a physical-cell mask to a primitive state and its conserved form."""
        self.W[:, 1:-1, 1:-1][:, mask] = W_state[:, None]
        U_state = prim_to_cons(W_state, self.gamma)
        self.U[:, 1:-1, 1:-1][:, mask] = U_state[:, None]

    def _initial_conditions(self) -> None:
        """Populate physical cells for the selected test problem."""
        ambient = np.array([0.01, 0.0, 0.0, 0.01])
        self.W[:, 1:-1, 1:-1] = ambient[:, None, None]
        self.U[:, 1:-1, 1:-1] = prim_to_cons(ambient, self.gamma)[:, None, None]
        X, Y = self.xx[1:-1, 1:-1], self.yy[1:-1, 1:-1]

        if self.initial_condition in ("two_clouds", "presentation_parameters"):
            if self.initial_condition == "two_clouds":
                left = np.array([2.5, 2.0, 0.0, 0.01]); right = np.array([0.5, -1.0, 0.0, 0.01])
                r_left, r_right = 0.4, 1.0
            else:
                left = np.array([1.0, 1.0, 0.0, 1.0]); right = np.array([0.125, -1.0, 0.0, 0.1])
                r_left, r_right = 1.5, 1.5
            maskL = (X-2.5)**2 + (Y-2.5)**2 < r_left**2
            maskR = (X-7.5)**2 + (Y-2.5)**2 < r_right**2
            self._set_region(maskL, left)
            self._set_region(maskR, right)
        elif self.initial_condition == "riemann":
            left = np.array([1.0, 2.0, 0.0, 1.0])
            right = np.array([0.125, 0.0, 0.0, 0.1])
            mid_x = self.xmin + 0.2*(self.xmax-self.xmin)
            self._set_region(X < mid_x, left)
            self._set_region(X >= mid_x, right)
        else:
            raise ValueError("initial_condition must be 'two_clouds', 'presentation_parameters', or 'riemann'")

    def apply_bc(self, U: Array) -> Array:
        """Fill the one-cell ghost layer in-place and return U."""
        if self.bc == "outflow":
            U[:, 0, 1:-1] = U[:, 1, 1:-1]
            U[:, -1, 1:-1] = U[:, -2, 1:-1]
            U[:, 1:-1, 0] = U[:, 1:-1, 1]
            U[:, 1:-1, -1] = U[:, 1:-1, -2]
        else:  # periodic
            U[:, 0, 1:-1] = U[:, -2, 1:-1]
            U[:, -1, 1:-1] = U[:, 1, 1:-1]
            U[:, 1:-1, 0] = U[:, 1:-1, -2]
            U[:, 1:-1, -1] = U[:, 1:-1, 1]
        # Fill corners consistently from the nearest already-filled ghost edges.
        U[:, 0, 0] = U[:, 0, 1]; U[:, 0, -1] = U[:, 0, -2]
        U[:, -1, 0] = U[:, -1, 1]; U[:, -1, -1] = U[:, -1, -2]
        return U

    def compute_dt(self, CFL: float = 0.45) -> float:
        """Compute dt = CFL*min(dx/max(|u|+c), dy/max(|v|+c))."""
        W = cons_to_prim(self.U, self.gamma)
        rho, u, v, p = W[:, 1:-1, 1:-1]
        c = sound_speed(p, rho, self.gamma)
        sx = float(np.max(np.abs(u)+c))
        sy = float(np.max(np.abs(v)+c))
        if sx <= 0.0 or sy <= 0.0:
            raise RuntimeError("Non-positive characteristic speed; cannot compute a CFL time step")
        return CFL * min(self.dx/sx, self.dy/sy)

    def get_godunov_update(self, U_in: Array) -> Array:
        """Return RHS(U) = -div(F,G), with shape identical to U_in.

        Fluxes are computed on reconstructed left/right states at faces.
        The divergence is evaluated only in physical cells; ghost-cell RHS
        values remain zero and are refilled by apply_bc after each update.
        """
        U = self.apply_bc(U_in.copy())
        W = cons_to_prim(U, self.gamma)
        WLx, WRx = reconstruct(W, idir=0)
        ULx, URx = prim_to_cons(WLx, self.gamma), prim_to_cons(WRx, self.gamma)
        flux_solver = hll_flux if self.scheme == "hll" else hllc_flux
        Fx = flux_solver(ULx, URx, idir=0, gamma=self.gamma)
        WLy, WRy = reconstruct(W, idir=1)
        ULy, URy = prim_to_cons(WLy, self.gamma), prim_to_cons(WRy, self.gamma)
        Fy = flux_solver(ULy, URy, idir=1, gamma=self.gamma)

        rhs = np.zeros_like(U)
        rhs[:, 1:-1, 1:-1] = -(
            (Fx[:, 1:-1, 1:-1]-Fx[:, 0:-2, 1:-1])/self.dx
            + (Fy[:, 1:-1, 1:-1]-Fy[:, 1:-1, 0:-2])/self.dy
        )
        return rhs

    def advance(self, dt: float) -> Array:
        """Advance U by dt using the two-stage SSP-RK2 (Heun) update."""
        rhs_n = self.get_godunov_update(self.U)
        U_stage = self.U + dt*rhs_n
        U_stage = self.apply_bc(U_stage)
        rhs_stage = self.get_godunov_update(U_stage)
        U_new = self.U + 0.5*dt*(rhs_n+rhs_stage)
        return self.apply_bc(U_new)

    def step(self, dt: float) -> None:
        """Perform one complete time step and synchronize U, W and t."""
        if dt <= 0.0:
            raise ValueError("dt must be positive")
        self.U = self.advance(dt)
        self.U = self.apply_bc(self.U)
        self.W = cons_to_prim(self.U, self.gamma)
        self.t += dt

    def plot_overview(self, output_dir: str = "Godunov") -> str:
        """Save density, pressure and speed maps; return the image filename."""
        if plt is None:
            raise ImportError("matplotlib is required for plotting")
        out = Path(output_dir); out.mkdir(parents=True, exist_ok=True)
        W = self.W[:, 1:-1, 1:-1]
        U = self.U[:, 1:-1, 1:-1]
        fields = (U[0], W[3], np.sqrt(W[1]**2+W[2]**2))
        titles = ("Density", "Pressure", "Velocity magnitude")
        fig, axes = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)
        for ax, field, title in zip(axes, fields, titles):
            im = ax.imshow(field.T, origin="lower", extent=(self.xmin,self.xmax,self.ymin,self.ymax), aspect="auto")
            ax.set_title(title); ax.set_xlabel("x"); ax.set_ylabel("y")
            fig.colorbar(im, ax=ax)
        fig.suptitle(f"t={self.t:.5e} | {self.scheme.upper()} | {self.bc}")
        filename = out / f"overview_{self.initial_condition}_{self.scheme}_{self.bc}_{len(self.fignames):05d}.png"
        fig.savefig(filename, dpi=140); plt.close(fig)
        self.fignames.append(str(filename))
        return str(filename)


def run(config: RunConfig) -> Godunov2D:
    """Run one simulation, save periodic plots, and return the final solver."""
    sim = Godunov2D(config.nx, config.ny, config.xmin, config.xmax,
                    config.ymin, config.ymax, config.gamma,
                    config.initial_condition, config.scheme, config.boundary)
    plot_interval = config.final_time/50.0 if config.final_time > 0 else 1.0
    next_plot = 0.0
    sim.plot_overview()
    while sim.t < config.final_time:
        dt = sim.compute_dt(config.cfl)
        dt = min(dt, config.final_time-sim.t)  # land exactly on final_time
        if dt < 1.0e-12:
            raise RuntimeError(f"Time step became too small ({dt:.3e}) at t={sim.t:.6e}")
        sim.step(dt)
        if sim.t >= next_plot + plot_interval - 1.0e-14:
            sim.plot_overview()
            next_plot += plot_interval
    sim.plot_overview()  # always save the final state
    return sim


if __name__ == "__main__":
    # Edit this block to define the test.  For a first run, use the Riemann test.
    cfg = RunConfig(nx=200, ny=100, final_time=1.0, cfl=0.4,
                    initial_condition="riemann", scheme="hllc", boundary="outflow")
    simulation = run(cfg)
    print(f"Finished at t={simulation.t:.6e}; saved {len(simulation.fignames)} plots.")
