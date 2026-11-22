# HPC Master Class — Projects

Two projects from the course.

| Folder | Topic |
| --- | --- |
| [`PIC/`](PIC) | 1D hybrid particle-in-cell plasma code (C++, OpenMP) |
| [`Hydro/`](Hydro) | 2D Euler solver (Godunov finite volumes): collision of two gas clouds (Python) |

The course slides and LaTeX notes are kept locally and are not in the repository.

## PIC — Hybrid particle-in-cell code

A 1D hybrid-kinetic plasma code: the ions are macro-particles, the electrons a massless isothermal
fluid. Positions are 1D and velocities 3D (1D3V), the box is periodic, and the units are normalised
(`μ0 = e = m_i = 1`).

One time step: the ions are pushed (Boris) in the fields gathered from the grid, deposited on the
grid (density, bulk velocity), then Ampère gives the current, the generalised Ohm's law gives `E`
(`-u×B`, Hall, electron pressure, resistivity, hyper-resistivity) and Faraday advances `B`, on a
staggered Yee grid. The step is done by iterated Crank–Nicolson: two predictions, one correction.
`push` and `deposit` run in parallel with OpenMP (one copy of the deposited sums per thread, then a
reduction).

```sh
cd PIC
cmake -S . -B build && cmake --build build
OMP_NUM_THREADS=4 ./build/hybrid sound     # or uniform, decay, timing
```

Each run writes `<case>.csv` (time, `By`, `Bz`, density) and prints where the time went (particles,
copies, fields) with the limit of Amdahl's law.

| Case | What happens |
| --- | --- |
| `uniform` | cold plasma at rest in a uniform field: nothing changes |
| `sound` | a small density ripple travels at the sound speed `sqrt(T_e)` |
| `decay` | a magnetic field fades at the rate `η k² + ν k⁴` |
| `timing` | 128 000 warm protons, only to measure the speed |

Layout of `PIC/src/`: `vec3.hpp` (3D vector), `grid` (periodic Yee grid, derivatives, shape
function), `population` (the ions: loader, push, deposit), `field_solver` (Ampère, Ohm, Faraday),
`hybrid` (one time step), `initial_conditions` and `main.cpp` (the cases and the program),
`stopwatch.hpp` (timing). The tests (Python, through pybind11) are kept locally.

## Hydro — Two-cloud collision

A 2D finite-volume solver of the Euler equations: MUSCL reconstruction with the minmod limiter,
HLL or HLLC Riemann solver, SSP-RK2 time stepping, CFL time step, outflow or periodic boundaries.
Two cold, dense clouds in pressure balance with a hot medium fly at each other.

```sh
cd Hydro
python3 gudanov2d.py        # the collision, saves collision_*.png
python3 sod_test.py         # Sod shock tube against the exact solution, saves sod_test.png
```

- `gudanov2d.py` — the solver and the collision run.
- `exact_riemann.py` — exact solution of the 1D Riemann problem, for the test.
- `sod_test.py` — the test: Sod shock tube with HLL and HLLC against the exact solution.
- `gudanov2d_professor.py` — the professor's reference code, kept for comparison.

Needs `numpy` and `matplotlib`.
