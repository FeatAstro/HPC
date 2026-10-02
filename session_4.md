## Session 3
shape: which order? (first one)
unit test for the density
Larmor radius with fields (periodic boundary conditions -- test already done) --> maybe think about simple tests
test: we write the expectation
tests: continuously with the new implementations (refactor should not break the tests)
need very strict tests as much as possible

## Session 4
vectorization: one cpu, one "core" (SIMD)
particles: all x's, all y's, ...
SOA vs. AOS (structures of arrays vs. arrays of structures)

CPU threads: multiple cores with same memory (does no have to go to the network)
one thread per core (software abstraction)
particles same cell,  different threads --> no control on the results (no problem on reading, just on writing) (different floating points errors for a+b and b+a)  
Race condition (deposit particles) --> either u use different memory grids and then the reduction (NP operation) or atomic access (sequential order, not parallel when there is a concurrent access)  
threads same cache --> destroys performance (use cache alignement, padding,...) -- false sharing

order of execution (cannot control in parallelization)
if parallelization: unit tests in sequential code but also tests on parallelization


debugging: ddt, totalview or prints thread ID

noise same in threads if same seed (so different seeds on different processes/threads)

always make sure that the bug only happens in parallelization

load balancing/load imbalancing  
we cannot control how the load process  
distributed parallelism: multi processing  
local communication  
ghost node (process 1 wants a node in process 2) --> overlap precessors (exchange neighbours)  
global communication --> do not scale  
one processor gives all communication or receive all communication from the other processes  
implement these communications

MPI --> standard library (OpenMPI, MPICH, IntelMPI,...)

Scalable parallelism --> faster the more threads/processes it uses

Super important --> depressing Amdhal's law (speed is limited by the sequential code -- communication for example)

Gustanfson's law --> parallel time = 1, parallel to increase workload

strong and weak scalability (Amdhal Gustafson)

MPI barrier




recommend for parall: field equation and all the others on a single process, cut the particles into processes (deposit-->reduction on process 0--grid)

deposit on fields (reduction on process 0) --> this fraction will be sequential so scaling will not be as high as wanted

grid decomposition problem-->dealing with particles going to other processes (synchronize things)

tests: Adham's law (runtime)