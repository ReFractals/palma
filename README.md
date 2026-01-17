<p align="center">
  <img src="https://raw.githubusercontent.com/ReFractals/palma/main/assets/palma_logo.svg" alt="PALMA Logo" width="200">
</p>
<h1 align="center">PALMA</h1>
<h3 align="center">Parallel Algebra Library for Max-plus Applications</h3>

<p align="center">
  <strong>A lightweight, high-performance tropical algebra library for ARM-based embedded systems</strong>
</p>

<p align="center">
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-green.svg" alt="License: MIT"></a>
  <a href="https://github.com/axiom-research/palma/releases"><img src="https://img.shields.io/badge/version-1.0.0-blue.svg" alt="Version"></a>
  <img src="https://img.shields.io/badge/platform-ARM%20%7C%20Raspberry%20Pi-red.svg" alt="Platform">
  <img src="https://img.shields.io/badge/C-C99-orange.svg" alt="C99">
</p>

---

## Overview

**PALMA** is a dependency-free C99 library for tropical (max-plus/min-plus) linear algebra, optimized for ARM NEON SIMD. It provides efficient implementations of semiring matrix operations, graph algorithms, eigenvalue computation, and real-time scheduling, all in a single lightweight package.

### Key Features

- **Five Semirings**: Max-plus, min-plus, max-min (bottleneck), min-max, Boolean
- **NEON Optimized**: Up to 1.8× speedup with ARM SIMD vectorization
- **Sparse Support**: CSR format for memory-efficient large graphs
- **Zero Dependencies**: Pure C99, no external libraries required
- **Real-Time Ready**: Deterministic execution, suitable for control systems

## Quick Start

```bash
# Clone and build
git clone https://github.com/axiom-research/palma.git
cd palma
make all

# Run examples
./build/bin/example_graphs
./build/bin/example_scheduling
./build/bin/benchmark
```

## What is Tropical Algebra?

Tropical algebra redefines arithmetic to "linearize" optimization problems:

| Semiring | ⊕ (add) | ⊗ (mult) | Zero | One | Application |
|----------|---------|----------|------|-----|-------------|
| **Max-Plus** | max | + | −∞ | 0 | Scheduling, critical paths |
| **Min-Plus** | min | + | +∞ | 0 | Shortest paths |
| **Max-Min** | max | min | −∞ | +∞ | Bandwidth/bottleneck |
| **Boolean** | OR | AND | 0 | 1 | Reachability |

**Key insight**: Floyd-Warshall, Bellman-Ford, and scheduling algorithms are all tropical matrix operations with different semirings.

## Example: Shortest Paths

```c
#include "palma.h"

int main() {
    // Create graph adjacency matrix
    palma_matrix_t *adj = palma_matrix_create_zero(4, 4, PALMA_MINPLUS);
    
    palma_matrix_set(adj, 0, 1, 5);   // Edge 0→1, weight 5
    palma_matrix_set(adj, 1, 2, 3);   // Edge 1→2, weight 3
    palma_matrix_set(adj, 2, 3, 2);   // Edge 2→3, weight 2
    
    // Compute ALL shortest paths (tropical closure = Floyd-Warshall)
    palma_matrix_t *dist = palma_matrix_closure(adj, PALMA_MINPLUS);
    
    printf("Shortest 0→3: %d\n", palma_matrix_get(dist, 0, 3)); // Output: 10
    
    palma_matrix_destroy(adj);
    palma_matrix_destroy(dist);
    return 0;
}
```

## Example: Task Scheduling

```c
#include "palma.h"

int main() {
    palma_scheduler_t *sched = palma_scheduler_create(4, true);
    
    // Task 0 (10ms) → Task 1 (20ms) → Task 2 (15ms) → Task 3
    palma_scheduler_add_constraint(sched, 0, 1, 10);
    palma_scheduler_add_constraint(sched, 1, 2, 20);
    palma_scheduler_add_constraint(sched, 2, 3, 15);
    
    palma_scheduler_solve(sched, 0);
    
    printf("Makespan: %d\n", palma_scheduler_makespan(sched));
    printf("Cycle time: %d\n", palma_scheduler_cycle_time(sched));
    
    palma_scheduler_destroy(sched);
    return 0;
}
```

## Performance (Raspberry Pi 4)

| Matrix Size | MatMul | MatVec | Closure | MOPS |
|-------------|--------|--------|---------|------|
| 32×32 | 44 μs | 3 μs | 138 μs | 1,495 |
| 64×64 | 314 μs | 12 μs | 1,001 μs | 1,668 |
| 128×128 | 10.5 ms | 48 μs | 7.8 ms | 399 |

- **NEON speedup**: Up to 1.8× over scalar
- **SSSP speedup**: Up to 11.9× vs Bellman-Ford
- **Sparse efficiency**: 3.5× faster at 50% sparsity

## API Overview

### Matrix Operations
```c
palma_matrix_t* palma_matrix_create(rows, cols);
palma_matrix_t* palma_matrix_mul(A, B, semiring);
palma_matrix_t* palma_matrix_closure(A, semiring);  // Kleene star
void palma_matvec(A, x, y, semiring);               // y = A ⊗ x
```

### Graph Algorithms
```c
void palma_single_source_paths(A, source, dist, semiring);
palma_matrix_t* palma_all_pairs_paths(A, semiring);
palma_matrix_t* palma_bottleneck_paths(A);
palma_matrix_t* palma_reachability(A);
```

### Eigenvalue (Cycle Time)
```c
double palma_eigenvalue(A, semiring);  // Maximum cycle mean
void palma_eigenvector(A, v, &lambda, semiring, max_iter);
```

### Scheduling
```c
palma_scheduler_t* palma_scheduler_create(n_tasks, cyclic);
void palma_scheduler_add_constraint(sched, from, to, weight);
int palma_scheduler_solve(sched, start_time);
palma_val_t palma_scheduler_makespan(sched);
palma_val_t palma_scheduler_cycle_time(sched);
```

## Build Options

```bash
make all        # Default build (auto-detects NEON)
make release    # Optimized with LTO
make debug      # With symbols and sanitizers
make scalar     # Without NEON (for comparison)
make openmp     # Multi-threaded
make install    # Install to /usr/local
```

## Documentation

- [API Reference](docs/API.md)
- [Mathematical Background](docs/MATH.md)
- [Performance Guide](docs/PERFORMANCE.md)
- [Examples](examples/)

## Applications

- **Real-time systems**: Task scheduling, control loop timing
- **Network routing**: Shortest paths, maximum bandwidth
- **Manufacturing**: Production line throughput, bottleneck analysis
- **Discrete event systems**: Petri nets, timed automata

## Author

**Gnankan Landry Regis N'guessan**

- Email: [rnguessan@aimsric.org](mailto:rnguessan@aimsric.org)
- Affiliations:
  - Axiom Research Group
  - NM-AIST, Arusha, Tanzania
  - AIMS Research and Innovation Centre, Kigali, Rwanda

## Citation

```bibtex
@software{palma2026,
  author = {N'guessan, Gnankan Landry Regis},
  title = {{PALMA}: Parallel Algebra Library for Max-plus Applications},
  year = {2026},
  url = {https://github.com/axiom-research/palma},
  version = {1.0.0}
}
```

## License

MIT License — see [LICENSE](LICENSE) for details.

## Contributing

Contributions welcome! See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

---

<p align="center">
  <sub>PALMA — Bringing tropical mathematics to embedded systems</sub>
</p>
