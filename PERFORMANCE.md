# Performance Guide

This document provides guidance on achieving optimal performance with PALMA.

## Table of Contents

- [Hardware Considerations](#hardware-considerations)
- [Build Optimization](#build-optimization)
- [Algorithm Selection](#algorithm-selection)
- [Memory Optimization](#memory-optimization)
- [NEON SIMD](#neon-simd)
- [Benchmarking](#benchmarking)
- [Real-Time Guidelines](#real-time-guidelines)

---

## Hardware Considerations

### Recommended Platforms

| Platform | NEON | Performance | Notes |
|----------|------|-------------|-------|
| Raspberry Pi 4 | ✓ | Excellent | Primary development target |
| Raspberry Pi 5 | ✓ | Better | Faster CPU, larger cache |
| Raspberry Pi 3 | ✓ | Good | 32-bit or 64-bit |
| Pi Zero 2 W | ✓ | Moderate | Low power option |
| x86_64 | ✗ | Good | Scalar fallback |

### Cache Hierarchy Impact

PALMA performance is often **memory-bound**, not compute-bound. Cache utilization is critical:

| Pi 4 Cache | Size | Impact |
|------------|------|--------|
| L1 Data | 32 KB | ~1 cycle access |
| L1 Instruction | 48 KB | Code must fit |
| L2 | 1 MB | ~10 cycle access |
| RAM | 8 GB | ~100 cycle access |

**Rule of thumb**: Keep working matrices in L2 cache (~1 MB) for best performance.

---

## Build Optimization

### Release Build

For production use:
```bash
make release
```

This enables:
- `-O3` optimization
- Link-time optimization (LTO)
- Architecture-specific tuning
- NEON auto-detection

### Build Variants

```bash
make all       # Standard build with NEON if available
make scalar    # Disable NEON (for comparison)
make debug     # Debug symbols, no optimization
make openmp    # Multi-threaded (experimental)
```

### Compiler Flags

The Makefile automatically sets optimal flags for ARM:
```
-mcpu=cortex-a72 -mtune=cortex-a72 -mfpu=neon-fp-armv8
```

---

## Algorithm Selection

### Operation Costs

| Operation | Complexity | 64×64 Time | 128×128 Time |
|-----------|------------|------------|--------------|
| Matrix-vector | O(n²) | 12 μs | 48 μs |
| Matrix multiply | O(n³) | 314 μs | 10.5 ms |
| Closure | O(n³) | 1 ms | 7.8 ms |
| Eigenvalue | O(n³) | ~1 ms | ~8 ms |

### When to Use Sparse

Use sparse matrices when:
- Sparsity > 50% (>50% zeros)
- Matrix size > 64×64
- Memory is constrained

```c
// Convert if beneficial
if (palma_matrix_sparsity(dense) > 0.5) {
    palma_sparse_t *sparse = palma_sparse_from_dense(dense, PALMA_MINPLUS);
    // Use sparse operations
}
```

### SSSP vs Closure

For single-source shortest paths:

| Method | Complexity | When to Use |
|--------|------------|-------------|
| SSSP (iterate) | O(n² · iterations) | Single source, dense |
| Closure | O(n³) | All pairs needed |
| Sparse SSSP | O(nnz · iterations) | Sparse graphs |

---

## Memory Optimization

### Dense Matrix Memory

```
Memory = rows × cols × 4 bytes
```

| Size | Memory |
|------|--------|
| 32×32 | 4 KB |
| 64×64 | 16 KB |
| 128×128 | 64 KB |
| 256×256 | 256 KB |
| 512×512 | 1 MB |

### Sparse Matrix Memory

```
Memory ≈ nnz × 12 bytes + rows × 4 bytes
```

At 30% density:
| Size | Dense | Sparse | Savings |
|------|-------|--------|---------|
| 64×64 | 16 KB | ~7 KB | 56% |
| 128×128 | 64 KB | ~28 KB | 56% |
| 256×256 | 256 KB | ~110 KB | 57% |

### Memory Tips

1. **Reuse matrices**: Avoid repeated allocation
2. **In-place operations**: Use when possible
3. **Destroy promptly**: Free memory after use
4. **Static allocation**: For real-time systems

---

## NEON SIMD

### How It Helps

NEON processes 4 × 32-bit integers per instruction:

```
Standard:  a[0]+b[0], a[1]+b[1], a[2]+b[2], a[3]+b[3]  (4 instructions)
NEON:      vaddq_s32(a, b)                              (1 instruction)
```

### Measured Speedups

| Operation | Size | Scalar | NEON | Speedup |
|-----------|------|--------|------|---------|
| MatVec | 64 | 22 μs | 12 μs | 1.8× |
| MatMul | 64 | 510 μs | 314 μs | 1.6× |
| Closure | 64 | 1.6 ms | 1.0 ms | 1.6× |

### Why Not 4×?

Theoretical 4× speedup is reduced by:
- Memory bandwidth limits
- Loop overhead
- Non-vectorizable code
- Cache misses

### Checking NEON Status

```c
#ifdef PALMA_USE_NEON
    printf("NEON enabled\n");
#else
    printf("NEON disabled (scalar mode)\n");
#endif
```

---

## Benchmarking

### Built-in Benchmark

```bash
make run-benchmark
```

This measures:
- Matrix multiplication (various sizes)
- Matrix-vector products
- Closure computation
- Sparse operations
- All semirings

### Custom Benchmarking

```c
#include "palma.h"
#include <time.h>

// Warm up cache
for (int i = 0; i < 10; i++) {
    palma_matrix_mul(A, B, PALMA_MAXPLUS);
}

// Measure
clock_t start = clock();
for (int i = 0; i < 100; i++) {
    palma_matrix_t *C = palma_matrix_mul(A, B, PALMA_MAXPLUS);
    palma_matrix_destroy(C);
}
clock_t end = clock();

double avg_us = (double)(end - start) / CLOCKS_PER_SEC * 1e6 / 100;
printf("Average: %.1f μs\n", avg_us);
```

### MOPS Calculation

```
MOPS = (operations per call) / (time in μs)
     = (n³ for matmul) / (time_μs)
```

---

## Real-Time Guidelines

### Worst-Case Execution Time

For real-time systems, use measured worst-case times:

| Operation | 32×32 | 64×64 | 128×128 |
|-----------|-------|-------|---------|
| MatVec | 5 μs | 20 μs | 80 μs |
| MatMul | 80 μs | 500 μs | 15 ms |
| Closure | 200 μs | 1.5 ms | 12 ms |

### Deadline Mapping

| Deadline | Max Matrix Size |
|----------|-----------------|
| 100 μs | 16×16 |
| 1 ms | 32×32 |
| 10 ms | 64×64 |
| 100 ms | 128×128 |

### Tips for Real-Time

1. **Pre-allocate**: No malloc in hot paths
2. **Bound iterations**: Set max_iter for eigenvalue
3. **Use sparse**: Faster for large, sparse graphs
4. **Profile**: Measure on target hardware
5. **Cache warming**: First call may be slower

### Static Allocation Mode

For safety-critical systems:
```c
// Pre-allocate workspace
static palma_val_t workspace[MAX_N * MAX_N];

// Use static buffer (not yet in public API)
palma_matrix_t mat = {
    .rows = n,
    .cols = n,
    .data = workspace
};
```

---

## Optimization Checklist

- [ ] Use release build (`make release`)
- [ ] Match matrix size to cache (≤256×256 for L2)
- [ ] Use sparse for >50% sparsity
- [ ] Reuse matrices to avoid allocation
- [ ] Warm up cache before timing
- [ ] Profile on target hardware
- [ ] Consider SSSP vs closure based on needs
