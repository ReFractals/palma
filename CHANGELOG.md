# Changelog

All notable changes to PALMA will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-01-17

### Added
- Initial public release
- Five tropical semirings: max-plus, min-plus, max-min, min-max, Boolean
- Dense matrix operations with ARM NEON SIMD optimization
- Sparse matrix support (CSR format)
- Tropical eigenvalue computation via Karp's algorithm
- Eigenvector computation via power iteration
- High-level scheduling API
- Graph algorithms: SSSP, APSP, bottleneck paths, reachability
- File I/O: CSV, binary, GraphViz DOT export
- Comprehensive examples: graphs, scheduling, eigenvalues, benchmarks
- Full API documentation
- MIT license

### Performance
- Up to 1.8× speedup with NEON over scalar
- Up to 11.9× faster than classical Bellman-Ford
- Peak 2,274 MOPS on Raspberry Pi 4 (64×64 max-min)

### Platforms
- Primary: Raspberry Pi 4 (ARM Cortex-A72)
- Tested: Raspberry Pi 3, Pi Zero 2 W
- Portable scalar fallback for non-ARM platforms

## [Unreleased]

### Planned
- OpenMP multi-threading support
- Python bindings
- Additional examples (Petri nets, control systems)
- RISC-V vector extension support
