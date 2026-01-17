# PALMA API Reference

This document provides a complete reference for the PALMA public API.

## Table of Contents

- [Core Types](#core-types)
- [Semiring Constants](#semiring-constants)
- [Matrix Operations](#matrix-operations)
- [Sparse Matrix Operations](#sparse-matrix-operations)
- [Vector Operations](#vector-operations)
- [Graph Algorithms](#graph-algorithms)
- [Eigenvalue/Eigenvector](#eigenvalueeigenvector)
- [Scheduling](#scheduling)
- [File I/O](#file-io)
- [Utility Functions](#utility-functions)

---

## Core Types

### `palma_val_t`
```c
typedef int32_t palma_val_t;
```
The fundamental value type for tropical algebra operations. Uses 32-bit signed integers.

### `palma_idx_t`
```c
typedef uint32_t palma_idx_t;
```
Index type for sparse matrix operations.

### `palma_matrix_t`
```c
typedef struct {
    size_t rows;
    size_t cols;
    palma_val_t *data;
} palma_matrix_t;
```
Dense matrix structure. Data is stored in row-major order.

### `palma_sparse_t`
```c
typedef struct {
    size_t rows;
    size_t cols;
    size_t nnz;           // Number of non-zero elements
    palma_val_t *values;  // Non-zero values
    palma_idx_t *col_idx; // Column indices
    palma_idx_t *row_ptr; // Row pointers (CSR format)
    palma_semiring_t semiring;
} palma_sparse_t;
```
Sparse matrix in Compressed Sparse Row (CSR) format.

---

## Semiring Constants

```c
typedef enum {
    PALMA_MAXPLUS = 0,  // ⊕ = max, ⊗ = +, ε = -∞, e = 0
    PALMA_MINPLUS = 1,  // ⊕ = min, ⊗ = +, ε = +∞, e = 0
    PALMA_MAXMIN  = 2,  // ⊕ = max, ⊗ = min, ε = -∞, e = +∞
    PALMA_MINMAX  = 3,  // ⊕ = min, ⊗ = max, ε = +∞, e = -∞
    PALMA_BOOLEAN = 4   // ⊕ = OR, ⊗ = AND, ε = 0, e = 1
} palma_semiring_t;
```

### Infinity Values
```c
#define PALMA_INF_POS  INT32_MAX  // +∞
#define PALMA_INF_NEG  INT32_MIN  // -∞
```

---

## Matrix Operations

### Creation

#### `palma_matrix_create`
```c
palma_matrix_t* palma_matrix_create(size_t rows, size_t cols);
```
Creates an uninitialized matrix. Returns `NULL` on allocation failure.

#### `palma_matrix_create_zero`
```c
palma_matrix_t* palma_matrix_create_zero(size_t rows, size_t cols, palma_semiring_t s);
```
Creates a matrix filled with the semiring's zero element (ε).

#### `palma_matrix_create_identity`
```c
palma_matrix_t* palma_matrix_create_identity(size_t n, palma_semiring_t s);
```
Creates an n×n identity matrix (diagonal = e, off-diagonal = ε).

#### `palma_matrix_clone`
```c
palma_matrix_t* palma_matrix_clone(const palma_matrix_t *A);
```
Creates a deep copy of matrix A.

### Destruction

#### `palma_matrix_destroy`
```c
void palma_matrix_destroy(palma_matrix_t *mat);
```
Frees all memory associated with the matrix.

### Element Access

#### `palma_matrix_get`
```c
palma_val_t palma_matrix_get(const palma_matrix_t *A, size_t i, size_t j);
```
Returns element A[i,j]. No bounds checking in release builds.

#### `palma_matrix_set`
```c
void palma_matrix_set(palma_matrix_t *A, size_t i, size_t j, palma_val_t v);
```
Sets element A[i,j] = v.

### Arithmetic

#### `palma_matrix_add`
```c
palma_matrix_t* palma_matrix_add(const palma_matrix_t *A, 
                                  const palma_matrix_t *B, 
                                  palma_semiring_t s);
```
Computes C = A ⊕ B (element-wise tropical addition).

#### `palma_matrix_mul`
```c
palma_matrix_t* palma_matrix_mul(const palma_matrix_t *A,
                                  const palma_matrix_t *B,
                                  palma_semiring_t s);
```
Computes C = A ⊗ B (tropical matrix multiplication).

#### `palma_matrix_power`
```c
palma_matrix_t* palma_matrix_power(const palma_matrix_t *A, 
                                    int k, 
                                    palma_semiring_t s);
```
Computes A^k using binary exponentiation.

#### `palma_matrix_closure`
```c
palma_matrix_t* palma_matrix_closure(const palma_matrix_t *A, 
                                      palma_semiring_t s);
```
Computes Kleene star: A* = I ⊕ A ⊕ A² ⊕ ... ⊕ A^(n-1).

---

## Vector Operations

#### `palma_matvec`
```c
void palma_matvec(const palma_matrix_t *A, 
                  const palma_val_t *x,
                  palma_val_t *y, 
                  palma_semiring_t s);
```
Computes y = A ⊗ x (tropical matrix-vector product).

#### `palma_matvec_neon`
```c
void palma_matvec_neon(const palma_matrix_t *A, 
                       const palma_val_t *x,
                       palma_val_t *y, 
                       palma_semiring_t s);
```
NEON-accelerated version (ARM only).

---

## Graph Algorithms

#### `palma_single_source_paths`
```c
void palma_single_source_paths(const palma_matrix_t *A, 
                                size_t src,
                                palma_val_t *dist, 
                                palma_semiring_t s);
```
Computes single-source paths from vertex `src`.
- Min-plus: shortest paths
- Max-plus: longest paths

#### `palma_all_pairs_paths`
```c
palma_matrix_t* palma_all_pairs_paths(const palma_matrix_t *A, 
                                       palma_semiring_t s);
```
Computes all-pairs paths (equivalent to `palma_matrix_closure`).

#### `palma_reachability`
```c
palma_matrix_t* palma_reachability(const palma_matrix_t *A);
```
Computes transitive closure (Boolean semiring).

#### `palma_bottleneck_paths`
```c
palma_matrix_t* palma_bottleneck_paths(const palma_matrix_t *A);
```
Computes maximum-capacity paths (max-min semiring).

---

## Eigenvalue/Eigenvector

#### `palma_eigenvalue`
```c
palma_val_t palma_eigenvalue(const palma_matrix_t *A, palma_semiring_t s);
```
Computes the tropical eigenvalue (maximum cycle mean) using Karp's algorithm.

#### `palma_eigenvector`
```c
void palma_eigenvector(const palma_matrix_t *A, 
                       palma_val_t *v,
                       palma_val_t *lambda, 
                       palma_semiring_t s, 
                       int max_iter);
```
Computes eigenvector via power iteration. Stores eigenvalue in `*lambda`.

#### `palma_critical_nodes`
```c
int palma_critical_nodes(const palma_matrix_t *A, 
                         int *critical, 
                         palma_semiring_t s);
```
Identifies nodes on critical cycles. Returns count of critical nodes.

---

## Scheduling

#### `palma_scheduler_create`
```c
palma_scheduler_t* palma_scheduler_create(size_t n_tasks, bool cyclic);
```
Creates a scheduler for n_tasks. Set `cyclic=true` for periodic systems.

#### `palma_scheduler_add_constraint`
```c
void palma_scheduler_add_constraint(palma_scheduler_t *sched,
                                    size_t from, 
                                    size_t to, 
                                    palma_val_t weight);
```
Adds precedence constraint: task `to` cannot start until `weight` time units after `from` starts.

#### `palma_scheduler_solve`
```c
int palma_scheduler_solve(palma_scheduler_t *sched, palma_val_t start_time);
```
Computes the schedule. Returns 0 on success.

#### `palma_scheduler_makespan`
```c
palma_val_t palma_scheduler_makespan(palma_scheduler_t *sched);
```
Returns the total schedule length (completion time of last task).

#### `palma_scheduler_cycle_time`
```c
palma_val_t palma_scheduler_cycle_time(palma_scheduler_t *sched);
```
Returns the cycle time for periodic systems (tropical eigenvalue).

#### `palma_scheduler_throughput`
```c
double palma_scheduler_throughput(palma_scheduler_t *sched);
```
Returns throughput = 1.0 / cycle_time.

#### `palma_scheduler_destroy`
```c
void palma_scheduler_destroy(palma_scheduler_t *sched);
```
Frees scheduler resources.

---

## File I/O

#### `palma_matrix_save_csv`
```c
int palma_matrix_save_csv(const palma_matrix_t *mat, 
                          const char *filename,
                          palma_semiring_t s);
```
Exports matrix to CSV file. Infinity values are written as "inf" or "-inf".

#### `palma_matrix_load_csv`
```c
palma_matrix_t* palma_matrix_load_csv(const char *filename, 
                                       palma_semiring_t s);
```
Imports matrix from CSV file.

#### `palma_matrix_export_dot`
```c
int palma_matrix_export_dot(const palma_matrix_t *mat,
                            const char *filename,
                            palma_semiring_t s,
                            const char **names);
```
Exports matrix as GraphViz DOT file. `names` is optional array of node labels.

---

## Utility Functions

#### `palma_semiring_zero`
```c
palma_val_t palma_semiring_zero(palma_semiring_t s);
```
Returns the zero element (ε) for the given semiring.

#### `palma_semiring_one`
```c
palma_val_t palma_semiring_one(palma_semiring_t s);
```
Returns the identity element (e) for the given semiring.

#### `palma_tropical_add`
```c
palma_val_t palma_tropical_add(palma_val_t a, palma_val_t b, palma_semiring_t s);
```
Computes a ⊕ b.

#### `palma_tropical_mul`
```c
palma_val_t palma_tropical_mul(palma_val_t a, palma_val_t b, palma_semiring_t s);
```
Computes a ⊗ b.
