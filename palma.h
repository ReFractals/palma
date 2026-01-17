/**
 * @file palma.h
 * @brief PALMA - Parallel Algebra Library for Max-plus Applications
 * 
 * A lightweight, high-performance tropical (max-plus/min-plus) algebra library
 * optimized for ARM-based embedded systems like Raspberry Pi.
 * 
 * Tropical algebra replaces standard arithmetic:
 *   - Max-plus semiring: a ⊕ b = max(a,b), a ⊗ b = a + b
 *   - Min-plus semiring: a ⊕ b = min(a,b), a ⊗ b = a + b
 *   - Bottleneck semiring: a ⊕ b = max(a,b), a ⊗ b = min(a,b)
 *   - Boolean semiring: a ⊕ b = a OR b, a ⊗ b = a AND b
 * 
 * This "linearizes" many optimization problems into matrix equations.
 * 
 * @author Gnankan Landry Regis N'guessan
 *         Axiom Research Group
 *         Department of Applied Mathematics and Computational Science,
 *         The Nelson Mandela African Institution of Science and Technology (NM-AIST),
 *         Arusha, Tanzania
 *         African Institute for Mathematical Sciences (AIMS),
 *         Research and Innovation Centre (RIC), Kigali, Rwanda
 * @email  rnguessan@aimsric.org
 * 
 * @version 1.0.0
 * @date    2024
 * @license MIT
 * 
 * @copyright Copyright (c) 2024 Gnankan Landry Regis N'guessan
 *            All rights reserved.
 */

#ifndef PALMA_H
#define PALMA_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <limits.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/*============================================================================
 * VERSION INFORMATION
 *============================================================================*/

#define PALMA_VERSION_MAJOR  1
#define PALMA_VERSION_MINOR  0
#define PALMA_VERSION_PATCH  0
#define PALMA_VERSION_STRING "1.0.0"

/*============================================================================
 * CONFIGURATION
 *============================================================================*/

/** Use NEON SIMD optimizations on ARM (auto-detected, can override) */
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    #define PALMA_USE_NEON 1
#else
    #define PALMA_USE_NEON 0
#endif

/** Use OpenMP for multi-core parallelization */
#ifdef _OPENMP
    #define PALMA_USE_OPENMP 1
#else
    #define PALMA_USE_OPENMP 0
#endif

/** Integer type for tropical values (32-bit for NEON alignment) */
typedef int32_t palma_val_t;

/** Index type for sparse matrices */
typedef uint32_t palma_idx_t;

/*============================================================================
 * CONSTANTS
 *============================================================================*/

/** Special values */
#define PALMA_NEG_INF       INT32_MIN   /**< Negative infinity (-∞) */
#define PALMA_POS_INF       INT32_MAX   /**< Positive infinity (+∞) */
#define PALMA_ZERO          0           /**< Multiplicative identity */

/** Default tolerances */
#define PALMA_DEFAULT_MAX_ITER  1000    /**< Default max iterations */
#define PALMA_DEFAULT_TOL       1       /**< Default tolerance for convergence */

/*============================================================================
 * ERROR HANDLING
 *============================================================================*/

/**
 * @brief Error codes returned by PALMA functions
 */
typedef enum {
    PALMA_SUCCESS = 0,              /**< Operation completed successfully */
    PALMA_ERR_NULL_PTR = -1,        /**< NULL pointer argument */
    PALMA_ERR_INVALID_DIM = -2,     /**< Invalid matrix dimensions */
    PALMA_ERR_OUT_OF_MEMORY = -3,   /**< Memory allocation failed */
    PALMA_ERR_INVALID_ARG = -4,     /**< Invalid argument value */
    PALMA_ERR_NOT_SQUARE = -5,      /**< Matrix must be square */
    PALMA_ERR_NOT_CONVERGED = -6,   /**< Algorithm did not converge */
    PALMA_ERR_FILE_OPEN = -7,       /**< Failed to open file */
    PALMA_ERR_FILE_READ = -8,       /**< Failed to read from file */
    PALMA_ERR_FILE_WRITE = -9,      /**< Failed to write to file */
    PALMA_ERR_FILE_FORMAT = -10,    /**< Invalid file format */
    PALMA_ERR_INDEX_BOUNDS = -11,   /**< Index out of bounds */
    PALMA_ERR_SPARSE_FORMAT = -12,  /**< Invalid sparse matrix format */
    PALMA_ERR_UNSUPPORTED = -13     /**< Unsupported operation */
} palma_error_t;

/**
 * @brief Get human-readable error message
 * @param err Error code
 * @return Descriptive string for the error
 */
const char* palma_strerror(palma_error_t err);

/**
 * @brief Last error code (thread-local where supported)
 */
palma_error_t palma_get_last_error(void);

/**
 * @brief Set last error code
 */
void palma_set_last_error(palma_error_t err);

/**
 * @brief Clear last error
 */
void palma_clear_error(void);

/*============================================================================
 * SEMIRING TYPES
 *============================================================================*/

/**
 * @brief Semiring type selector
 * 
 * Each semiring defines different ⊕ (add) and ⊗ (mul) operations:
 * 
 * | Semiring    | a ⊕ b      | a ⊗ b      | Zero (ε) | One (e) | Use Case          |
 * |-------------|------------|------------|----------|---------|-------------------|
 * | MAXPLUS     | max(a,b)   | a + b      | -∞       | 0       | Longest paths     |
 * | MINPLUS     | min(a,b)   | a + b      | +∞       | 0       | Shortest paths    |
 * | MAXMIN      | max(a,b)   | min(a,b)   | -∞       | +∞      | Bottleneck/bandwidth |
 * | MINMAX      | min(a,b)   | max(a,b)   | +∞       | -∞      | Reliability       |
 * | BOOLEAN     | a OR b     | a AND b    | 0        | 1       | Reachability      |
 */
typedef enum {
    PALMA_MAXPLUS = 0,   /**< Max-plus: (max, +) for longest paths, scheduling */
    PALMA_MINPLUS = 1,   /**< Min-plus: (min, +) for shortest paths */
    PALMA_MAXMIN  = 2,   /**< Max-min (bottleneck): for bandwidth, capacity */
    PALMA_MINMAX  = 3,   /**< Min-max: for reliability paths */
    PALMA_BOOLEAN = 4    /**< Boolean: (OR, AND) for reachability analysis */
} palma_semiring_t;

/**
 * @brief Get the zero element (additive identity) for a semiring
 * @param semiring Semiring type
 * @return Zero element ε such that a ⊕ ε = a
 */
palma_val_t palma_zero(palma_semiring_t semiring);

/**
 * @brief Get the one element (multiplicative identity) for a semiring
 * @param semiring Semiring type
 * @return One element e such that a ⊗ e = a
 */
palma_val_t palma_one(palma_semiring_t semiring);

/**
 * @brief Tropical addition: a ⊕ b
 * @param a First operand
 * @param b Second operand
 * @param semiring Semiring type
 * @return a ⊕ b according to semiring
 */
palma_val_t palma_add(palma_val_t a, palma_val_t b, palma_semiring_t semiring);

/**
 * @brief Tropical multiplication: a ⊗ b
 * @param a First operand
 * @param b Second operand
 * @param semiring Semiring type
 * @return a ⊗ b according to semiring
 */
palma_val_t palma_mul(palma_val_t a, palma_val_t b, palma_semiring_t semiring);

/**
 * @brief Check if value is the zero element
 * @param a Value to check
 * @param semiring Semiring type
 * @return true if a equals the semiring's zero element
 */
bool palma_is_zero(palma_val_t a, palma_semiring_t semiring);

/**
 * @brief Get semiring name as string
 * @param semiring Semiring type
 * @return Human-readable name
 */
const char* palma_semiring_name(palma_semiring_t semiring);

/*============================================================================
 * DENSE MATRIX STRUCTURE
 *============================================================================*/

/**
 * @brief Dense tropical matrix structure
 * 
 * Row-major storage for cache efficiency.
 * Rows are 16-byte aligned for NEON optimization.
 */
typedef struct {
    palma_val_t *data;      /**< Matrix data (row-major, aligned) */
    size_t rows;            /**< Number of rows */
    size_t cols;            /**< Number of columns */
    size_t stride;          /**< Row stride (>= cols, aligned) */
    bool owns_data;         /**< Whether to free data on destroy */
} palma_matrix_t;

/*============================================================================
 * SPARSE MATRIX STRUCTURE (CSR FORMAT)
 *============================================================================*/

/**
 * @brief Sparse tropical matrix in Compressed Sparse Row (CSR) format
 * 
 * CSR format stores only non-zero (non-infinite) elements:
 *   - values[]:    Non-zero values
 *   - col_idx[]:   Column index for each value
 *   - row_ptr[]:   Index into values/col_idx for start of each row
 * 
 * For row i, elements are in values[row_ptr[i]] to values[row_ptr[i+1]-1]
 * 
 * Memory: O(nnz + n) instead of O(n²) for dense
 */
typedef struct {
    palma_val_t *values;    /**< Non-zero values (length = nnz) */
    palma_idx_t *col_idx;   /**< Column indices (length = nnz) */
    palma_idx_t *row_ptr;   /**< Row pointers (length = rows + 1) */
    size_t rows;            /**< Number of rows */
    size_t cols;            /**< Number of columns */
    size_t nnz;             /**< Number of non-zero entries */
    size_t capacity;        /**< Allocated capacity for values/col_idx */
    palma_semiring_t semiring; /**< Semiring (determines what "zero" means) */
} palma_sparse_t;

/*============================================================================
 * DENSE MATRIX LIFECYCLE
 *============================================================================*/

/**
 * @brief Create a new dense matrix (uninitialized values)
 * @param rows Number of rows
 * @param cols Number of columns
 * @return Pointer to new matrix, or NULL on failure
 */
palma_matrix_t* palma_matrix_create(size_t rows, size_t cols);

/**
 * @brief Create a dense matrix initialized to semiring zero
 * @param rows Number of rows
 * @param cols Number of columns
 * @param semiring Semiring type (determines zero value)
 * @return Pointer to new matrix, or NULL on failure
 */
palma_matrix_t* palma_matrix_create_zero(size_t rows, size_t cols, palma_semiring_t semiring);

/**
 * @brief Create a tropical identity matrix
 * @param n Size of square matrix
 * @param semiring Semiring type
 * @return Pointer to new identity matrix, or NULL on failure
 */
palma_matrix_t* palma_matrix_create_identity(size_t n, palma_semiring_t semiring);

/**
 * @brief Create a matrix wrapping existing data (no copy)
 * @param data Pointer to existing data (must remain valid)
 * @param rows Number of rows
 * @param cols Number of columns
 * @param stride Row stride
 * @return Pointer to matrix wrapper, or NULL on failure
 */
palma_matrix_t* palma_matrix_wrap(palma_val_t *data, size_t rows, size_t cols, size_t stride);

/**
 * @brief Clone (deep copy) a matrix
 * @param src Source matrix
 * @return Pointer to new copy, or NULL on failure
 */
palma_matrix_t* palma_matrix_clone(const palma_matrix_t *src);

/**
 * @brief Destroy a matrix and free resources
 * @param mat Matrix to destroy (NULL-safe)
 */
void palma_matrix_destroy(palma_matrix_t *mat);

/*============================================================================
 * DENSE MATRIX ACCESS
 *============================================================================*/

/**
 * @brief Get element at (row, col)
 */
static inline palma_val_t palma_matrix_get(const palma_matrix_t *mat, size_t row, size_t col) {
    return mat->data[row * mat->stride + col];
}

/**
 * @brief Set element at (row, col)
 */
static inline void palma_matrix_set(palma_matrix_t *mat, size_t row, size_t col, palma_val_t val) {
    mat->data[row * mat->stride + col] = val;
}

/**
 * @brief Get pointer to start of a row
 */
static inline palma_val_t* palma_matrix_row(palma_matrix_t *mat, size_t row) {
    return &mat->data[row * mat->stride];
}

/**
 * @brief Safe element access with bounds checking
 * @param mat Matrix
 * @param row Row index
 * @param col Column index
 * @param out Output value pointer
 * @return PALMA_SUCCESS or error code
 */
palma_error_t palma_matrix_get_safe(const palma_matrix_t *mat, size_t row, size_t col, palma_val_t *out);

/**
 * @brief Safe element set with bounds checking
 * @param mat Matrix
 * @param row Row index
 * @param col Column index
 * @param val Value to set
 * @return PALMA_SUCCESS or error code
 */
palma_error_t palma_matrix_set_safe(palma_matrix_t *mat, size_t row, size_t col, palma_val_t val);

/*============================================================================
 * SPARSE MATRIX LIFECYCLE
 *============================================================================*/

/**
 * @brief Create an empty sparse matrix with given capacity
 * @param rows Number of rows
 * @param cols Number of columns
 * @param capacity Initial capacity for non-zero elements
 * @param semiring Semiring type
 * @return Pointer to new sparse matrix, or NULL on failure
 */
palma_sparse_t* palma_sparse_create(size_t rows, size_t cols, size_t capacity, palma_semiring_t semiring);

/**
 * @brief Create sparse matrix from dense matrix
 * @param dense Source dense matrix
 * @param semiring Semiring type (determines which values are "zero")
 * @return Pointer to new sparse matrix, or NULL on failure
 */
palma_sparse_t* palma_sparse_from_dense(const palma_matrix_t *dense, palma_semiring_t semiring);

/**
 * @brief Convert sparse matrix to dense
 * @param sparse Source sparse matrix
 * @return Pointer to new dense matrix, or NULL on failure
 */
palma_matrix_t* palma_sparse_to_dense(const palma_sparse_t *sparse);

/**
 * @brief Clone a sparse matrix
 * @param src Source sparse matrix
 * @return Pointer to new copy, or NULL on failure
 */
palma_sparse_t* palma_sparse_clone(const palma_sparse_t *src);

/**
 * @brief Destroy a sparse matrix and free resources
 * @param sp Sparse matrix to destroy (NULL-safe)
 */
void palma_sparse_destroy(palma_sparse_t *sp);

/*============================================================================
 * SPARSE MATRIX ACCESS & MODIFICATION
 *============================================================================*/

/**
 * @brief Get element from sparse matrix
 * @param sp Sparse matrix
 * @param row Row index
 * @param col Column index
 * @return Element value (returns semiring zero if not stored)
 */
palma_val_t palma_sparse_get(const palma_sparse_t *sp, size_t row, size_t col);

/**
 * @brief Set element in sparse matrix
 * 
 * Note: Setting a value to semiring zero will NOT remove the entry.
 * Use palma_sparse_compress() to remove zeros after bulk modifications.
 * 
 * @param sp Sparse matrix
 * @param row Row index
 * @param col Column index
 * @param val Value to set
 * @return PALMA_SUCCESS or error code
 */
palma_error_t palma_sparse_set(palma_sparse_t *sp, size_t row, size_t col, palma_val_t val);

/**
 * @brief Add a new entry (assumes entry doesn't exist - faster than set)
 * @param sp Sparse matrix
 * @param row Row index
 * @param col Column index
 * @param val Value to add
 * @return PALMA_SUCCESS or error code
 */
palma_error_t palma_sparse_add_entry(palma_sparse_t *sp, size_t row, size_t col, palma_val_t val);

/**
 * @brief Remove zero entries and sort columns within each row
 * @param sp Sparse matrix
 * @return PALMA_SUCCESS or error code
 */
palma_error_t palma_sparse_compress(palma_sparse_t *sp);

/**
 * @brief Get number of non-zeros in a row
 * @param sp Sparse matrix
 * @param row Row index
 * @return Number of non-zero entries in row
 */
size_t palma_sparse_row_nnz(const palma_sparse_t *sp, size_t row);

/**
 * @brief Get sparsity ratio (fraction of zeros)
 * @param sp Sparse matrix
 * @return Sparsity as a value between 0.0 and 1.0
 */
double palma_sparse_sparsity(const palma_sparse_t *sp);

/*============================================================================
 * DENSE MATRIX OPERATIONS
 *============================================================================*/

/**
 * @brief Tropical matrix multiplication: C = A ⊗ B
 * 
 * C[i,j] = ⊕_k (A[i,k] ⊗ B[k,j])
 * 
 * @param A Left matrix (m × n)
 * @param B Right matrix (n × p)
 * @param semiring Semiring type
 * @return Result matrix (m × p), or NULL on failure
 */
palma_matrix_t* palma_matrix_mul(const palma_matrix_t *A, const palma_matrix_t *B, 
                                  palma_semiring_t semiring);

/**
 * @brief In-place tropical matrix multiplication: C = A ⊗ B
 * @param C Pre-allocated result matrix (m × p)
 * @param A Left matrix (m × n)
 * @param B Right matrix (n × p)
 * @param semiring Semiring type
 * @return PALMA_SUCCESS or error code
 */
palma_error_t palma_matrix_mul_into(palma_matrix_t *C, const palma_matrix_t *A, 
                                     const palma_matrix_t *B, palma_semiring_t semiring);

/**
 * @brief Tropical matrix addition: C = A ⊕ B (element-wise)
 * @param A First matrix
 * @param B Second matrix (same dimensions)
 * @param semiring Semiring type
 * @return Result matrix, or NULL on failure
 */
palma_matrix_t* palma_matrix_add(const palma_matrix_t *A, const palma_matrix_t *B,
                                  palma_semiring_t semiring);

/**
 * @brief Tropical matrix power: A^n
 * 
 * Uses binary exponentiation: O(log n) matrix multiplications.
 * 
 * @param A Square matrix
 * @param n Power (n >= 0)
 * @param semiring Semiring type
 * @return A^n, or NULL on failure
 */
palma_matrix_t* palma_matrix_power(const palma_matrix_t *A, unsigned int n,
                                    palma_semiring_t semiring);

/**
 * @brief Tropical closure (Kleene star): A* = I ⊕ A ⊕ A² ⊕ A³ ⊕ ...
 * 
 * Computes the reflexive-transitive closure.
 * For path problems: A*[i,j] = optimal path from i to j.
 * 
 * @param A Square matrix
 * @param semiring Semiring type
 * @return A*, or NULL on failure
 */
palma_matrix_t* palma_matrix_closure(const palma_matrix_t *A, palma_semiring_t semiring);

/**
 * @brief Transitive closure: A+ = A ⊕ A² ⊕ A³ ⊕ ...
 * 
 * Like closure but without identity (requires at least one step).
 * 
 * @param A Square matrix
 * @param semiring Semiring type
 * @return A+, or NULL on failure
 */
palma_matrix_t* palma_matrix_transitive_closure(const palma_matrix_t *A, palma_semiring_t semiring);

/*============================================================================
 * SPARSE MATRIX OPERATIONS
 *============================================================================*/

/**
 * @brief Sparse tropical matrix multiplication: C = A ⊗ B
 * @param A Left sparse matrix
 * @param B Right sparse matrix
 * @return Result sparse matrix, or NULL on failure
 */
palma_sparse_t* palma_sparse_mul(const palma_sparse_t *A, const palma_sparse_t *B);

/**
 * @brief Sparse matrix-vector multiplication: y = A ⊗ x
 * @param A Sparse matrix (m × n)
 * @param x Input vector (length n)
 * @param y Output vector (length m, pre-allocated)
 * @return PALMA_SUCCESS or error code
 */
palma_error_t palma_sparse_matvec(const palma_sparse_t *A, const palma_val_t *x, 
                                   palma_val_t *y);

/**
 * @brief Sparse matrix closure: A*
 * @param A Square sparse matrix
 * @return Sparse closure, or NULL on failure
 */
palma_sparse_t* palma_sparse_closure(const palma_sparse_t *A);

/*============================================================================
 * VECTOR OPERATIONS
 *============================================================================*/

/**
 * @brief Tropical matrix-vector multiplication: y = A ⊗ x
 * @param A Dense matrix (m × n)
 * @param x Input vector (length n)
 * @param y Output vector (length m, pre-allocated)
 * @param semiring Semiring type
 * @return PALMA_SUCCESS or error code
 */
palma_error_t palma_matvec(const palma_matrix_t *A, const palma_val_t *x, 
                            palma_val_t *y, palma_semiring_t semiring);

/**
 * @brief Iterate system: x(k+1) = A ⊗ x(k)
 * @param A System matrix (square)
 * @param x State vector (modified in place)
 * @param n Number of iterations
 * @param semiring Semiring type
 * @return PALMA_SUCCESS or error code
 */
palma_error_t palma_iterate(const palma_matrix_t *A, palma_val_t *x, 
                             unsigned int n, palma_semiring_t semiring);

/**
 * @brief Tropical vector dot product: x · y = ⊕_i (x[i] ⊗ y[i])
 * @param x First vector
 * @param y Second vector
 * @param len Vector length
 * @param semiring Semiring type
 * @return Dot product result
 */
palma_val_t palma_dot(const palma_val_t *x, const palma_val_t *y, 
                       size_t len, palma_semiring_t semiring);

/*============================================================================
 * EIGENVALUE & EIGENVECTOR COMPUTATION
 *============================================================================*/

/**
 * @brief Compute the maximum cycle mean (tropical eigenvalue)
 * 
 * For max-plus, this is the growth rate of A^k as k → ∞.
 * Equivalent to the largest average weight over all cycles.
 * Uses Karp's algorithm: O(n³) time.
 * 
 * @param A Square matrix
 * @param semiring Semiring type
 * @return Tropical eigenvalue λ, or PALMA_NEG_INF if acyclic
 */
palma_val_t palma_eigenvalue(const palma_matrix_t *A, palma_semiring_t semiring);

/**
 * @brief Compute tropical eigenvector
 * 
 * Finds v such that A ⊗ v = λ ⊗ v (where λ is the eigenvalue).
 * Uses power iteration with normalization.
 * 
 * @param A Square matrix
 * @param eigenvector Output vector (length n, pre-allocated)
 * @param eigenvalue Output eigenvalue (optional, can be NULL)
 * @param semiring Semiring type
 * @param max_iter Maximum iterations (0 for default)
 * @return PALMA_SUCCESS or error code (including PALMA_ERR_NOT_CONVERGED)
 */
palma_error_t palma_eigenvector(const palma_matrix_t *A, palma_val_t *eigenvector,
                                 palma_val_t *eigenvalue, palma_semiring_t semiring,
                                 unsigned int max_iter);

/**
 * @brief Compute all nodes on critical cycles
 * 
 * Returns indices of nodes that participate in cycles achieving
 * the maximum cycle mean.
 * 
 * @param A Square matrix
 * @param critical_nodes Output array (length n, set to 1 if on critical cycle)
 * @param semiring Semiring type
 * @return Number of critical nodes, or negative error code
 */
int palma_critical_nodes(const palma_matrix_t *A, int *critical_nodes,
                          palma_semiring_t semiring);

/*============================================================================
 * GRAPH ALGORITHMS
 *============================================================================*/

/**
 * @brief All-pairs optimal paths
 * 
 * Computes optimal path weights between all pairs of nodes.
 * Uses tropical matrix closure.
 * 
 * @param adj Adjacency matrix with edge weights
 * @param semiring PALMA_MINPLUS for shortest, PALMA_MAXPLUS for longest
 * @return Distance matrix, or NULL on failure
 */
palma_matrix_t* palma_all_pairs_paths(const palma_matrix_t *adj, palma_semiring_t semiring);

/**
 * @brief Single-source optimal paths (Bellman-Ford style)
 * @param adj Adjacency matrix
 * @param source Source vertex index
 * @param dist Output distance vector (length n, pre-allocated)
 * @param semiring PALMA_MINPLUS for shortest, PALMA_MAXPLUS for longest
 * @return PALMA_SUCCESS or error code
 */
palma_error_t palma_single_source_paths(const palma_matrix_t *adj, size_t source,
                                         palma_val_t *dist, palma_semiring_t semiring);

/**
 * @brief Reachability analysis using Boolean semiring
 * @param adj Adjacency matrix (non-zero = edge exists)
 * @return Reachability matrix (1 if path exists, 0 otherwise), or NULL on failure
 */
palma_matrix_t* palma_reachability(const palma_matrix_t *adj);

/**
 * @brief Bottleneck paths (maximum capacity paths)
 * 
 * Finds paths maximizing the minimum edge weight (bandwidth/capacity).
 * Uses max-min semiring.
 * 
 * @param adj Adjacency matrix with edge capacities
 * @return Capacity matrix, or NULL on failure
 */
palma_matrix_t* palma_bottleneck_paths(const palma_matrix_t *adj);

/*============================================================================
 * SCHEDULING APPLICATIONS
 *============================================================================*/

/**
 * @brief Task scheduling system
 * 
 * Models precedence-constrained scheduling:
 * x(k+1) = A ⊗ x(k) ⊕ b
 * 
 * A[i,j] = duration of j if j precedes i
 * x[i] = completion time of task i
 * b[i] = external input (ready time)
 */
typedef struct {
    palma_matrix_t *system;     /**< System matrix A */
    palma_val_t *state;         /**< Current state x */
    palma_val_t *input;         /**< External input b */
    size_t n_tasks;             /**< Number of tasks */
    palma_semiring_t semiring;  /**< Semiring (usually MAXPLUS) */
    char **task_names;          /**< Optional task names */
} palma_scheduler_t;

/**
 * @brief Create a scheduler
 * @param n_tasks Number of tasks
 * @param use_maxplus If true, compute latest times; if false, earliest times
 * @return Scheduler instance, or NULL on failure
 */
palma_scheduler_t* palma_scheduler_create(size_t n_tasks, bool use_maxplus);

/**
 * @brief Destroy a scheduler
 * @param sched Scheduler to destroy (NULL-safe)
 */
void palma_scheduler_destroy(palma_scheduler_t *sched);

/**
 * @brief Set task name (optional, for display)
 * @param sched Scheduler
 * @param task Task index
 * @param name Task name (copied)
 * @return PALMA_SUCCESS or error code
 */
palma_error_t palma_scheduler_set_name(palma_scheduler_t *sched, size_t task, const char *name);

/**
 * @brief Add precedence constraint: 'from' must complete before 'to' starts
 * @param sched Scheduler
 * @param from Predecessor task index
 * @param to Successor task index
 * @param duration Processing time of 'from'
 * @return PALMA_SUCCESS or error code
 */
palma_error_t palma_scheduler_add_constraint(palma_scheduler_t *sched, size_t from,
                                              size_t to, palma_val_t duration);

/**
 * @brief Set ready time for a task
 * @param sched Scheduler
 * @param task Task index
 * @param ready_time When task becomes available
 * @return PALMA_SUCCESS or error code
 */
palma_error_t palma_scheduler_set_ready_time(palma_scheduler_t *sched, size_t task,
                                              palma_val_t ready_time);

/**
 * @brief Solve the schedule
 * @param sched Scheduler
 * @param max_iter Maximum iterations (0 for default = n_tasks)
 * @return Number of iterations used, or negative error code
 */
int palma_scheduler_solve(palma_scheduler_t *sched, unsigned int max_iter);

/**
 * @brief Get completion time for a task
 * @param sched Scheduler (after solve)
 * @param task Task index
 * @return Completion time
 */
palma_val_t palma_scheduler_get_completion(const palma_scheduler_t *sched, size_t task);

/**
 * @brief Compute cycle time (maximum eigenvalue) for periodic schedules
 * @param sched Scheduler with cyclic dependencies
 * @return Cycle time, or PALMA_NEG_INF if acyclic
 */
palma_val_t palma_scheduler_cycle_time(const palma_scheduler_t *sched);

/**
 * @brief Compute throughput (1/cycle_time) for periodic schedules
 * @param sched Scheduler
 * @return Throughput in tasks per time unit, or 0 if not cyclic
 */
double palma_scheduler_throughput(const palma_scheduler_t *sched);

/**
 * @brief Find critical path
 * @param sched Scheduler (after solve)
 * @param path Output array for task indices on critical path
 * @param max_len Maximum length of path array
 * @return Length of critical path, or negative error code
 */
int palma_scheduler_critical_path(const palma_scheduler_t *sched, size_t *path, size_t max_len);

/*============================================================================
 * FILE I/O
 *============================================================================*/

/**
 * @brief Save dense matrix to CSV file
 * @param mat Matrix to save
 * @param filename Output file path
 * @param semiring Semiring (for formatting infinities)
 * @return PALMA_SUCCESS or error code
 */
palma_error_t palma_matrix_save_csv(const palma_matrix_t *mat, const char *filename,
                                     palma_semiring_t semiring);

/**
 * @brief Load dense matrix from CSV file
 * @param filename Input file path
 * @param semiring Semiring (for parsing infinities)
 * @return Loaded matrix, or NULL on failure
 */
palma_matrix_t* palma_matrix_load_csv(const char *filename, palma_semiring_t semiring);

/**
 * @brief Save dense matrix to binary file
 * @param mat Matrix to save
 * @param filename Output file path
 * @return PALMA_SUCCESS or error code
 */
palma_error_t palma_matrix_save_binary(const palma_matrix_t *mat, const char *filename);

/**
 * @brief Load dense matrix from binary file
 * @param filename Input file path
 * @return Loaded matrix, or NULL on failure
 */
palma_matrix_t* palma_matrix_load_binary(const char *filename);

/**
 * @brief Save sparse matrix to CSV file (COO format: row,col,value)
 * @param sp Sparse matrix to save
 * @param filename Output file path
 * @return PALMA_SUCCESS or error code
 */
palma_error_t palma_sparse_save_csv(const palma_sparse_t *sp, const char *filename);

/**
 * @brief Load sparse matrix from CSV file (COO format)
 * @param filename Input file path
 * @param semiring Semiring type
 * @return Loaded sparse matrix, or NULL on failure
 */
palma_sparse_t* palma_sparse_load_csv(const char *filename, palma_semiring_t semiring);

/**
 * @brief Export matrix to GraphViz DOT format
 * @param mat Matrix (adjacency matrix)
 * @param filename Output file path
 * @param semiring Semiring (for labeling edges)
 * @param node_names Optional array of node names (can be NULL)
 * @return PALMA_SUCCESS or error code
 */
palma_error_t palma_matrix_export_dot(const palma_matrix_t *mat, const char *filename,
                                       palma_semiring_t semiring, const char **node_names);

/*============================================================================
 * NEON-OPTIMIZED OPERATIONS (ARM only)
 *============================================================================*/

#if PALMA_USE_NEON

/**
 * @brief NEON-optimized matrix multiplication
 */
palma_error_t palma_matrix_mul_neon(palma_matrix_t *C, const palma_matrix_t *A,
                                     const palma_matrix_t *B, palma_semiring_t semiring);

/**
 * @brief NEON-optimized matrix-vector multiplication
 */
palma_error_t palma_matvec_neon(const palma_matrix_t *A, const palma_val_t *x,
                                 palma_val_t *y, palma_semiring_t semiring);

#endif /* PALMA_USE_NEON */

/*============================================================================
 * UTILITY FUNCTIONS
 *============================================================================*/

/**
 * @brief Print dense matrix to file/stdout
 * @param mat Matrix to print
 * @param name Optional name string
 * @param semiring Semiring (for formatting)
 * @param fp Output file (use stdout for console)
 */
void palma_matrix_print(const palma_matrix_t *mat, const char *name,
                         palma_semiring_t semiring, FILE *fp);

/**
 * @brief Print sparse matrix to file/stdout
 * @param sp Sparse matrix to print
 * @param name Optional name string
 * @param fp Output file (use stdout for console)
 */
void palma_sparse_print(const palma_sparse_t *sp, const char *name, FILE *fp);

/**
 * @brief Print vector to file/stdout
 * @param vec Vector to print
 * @param len Vector length
 * @param name Optional name string
 * @param semiring Semiring (for formatting)
 * @param fp Output file
 */
void palma_vector_print(const palma_val_t *vec, size_t len, const char *name,
                         palma_semiring_t semiring, FILE *fp);

/**
 * @brief Get library version string
 * @return Version string "major.minor.patch"
 */
const char* palma_version(void);

/**
 * @brief Get library version components
 * @param major Output major version
 * @param minor Output minor version
 * @param patch Output patch version
 */
void palma_version_components(int *major, int *minor, int *patch);

/**
 * @brief Check if NEON optimizations are available
 * @return true if NEON is available and enabled
 */
bool palma_has_neon(void);

/**
 * @brief Check if OpenMP is available
 * @return true if OpenMP is available and enabled
 */
bool palma_has_openmp(void);

/**
 * @brief Get build configuration string
 * @return Description of compile-time options
 */
const char* palma_build_config(void);

#ifdef __cplusplus
}
#endif

#endif /* PALMA_H */
