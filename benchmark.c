/**
 * @file benchmark.c
 * @brief Performance Benchmark for PALMA Library
 * 
 * Measures execution time for tropical matrix operations across different
 * sizes and semirings. Essential for profiling on Raspberry Pi.
 * 
 * @author Gnankan Landry Regis N'guessan
 *         Axiom Research Group
 *         NM-AIST / AIMS-RIC
 * @email  rnguessan@aimsric.org
 */

#define _POSIX_C_SOURCE 199309L

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "palma.h"

/* Timing function */
static double get_time_us(void) {
#if defined(_POSIX_C_SOURCE) && _POSIX_C_SOURCE >= 199309L
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1e6 + ts.tv_nsec / 1e3;
#else
    return (double)clock() / CLOCKS_PER_SEC * 1e6;
#endif
}

/* Fill matrix with random values */
static void fill_random(palma_matrix_t *mat, palma_semiring_t semiring) {
    palma_val_t zero = palma_zero(semiring);
    
    for (size_t i = 0; i < mat->rows; i++) {
        for (size_t j = 0; j < mat->cols; j++) {
            /* 30% sparse */
            if (rand() % 100 < 30) {
                palma_matrix_set(mat, i, j, zero);
            } else {
                palma_matrix_set(mat, i, j, (rand() % 100) + 1);
            }
        }
    }
}

typedef struct {
    size_t size;
    double dense_mul_us;
    double sparse_mul_us;
    double matvec_us;
    double closure_us;
    double eigenvalue_us;
    int iterations;
} benchmark_result_t;

static benchmark_result_t run_benchmark(size_t n, int iterations, palma_semiring_t semiring) {
    benchmark_result_t result = {0};
    result.size = n;
    result.iterations = iterations;
    
    /* Create test matrices */
    palma_matrix_t *A = palma_matrix_create(n, n);
    palma_matrix_t *B = palma_matrix_create(n, n);
    palma_val_t *x = (palma_val_t*)malloc(n * sizeof(palma_val_t));
    palma_val_t *y = (palma_val_t*)malloc(n * sizeof(palma_val_t));
    
    if (!A || !B || !x || !y) {
        fprintf(stderr, "Memory allocation failed for size %zu\n", n);
        goto cleanup;
    }
    
    fill_random(A, semiring);
    fill_random(B, semiring);
    for (size_t i = 0; i < n; i++) {
        x[i] = rand() % 100;
    }
    
    /* Warm up */
    palma_matrix_t *C = palma_matrix_mul(A, B, semiring);
    palma_matrix_destroy(C);
    
    /* Benchmark dense matrix multiplication */
    double start = get_time_us();
    for (int i = 0; i < iterations; i++) {
        C = palma_matrix_mul(A, B, semiring);
        palma_matrix_destroy(C);
    }
    result.dense_mul_us = (get_time_us() - start) / iterations;
    
    /* Benchmark sparse multiplication */
    palma_sparse_t *spA = palma_sparse_from_dense(A, semiring);
    palma_sparse_t *spB = palma_sparse_from_dense(B, semiring);
    
    start = get_time_us();
    for (int i = 0; i < iterations; i++) {
        palma_sparse_t *spC = palma_sparse_mul(spA, spB);
        palma_sparse_destroy(spC);
    }
    result.sparse_mul_us = (get_time_us() - start) / iterations;
    
    palma_sparse_destroy(spA);
    palma_sparse_destroy(spB);
    
    /* Benchmark matrix-vector multiplication */
    int matvec_iters = iterations * 10;
    start = get_time_us();
    for (int i = 0; i < matvec_iters; i++) {
        palma_matvec(A, x, y, semiring);
    }
    result.matvec_us = (get_time_us() - start) / matvec_iters;
    
    /* Benchmark closure (smaller sizes only) */
    if (n <= 128) {
        int closure_iters = (n <= 32) ? iterations : (iterations / 4);
        if (closure_iters < 1) closure_iters = 1;
        
        start = get_time_us();
        for (int i = 0; i < closure_iters; i++) {
            palma_matrix_t *closure = palma_matrix_closure(A, semiring);
            palma_matrix_destroy(closure);
        }
        result.closure_us = (get_time_us() - start) / closure_iters;
    }
    
    /* Benchmark eigenvalue computation (smaller sizes only) */
    if (n <= 64) {
        int eigen_iters = (n <= 16) ? iterations : (iterations / 10);
        if (eigen_iters < 1) eigen_iters = 1;
        
        start = get_time_us();
        for (int i = 0; i < eigen_iters; i++) {
            palma_eigenvalue(A, semiring);
        }
        result.eigenvalue_us = (get_time_us() - start) / eigen_iters;
    }
    
cleanup:
    palma_matrix_destroy(A);
    palma_matrix_destroy(B);
    free(x);
    free(y);
    
    return result;
}

static void print_results(benchmark_result_t *results, int count, const char *title) {
    printf("\n%s\n", title);
    printf("%-6s | %10s | %10s | %10s | %10s | %10s | %10s\n",
           "Size", "Dense Mul", "Sparse Mul", "MatVec", "Closure", "Eigenval", "MOPS");
    printf("-------+------------+------------+------------+------------+------------+-----------\n");
    
    for (int i = 0; i < count; i++) {
        size_t n = results[i].size;
        double ops = 2.0 * n * n * n;
        double mops = ops / results[i].dense_mul_us;
        
        printf("%-6zu | %8.1f us | %8.1f us | %8.1f us |",
               n, results[i].dense_mul_us, results[i].sparse_mul_us, results[i].matvec_us);
        
        if (results[i].closure_us > 0) {
            printf(" %8.1f us |", results[i].closure_us);
        } else {
            printf(" %10s |", "N/A");
        }
        
        if (results[i].eigenvalue_us > 0) {
            printf(" %8.1f us |", results[i].eigenvalue_us);
        } else {
            printf(" %10s |", "N/A");
        }
        
        printf(" %8.1f\n", mops);
    }
    
    printf("\nMOPS = Million operations per second (2n^3 ops for nxn matrix multiply)\n");
}

int main(int argc, char *argv[]) {
    (void)argc; (void)argv;
    
    printf("==================================================================\n");
    printf("  PALMA - Parallel Algebra Library for Max-plus Applications\n");
    printf("  Performance Benchmark\n");
    printf("==================================================================\n\n");
    
    printf("Library: %s\n", palma_build_config());
    printf("Author: Gnankan Landry Regis N'guessan\n\n");
    
    /* Platform info */
#if defined(__aarch64__)
    printf("Architecture: ARM64 (aarch64)\n");
#elif defined(__arm__)
    printf("Architecture: ARM32\n");
#elif defined(__x86_64__)
    printf("Architecture: x86-64\n");
#else
    printf("Architecture: Unknown\n");
#endif
    
    printf("NEON SIMD: %s\n", palma_has_neon() ? "ENABLED" : "DISABLED");
    printf("OpenMP: %s\n\n", palma_has_openmp() ? "ENABLED" : "DISABLED");
    
    srand((unsigned int)time(NULL));
    
    /* Matrix sizes to test */
    size_t sizes[] = {8, 16, 32, 64, 128, 256, 512};
    int num_sizes = sizeof(sizes) / sizeof(sizes[0]);
    int base_iterations = 100;
    
    benchmark_result_t results_maxplus[7];
    benchmark_result_t results_minplus[7];
    benchmark_result_t results_maxmin[7];
    
    /* Run max-plus benchmarks */
    printf("Running max-plus benchmarks");
    fflush(stdout);
    for (int i = 0; i < num_sizes; i++) {
        int iters = base_iterations;
        if (sizes[i] >= 256) iters = 10;
        if (sizes[i] >= 512) iters = 3;
        
        results_maxplus[i] = run_benchmark(sizes[i], iters, PALMA_MAXPLUS);
        printf(".");
        fflush(stdout);
    }
    printf(" done\n");
    
    /* Run min-plus benchmarks */
    printf("Running min-plus benchmarks");
    fflush(stdout);
    for (int i = 0; i < num_sizes; i++) {
        int iters = base_iterations;
        if (sizes[i] >= 256) iters = 10;
        if (sizes[i] >= 512) iters = 3;
        
        results_minplus[i] = run_benchmark(sizes[i], iters, PALMA_MINPLUS);
        printf(".");
        fflush(stdout);
    }
    printf(" done\n");
    
    /* Run max-min (bottleneck) benchmarks */
    printf("Running max-min benchmarks");
    fflush(stdout);
    for (int i = 0; i < num_sizes; i++) {
        int iters = base_iterations;
        if (sizes[i] >= 256) iters = 10;
        if (sizes[i] >= 512) iters = 3;
        
        results_maxmin[i] = run_benchmark(sizes[i], iters, PALMA_MAXMIN);
        printf(".");
        fflush(stdout);
    }
    printf(" done\n");
    
    /* Print results */
    print_results(results_maxplus, num_sizes, 
                  "=== Max-Plus Semiring (scheduling, longest paths) ===");
    print_results(results_minplus, num_sizes, 
                  "=== Min-Plus Semiring (shortest paths) ===");
    print_results(results_maxmin, num_sizes, 
                  "=== Max-Min Semiring (bottleneck/bandwidth) ===");
    
    /* Memory usage */
    printf("\n=== Memory Usage ===\n");
    printf("%-6s | %12s | %s\n", "Size", "Dense (KB)", "Typical Use Case");
    printf("-------+--------------+-----------------------------\n");
    
    for (int i = 0; i < num_sizes; i++) {
        size_t n = sizes[i];
        double kb = (n * n * sizeof(palma_val_t)) / 1024.0;
        
        const char *use;
        if (n <= 16) use = "Small embedded MCU";
        else if (n <= 64) use = "Real-time scheduling";
        else if (n <= 256) use = "Network routing";
        else use = "Large-scale analysis";
        
        printf("%-6zu | %10.1f KB | %s\n", n, kb, use);
    }
    
    /* Raspberry Pi recommendations */
    printf("\n=== Raspberry Pi Recommendations ===\n");
    printf("Pi Zero:  Up to 256x256 dense, use sparse for larger\n");
    printf("Pi 3B+:   Up to 512x512 dense, NEON gives 2-3x speedup\n");
    printf("Pi 4:     Up to 1024x1024, enable OpenMP for 4x speedup\n");
    printf("Pi 5:     Up to 2048x2048, best with NEON+OpenMP\n");
    
    printf("\n=== Real-Time Constraints ===\n");
    printf("For 1ms deadline:   Use <=32x32 matrices\n");
    printf("For 10ms deadline:  Use <=128x128 matrices (with NEON)\n");
    printf("For 100ms deadline: Use <=512x512 matrices\n");
    
    printf("\n=== Benchmark Complete ===\n");
    return 0;
}
