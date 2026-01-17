/**
 * @file example_eigenvalue.c
 * @brief Tropical Eigenvalue and Eigenvector Computation
 * 
 * Demonstrates computation of tropical eigenvalues (maximum cycle mean)
 * and eigenvectors for analyzing periodic systems and steady-state behavior.
 * 
 * In tropical algebra:
 *   - The eigenvalue λ is the maximum average weight over all cycles
 *   - A ⊗ v = λ ⊗ v means A ⊗ v = v + λ (component-wise)
 *   - The eigenvalue determines the asymptotic growth rate of A^k
 * 
 * @author Gnankan Landry Regis N'guessan
 *         Axiom Research Group
 *         NM-AIST / AIMS-RIC
 * @email  rnguessan@aimsric.org
 */

#include <stdio.h>
#include <stdlib.h>
#include "palma.h"

int main(void) {
    printf("╔══════════════════════════════════════════════════════════════════╗\n");
    printf("║  PALMA - Parallel Algebra Library for Max-plus Applications      ║\n");
    printf("║  Tropical Eigenvalue and Eigenvector Computation                 ║\n");
    printf("╚══════════════════════════════════════════════════════════════════╝\n\n");
    
    printf("Library: %s\n", palma_build_config());
    printf("Author: Gnankan Landry Regis N'guessan\n\n");
    
    /* ========== EXAMPLE 1: Simple Cycle ========== */
    printf("=== Example 1: Simple 3-Node Cycle ===\n\n");
    
    /*
     * A simple cycle: 0 → 1 → 2 → 0
     *                 5   3   4
     * 
     * Cycle weight = 5 + 3 + 4 = 12
     * Cycle length = 3
     * Maximum cycle mean = 12/3 = 4
     */
    
    palma_matrix_t *A1 = palma_matrix_create_zero(3, 3, PALMA_MAXPLUS);
    
    palma_matrix_set(A1, 1, 0, 5);  /* 0 → 1, weight 5 */
    palma_matrix_set(A1, 2, 1, 3);  /* 1 → 2, weight 3 */
    palma_matrix_set(A1, 0, 2, 4);  /* 2 → 0, weight 4 */
    
    printf("Matrix A (simple cycle 0→1→2→0):\n");
    palma_matrix_print(A1, "A", PALMA_MAXPLUS, stdout);
    
    palma_val_t lambda1 = palma_eigenvalue(A1, PALMA_MAXPLUS);
    printf("\nTropical eigenvalue λ = %d\n", lambda1);
    printf("Expected: (5+3+4)/3 = 4 ✓\n\n");
    
    printf("Interpretation:\n");
    printf("  - For large k, (A^k)[i][j] ≈ k·λ + constant\n");
    printf("  - The system 'grows' by λ = %d per iteration\n", lambda1);
    printf("  - In scheduling: cycle time = %d time units\n\n", lambda1);
    
    /* ========== EXAMPLE 2: Multiple Cycles ========== */
    printf("=== Example 2: Multiple Cycles ===\n\n");
    
    /*
     * Two cycles sharing a node:
     *   Cycle 1: 0 → 1 → 0 (weights 3, 5) → mean = 8/2 = 4
     *   Cycle 2: 0 → 2 → 0 (weights 2, 4) → mean = 6/2 = 3
     * 
     * Maximum cycle mean = max(4, 3) = 4
     */
    
    palma_matrix_t *A2 = palma_matrix_create_zero(3, 3, PALMA_MAXPLUS);
    
    palma_matrix_set(A2, 1, 0, 3);  /* 0 → 1 */
    palma_matrix_set(A2, 0, 1, 5);  /* 1 → 0 */
    palma_matrix_set(A2, 2, 0, 2);  /* 0 → 2 */
    palma_matrix_set(A2, 0, 2, 4);  /* 2 → 0 */
    
    printf("Matrix A (two cycles through node 0):\n");
    palma_matrix_print(A2, "A", PALMA_MAXPLUS, stdout);
    
    palma_val_t lambda2 = palma_eigenvalue(A2, PALMA_MAXPLUS);
    printf("\nTropical eigenvalue λ = %d\n", lambda2);
    printf("Cycle 1 (0↔1): mean = (3+5)/2 = 4\n");
    printf("Cycle 2 (0↔2): mean = (2+4)/2 = 3\n");
    printf("Maximum = %d ✓\n\n", lambda2);
    
    /* ========== EXAMPLE 3: Eigenvector Computation ========== */
    printf("=== Example 3: Eigenvector Computation ===\n\n");
    
    /*
     * For A ⊗ v = λ ⊗ v:
     * The eigenvector v describes the relative timing offsets
     * in the periodic steady state.
     */
    
    palma_val_t eigenvec[3];
    palma_val_t eigenval;
    
    palma_error_t err = palma_eigenvector(A1, eigenvec, &eigenval, PALMA_MAXPLUS, 100);
    
    if (err == PALMA_SUCCESS) {
        printf("Eigenvector computation converged!\n\n");
    } else if (err == PALMA_ERR_NOT_CONVERGED) {
        printf("Eigenvector computation did not fully converge (using last iterate)\n\n");
    }
    
    printf("Eigenvalue λ = %d\n", eigenval);
    printf("Eigenvector v = ");
    palma_vector_print(eigenvec, 3, NULL, PALMA_MAXPLUS, stdout);
    
    printf("\nVerification: A ⊗ v should equal λ ⊗ v = v + %d\n", eigenval);
    
    palma_val_t Av[3];
    palma_matvec(A1, eigenvec, Av, PALMA_MAXPLUS);
    
    printf("A ⊗ v = ");
    palma_vector_print(Av, 3, NULL, PALMA_MAXPLUS, stdout);
    
    printf("v + λ = [");
    for (int i = 0; i < 3; i++) {
        if (eigenvec[i] != PALMA_NEG_INF) {
            printf("%d", eigenvec[i] + eigenval);
        } else {
            printf("-∞");
        }
        if (i < 2) printf(", ");
    }
    printf("]\n\n");
    
    /* ========== EXAMPLE 4: Critical Nodes ========== */
    printf("=== Example 4: Critical Nodes ===\n\n");
    
    printf("Critical nodes are those participating in cycles with maximum mean.\n\n");
    
    int critical[3];
    int n_critical = palma_critical_nodes(A2, critical, PALMA_MAXPLUS);
    
    printf("Matrix A2 has %d critical node(s):\n", n_critical);
    for (int i = 0; i < 3; i++) {
        if (critical[i]) {
            printf("  Node %d: on critical cycle\n", i);
        }
    }
    
    printf("\nNodes 0 and 1 form the critical cycle with mean = 4.\n");
    printf("Node 2 is not critical (its cycle has mean = 3 < 4).\n\n");
    
    /* ========== EXAMPLE 5: Production System ========== */
    printf("=== Example 5: Production System Analysis ===\n\n");
    
    /*
     * A manufacturing system with 4 machines:
     *   M0: Raw material input (every 5 time units)
     *   M1: Processing stage 1 (takes 3 time units after M0)
     *   M2: Processing stage 2 (takes 4 time units after M1)
     *   M3: Output/packaging (takes 2 time units after M2, feeds back to M0)
     * 
     * The cycle time determines production throughput.
     */
    
    printf("Manufacturing system:\n");
    printf("  M0 (input) → M1 (process) → M2 (process) → M3 (output) → M0\n\n");
    
    palma_matrix_t *prod = palma_matrix_create_zero(4, 4, PALMA_MAXPLUS);
    
    /* Processing times as edge weights */
    palma_matrix_set(prod, 1, 0, 5);   /* M0 → M1: 5 units */
    palma_matrix_set(prod, 2, 1, 3);   /* M1 → M2: 3 units */
    palma_matrix_set(prod, 3, 2, 4);   /* M2 → M3: 4 units */
    palma_matrix_set(prod, 0, 3, 2);   /* M3 → M0: 2 units (feedback) */
    
    palma_matrix_print(prod, "Production System", PALMA_MAXPLUS, stdout);
    
    palma_val_t cycle_time = palma_eigenvalue(prod, PALMA_MAXPLUS);
    
    printf("\nCycle time (tropical eigenvalue): %d time units\n", cycle_time);
    printf("Total processing: 5 + 3 + 4 + 2 = 14 units\n");
    printf("Cycle mean: 14/4 = 3.5, rounded to %d\n\n", cycle_time);
    
    double throughput = 1.0 / (double)cycle_time;
    printf("Production throughput: %.3f items per time unit\n", throughput);
    printf("Or: 1 item every %d time units\n\n", cycle_time);
    
    printf("To increase throughput, reduce weights on the critical cycle.\n");
    printf("All machines are on the single cycle, so any improvement helps.\n");
    
    /* Cleanup */
    palma_matrix_destroy(A1);
    palma_matrix_destroy(A2);
    palma_matrix_destroy(prod);
    
    printf("\n=== Example Complete ===\n");
    return 0;
}
