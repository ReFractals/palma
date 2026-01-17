/**
 * @file example_graphs.c
 * @brief Graph Algorithms using Multiple Tropical Semirings
 * 
 * Demonstrates how different semirings solve different graph problems:
 *   - Min-plus: Shortest paths
 *   - Max-plus: Longest paths (critical paths)
 *   - Max-min: Bottleneck paths (maximum bandwidth)
 *   - Boolean: Reachability analysis
 * 
 * @author Gnankan Landry Regis N'guessan
 *         Axiom Research Group
 *         NM-AIST / AIMS-RIC
 * @email  rnguessan@aimsric.org
 */

#include <stdio.h>
#include <stdlib.h>
#include "palma.h"

/* Network nodes */
enum {
    NODE_SERVER = 0,
    NODE_ROUTER_A,
    NODE_ROUTER_B,
    NODE_CLIENT_1,
    NODE_CLIENT_2,
    NODE_CLIENT_3,
    NUM_NODES
};

static const char *node_names[] = {
    "Server", "Router_A", "Router_B", "Client_1", "Client_2", "Client_3"
};

static void print_distance_table(const palma_matrix_t *dist, const char *title, 
                                  palma_semiring_t semiring) {
    printf("\n%s:\n", title);
    printf("%-10s", "From\\To");
    for (int j = 0; j < NUM_NODES; j++) {
        printf(" %9s", node_names[j]);
    }
    printf("\n");
    
    for (int i = 0; i < 10 + NUM_NODES * 10; i++) printf("-");
    printf("\n");
    
    palma_val_t zero = palma_zero(semiring);
    
    for (int i = 0; i < NUM_NODES; i++) {
        printf("%-10s", node_names[i]);
        for (int j = 0; j < NUM_NODES; j++) {
            palma_val_t d = palma_matrix_get(dist, i, j);
            if (d == zero) {
                printf(" %9s", semiring == PALMA_MINPLUS ? "∞" : "-∞");
            } else if (d == PALMA_POS_INF) {
                printf(" %9s", "∞");
            } else if (d == PALMA_NEG_INF) {
                printf(" %9s", "-∞");
            } else {
                printf(" %9d", d);
            }
        }
        printf("\n");
    }
}

int main(void) {
    printf("╔══════════════════════════════════════════════════════════════════╗\n");
    printf("║  PALMA - Parallel Algebra Library for Max-plus Applications      ║\n");
    printf("║  Graph Algorithms with Multiple Semirings                        ║\n");
    printf("╚══════════════════════════════════════════════════════════════════╝\n\n");
    
    printf("Library: %s\n", palma_build_config());
    printf("Author: Gnankan Landry Regis N'guessan\n\n");
    
    /*
     * Network topology:
     *                    ┌─────────────┐
     *          5ms       │  Router A   │      3ms
     *    ┌───────────────┤  (bw: 100)  ├───────────────┐
     *    │               └──────┬──────┘               │
     *    │                      │ 2ms                  │
     *    ▼                      ▼                      ▼
     * ┌──────┐             ┌─────────────┐         ┌──────────┐
     * │Server│◄────────────┤  Router B   ├────────►│ Client 1 │
     * └──┬───┘    8ms      │  (bw: 50)   │  4ms    └──────────┘
     *    │                 └──────┬──────┘
     *    │                        │ 6ms
     *    │    10ms           ┌────┴────┐
     *    └───────────────────┤Client 2 │
     *                        └─────────┘
     *                             │ 7ms
     *                        ┌────┴────┐
     *                        │Client 3 │
     *                        └─────────┘
     */
    
    printf("=== Network Topology ===\n");
    printf("Modeling a network with latencies and bandwidths.\n\n");
    
    /* Create adjacency matrix with latencies */
    palma_matrix_t *latency = palma_matrix_create_zero(NUM_NODES, NUM_NODES, PALMA_MINPLUS);
    palma_matrix_t *bandwidth = palma_matrix_create_zero(NUM_NODES, NUM_NODES, PALMA_MAXMIN);
    
    #define ADD_EDGE(from, to, lat, bw) do { \
        palma_matrix_set(latency, from, to, lat); \
        palma_matrix_set(latency, to, from, lat); \
        palma_matrix_set(bandwidth, from, to, bw); \
        palma_matrix_set(bandwidth, to, from, bw); \
    } while(0)
    
    /* Define network edges */
    ADD_EDGE(NODE_SERVER, NODE_ROUTER_A, 5, 100);
    ADD_EDGE(NODE_SERVER, NODE_ROUTER_B, 8, 50);
    ADD_EDGE(NODE_SERVER, NODE_CLIENT_2, 10, 30);
    ADD_EDGE(NODE_ROUTER_A, NODE_ROUTER_B, 2, 80);
    ADD_EDGE(NODE_ROUTER_A, NODE_CLIENT_1, 3, 100);
    ADD_EDGE(NODE_ROUTER_B, NODE_CLIENT_1, 4, 60);
    ADD_EDGE(NODE_ROUTER_B, NODE_CLIENT_2, 6, 40);
    ADD_EDGE(NODE_CLIENT_2, NODE_CLIENT_3, 7, 20);
    
    /* Self-loops with identity */
    for (int i = 0; i < NUM_NODES; i++) {
        palma_matrix_set(latency, i, i, 0);
        palma_matrix_set(bandwidth, i, i, PALMA_POS_INF);
    }
    
    printf("Latency matrix (edge weights in ms):\n");
    palma_matrix_print(latency, "L", PALMA_MINPLUS, stdout);
    
    /* ========== SHORTEST PATHS (Min-Plus) ========== */
    printf("\n=== 1. SHORTEST PATHS (Min-Plus Semiring) ===\n");
    printf("Semiring: (min, +) with zero = +∞, one = 0\n");
    printf("Finds minimum total latency paths.\n");
    
    palma_matrix_t *shortest = palma_matrix_closure(latency, PALMA_MINPLUS);
    print_distance_table(shortest, "Shortest Path Latencies (ms)", PALMA_MINPLUS);
    
    /* Single-source example */
    palma_val_t dist[NUM_NODES];
    palma_single_source_paths(latency, NODE_SERVER, dist, PALMA_MINPLUS);
    printf("\nFrom Server to all nodes:\n");
    for (int i = 0; i < NUM_NODES; i++) {
        printf("  → %s: %dms\n", node_names[i], dist[i]);
    }
    
    /* ========== BOTTLENECK PATHS (Max-Min) ========== */
    printf("\n=== 2. BOTTLENECK PATHS (Max-Min Semiring) ===\n");
    printf("Semiring: (max, min) with zero = -∞, one = +∞\n");
    printf("Finds maximum bandwidth paths (limited by smallest edge).\n");
    
    palma_matrix_t *bottleneck = palma_bottleneck_paths(bandwidth);
    print_distance_table(bottleneck, "Maximum Bandwidth Paths (Mbps)", PALMA_MAXMIN);
    
    printf("\nInterpretation: Server→Client_3 max bandwidth is %d Mbps\n",
           palma_matrix_get(bottleneck, NODE_SERVER, NODE_CLIENT_3));
    printf("  (Limited by the Client_2→Client_3 link at 20 Mbps)\n");
    
    /* ========== REACHABILITY (Boolean) ========== */
    printf("\n=== 3. REACHABILITY (Boolean Semiring) ===\n");
    printf("Semiring: (OR, AND) with zero = 0, one = 1\n");
    printf("Determines which nodes can reach which.\n");
    
    palma_matrix_t *reach = palma_reachability(latency);
    
    printf("\nReachability matrix (1 = path exists):\n");
    printf("%-10s", "From\\To");
    for (int j = 0; j < NUM_NODES; j++) {
        printf(" %5s", node_names[j]);
    }
    printf("\n");
    for (int i = 0; i < 10 + NUM_NODES * 6; i++) printf("-");
    printf("\n");
    
    for (int i = 0; i < NUM_NODES; i++) {
        printf("%-10s", node_names[i]);
        for (int j = 0; j < NUM_NODES; j++) {
            palma_val_t r = palma_matrix_get(reach, i, j);
            printf(" %5s", r ? "✓" : "✗");
        }
        printf("\n");
    }
    
    printf("\nAll nodes are mutually reachable (fully connected network).\n");
    
    /* ========== MATRIX POWERS ========== */
    printf("\n=== 4. MATRIX POWERS: k-Hop Paths ===\n");
    printf("L^k gives shortest paths using exactly k hops.\n\n");
    
    palma_matrix_t *L2 = palma_matrix_power(latency, 2, PALMA_MINPLUS);
    printf("L² (shortest 2-hop paths):\n");
    palma_matrix_print(L2, "L²", PALMA_MINPLUS, stdout);
    
    printf("\nExample: Server→Client_1 in exactly 2 hops:\n");
    printf("  Server → Router_A → Client_1: 5 + 3 = 8ms\n");
    printf("  Server → Router_B → Client_1: 8 + 4 = 12ms\n");
    printf("  L²[Server][Client_1] = min(8, 12) = %d ms ✓\n",
           palma_matrix_get(L2, NODE_SERVER, NODE_CLIENT_1));
    
    /* ========== SPARSE MATRIX DEMO ========== */
    printf("\n=== 5. SPARSE MATRIX OPERATIONS ===\n");
    
    palma_sparse_t *sp_latency = palma_sparse_from_dense(latency, PALMA_MINPLUS);
    
    printf("Converted to sparse format:\n");
    printf("  Original size: %zu × %zu = %zu elements\n", 
           latency->rows, latency->cols, latency->rows * latency->cols);
    printf("  Sparse: %zu non-zeros (%.1f%% sparsity)\n",
           sp_latency->nnz, palma_sparse_sparsity(sp_latency) * 100);
    
    palma_sparse_print(sp_latency, "Sparse L", stdout);
    
    /* Sparse matrix-vector multiplication */
    palma_val_t x[NUM_NODES] = {0, PALMA_POS_INF, PALMA_POS_INF, 
                                 PALMA_POS_INF, PALMA_POS_INF, PALMA_POS_INF};
    palma_val_t y[NUM_NODES];
    
    palma_sparse_matvec(sp_latency, x, y);
    
    printf("\nSparse matrix-vector: L ⊗ [0, ∞, ∞, ∞, ∞, ∞]ᵀ\n");
    printf("Result (1-hop distances from Server): ");
    palma_vector_print(y, NUM_NODES, NULL, PALMA_MINPLUS, stdout);
    
    /* ========== FILE I/O ========== */
    printf("\n=== 6. FILE I/O ===\n");
    
    palma_matrix_save_csv(latency, "network_latency.csv", PALMA_MINPLUS);
    printf("Saved latency matrix to 'network_latency.csv'\n");
    
    palma_sparse_save_csv(sp_latency, "network_sparse.csv");
    printf("Saved sparse matrix to 'network_sparse.csv'\n");
    
    palma_matrix_save_binary(shortest, "shortest_paths.bin");
    printf("Saved shortest paths to 'shortest_paths.bin' (binary)\n");
    
    palma_matrix_export_dot(latency, "network.dot", PALMA_MINPLUS, node_names);
    printf("Exported GraphViz DOT to 'network.dot'\n");
    printf("  Visualize: dot -Tpng network.dot -o network.png\n");
    
    /* Cleanup */
    palma_matrix_destroy(latency);
    palma_matrix_destroy(bandwidth);
    palma_matrix_destroy(shortest);
    palma_matrix_destroy(bottleneck);
    palma_matrix_destroy(reach);
    palma_matrix_destroy(L2);
    palma_sparse_destroy(sp_latency);
    
    printf("\n=== Example Complete ===\n");
    return 0;
}
