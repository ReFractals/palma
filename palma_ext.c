/**
 * @file palma_ext.c
 * @brief PALMA Extended Features - Eigenvalues, Scheduling, File I/O
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

#define _POSIX_C_SOURCE 200809L
#define _DEFAULT_SOURCE

#include "palma.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>

#if PALMA_USE_NEON
#include <arm_neon.h>
#endif

/* Binary file magic */
#define PALMA_BINARY_MAGIC 0x504C4D41
#define PALMA_BINARY_VERSION 1

/*============================================================================
 * EIGENVALUE & EIGENVECTOR COMPUTATION
 *============================================================================*/

palma_val_t palma_eigenvalue(const palma_matrix_t *A, palma_semiring_t semiring) {
    if (!A || A->rows != A->cols) {
        palma_set_last_error(A ? PALMA_ERR_NOT_SQUARE : PALMA_ERR_NULL_PTR);
        return PALMA_NEG_INF;
    }
    
    size_t n = A->rows;
    
    /* Karp's algorithm for maximum cycle mean */
    /* D[k][v] = optimal weight of path of exactly k edges ending at v */
    palma_val_t **D = (palma_val_t**)malloc((n + 1) * sizeof(palma_val_t*));
    if (!D) {
        palma_set_last_error(PALMA_ERR_OUT_OF_MEMORY);
        return PALMA_NEG_INF;
    }
    
    palma_val_t zero = palma_zero(semiring);
    palma_val_t one = palma_one(semiring);
    
    for (size_t k = 0; k <= n; k++) {
        D[k] = (palma_val_t*)malloc(n * sizeof(palma_val_t));
        if (!D[k]) {
            for (size_t j = 0; j < k; j++) free(D[j]);
            free(D);
            palma_set_last_error(PALMA_ERR_OUT_OF_MEMORY);
            return PALMA_NEG_INF;
        }
        
        for (size_t v = 0; v < n; v++) {
            D[k][v] = (k == 0) ? one : zero;
        }
    }
    
    /* Dynamic programming: D[k][v] = ⊕_u (D[k-1][u] ⊗ A[u][v]) */
    for (size_t k = 1; k <= n; k++) {
        for (size_t v = 0; v < n; v++) {
            palma_val_t best = zero;
            for (size_t u = 0; u < n; u++) {
                palma_val_t edge = palma_matrix_get(A, u, v);
                if (edge != zero && D[k-1][u] != zero) {
                    palma_val_t path = palma_mul(D[k-1][u], edge, semiring);
                    best = palma_add(best, path, semiring);
                }
            }
            D[k][v] = best;
        }
    }
    
    /* Find maximum cycle mean */
    palma_val_t max_mean = PALMA_NEG_INF;
    
    for (size_t v = 0; v < n; v++) {
        if (D[n][v] == zero) continue;
        
        palma_val_t min_for_v = PALMA_POS_INF;
        for (size_t k = 0; k < n; k++) {
            if (D[k][v] != zero) {
                /* (D[n][v] - D[k][v]) / (n - k) */
                int64_t diff;
                if (semiring == PALMA_MAXPLUS || semiring == PALMA_MINPLUS) {
                    diff = (int64_t)D[n][v] - (int64_t)D[k][v];
                } else {
                    /* For other semirings, eigenvalue computation may not apply */
                    diff = 0;
                }
                palma_val_t mean = (palma_val_t)(diff / (int64_t)(n - k));
                if (mean < min_for_v) min_for_v = mean;
            }
        }
        
        if (min_for_v != PALMA_POS_INF && min_for_v > max_mean) {
            max_mean = min_for_v;
        }
    }
    
    /* Cleanup */
    for (size_t k = 0; k <= n; k++) free(D[k]);
    free(D);
    
    palma_clear_error();
    return max_mean;
}

palma_error_t palma_eigenvector(const palma_matrix_t *A, palma_val_t *eigenvector,
                                 palma_val_t *eigenvalue, palma_semiring_t semiring,
                                 unsigned int max_iter) {
    if (!A || !eigenvector) return PALMA_ERR_NULL_PTR;
    if (A->rows != A->cols) return PALMA_ERR_NOT_SQUARE;
    
    size_t n = A->rows;
    if (max_iter == 0) max_iter = PALMA_DEFAULT_MAX_ITER;
    
    /* Compute eigenvalue first */
    palma_val_t lambda = palma_eigenvalue(A, semiring);
    if (eigenvalue) *eigenvalue = lambda;
    
    if (lambda == PALMA_NEG_INF) {
        /* Acyclic - no eigenvector in the usual sense */
        for (size_t i = 0; i < n; i++) {
            eigenvector[i] = palma_zero(semiring);
        }
        return PALMA_ERR_NOT_CONVERGED;
    }
    
    /* Power iteration with normalization */
    palma_val_t *x = (palma_val_t*)malloc(n * sizeof(palma_val_t));
    palma_val_t *y = (palma_val_t*)malloc(n * sizeof(palma_val_t));
    
    if (!x || !y) {
        free(x);
        free(y);
        return PALMA_ERR_OUT_OF_MEMORY;
    }
    
    /* Initialize with ones */
    palma_val_t one = palma_one(semiring);
    for (size_t i = 0; i < n; i++) {
        x[i] = one;
    }
    
    /* Iterate: x ← (A ⊗ x) / λ (in tropical: subtract λ from each component) */
    for (unsigned int iter = 0; iter < max_iter; iter++) {
        palma_matvec(A, x, y, semiring);
        
        /* Normalize by subtracting λ (for max-plus) */
        if (semiring == PALMA_MAXPLUS || semiring == PALMA_MINPLUS) {
            for (size_t i = 0; i < n; i++) {
                if (y[i] != palma_zero(semiring)) {
                    y[i] -= lambda;
                }
            }
        }
        
        /* Check convergence */
        bool converged = true;
        for (size_t i = 0; i < n; i++) {
            if (y[i] != x[i]) {
                converged = false;
                break;
            }
        }
        
        memcpy(x, y, n * sizeof(palma_val_t));
        
        if (converged) {
            memcpy(eigenvector, x, n * sizeof(palma_val_t));
            free(x);
            free(y);
            return PALMA_SUCCESS;
        }
    }
    
    /* Return last iterate even if not converged */
    memcpy(eigenvector, x, n * sizeof(palma_val_t));
    free(x);
    free(y);
    
    return PALMA_ERR_NOT_CONVERGED;
}

int palma_critical_nodes(const palma_matrix_t *A, int *critical_nodes,
                          palma_semiring_t semiring) {
    if (!A || !critical_nodes) {
        palma_set_last_error(PALMA_ERR_NULL_PTR);
        return -1;
    }
    if (A->rows != A->cols) {
        palma_set_last_error(PALMA_ERR_NOT_SQUARE);
        return -1;
    }
    
    size_t n = A->rows;
    palma_val_t lambda = palma_eigenvalue(A, semiring);
    
    if (lambda == PALMA_NEG_INF) {
        /* Acyclic - no critical cycles */
        memset(critical_nodes, 0, n * sizeof(int));
        return 0;
    }
    
    /* A node is critical if it's on a cycle with mean equal to λ */
    /* Compute A - λI and find strongly connected components with zero cycle mean */
    
    int count = 0;
    memset(critical_nodes, 0, n * sizeof(int));
    
    /* Simple approach: check each node's participation in critical cycles */
    for (size_t i = 0; i < n; i++) {
        /* Check if node i is on a cycle achieving λ */
        for (size_t j = 0; j < n; j++) {
            palma_val_t a_ij = palma_matrix_get(A, i, j);
            palma_val_t a_ji = palma_matrix_get(A, j, i);
            
            if (a_ij != palma_zero(semiring) && a_ji != palma_zero(semiring)) {
                /* 2-cycle mean */
                palma_val_t cycle_weight = palma_mul(a_ij, a_ji, semiring);
                if (semiring == PALMA_MAXPLUS || semiring == PALMA_MINPLUS) {
                    palma_val_t mean = cycle_weight / 2;
                    if (mean >= lambda - PALMA_DEFAULT_TOL) {
                        critical_nodes[i] = 1;
                        critical_nodes[j] = 1;
                    }
                }
            }
        }
        
        /* Self-loop */
        palma_val_t diag = palma_matrix_get(A, i, i);
        if (diag != palma_zero(semiring) && diag >= lambda - PALMA_DEFAULT_TOL) {
            critical_nodes[i] = 1;
        }
    }
    
    for (size_t i = 0; i < n; i++) {
        if (critical_nodes[i]) count++;
    }
    
    palma_clear_error();
    return count;
}

/*============================================================================
 * GRAPH ALGORITHMS
 *============================================================================*/

palma_matrix_t* palma_all_pairs_paths(const palma_matrix_t *adj, palma_semiring_t semiring) {
    return palma_matrix_closure(adj, semiring);
}

palma_error_t palma_single_source_paths(const palma_matrix_t *adj, size_t source,
                                         palma_val_t *dist, palma_semiring_t semiring) {
    if (!adj || !dist) return PALMA_ERR_NULL_PTR;
    if (source >= adj->rows) return PALMA_ERR_INDEX_BOUNDS;
    
    palma_val_t zero = palma_zero(semiring);
    palma_val_t one = palma_one(semiring);
    
    /* Initialize distances */
    for (size_t i = 0; i < adj->rows; i++) {
        dist[i] = zero;
    }
    dist[source] = one;
    
    /* Bellman-Ford via tropical iteration */
    palma_iterate(adj, dist, (unsigned int)adj->rows, semiring);
    
    return PALMA_SUCCESS;
}

palma_matrix_t* palma_reachability(const palma_matrix_t *adj) {
    if (!adj) {
        palma_set_last_error(PALMA_ERR_NULL_PTR);
        return NULL;
    }
    
    /* Convert to Boolean: any non-zero value becomes 1 */
    palma_matrix_t *bool_adj = palma_matrix_create(adj->rows, adj->cols);
    if (!bool_adj) return NULL;
    
    for (size_t i = 0; i < adj->rows; i++) {
        for (size_t j = 0; j < adj->cols; j++) {
            palma_val_t val = palma_matrix_get(adj, i, j);
            /* Consider anything that's not the typical "no edge" as an edge */
            bool is_edge = (val != PALMA_NEG_INF && val != PALMA_POS_INF) || 
                           (i == j);  /* Self-reachable */
            palma_matrix_set(bool_adj, i, j, is_edge ? 1 : 0);
        }
    }
    
    /* Compute transitive closure with Boolean semiring */
    palma_matrix_t *reach = palma_matrix_closure(bool_adj, PALMA_BOOLEAN);
    palma_matrix_destroy(bool_adj);
    
    return reach;
}

palma_matrix_t* palma_bottleneck_paths(const palma_matrix_t *adj) {
    if (!adj) {
        palma_set_last_error(PALMA_ERR_NULL_PTR);
        return NULL;
    }
    
    /* Use max-min semiring for bottleneck/capacity paths */
    return palma_matrix_closure(adj, PALMA_MAXMIN);
}

/*============================================================================
 * SCHEDULING
 *============================================================================*/

palma_scheduler_t* palma_scheduler_create(size_t n_tasks, bool use_maxplus) {
    palma_scheduler_t *sched = (palma_scheduler_t*)malloc(sizeof(palma_scheduler_t));
    if (!sched) {
        palma_set_last_error(PALMA_ERR_OUT_OF_MEMORY);
        return NULL;
    }
    
    sched->n_tasks = n_tasks;
    sched->semiring = use_maxplus ? PALMA_MAXPLUS : PALMA_MINPLUS;
    sched->task_names = NULL;
    
    sched->system = palma_matrix_create_zero(n_tasks, n_tasks, sched->semiring);
    sched->state = (palma_val_t*)malloc(n_tasks * sizeof(palma_val_t));
    sched->input = (palma_val_t*)malloc(n_tasks * sizeof(palma_val_t));
    
    if (!sched->system || !sched->state || !sched->input) {
        palma_scheduler_destroy(sched);
        palma_set_last_error(PALMA_ERR_OUT_OF_MEMORY);
        return NULL;
    }
    
    palma_val_t zero = palma_zero(sched->semiring);
    for (size_t i = 0; i < n_tasks; i++) {
        sched->state[i] = zero;
        sched->input[i] = zero;
    }
    
    palma_clear_error();
    return sched;
}

void palma_scheduler_destroy(palma_scheduler_t *sched) {
    if (!sched) return;
    
    palma_matrix_destroy(sched->system);
    free(sched->state);
    free(sched->input);
    
    if (sched->task_names) {
        for (size_t i = 0; i < sched->n_tasks; i++) {
            free(sched->task_names[i]);
        }
        free(sched->task_names);
    }
    
    free(sched);
}

palma_error_t palma_scheduler_set_name(palma_scheduler_t *sched, size_t task, const char *name) {
    if (!sched || !name) return PALMA_ERR_NULL_PTR;
    if (task >= sched->n_tasks) return PALMA_ERR_INDEX_BOUNDS;
    
    if (!sched->task_names) {
        sched->task_names = (char**)calloc(sched->n_tasks, sizeof(char*));
        if (!sched->task_names) return PALMA_ERR_OUT_OF_MEMORY;
    }
    
    free(sched->task_names[task]);
    sched->task_names[task] = strdup(name);
    if (!sched->task_names[task]) return PALMA_ERR_OUT_OF_MEMORY;
    
    return PALMA_SUCCESS;
}

palma_error_t palma_scheduler_add_constraint(palma_scheduler_t *sched, size_t from,
                                              size_t to, palma_val_t duration) {
    if (!sched) return PALMA_ERR_NULL_PTR;
    if (from >= sched->n_tasks || to >= sched->n_tasks) return PALMA_ERR_INDEX_BOUNDS;
    
    palma_val_t current = palma_matrix_get(sched->system, to, from);
    palma_val_t updated = palma_add(current, duration, sched->semiring);
    palma_matrix_set(sched->system, to, from, updated);
    
    return PALMA_SUCCESS;
}

palma_error_t palma_scheduler_set_ready_time(palma_scheduler_t *sched, size_t task,
                                              palma_val_t ready_time) {
    if (!sched) return PALMA_ERR_NULL_PTR;
    if (task >= sched->n_tasks) return PALMA_ERR_INDEX_BOUNDS;
    
    sched->input[task] = palma_add(sched->input[task], ready_time, sched->semiring);
    sched->state[task] = palma_add(sched->state[task], ready_time, sched->semiring);
    
    return PALMA_SUCCESS;
}

int palma_scheduler_solve(palma_scheduler_t *sched, unsigned int max_iter) {
    if (!sched) {
        palma_set_last_error(PALMA_ERR_NULL_PTR);
        return -1;
    }
    
    if (max_iter == 0) max_iter = (unsigned int)sched->n_tasks;
    
    palma_val_t *prev = (palma_val_t*)malloc(sched->n_tasks * sizeof(palma_val_t));
    palma_val_t *temp = (palma_val_t*)malloc(sched->n_tasks * sizeof(palma_val_t));
    
    if (!prev || !temp) {
        free(prev);
        free(temp);
        palma_set_last_error(PALMA_ERR_OUT_OF_MEMORY);
        return -1;
    }
    
    unsigned int iter;
    for (iter = 0; iter < max_iter; iter++) {
        memcpy(prev, sched->state, sched->n_tasks * sizeof(palma_val_t));
        
        /* x = A ⊗ x ⊕ b */
        palma_matvec(sched->system, prev, temp, sched->semiring);
        
        for (size_t i = 0; i < sched->n_tasks; i++) {
            sched->state[i] = palma_add(temp[i], sched->input[i], sched->semiring);
            /* Also ensure monotonicity */
            sched->state[i] = palma_add(sched->state[i], prev[i], sched->semiring);
        }
        
        /* Check convergence */
        bool converged = true;
        for (size_t i = 0; i < sched->n_tasks; i++) {
            if (sched->state[i] != prev[i]) {
                converged = false;
                break;
            }
        }
        
        if (converged) {
            free(prev);
            free(temp);
            palma_clear_error();
            return (int)iter + 1;
        }
    }
    
    free(prev);
    free(temp);
    palma_clear_error();
    return (int)max_iter;
}

palma_val_t palma_scheduler_get_completion(const palma_scheduler_t *sched, size_t task) {
    if (!sched || task >= sched->n_tasks) {
        return PALMA_NEG_INF;
    }
    return sched->state[task];
}

palma_val_t palma_scheduler_cycle_time(const palma_scheduler_t *sched) {
    if (!sched) {
        palma_set_last_error(PALMA_ERR_NULL_PTR);
        return PALMA_NEG_INF;
    }
    
    return palma_eigenvalue(sched->system, sched->semiring);
}

double palma_scheduler_throughput(const palma_scheduler_t *sched) {
    if (!sched) return 0.0;
    
    palma_val_t cycle = palma_scheduler_cycle_time(sched);
    if (cycle == PALMA_NEG_INF || cycle == 0) return 0.0;
    
    return 1.0 / (double)cycle;
}

int palma_scheduler_critical_path(const palma_scheduler_t *sched, size_t *path, size_t max_len) {
    if (!sched || !path || max_len == 0) {
        palma_set_last_error(PALMA_ERR_NULL_PTR);
        return -1;
    }
    
    /* Find task with maximum completion time */
    palma_val_t max_time = PALMA_NEG_INF;
    size_t end_task = 0;
    
    for (size_t i = 0; i < sched->n_tasks; i++) {
        if (sched->state[i] > max_time) {
            max_time = sched->state[i];
            end_task = i;
        }
    }
    
    /* Backtrack to find critical path */
    size_t len = 0;
    size_t current = end_task;
    
    /* Temporary array to build path in reverse */
    size_t *temp_path = (size_t*)malloc(sched->n_tasks * sizeof(size_t));
    if (!temp_path) {
        palma_set_last_error(PALMA_ERR_OUT_OF_MEMORY);
        return -1;
    }
    
    temp_path[len++] = current;
    
    while (len < sched->n_tasks) {
        palma_val_t current_time = sched->state[current];
        size_t predecessor = current;
        bool found = false;
        
        for (size_t j = 0; j < sched->n_tasks; j++) {
            palma_val_t edge = palma_matrix_get(sched->system, current, j);
            if (edge != palma_zero(sched->semiring)) {
                palma_val_t pred_time = sched->state[j];
                palma_val_t expected = palma_mul(pred_time, edge, sched->semiring);
                
                if (expected == current_time) {
                    predecessor = j;
                    found = true;
                    break;
                }
            }
        }
        
        if (!found || predecessor == current) break;
        
        current = predecessor;
        temp_path[len++] = current;
    }
    
    /* Reverse into output */
    size_t out_len = (len < max_len) ? len : max_len;
    for (size_t i = 0; i < out_len; i++) {
        path[i] = temp_path[len - 1 - i];
    }
    
    free(temp_path);
    palma_clear_error();
    return (int)out_len;
}

/*============================================================================
 * FILE I/O
 *============================================================================*/

palma_error_t palma_matrix_save_csv(const palma_matrix_t *mat, const char *filename,
                                     palma_semiring_t semiring) {
    if (!mat || !filename) return PALMA_ERR_NULL_PTR;
    
    FILE *fp = fopen(filename, "w");
    if (!fp) return PALMA_ERR_FILE_OPEN;
    
    palma_val_t neg_inf_val = PALMA_NEG_INF;
    palma_val_t pos_inf_val = PALMA_POS_INF;
    
    /* Write header comment */
    fprintf(fp, "# PALMA matrix %zux%zu, semiring=%s\n", mat->rows, mat->cols, 
            palma_semiring_name(semiring));
    
    for (size_t i = 0; i < mat->rows; i++) {
        for (size_t j = 0; j < mat->cols; j++) {
            palma_val_t val = palma_matrix_get(mat, i, j);
            
            if (val == neg_inf_val) {
                fprintf(fp, "-inf");
            } else if (val == pos_inf_val) {
                fprintf(fp, "inf");
            } else {
                fprintf(fp, "%d", val);
            }
            
            if (j < mat->cols - 1) fprintf(fp, ",");
        }
        fprintf(fp, "\n");
    }
    
    fclose(fp);
    return PALMA_SUCCESS;
}

palma_matrix_t* palma_matrix_load_csv(const char *filename, palma_semiring_t semiring) {
    (void)semiring;  /* Used only for parsing special values */
    
    if (!filename) {
        palma_set_last_error(PALMA_ERR_NULL_PTR);
        return NULL;
    }
    
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        palma_set_last_error(PALMA_ERR_FILE_OPEN);
        return NULL;
    }
    
    /* First pass: count rows and columns */
    size_t rows = 0, cols = 0;
    char line[65536];
    
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        
        if (rows == 0) {
            /* Count columns in first data row */
            cols = 1;
            for (char *p = line; *p; p++) {
                if (*p == ',') cols++;
            }
        }
        rows++;
    }
    
    if (rows == 0 || cols == 0) {
        fclose(fp);
        palma_set_last_error(PALMA_ERR_FILE_FORMAT);
        return NULL;
    }
    
    /* Create matrix */
    palma_matrix_t *mat = palma_matrix_create(rows, cols);
    if (!mat) {
        fclose(fp);
        return NULL;
    }
    
    /* Second pass: read values */
    rewind(fp);
    size_t row = 0;
    
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        
        char *token = strtok(line, ",\n\r");
        size_t col = 0;
        
        while (token && col < cols) {
            /* Skip whitespace */
            while (*token && isspace(*token)) token++;
            
            palma_val_t val;
            if (strncmp(token, "-inf", 4) == 0 || strncmp(token, "-Inf", 4) == 0) {
                val = PALMA_NEG_INF;
            } else if (strncmp(token, "inf", 3) == 0 || strncmp(token, "Inf", 3) == 0) {
                val = PALMA_POS_INF;
            } else {
                val = (palma_val_t)atoi(token);
            }
            
            palma_matrix_set(mat, row, col, val);
            token = strtok(NULL, ",\n\r");
            col++;
        }
        row++;
    }
    
    fclose(fp);
    palma_clear_error();
    return mat;
}

palma_error_t palma_matrix_save_binary(const palma_matrix_t *mat, const char *filename) {
    if (!mat || !filename) return PALMA_ERR_NULL_PTR;
    
    FILE *fp = fopen(filename, "wb");
    if (!fp) return PALMA_ERR_FILE_OPEN;
    
    /* Write header */
    uint32_t magic = PALMA_BINARY_MAGIC;
    uint32_t version = PALMA_BINARY_VERSION;
    uint32_t rows = (uint32_t)mat->rows;
    uint32_t cols = (uint32_t)mat->cols;
    
    fwrite(&magic, sizeof(magic), 1, fp);
    fwrite(&version, sizeof(version), 1, fp);
    fwrite(&rows, sizeof(rows), 1, fp);
    fwrite(&cols, sizeof(cols), 1, fp);
    
    /* Write data row by row (to handle stride) */
    for (size_t i = 0; i < mat->rows; i++) {
        fwrite(&mat->data[i * mat->stride], sizeof(palma_val_t), mat->cols, fp);
    }
    
    fclose(fp);
    return PALMA_SUCCESS;
}

palma_matrix_t* palma_matrix_load_binary(const char *filename) {
    if (!filename) {
        palma_set_last_error(PALMA_ERR_NULL_PTR);
        return NULL;
    }
    
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        palma_set_last_error(PALMA_ERR_FILE_OPEN);
        return NULL;
    }
    
    /* Read header */
    uint32_t magic, version, rows, cols;
    
    if (fread(&magic, sizeof(magic), 1, fp) != 1 ||
        fread(&version, sizeof(version), 1, fp) != 1 ||
        fread(&rows, sizeof(rows), 1, fp) != 1 ||
        fread(&cols, sizeof(cols), 1, fp) != 1) {
        fclose(fp);
        palma_set_last_error(PALMA_ERR_FILE_READ);
        return NULL;
    }
    
    if (magic != PALMA_BINARY_MAGIC) {
        fclose(fp);
        palma_set_last_error(PALMA_ERR_FILE_FORMAT);
        return NULL;
    }
    
    palma_matrix_t *mat = palma_matrix_create(rows, cols);
    if (!mat) {
        fclose(fp);
        return NULL;
    }
    
    /* Read data */
    for (size_t i = 0; i < rows; i++) {
        if (fread(&mat->data[i * mat->stride], sizeof(palma_val_t), cols, fp) != cols) {
            palma_matrix_destroy(mat);
            fclose(fp);
            palma_set_last_error(PALMA_ERR_FILE_READ);
            return NULL;
        }
    }
    
    fclose(fp);
    palma_clear_error();
    return mat;
}

palma_error_t palma_sparse_save_csv(const palma_sparse_t *sp, const char *filename) {
    if (!sp || !filename) return PALMA_ERR_NULL_PTR;
    
    FILE *fp = fopen(filename, "w");
    if (!fp) return PALMA_ERR_FILE_OPEN;
    
    /* Write header */
    fprintf(fp, "# PALMA sparse matrix %zux%zu, nnz=%zu, semiring=%s\n",
            sp->rows, sp->cols, sp->nnz, palma_semiring_name(sp->semiring));
    fprintf(fp, "# Format: row,col,value (COO format)\n");
    fprintf(fp, "%zu,%zu,%zu\n", sp->rows, sp->cols, sp->nnz);
    
    /* Write entries */
    for (size_t i = 0; i < sp->rows; i++) {
        for (palma_idx_t k = sp->row_ptr[i]; k < sp->row_ptr[i + 1]; k++) {
            fprintf(fp, "%zu,%u,%d\n", i, sp->col_idx[k], sp->values[k]);
        }
    }
    
    fclose(fp);
    return PALMA_SUCCESS;
}

palma_sparse_t* palma_sparse_load_csv(const char *filename, palma_semiring_t semiring) {
    if (!filename) {
        palma_set_last_error(PALMA_ERR_NULL_PTR);
        return NULL;
    }
    
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        palma_set_last_error(PALMA_ERR_FILE_OPEN);
        return NULL;
    }
    
    char line[1024];
    size_t rows = 0, cols = 0, nnz = 0;
    
    /* Read past comments and get dimensions */
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#') continue;
        if (sscanf(line, "%zu,%zu,%zu", &rows, &cols, &nnz) == 3) break;
    }
    
    if (rows == 0 || cols == 0) {
        fclose(fp);
        palma_set_last_error(PALMA_ERR_FILE_FORMAT);
        return NULL;
    }
    
    palma_sparse_t *sp = palma_sparse_create(rows, cols, nnz, semiring);
    if (!sp) {
        fclose(fp);
        return NULL;
    }
    
    /* Read entries */
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#') continue;
        
        size_t row;
        unsigned int col;
        int val;
        
        if (sscanf(line, "%zu,%u,%d", &row, &col, &val) == 3) {
            palma_sparse_set(sp, row, col, val);
        }
    }
    
    fclose(fp);
    palma_clear_error();
    return sp;
}

palma_error_t palma_matrix_export_dot(const palma_matrix_t *mat, const char *filename,
                                       palma_semiring_t semiring, const char **node_names) {
    if (!mat || !filename) return PALMA_ERR_NULL_PTR;
    
    FILE *fp = fopen(filename, "w");
    if (!fp) return PALMA_ERR_FILE_OPEN;
    
    palma_val_t zero = palma_zero(semiring);
    
    fprintf(fp, "digraph PALMA {\n");
    fprintf(fp, "  // Generated by PALMA - Parallel Algebra Library for Max-plus Applications\n");
    fprintf(fp, "  // Author: Gnankan Landry Regis N'guessan\n");
    fprintf(fp, "  rankdir=LR;\n");
    fprintf(fp, "  node [shape=circle];\n\n");
    
    /* Define nodes */
    for (size_t i = 0; i < mat->rows; i++) {
        if (node_names && node_names[i]) {
            fprintf(fp, "  %zu [label=\"%s\"];\n", i, node_names[i]);
        } else {
            fprintf(fp, "  %zu;\n", i);
        }
    }
    fprintf(fp, "\n");
    
    /* Define edges */
    for (size_t i = 0; i < mat->rows; i++) {
        for (size_t j = 0; j < mat->cols; j++) {
            palma_val_t val = palma_matrix_get(mat, i, j);
            if (val != zero && i != j) {
                if (val == PALMA_NEG_INF) {
                    fprintf(fp, "  %zu -> %zu [label=\"-∞\"];\n", j, i);
                } else if (val == PALMA_POS_INF) {
                    fprintf(fp, "  %zu -> %zu [label=\"∞\"];\n", j, i);
                } else {
                    fprintf(fp, "  %zu -> %zu [label=\"%d\"];\n", j, i, val);
                }
            }
        }
    }
    
    fprintf(fp, "}\n");
    fclose(fp);
    
    return PALMA_SUCCESS;
}

/*============================================================================
 * NEON OPTIMIZATIONS
 *============================================================================*/

#if PALMA_USE_NEON

palma_error_t palma_matrix_mul_neon(palma_matrix_t *C, const palma_matrix_t *A,
                                     const palma_matrix_t *B, palma_semiring_t semiring) {
    if (!C || !A || !B) return PALMA_ERR_NULL_PTR;
    if (A->cols != B->rows) return PALMA_ERR_INVALID_DIM;
    if (C->rows != A->rows || C->cols != B->cols) return PALMA_ERR_INVALID_DIM;
    
    /* Only support max-plus and min-plus with NEON */
    if (semiring != PALMA_MAXPLUS && semiring != PALMA_MINPLUS) {
        return palma_matrix_mul_into(C, A, B, semiring);
    }
    
    palma_val_t zero = palma_zero(semiring);
    int32x4_t zero_vec = vdupq_n_s32(zero);
    
    for (size_t i = 0; i < A->rows; i++) {
        for (size_t j = 0; j < B->cols; j++) {
            int32x4_t acc = zero_vec;
            
            size_t k = 0;
            for (; k + 4 <= A->cols; k += 4) {
                int32x4_t a_vec = vld1q_s32(&A->data[i * A->stride + k]);
                
                int32_t b_vals[4] = {
                    palma_matrix_get(B, k, j),
                    palma_matrix_get(B, k+1, j),
                    palma_matrix_get(B, k+2, j),
                    palma_matrix_get(B, k+3, j)
                };
                int32x4_t b_vec = vld1q_s32(b_vals);
                
                int32x4_t prod = vaddq_s32(a_vec, b_vec);
                
                if (semiring == PALMA_MAXPLUS) {
                    acc = vmaxq_s32(acc, prod);
                } else {
                    acc = vminq_s32(acc, prod);
                }
            }
            
            int32_t result_arr[4];
            vst1q_s32(result_arr, acc);
            palma_val_t result = result_arr[0];
            for (int p = 1; p < 4; p++) {
                result = palma_add(result, result_arr[p], semiring);
            }
            
            for (; k < A->cols; k++) {
                palma_val_t a_ik = palma_matrix_get(A, i, k);
                palma_val_t b_kj = palma_matrix_get(B, k, j);
                palma_val_t prod = palma_mul(a_ik, b_kj, semiring);
                result = palma_add(result, prod, semiring);
            }
            
            palma_matrix_set(C, i, j, result);
        }
    }
    
    return PALMA_SUCCESS;
}

palma_error_t palma_matvec_neon(const palma_matrix_t *A, const palma_val_t *x,
                                 palma_val_t *y, palma_semiring_t semiring) {
    if (!A || !x || !y) return PALMA_ERR_NULL_PTR;
    
    if (semiring != PALMA_MAXPLUS && semiring != PALMA_MINPLUS) {
        palma_val_t zero = palma_zero(semiring);
        for (size_t i = 0; i < A->rows; i++) {
            palma_val_t sum = zero;
            for (size_t j = 0; j < A->cols; j++) {
                palma_val_t prod = palma_mul(palma_matrix_get(A, i, j), x[j], semiring);
                sum = palma_add(sum, prod, semiring);
            }
            y[i] = sum;
        }
        return PALMA_SUCCESS;
    }
    
    palma_val_t zero = palma_zero(semiring);
    
    for (size_t i = 0; i < A->rows; i++) {
        int32x4_t acc = vdupq_n_s32(zero);
        
        size_t j = 0;
        for (; j + 4 <= A->cols; j += 4) {
            int32x4_t a_vec = vld1q_s32(&A->data[i * A->stride + j]);
            int32x4_t x_vec = vld1q_s32(&x[j]);
            
            int32x4_t prod = vaddq_s32(a_vec, x_vec);
            
            if (semiring == PALMA_MAXPLUS) {
                acc = vmaxq_s32(acc, prod);
            } else {
                acc = vminq_s32(acc, prod);
            }
        }
        
        int32_t result_arr[4];
        vst1q_s32(result_arr, acc);
        palma_val_t result = result_arr[0];
        for (int p = 1; p < 4; p++) {
            result = palma_add(result, result_arr[p], semiring);
        }
        
        for (; j < A->cols; j++) {
            palma_val_t prod = palma_mul(palma_matrix_get(A, i, j), x[j], semiring);
            result = palma_add(result, prod, semiring);
        }
        
        y[i] = result;
    }
    
    return PALMA_SUCCESS;
}

#endif /* PALMA_USE_NEON */

/*============================================================================
 * UTILITY FUNCTIONS
 *============================================================================*/

void palma_matrix_print(const palma_matrix_t *mat, const char *name,
                         palma_semiring_t semiring, FILE *fp) {
    if (!fp) fp = stdout;
    
    if (!mat) {
        fprintf(fp, "%s: NULL\n", name ? name : "Matrix");
        return;
    }
    
    palma_val_t zero = palma_zero(semiring);
    const char *zero_str;
    
    switch (semiring) {
        case PALMA_MAXPLUS:
        case PALMA_MAXMIN:
            zero_str = "-∞";
            break;
        case PALMA_MINPLUS:
        case PALMA_MINMAX:
            zero_str = "+∞";
            break;
        default:
            zero_str = "0";
    }
    
    fprintf(fp, "%s (%zu × %zu):\n", name ? name : "Matrix", mat->rows, mat->cols);
    
    for (size_t i = 0; i < mat->rows; i++) {
        fprintf(fp, "  [");
        for (size_t j = 0; j < mat->cols; j++) {
            palma_val_t val = palma_matrix_get(mat, i, j);
            if (val == zero) {
                fprintf(fp, "%6s", zero_str);
            } else if (val == PALMA_POS_INF && zero != PALMA_POS_INF) {
                fprintf(fp, "%6s", "+∞");
            } else if (val == PALMA_NEG_INF && zero != PALMA_NEG_INF) {
                fprintf(fp, "%6s", "-∞");
            } else {
                fprintf(fp, "%6d", val);
            }
            if (j < mat->cols - 1) fprintf(fp, ", ");
        }
        fprintf(fp, "]\n");
    }
}

void palma_sparse_print(const palma_sparse_t *sp, const char *name, FILE *fp) {
    if (!fp) fp = stdout;
    
    if (!sp) {
        fprintf(fp, "%s: NULL\n", name ? name : "Sparse Matrix");
        return;
    }
    
    fprintf(fp, "%s (%zu × %zu, nnz=%zu, sparsity=%.1f%%):\n",
            name ? name : "Sparse Matrix", sp->rows, sp->cols, sp->nnz,
            palma_sparse_sparsity(sp) * 100.0);
    
    for (size_t i = 0; i < sp->rows; i++) {
        if (sp->row_ptr[i] < sp->row_ptr[i + 1]) {
            fprintf(fp, "  Row %zu:", i);
            for (palma_idx_t k = sp->row_ptr[i]; k < sp->row_ptr[i + 1]; k++) {
                fprintf(fp, " [%u]=%d", sp->col_idx[k], sp->values[k]);
            }
            fprintf(fp, "\n");
        }
    }
}

void palma_vector_print(const palma_val_t *vec, size_t len, const char *name,
                         palma_semiring_t semiring, FILE *fp) {
    if (!fp) fp = stdout;
    
    if (!vec) {
        fprintf(fp, "%s: NULL\n", name ? name : "Vector");
        return;
    }
    
    palma_val_t zero = palma_zero(semiring);
    const char *zero_str = (semiring == PALMA_MAXPLUS || semiring == PALMA_MAXMIN) ? "-∞" : 
                           (semiring == PALMA_MINPLUS || semiring == PALMA_MINMAX) ? "+∞" : "0";
    
    fprintf(fp, "%s (%zu): [", name ? name : "Vector", len);
    
    for (size_t i = 0; i < len; i++) {
        if (vec[i] == zero) {
            fprintf(fp, "%s", zero_str);
        } else if (vec[i] == PALMA_POS_INF && zero != PALMA_POS_INF) {
            fprintf(fp, "+∞");
        } else if (vec[i] == PALMA_NEG_INF && zero != PALMA_NEG_INF) {
            fprintf(fp, "-∞");
        } else {
            fprintf(fp, "%d", vec[i]);
        }
        if (i < len - 1) fprintf(fp, ", ");
    }
    fprintf(fp, "]\n");
}

const char* palma_version(void) {
    return PALMA_VERSION_STRING;
}

void palma_version_components(int *major, int *minor, int *patch) {
    if (major) *major = PALMA_VERSION_MAJOR;
    if (minor) *minor = PALMA_VERSION_MINOR;
    if (patch) *patch = PALMA_VERSION_PATCH;
}

bool palma_has_neon(void) {
#if PALMA_USE_NEON
    return true;
#else
    return false;
#endif
}

bool palma_has_openmp(void) {
#if PALMA_USE_OPENMP
    return true;
#else
    return false;
#endif
}

const char* palma_build_config(void) {
    static char config[256];
    snprintf(config, sizeof(config),
             "PALMA v%s [NEON:%s, OpenMP:%s]",
             PALMA_VERSION_STRING,
             palma_has_neon() ? "ON" : "OFF",
             palma_has_openmp() ? "ON" : "OFF");
    return config;
}
