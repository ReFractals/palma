/**
 * @file palma.c
 * @brief PALMA - Parallel Algebra Library for Max-plus Applications
 * 
 * Implementation of tropical algebra operations for embedded systems.
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

#define _POSIX_C_SOURCE 200112L

#include "palma.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>

#if PALMA_USE_NEON
#include <arm_neon.h>
#endif

#if PALMA_USE_OPENMP
#include <omp.h>
#endif

/*============================================================================
 * INTERNAL CONSTANTS
 *============================================================================*/

#define ALIGN_SIZE 16
#define ALIGN_STRIDE(cols) (((cols) + 3) & ~3)

/* Binary file magic number */
#define PALMA_BINARY_MAGIC 0x504C4D41  /* "PLMA" */
#define PALMA_BINARY_VERSION 1

/*============================================================================
 * ERROR HANDLING
 *============================================================================*/

/* Thread-local error storage (fallback to global if not supported) */
#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L && !defined(__STDC_NO_THREADS__)
    #include <threads.h>
    static _Thread_local palma_error_t g_last_error = PALMA_SUCCESS;
#elif defined(__GNUC__) || defined(__clang__)
    static __thread palma_error_t g_last_error = PALMA_SUCCESS;
#else
    static palma_error_t g_last_error = PALMA_SUCCESS;
#endif

static const char* error_messages[] = {
    "Success",
    "NULL pointer argument",
    "Invalid matrix dimensions",
    "Out of memory",
    "Invalid argument value",
    "Matrix must be square",
    "Algorithm did not converge",
    "Failed to open file",
    "Failed to read from file",
    "Failed to write to file",
    "Invalid file format",
    "Index out of bounds",
    "Invalid sparse matrix format",
    "Unsupported operation"
};

const char* palma_strerror(palma_error_t err) {
    int idx = (err <= 0) ? -err : 0;
    if (idx >= (int)(sizeof(error_messages) / sizeof(error_messages[0]))) {
        return "Unknown error";
    }
    return error_messages[idx];
}

palma_error_t palma_get_last_error(void) {
    return g_last_error;
}

void palma_set_last_error(palma_error_t err) {
    g_last_error = err;
}

void palma_clear_error(void) {
    g_last_error = PALMA_SUCCESS;
}

/* Internal macro for setting error and returning */
#define PALMA_RETURN_ERROR(err) do { \
    palma_set_last_error(err); \
    return (err); \
} while(0)

#define PALMA_RETURN_NULL(err) do { \
    palma_set_last_error(err); \
    return NULL; \
} while(0)

/*============================================================================
 * MEMORY HELPERS
 *============================================================================*/

static void* aligned_alloc_wrapper(size_t alignment, size_t size) {
#if defined(_POSIX_C_SOURCE) && _POSIX_C_SOURCE >= 200112L
    void *ptr = NULL;
    if (posix_memalign(&ptr, alignment, size) != 0) return NULL;
    return ptr;
#else
    (void)alignment;
    return malloc(size);
#endif
}

/*============================================================================
 * SEMIRING OPERATIONS
 *============================================================================*/

palma_val_t palma_zero(palma_semiring_t semiring) {
    switch (semiring) {
        case PALMA_MAXPLUS:  return PALMA_NEG_INF;
        case PALMA_MINPLUS:  return PALMA_POS_INF;
        case PALMA_MAXMIN:   return PALMA_NEG_INF;
        case PALMA_MINMAX:   return PALMA_POS_INF;
        case PALMA_BOOLEAN:  return 0;
        default:             return PALMA_NEG_INF;
    }
}

palma_val_t palma_one(palma_semiring_t semiring) {
    switch (semiring) {
        case PALMA_MAXPLUS:  return 0;
        case PALMA_MINPLUS:  return 0;
        case PALMA_MAXMIN:   return PALMA_POS_INF;
        case PALMA_MINMAX:   return PALMA_NEG_INF;
        case PALMA_BOOLEAN:  return 1;
        default:             return 0;
    }
}

palma_val_t palma_add(palma_val_t a, palma_val_t b, palma_semiring_t semiring) {
    switch (semiring) {
        case PALMA_MAXPLUS:
        case PALMA_MAXMIN:
            return (a > b) ? a : b;
        case PALMA_MINPLUS:
        case PALMA_MINMAX:
            return (a < b) ? a : b;
        case PALMA_BOOLEAN:
            return (a || b) ? 1 : 0;
        default:
            return (a > b) ? a : b;
    }
}

palma_val_t palma_mul(palma_val_t a, palma_val_t b, palma_semiring_t semiring) {
    switch (semiring) {
        case PALMA_MAXPLUS:
        case PALMA_MINPLUS: {
            /* Standard addition with infinity handling */
            if (a == PALMA_NEG_INF || b == PALMA_NEG_INF) return PALMA_NEG_INF;
            if (a == PALMA_POS_INF || b == PALMA_POS_INF) return PALMA_POS_INF;
            int64_t sum = (int64_t)a + (int64_t)b;
            if (sum > INT32_MAX) return PALMA_POS_INF;
            if (sum < INT32_MIN) return PALMA_NEG_INF;
            return (palma_val_t)sum;
        }
        case PALMA_MAXMIN:
            return (a < b) ? a : b;  /* min */
        case PALMA_MINMAX:
            return (a > b) ? a : b;  /* max */
        case PALMA_BOOLEAN:
            return (a && b) ? 1 : 0;
        default:
            return a + b;
    }
}

bool palma_is_zero(palma_val_t a, palma_semiring_t semiring) {
    return a == palma_zero(semiring);
}

const char* palma_semiring_name(palma_semiring_t semiring) {
    switch (semiring) {
        case PALMA_MAXPLUS:  return "max-plus";
        case PALMA_MINPLUS:  return "min-plus";
        case PALMA_MAXMIN:   return "max-min (bottleneck)";
        case PALMA_MINMAX:   return "min-max";
        case PALMA_BOOLEAN:  return "Boolean";
        default:             return "unknown";
    }
}

/*============================================================================
 * DENSE MATRIX LIFECYCLE
 *============================================================================*/

palma_matrix_t* palma_matrix_create(size_t rows, size_t cols) {
    if (rows == 0 || cols == 0) {
        PALMA_RETURN_NULL(PALMA_ERR_INVALID_DIM);
    }
    
    palma_matrix_t *mat = (palma_matrix_t*)malloc(sizeof(palma_matrix_t));
    if (!mat) {
        PALMA_RETURN_NULL(PALMA_ERR_OUT_OF_MEMORY);
    }
    
    mat->rows = rows;
    mat->cols = cols;
    mat->stride = ALIGN_STRIDE(cols);
    mat->owns_data = true;
    
    size_t data_size = mat->rows * mat->stride * sizeof(palma_val_t);
    mat->data = (palma_val_t*)aligned_alloc_wrapper(ALIGN_SIZE, data_size);
    
    if (!mat->data) {
        free(mat);
        PALMA_RETURN_NULL(PALMA_ERR_OUT_OF_MEMORY);
    }
    
    palma_clear_error();
    return mat;
}

palma_matrix_t* palma_matrix_create_zero(size_t rows, size_t cols, palma_semiring_t semiring) {
    palma_matrix_t *mat = palma_matrix_create(rows, cols);
    if (!mat) return NULL;
    
    palma_val_t zero = palma_zero(semiring);
    size_t total = mat->rows * mat->stride;
    
    for (size_t i = 0; i < total; i++) {
        mat->data[i] = zero;
    }
    
    return mat;
}

palma_matrix_t* palma_matrix_create_identity(size_t n, palma_semiring_t semiring) {
    palma_matrix_t *mat = palma_matrix_create_zero(n, n, semiring);
    if (!mat) return NULL;
    
    palma_val_t one = palma_one(semiring);
    for (size_t i = 0; i < n; i++) {
        palma_matrix_set(mat, i, i, one);
    }
    
    return mat;
}

palma_matrix_t* palma_matrix_wrap(palma_val_t *data, size_t rows, size_t cols, size_t stride) {
    if (!data || rows == 0 || cols == 0) {
        PALMA_RETURN_NULL(PALMA_ERR_INVALID_ARG);
    }
    
    palma_matrix_t *mat = (palma_matrix_t*)malloc(sizeof(palma_matrix_t));
    if (!mat) {
        PALMA_RETURN_NULL(PALMA_ERR_OUT_OF_MEMORY);
    }
    
    mat->data = data;
    mat->rows = rows;
    mat->cols = cols;
    mat->stride = stride;
    mat->owns_data = false;
    
    palma_clear_error();
    return mat;
}

palma_matrix_t* palma_matrix_clone(const palma_matrix_t *src) {
    if (!src) {
        PALMA_RETURN_NULL(PALMA_ERR_NULL_PTR);
    }
    
    palma_matrix_t *dst = palma_matrix_create(src->rows, src->cols);
    if (!dst) return NULL;
    
    for (size_t i = 0; i < src->rows; i++) {
        memcpy(palma_matrix_row(dst, i),
               &src->data[i * src->stride],
               src->cols * sizeof(palma_val_t));
    }
    
    return dst;
}

void palma_matrix_destroy(palma_matrix_t *mat) {
    if (!mat) return;
    if (mat->owns_data && mat->data) {
        free(mat->data);
    }
    free(mat);
}

/*============================================================================
 * DENSE MATRIX ACCESS
 *============================================================================*/

palma_error_t palma_matrix_get_safe(const palma_matrix_t *mat, size_t row, size_t col, palma_val_t *out) {
    if (!mat || !out) PALMA_RETURN_ERROR(PALMA_ERR_NULL_PTR);
    if (row >= mat->rows || col >= mat->cols) PALMA_RETURN_ERROR(PALMA_ERR_INDEX_BOUNDS);
    
    *out = palma_matrix_get(mat, row, col);
    return PALMA_SUCCESS;
}

palma_error_t palma_matrix_set_safe(palma_matrix_t *mat, size_t row, size_t col, palma_val_t val) {
    if (!mat) PALMA_RETURN_ERROR(PALMA_ERR_NULL_PTR);
    if (row >= mat->rows || col >= mat->cols) PALMA_RETURN_ERROR(PALMA_ERR_INDEX_BOUNDS);
    
    palma_matrix_set(mat, row, col, val);
    return PALMA_SUCCESS;
}

/*============================================================================
 * SPARSE MATRIX LIFECYCLE
 *============================================================================*/

palma_sparse_t* palma_sparse_create(size_t rows, size_t cols, size_t capacity, palma_semiring_t semiring) {
    if (rows == 0 || cols == 0) {
        PALMA_RETURN_NULL(PALMA_ERR_INVALID_DIM);
    }
    
    palma_sparse_t *sp = (palma_sparse_t*)malloc(sizeof(palma_sparse_t));
    if (!sp) {
        PALMA_RETURN_NULL(PALMA_ERR_OUT_OF_MEMORY);
    }
    
    sp->rows = rows;
    sp->cols = cols;
    sp->nnz = 0;
    sp->capacity = (capacity > 0) ? capacity : 16;
    sp->semiring = semiring;
    
    sp->values = (palma_val_t*)malloc(sp->capacity * sizeof(palma_val_t));
    sp->col_idx = (palma_idx_t*)malloc(sp->capacity * sizeof(palma_idx_t));
    sp->row_ptr = (palma_idx_t*)calloc(rows + 1, sizeof(palma_idx_t));
    
    if (!sp->values || !sp->col_idx || !sp->row_ptr) {
        free(sp->values);
        free(sp->col_idx);
        free(sp->row_ptr);
        free(sp);
        PALMA_RETURN_NULL(PALMA_ERR_OUT_OF_MEMORY);
    }
    
    palma_clear_error();
    return sp;
}

palma_sparse_t* palma_sparse_from_dense(const palma_matrix_t *dense, palma_semiring_t semiring) {
    if (!dense) {
        PALMA_RETURN_NULL(PALMA_ERR_NULL_PTR);
    }
    
    /* Count non-zeros */
    palma_val_t zero = palma_zero(semiring);
    size_t nnz = 0;
    for (size_t i = 0; i < dense->rows; i++) {
        for (size_t j = 0; j < dense->cols; j++) {
            if (palma_matrix_get(dense, i, j) != zero) {
                nnz++;
            }
        }
    }
    
    palma_sparse_t *sp = palma_sparse_create(dense->rows, dense->cols, nnz, semiring);
    if (!sp) return NULL;
    
    /* Fill CSR structure */
    size_t idx = 0;
    for (size_t i = 0; i < dense->rows; i++) {
        sp->row_ptr[i] = (palma_idx_t)idx;
        for (size_t j = 0; j < dense->cols; j++) {
            palma_val_t val = palma_matrix_get(dense, i, j);
            if (val != zero) {
                sp->values[idx] = val;
                sp->col_idx[idx] = (palma_idx_t)j;
                idx++;
            }
        }
    }
    sp->row_ptr[dense->rows] = (palma_idx_t)idx;
    sp->nnz = idx;
    
    return sp;
}

palma_matrix_t* palma_sparse_to_dense(const palma_sparse_t *sparse) {
    if (!sparse) {
        PALMA_RETURN_NULL(PALMA_ERR_NULL_PTR);
    }
    
    palma_matrix_t *dense = palma_matrix_create_zero(sparse->rows, sparse->cols, sparse->semiring);
    if (!dense) return NULL;
    
    for (size_t i = 0; i < sparse->rows; i++) {
        for (palma_idx_t k = sparse->row_ptr[i]; k < sparse->row_ptr[i + 1]; k++) {
            palma_matrix_set(dense, i, sparse->col_idx[k], sparse->values[k]);
        }
    }
    
    return dense;
}

palma_sparse_t* palma_sparse_clone(const palma_sparse_t *src) {
    if (!src) {
        PALMA_RETURN_NULL(PALMA_ERR_NULL_PTR);
    }
    
    palma_sparse_t *dst = palma_sparse_create(src->rows, src->cols, src->nnz, src->semiring);
    if (!dst) return NULL;
    
    memcpy(dst->values, src->values, src->nnz * sizeof(palma_val_t));
    memcpy(dst->col_idx, src->col_idx, src->nnz * sizeof(palma_idx_t));
    memcpy(dst->row_ptr, src->row_ptr, (src->rows + 1) * sizeof(palma_idx_t));
    dst->nnz = src->nnz;
    
    return dst;
}

void palma_sparse_destroy(palma_sparse_t *sp) {
    if (!sp) return;
    free(sp->values);
    free(sp->col_idx);
    free(sp->row_ptr);
    free(sp);
}

/*============================================================================
 * SPARSE MATRIX ACCESS
 *============================================================================*/

palma_val_t palma_sparse_get(const palma_sparse_t *sp, size_t row, size_t col) {
    if (!sp || row >= sp->rows || col >= sp->cols) {
        return palma_zero(sp ? sp->semiring : PALMA_MAXPLUS);
    }
    
    /* Binary search in the row */
    palma_idx_t start = sp->row_ptr[row];
    palma_idx_t end = sp->row_ptr[row + 1];
    
    while (start < end) {
        palma_idx_t mid = start + (end - start) / 2;
        if (sp->col_idx[mid] == col) {
            return sp->values[mid];
        } else if (sp->col_idx[mid] < col) {
            start = mid + 1;
        } else {
            end = mid;
        }
    }
    
    return palma_zero(sp->semiring);
}

static palma_error_t sparse_ensure_capacity(palma_sparse_t *sp, size_t needed) {
    if (needed <= sp->capacity) return PALMA_SUCCESS;
    
    size_t new_cap = sp->capacity * 2;
    while (new_cap < needed) new_cap *= 2;
    
    palma_val_t *new_values = (palma_val_t*)realloc(sp->values, new_cap * sizeof(palma_val_t));
    palma_idx_t *new_col_idx = (palma_idx_t*)realloc(sp->col_idx, new_cap * sizeof(palma_idx_t));
    
    if (!new_values || !new_col_idx) {
        PALMA_RETURN_ERROR(PALMA_ERR_OUT_OF_MEMORY);
    }
    
    sp->values = new_values;
    sp->col_idx = new_col_idx;
    sp->capacity = new_cap;
    
    return PALMA_SUCCESS;
}

palma_error_t palma_sparse_set(palma_sparse_t *sp, size_t row, size_t col, palma_val_t val) {
    if (!sp) PALMA_RETURN_ERROR(PALMA_ERR_NULL_PTR);
    if (row >= sp->rows || col >= sp->cols) PALMA_RETURN_ERROR(PALMA_ERR_INDEX_BOUNDS);
    
    /* Find position in row */
    palma_idx_t start = sp->row_ptr[row];
    palma_idx_t end = sp->row_ptr[row + 1];
    palma_idx_t pos = start;
    
    while (pos < end && sp->col_idx[pos] < col) {
        pos++;
    }
    
    if (pos < end && sp->col_idx[pos] == col) {
        /* Update existing entry */
        sp->values[pos] = val;
        return PALMA_SUCCESS;
    }
    
    /* Insert new entry */
    palma_error_t err = sparse_ensure_capacity(sp, sp->nnz + 1);
    if (err != PALMA_SUCCESS) return err;
    
    /* Shift elements to make room */
    memmove(&sp->values[pos + 1], &sp->values[pos], (sp->nnz - pos) * sizeof(palma_val_t));
    memmove(&sp->col_idx[pos + 1], &sp->col_idx[pos], (sp->nnz - pos) * sizeof(palma_idx_t));
    
    sp->values[pos] = val;
    sp->col_idx[pos] = (palma_idx_t)col;
    
    /* Update row pointers */
    for (size_t r = row + 1; r <= sp->rows; r++) {
        sp->row_ptr[r]++;
    }
    sp->nnz++;
    
    return PALMA_SUCCESS;
}

palma_error_t palma_sparse_add_entry(palma_sparse_t *sp, size_t row, size_t col, palma_val_t val) {
    /* This is a simplified version that appends to the end - requires compress after */
    if (!sp) PALMA_RETURN_ERROR(PALMA_ERR_NULL_PTR);
    if (row >= sp->rows || col >= sp->cols) PALMA_RETURN_ERROR(PALMA_ERR_INDEX_BOUNDS);
    
    palma_error_t err = sparse_ensure_capacity(sp, sp->nnz + 1);
    if (err != PALMA_SUCCESS) return err;
    
    /* For now, use the slower palma_sparse_set */
    return palma_sparse_set(sp, row, col, val);
}

palma_error_t palma_sparse_compress(palma_sparse_t *sp) {
    if (!sp) PALMA_RETURN_ERROR(PALMA_ERR_NULL_PTR);
    
    palma_val_t zero = palma_zero(sp->semiring);
    size_t write_idx = 0;
    
    for (size_t row = 0; row < sp->rows; row++) {
        palma_idx_t new_row_start = (palma_idx_t)write_idx;
        
        for (palma_idx_t k = sp->row_ptr[row]; k < sp->row_ptr[row + 1]; k++) {
            if (sp->values[k] != zero) {
                sp->values[write_idx] = sp->values[k];
                sp->col_idx[write_idx] = sp->col_idx[k];
                write_idx++;
            }
        }
        
        sp->row_ptr[row] = new_row_start;
    }
    sp->row_ptr[sp->rows] = (palma_idx_t)write_idx;
    sp->nnz = write_idx;
    
    return PALMA_SUCCESS;
}

size_t palma_sparse_row_nnz(const palma_sparse_t *sp, size_t row) {
    if (!sp || row >= sp->rows) return 0;
    return sp->row_ptr[row + 1] - sp->row_ptr[row];
}

double palma_sparse_sparsity(const palma_sparse_t *sp) {
    if (!sp || sp->rows == 0 || sp->cols == 0) return 1.0;
    size_t total = sp->rows * sp->cols;
    return 1.0 - (double)sp->nnz / (double)total;
}

/*============================================================================
 * DENSE MATRIX OPERATIONS
 *============================================================================*/

palma_error_t palma_matrix_mul_into(palma_matrix_t *C, const palma_matrix_t *A,
                                     const palma_matrix_t *B, palma_semiring_t semiring) {
    if (!C || !A || !B) PALMA_RETURN_ERROR(PALMA_ERR_NULL_PTR);
    if (A->cols != B->rows) PALMA_RETURN_ERROR(PALMA_ERR_INVALID_DIM);
    if (C->rows != A->rows || C->cols != B->cols) PALMA_RETURN_ERROR(PALMA_ERR_INVALID_DIM);
    
    palma_val_t zero = palma_zero(semiring);
    
#if PALMA_USE_OPENMP
    #pragma omp parallel for collapse(2) if(A->rows * B->cols > 1000)
#endif
    for (size_t i = 0; i < A->rows; i++) {
        for (size_t j = 0; j < B->cols; j++) {
            palma_val_t sum = zero;
            
            for (size_t k = 0; k < A->cols; k++) {
                palma_val_t a_ik = palma_matrix_get(A, i, k);
                palma_val_t b_kj = palma_matrix_get(B, k, j);
                palma_val_t prod = palma_mul(a_ik, b_kj, semiring);
                sum = palma_add(sum, prod, semiring);
            }
            
            palma_matrix_set(C, i, j, sum);
        }
    }
    
    return PALMA_SUCCESS;
}

palma_matrix_t* palma_matrix_mul(const palma_matrix_t *A, const palma_matrix_t *B,
                                  palma_semiring_t semiring) {
    if (!A || !B) {
        PALMA_RETURN_NULL(PALMA_ERR_NULL_PTR);
    }
    if (A->cols != B->rows) {
        PALMA_RETURN_NULL(PALMA_ERR_INVALID_DIM);
    }
    
    palma_matrix_t *C = palma_matrix_create(A->rows, B->cols);
    if (!C) return NULL;
    
#if PALMA_USE_NEON
    if (palma_matrix_mul_neon(C, A, B, semiring) == PALMA_SUCCESS) {
        return C;
    }
#endif
    
    palma_error_t err = palma_matrix_mul_into(C, A, B, semiring);
    if (err != PALMA_SUCCESS) {
        palma_matrix_destroy(C);
        return NULL;
    }
    
    return C;
}

palma_matrix_t* palma_matrix_add(const palma_matrix_t *A, const palma_matrix_t *B,
                                  palma_semiring_t semiring) {
    if (!A || !B) {
        PALMA_RETURN_NULL(PALMA_ERR_NULL_PTR);
    }
    if (A->rows != B->rows || A->cols != B->cols) {
        PALMA_RETURN_NULL(PALMA_ERR_INVALID_DIM);
    }
    
    palma_matrix_t *C = palma_matrix_create(A->rows, A->cols);
    if (!C) return NULL;
    
    for (size_t i = 0; i < A->rows; i++) {
        for (size_t j = 0; j < A->cols; j++) {
            palma_val_t a = palma_matrix_get(A, i, j);
            palma_val_t b = palma_matrix_get(B, i, j);
            palma_matrix_set(C, i, j, palma_add(a, b, semiring));
        }
    }
    
    return C;
}

palma_matrix_t* palma_matrix_power(const palma_matrix_t *A, unsigned int n,
                                    palma_semiring_t semiring) {
    if (!A) {
        PALMA_RETURN_NULL(PALMA_ERR_NULL_PTR);
    }
    if (A->rows != A->cols) {
        PALMA_RETURN_NULL(PALMA_ERR_NOT_SQUARE);
    }
    
    if (n == 0) {
        return palma_matrix_create_identity(A->rows, semiring);
    }
    
    /* Binary exponentiation */
    palma_matrix_t *result = palma_matrix_create_identity(A->rows, semiring);
    palma_matrix_t *base = palma_matrix_clone(A);
    palma_matrix_t *temp = NULL;
    
    if (!result || !base) {
        palma_matrix_destroy(result);
        palma_matrix_destroy(base);
        PALMA_RETURN_NULL(PALMA_ERR_OUT_OF_MEMORY);
    }
    
    while (n > 0) {
        if (n & 1) {
            temp = palma_matrix_mul(result, base, semiring);
            palma_matrix_destroy(result);
            result = temp;
            if (!result) {
                palma_matrix_destroy(base);
                return NULL;
            }
        }
        
        n >>= 1;
        if (n > 0) {
            temp = palma_matrix_mul(base, base, semiring);
            palma_matrix_destroy(base);
            base = temp;
            if (!base) {
                palma_matrix_destroy(result);
                return NULL;
            }
        }
    }
    
    palma_matrix_destroy(base);
    return result;
}

palma_matrix_t* palma_matrix_closure(const palma_matrix_t *A, palma_semiring_t semiring) {
    if (!A) {
        PALMA_RETURN_NULL(PALMA_ERR_NULL_PTR);
    }
    if (A->rows != A->cols) {
        PALMA_RETURN_NULL(PALMA_ERR_NOT_SQUARE);
    }
    
    size_t n = A->rows;
    
    /* Floyd-Warshall style computation */
    palma_matrix_t *D = palma_matrix_clone(A);
    if (!D) return NULL;
    
    /* Add identity */
    palma_val_t one = palma_one(semiring);
    for (size_t i = 0; i < n; i++) {
        palma_val_t diag = palma_matrix_get(D, i, i);
        palma_matrix_set(D, i, i, palma_add(diag, one, semiring));
    }
    
    /* Floyd-Warshall iterations */
    for (size_t k = 0; k < n; k++) {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                palma_val_t d_ij = palma_matrix_get(D, i, j);
                palma_val_t d_ik = palma_matrix_get(D, i, k);
                palma_val_t d_kj = palma_matrix_get(D, k, j);
                palma_val_t via_k = palma_mul(d_ik, d_kj, semiring);
                palma_matrix_set(D, i, j, palma_add(d_ij, via_k, semiring));
            }
        }
    }
    
    return D;
}

palma_matrix_t* palma_matrix_transitive_closure(const palma_matrix_t *A, palma_semiring_t semiring) {
    if (!A) {
        PALMA_RETURN_NULL(PALMA_ERR_NULL_PTR);
    }
    if (A->rows != A->cols) {
        PALMA_RETURN_NULL(PALMA_ERR_NOT_SQUARE);
    }
    
    /* A+ = A* ⊗ A = A ⊗ A* */
    palma_matrix_t *star = palma_matrix_closure(A, semiring);
    if (!star) return NULL;
    
    palma_matrix_t *plus = palma_matrix_mul(star, A, semiring);
    palma_matrix_destroy(star);
    
    return plus;
}

/*============================================================================
 * SPARSE MATRIX OPERATIONS
 *============================================================================*/

palma_sparse_t* palma_sparse_mul(const palma_sparse_t *A, const palma_sparse_t *B) {
    if (!A || !B) {
        PALMA_RETURN_NULL(PALMA_ERR_NULL_PTR);
    }
    if (A->cols != B->rows || A->semiring != B->semiring) {
        PALMA_RETURN_NULL(PALMA_ERR_INVALID_DIM);
    }
    
    palma_semiring_t semiring = A->semiring;
    palma_val_t zero = palma_zero(semiring);
    
    /* Estimate capacity */
    size_t est_nnz = (A->nnz + B->nnz) * 2;
    if (est_nnz > A->rows * B->cols) est_nnz = A->rows * B->cols;
    
    palma_sparse_t *C = palma_sparse_create(A->rows, B->cols, est_nnz, semiring);
    if (!C) return NULL;
    
    /* Temporary array for accumulating row results */
    palma_val_t *row_vals = (palma_val_t*)malloc(B->cols * sizeof(palma_val_t));
    if (!row_vals) {
        palma_sparse_destroy(C);
        PALMA_RETURN_NULL(PALMA_ERR_OUT_OF_MEMORY);
    }
    
    for (size_t i = 0; i < A->rows; i++) {
        C->row_ptr[i] = (palma_idx_t)C->nnz;
        
        /* Initialize row accumulator to zero */
        for (size_t j = 0; j < B->cols; j++) {
            row_vals[j] = zero;
        }
        
        /* Compute row i of C */
        for (palma_idx_t ka = A->row_ptr[i]; ka < A->row_ptr[i + 1]; ka++) {
            palma_idx_t k = A->col_idx[ka];
            palma_val_t a_ik = A->values[ka];
            
            for (palma_idx_t kb = B->row_ptr[k]; kb < B->row_ptr[k + 1]; kb++) {
                palma_idx_t j = B->col_idx[kb];
                palma_val_t b_kj = B->values[kb];
                palma_val_t prod = palma_mul(a_ik, b_kj, semiring);
                row_vals[j] = palma_add(row_vals[j], prod, semiring);
            }
        }
        
        /* Store non-zeros */
        for (size_t j = 0; j < B->cols; j++) {
            if (row_vals[j] != zero) {
                sparse_ensure_capacity(C, C->nnz + 1);
                C->values[C->nnz] = row_vals[j];
                C->col_idx[C->nnz] = (palma_idx_t)j;
                C->nnz++;
            }
        }
    }
    C->row_ptr[A->rows] = (palma_idx_t)C->nnz;
    
    free(row_vals);
    return C;
}

palma_error_t palma_sparse_matvec(const palma_sparse_t *A, const palma_val_t *x,
                                   palma_val_t *y) {
    if (!A || !x || !y) PALMA_RETURN_ERROR(PALMA_ERR_NULL_PTR);
    
    palma_semiring_t semiring = A->semiring;
    palma_val_t zero = palma_zero(semiring);
    
    for (size_t i = 0; i < A->rows; i++) {
        palma_val_t sum = zero;
        
        for (palma_idx_t k = A->row_ptr[i]; k < A->row_ptr[i + 1]; k++) {
            palma_idx_t j = A->col_idx[k];
            palma_val_t prod = palma_mul(A->values[k], x[j], semiring);
            sum = palma_add(sum, prod, semiring);
        }
        
        y[i] = sum;
    }
    
    return PALMA_SUCCESS;
}

palma_sparse_t* palma_sparse_closure(const palma_sparse_t *A) {
    if (!A) {
        PALMA_RETURN_NULL(PALMA_ERR_NULL_PTR);
    }
    if (A->rows != A->cols) {
        PALMA_RETURN_NULL(PALMA_ERR_NOT_SQUARE);
    }
    
    /* Convert to dense, compute closure, convert back */
    /* For sparse implementation, this could be optimized with sparse Floyd-Warshall */
    palma_matrix_t *dense = palma_sparse_to_dense(A);
    if (!dense) return NULL;
    
    palma_matrix_t *closure = palma_matrix_closure(dense, A->semiring);
    palma_matrix_destroy(dense);
    if (!closure) return NULL;
    
    palma_sparse_t *result = palma_sparse_from_dense(closure, A->semiring);
    palma_matrix_destroy(closure);
    
    return result;
}

/*============================================================================
 * VECTOR OPERATIONS
 *============================================================================*/

palma_error_t palma_matvec(const palma_matrix_t *A, const palma_val_t *x,
                            palma_val_t *y, palma_semiring_t semiring) {
    if (!A || !x || !y) PALMA_RETURN_ERROR(PALMA_ERR_NULL_PTR);
    
#if PALMA_USE_NEON
    return palma_matvec_neon(A, x, y, semiring);
#endif
    
    palma_val_t zero = palma_zero(semiring);
    
    for (size_t i = 0; i < A->rows; i++) {
        palma_val_t sum = zero;
        for (size_t j = 0; j < A->cols; j++) {
            palma_val_t a_ij = palma_matrix_get(A, i, j);
            palma_val_t prod = palma_mul(a_ij, x[j], semiring);
            sum = palma_add(sum, prod, semiring);
        }
        y[i] = sum;
    }
    
    return PALMA_SUCCESS;
}

palma_error_t palma_iterate(const palma_matrix_t *A, palma_val_t *x,
                             unsigned int n, palma_semiring_t semiring) {
    if (!A || !x) PALMA_RETURN_ERROR(PALMA_ERR_NULL_PTR);
    if (A->rows != A->cols) PALMA_RETURN_ERROR(PALMA_ERR_NOT_SQUARE);
    
    palma_val_t *y = (palma_val_t*)malloc(A->rows * sizeof(palma_val_t));
    if (!y) PALMA_RETURN_ERROR(PALMA_ERR_OUT_OF_MEMORY);
    
    for (unsigned int iter = 0; iter < n; iter++) {
        palma_matvec(A, x, y, semiring);
        memcpy(x, y, A->rows * sizeof(palma_val_t));
    }
    
    free(y);
    return PALMA_SUCCESS;
}

palma_val_t palma_dot(const palma_val_t *x, const palma_val_t *y,
                       size_t len, palma_semiring_t semiring) {
    palma_val_t result = palma_zero(semiring);
    
    for (size_t i = 0; i < len; i++) {
        palma_val_t prod = palma_mul(x[i], y[i], semiring);
        result = palma_add(result, prod, semiring);
    }
    
    return result;
}
