/**
 * @file example_scheduling.c
 * @brief Task Scheduling Example using PALMA Tropical Algebra
 * 
 * Demonstrates how to model and solve precedence-constrained task scheduling
 * problems using max-plus algebra. This is one of the most powerful
 * applications of tropical algebra for embedded real-time systems.
 * 
 * @author Gnankan Landry Regis N'guessan
 *         Axiom Research Group
 *         NM-AIST / AIMS-RIC
 * @email  rnguessan@aimsric.org
 * 
 * @example
 * Compile: make example_scheduling
 * Run: ./build/bin/example_scheduling
 */

#include <stdio.h>
#include <stdlib.h>
#include "palma.h"

/* Task indices for an embedded boot sequence */
enum {
    TASK_HW_INIT = 0,
    TASK_KERNEL,
    TASK_DRIVERS,
    TASK_NETWORK,
    TASK_FILESYSTEM,
    TASK_SERVICES,
    NUM_TASKS
};

/* Task durations in milliseconds */
static const palma_val_t task_durations[] = {
    10,  /* HW_INIT */
    20,  /* KERNEL */
    15,  /* DRIVERS */
    25,  /* NETWORK */
    30,  /* FILESYSTEM */
    10   /* SERVICES */
};

int main(void) {
    printf("╔══════════════════════════════════════════════════════════════════╗\n");
    printf("║  PALMA - Parallel Algebra Library for Max-plus Applications      ║\n");
    printf("║  Task Scheduling Example                                         ║\n");
    printf("╚══════════════════════════════════════════════════════════════════╝\n\n");
    
    printf("Library: %s\n", palma_build_config());
    printf("Author: Gnankan Landry Regis N'guessan\n\n");
    
    /* Create scheduler */
    printf("=== Embedded System Boot Sequence ===\n\n");
    printf("Modeling task dependencies with max-plus algebra:\n");
    printf("  x(k+1) = A ⊗ x(k) ⊕ b\n\n");
    
    palma_scheduler_t *sched = palma_scheduler_create(NUM_TASKS, true);
    if (!sched) {
        fprintf(stderr, "Failed to create scheduler: %s\n", 
                palma_strerror(palma_get_last_error()));
        return 1;
    }
    
    /* Set task names */
    palma_scheduler_set_name(sched, TASK_HW_INIT, "Hardware Init");
    palma_scheduler_set_name(sched, TASK_KERNEL, "Load Kernel");
    palma_scheduler_set_name(sched, TASK_DRIVERS, "Init Drivers");
    palma_scheduler_set_name(sched, TASK_NETWORK, "Start Network");
    palma_scheduler_set_name(sched, TASK_FILESYSTEM, "Mount Filesystem");
    palma_scheduler_set_name(sched, TASK_SERVICES, "Start Services");
    
    /* Define precedence constraints */
    /* Task 0 (HW_INIT) has no predecessors */
    palma_scheduler_set_ready_time(sched, TASK_HW_INIT, 0);
    
    /* Task 1 depends on Task 0 */
    palma_scheduler_add_constraint(sched, TASK_HW_INIT, TASK_KERNEL, 
                                   task_durations[TASK_HW_INIT]);
    
    /* Tasks 2, 3, 4 depend on Task 1 (parallel execution) */
    palma_scheduler_add_constraint(sched, TASK_KERNEL, TASK_DRIVERS, 
                                   task_durations[TASK_KERNEL]);
    palma_scheduler_add_constraint(sched, TASK_KERNEL, TASK_NETWORK, 
                                   task_durations[TASK_KERNEL]);
    palma_scheduler_add_constraint(sched, TASK_KERNEL, TASK_FILESYSTEM, 
                                   task_durations[TASK_KERNEL]);
    
    /* Task 5 depends on Tasks 2, 3, 4 (sync point) */
    palma_scheduler_add_constraint(sched, TASK_DRIVERS, TASK_SERVICES, 
                                   task_durations[TASK_DRIVERS]);
    palma_scheduler_add_constraint(sched, TASK_NETWORK, TASK_SERVICES, 
                                   task_durations[TASK_NETWORK]);
    palma_scheduler_add_constraint(sched, TASK_FILESYSTEM, TASK_SERVICES, 
                                   task_durations[TASK_FILESYSTEM]);
    
    /* Print system matrix */
    printf("System matrix A (precedence with durations):\n");
    palma_matrix_print(sched->system, "A", PALMA_MAXPLUS, stdout);
    printf("\n");
    
    /* Solve */
    printf("Solving schedule...\n");
    int iters = palma_scheduler_solve(sched, 0);
    
    if (iters < 0) {
        fprintf(stderr, "Scheduler failed: %s\n", 
                palma_strerror(palma_get_last_error()));
        palma_scheduler_destroy(sched);
        return 1;
    }
    
    printf("Converged in %d iterations\n\n", iters);
    
    /* Print results */
    printf("Task Completion Schedule:\n");
    printf("┌────────────────────┬─────────┬───────────┐\n");
    printf("│ Task               │ Start   │ Complete  │\n");
    printf("├────────────────────┼─────────┼───────────┤\n");
    
    for (int i = 0; i < NUM_TASKS; i++) {
        palma_val_t completion = palma_scheduler_get_completion(sched, i);
        palma_val_t start = completion;
        
        /* Find actual start time by looking at predecessors */
        palma_val_t max_pred = 0;
        for (int j = 0; j < NUM_TASKS; j++) {
            palma_val_t edge = palma_matrix_get(sched->system, i, j);
            if (edge != PALMA_NEG_INF) {
                palma_val_t pred_end = palma_scheduler_get_completion(sched, j);
                if (pred_end > max_pred) max_pred = pred_end;
            }
        }
        start = max_pred;
        
        const char *name = sched->task_names[i] ? sched->task_names[i] : "Unknown";
        printf("│ %-18s │ %5dms │ %7dms │\n", 
               name, start, completion + task_durations[i]);
    }
    
    printf("└────────────────────┴─────────┴───────────┘\n\n");
    
    palma_val_t total = palma_scheduler_get_completion(sched, TASK_SERVICES) + 
                        task_durations[TASK_SERVICES];
    printf("Total boot time: %dms\n\n", total);
    
    /* Critical path */
    printf("Critical path (longest dependency chain):\n");
    printf("  HW_INIT → KERNEL → FILESYSTEM → SERVICES\n");
    printf("  10 + 20 + 30 + 10 = 70ms ✓\n\n");
    
    /* Export to DOT file for visualization */
    const char *node_names[] = {
        "HW_INIT", "KERNEL", "DRIVERS", "NETWORK", "FILESYSTEM", "SERVICES"
    };
    palma_matrix_export_dot(sched->system, "boot_sequence.dot", PALMA_MAXPLUS, node_names);
    printf("Exported dependency graph to 'boot_sequence.dot'\n");
    printf("Visualize with: dot -Tpng boot_sequence.dot -o boot_sequence.png\n\n");
    
    /* Cyclic scheduling example */
    printf("=== Cyclic Scheduling (Periodic Systems) ===\n\n");
    
    palma_scheduler_t *cyclic = palma_scheduler_create(3, true);
    
    /* 3-task cycle: A → B → C → A */
    palma_scheduler_set_name(cyclic, 0, "Task A");
    palma_scheduler_set_name(cyclic, 1, "Task B");
    palma_scheduler_set_name(cyclic, 2, "Task C");
    
    palma_scheduler_add_constraint(cyclic, 0, 1, 10);  /* A → B */
    palma_scheduler_add_constraint(cyclic, 1, 2, 15);  /* B → C */
    palma_scheduler_add_constraint(cyclic, 2, 0, 20);  /* C → A (feedback) */
    
    printf("Cyclic system matrix:\n");
    palma_matrix_print(cyclic->system, "A_cyclic", PALMA_MAXPLUS, stdout);
    
    palma_val_t cycle_time = palma_scheduler_cycle_time(cyclic);
    double throughput = palma_scheduler_throughput(cyclic);
    
    printf("\nCycle time (tropical eigenvalue λ): %dms\n", cycle_time);
    printf("Throughput: %.2f iterations/second\n", throughput * 1000);
    printf("\nInterpretation: The system can complete one full cycle every %dms.\n", cycle_time);
    
    palma_scheduler_destroy(cyclic);
    palma_scheduler_destroy(sched);
    
    printf("\n=== Example Complete ===\n");
    return 0;
}
