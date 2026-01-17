# Mathematical Background

This document provides the mathematical foundations of tropical algebra as implemented in PALMA.

## Table of Contents

- [Tropical Semirings](#tropical-semirings)
- [Tropical Linear Algebra](#tropical-linear-algebra)
- [Graph Interpretation](#graph-interpretation)
- [Tropical Eigenvalues](#tropical-eigenvalues)
- [Applications](#applications)
- [References](#references)

---

## Tropical Semirings

### Definition

A **semiring** (S, ⊕, ⊗, ε, e) consists of:
- A set S
- Addition ⊕ (commutative, associative, with identity ε)
- Multiplication ⊗ (associative, with identity e)
- Multiplication distributes over addition
- ε is absorbing: a ⊗ ε = ε ⊗ a = ε

### Tropical Semirings in PALMA

| Semiring | S | ⊕ | ⊗ | ε | e |
|----------|---|---|---|---|---|
| Max-Plus (ℝ_max) | ℝ ∪ {-∞} | max | + | -∞ | 0 |
| Min-Plus (ℝ_min) | ℝ ∪ {+∞} | min | + | +∞ | 0 |
| Max-Min | ℝ ∪ {-∞, +∞} | max | min | -∞ | +∞ |
| Min-Max | ℝ ∪ {-∞, +∞} | min | max | +∞ | -∞ |
| Boolean | {0, 1} | ∨ (OR) | ∧ (AND) | 0 | 1 |

### Key Property: Idempotence

All tropical semirings are **idempotent**: a ⊕ a = a.

This makes them fundamentally different from classical algebra and is the key to their optimization properties.

---

## Tropical Linear Algebra

### Matrix Operations

For matrices A, B over a tropical semiring:

**Addition** (element-wise):
```
(A ⊕ B)_ij = A_ij ⊕ B_ij
```

**Multiplication**:
```
(A ⊗ B)_ij = ⊕_k (A_ik ⊗ B_kj)
```

For max-plus, this becomes:
```
(A ⊗ B)_ij = max_k (A_ik + B_kj)
```

### Identity Matrix

The tropical identity matrix I has:
- I_ii = e (identity element)
- I_ij = ε (zero element) for i ≠ j

### Kleene Star (Closure)

For an n×n matrix A:
```
A* = I ⊕ A ⊕ A² ⊕ ... ⊕ A^(n-1)
```

This converges in finite steps for idempotent semirings.

---

## Graph Interpretation

### Adjacency Matrix

A weighted directed graph G = (V, E, w) is represented by its adjacency matrix A where:
- A_ij = w(i,j) if edge (i,j) exists
- A_ij = ε (semiring zero) otherwise

### Path Weights

The (i,j) entry of A^k gives the optimal k-hop path weight from i to j:
- **Min-plus**: minimum weight (shortest path)
- **Max-plus**: maximum weight (longest/critical path)
- **Max-min**: maximum bottleneck (widest path)
- **Boolean**: reachability (1 if path exists)

### All-Pairs Paths

The closure A* gives optimal paths of any length:
```
(A*)_ij = optimal path weight from i to j
```

For min-plus, this is the Floyd-Warshall algorithm.

---

## Tropical Eigenvalues

### Definition

A scalar λ is a **tropical eigenvalue** of matrix A if there exists a non-trivial vector v such that:
```
A ⊗ v = λ ⊗ v
```

where λ ⊗ v means adding λ to each component (in max-plus/min-plus).

### Maximum Cycle Mean

For max-plus, the eigenvalue equals the **maximum cycle mean**:
```
λ(A) = max_{C} (w(C) / |C|)
```

where w(C) is the total weight of cycle C and |C| is its length.

### Karp's Algorithm

PALMA computes eigenvalues using Karp's algorithm:
```
λ = max_i min_{0≤k<n} (A^n)_ii - (A^k)_ii) / (n - k)
```

This runs in O(n³) time and O(n²) space.

### Significance

The eigenvalue determines:
- **Asymptotic growth rate**: (A^k)_ij ≈ k·λ as k → ∞
- **Cycle time**: In periodic systems, λ is the minimum period
- **Throughput**: 1/λ items per time unit

---

## Applications

### 1. Shortest Paths (Min-Plus)

Standard shortest path problem:
```
d_ij = minimum total edge weight from i to j
```

**Solution**: d = A* (min-plus closure)

This unifies:
- Dijkstra's algorithm (single source)
- Bellman-Ford (negative weights)
- Floyd-Warshall (all pairs)

### 2. Critical Path (Max-Plus)

Project scheduling with task dependencies:
```
s_j = earliest start time of task j
```

If task i must complete before j, with duration d_i:
```
s_j ≥ s_i + d_i
```

**Solution**: s = A* ⊗ r where r is the ready time vector.

### 3. Maximum Bandwidth (Max-Min)

Network capacity problem:
```
b_ij = maximum bandwidth achievable from i to j
```

**Solution**: b = A* (max-min closure)

### 4. Reachability (Boolean)

Transitive closure:
```
R_ij = 1 if path exists from i to j
```

**Solution**: R = A* (Boolean closure)

### 5. Periodic Systems

Production lines, transportation schedules:
```
x(k+1) = A ⊗ x(k)
```

- **Cycle time**: λ(A) = tropical eigenvalue
- **Throughput**: 1/λ(A)
- **Bottleneck**: critical cycle achieving λ

---

## References

1. **Baccelli, F., Cohen, G., Olsder, G.J., Quadrat, J.P.** (1992). *Synchronization and Linearity: An Algebra for Discrete Event Systems*. Wiley.
   - [Free online version](https://www.rocq.inria.fr/metalau/cohen/SED/book-online.html)
   - The definitive reference for max-plus algebra

2. **Heidergott, B., Olsder, G.J., van der Woude, J.** (2006). *Max Plus at Work*. Princeton University Press.
   - Accessible introduction with applications

3. **Butkovič, P.** (2010). *Max-linear Systems: Theory and Algorithms*. Springer.
   - Comprehensive algorithmic treatment

4. **Karp, R.M.** (1978). A characterization of the minimum cycle mean in a digraph. *Discrete Mathematics*, 23(3), 309-311.
   - Original eigenvalue algorithm

5. **Mohri, M.** (2002). Semiring frameworks and algorithms for shortest-distance problems. *Journal of Automata, Languages and Combinatorics*, 7(3), 321-350.
   - Unifying semiring perspective
