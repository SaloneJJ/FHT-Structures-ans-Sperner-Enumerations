# FHT-Structures-and-Sperner-Enumerations
Python scripts for the exact enumeration, structural generation, symmetry analysis and metric properties of Fractal Hyper-Tree (FHT) structures.
This release introduces complete enumeration capabilities for Fractal Hyper-Trees (FHTs) alongside rigorous symmetry analysis, including automorphism group sizes, invariant (fixed) points, and vertex orbits.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22163671.svg)](https://doi.org/10.5281/zenodo.22163671)

## Overview
This repository provides computational frameworks to enumerate, classify, and analyze Finite Fractal Hyper-Trees (FHTs). 

# Definition of a Fractal Hyper-Tree Structure

Let $F$ denote a Fractal Hyper-Tree (FHT) defined over a finite set of vertices $V = \{v_1, v_2, \dots, v_n\}$. It is structured as a graded family of hypergraph sets $\mathcal{H}_n$ for $n = 0, \dots, N$, where $N$ is a fixed maximum depth. For each level $n$, the set of hypergraphs is written as $\mathcal{H}_n = \{H_{n,i} \mid i \in I_n\}$, where each element is a pair $H_{n,i} = (V_{n,i}, E_{n,i})$ with $V_{n,i} \subseteq V$.

The family satisfies the following axiomatic conditions:

### (i) Partition by Connected Components
For each level $n$, the set $\mathcal{H}_n$ consists of hypergraphs $H_{n,i}$ that correspond precisely to the disjoint connected components of the structure at level $n$.

### (ii) Top-Level Singularity and Connectedness
At the maximum depth level $N$, the set $\mathcal{H}_N$ consists of a single hypergraph $H_{N,1} = (V, E_{N,1})$ whose vertex set spans the entire space $V$.

### (iii) Base Level Structure
The base level $\mathcal{H}_0$ partitions the base elements into trivial singleton hypergraphs. Specifically, for each element $v \in V$, the union of all base vertex sets recovers $V$:
$$\bigcup_{i \in I_0} V_{0,i} = V$$

### (iv) Fractal Inclusion (Strict Hierarchy)
For each level $n$, for every $H \in \mathcal{H}_n$, for every $e \in E_H$, and for every lower-level $H'$, if the hyperedge $e$ intersects the vertex set $V_{H'}$, then $e$ must strictly contain the entire vertex set:
$$e \cap V_{H'} \neq \emptyset \implies V_{H'} \subsetneq e$$

### (v) Hierarchical Filiation
For any level $n$ and each $H \in \mathcal{H}_n$, there exists a lower-level hypergraph $H'$ and at least one hyperedge $e$ such that:
$$V_{H'} \subsetneq e$$

### (vi) Hyperedge Non-Triviality
To prevent degenerate unary branching, the set $C$ of maximal lower-level hypergraphs strictly contained in $e$ satisfies the non-triviality condition of branching into at least two components:
$$|C| \ge 2$$



It covers both **flat Sperner hypergraph topologies ($S(v)$)** and **hierarchical imbricated structures ($I(v)$)** through formal integer partition compositions.

## Key Features
- **Rigorous Isomorphism Testing:** Relies on canonical labeling via `pynauty` to filter out isomorphic duplicates.
- **Complete Symmetry Analysis:** Utilizes `pynauty.autgrp` to extract key algebraic metrics:
  - Automorphism group size ($|\operatorname{Aut}(H)|$)
  - Total number of vertex orbits ($k$)
  - Invariant points (fixed points under symmetry)
- **Automated Decomposition:** Computes exact recursive counts for structural hierarchies up to small-scale vertices.

## Dependencies
- Python 3.x
- `pynauty`
- `networkX`

You can install the required packages via pip:
```bash
!pip install pynauty networkx
