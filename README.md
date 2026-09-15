# FHT-Structures-and-Sperner-Enumerations
Python scripts for the exact enumeration, structural generation, symmetry analysis and metric properties of Fractal Hyper-Tree (FHT) structures.
This release introduces complete enumeration capabilities for Fractal Hyper-Trees (FHTs) alongside rigorous symmetry analysis, including automorphism group sizes, invariant (fixed) points, and vertex orbits.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22163671.svg)](https://doi.org/10.5281/zenodo.22163671)

## Overview
This repository provides computational frameworks to enumerate, classify, and analyze Finite Fractal Hyper-Trees (FHTs). 

# Definition of a Fractal Hyper-Tree Structure

Let $\mathcal{F}$ denote a Fractal Hyper-Tree (FHT) defined over a finite set of vertices $V = \{v_1, v_2, \dots, v_{|V|}\}$. It is structured as a graded family of hypergraph sets $(\mathcal{H}_n)_{n=0}^N$, where $N \in \mathbb{N}$ is a fixed maximum depth. For each level $n$, the set of hypergraphs is written as $\mathcal{H}_n = \{H_{n,i}\}_{i \in I_n}$, where each element is a pair $H_{n,i} = (V_{n,i}, E_{n,i})$ with $V_{n,i} \subseteq V$, and $I_n$ is a finite index set.

The family satisfies the following axiomatic conditions:

### (i) Partition by Connected Components
For each level $n \in \{0, \dots, N\}$, the set $\mathcal{H}_n = \{H_{n,i}\}_{i \in I_n}$ consists of hypergraphs $H_{n,i}$ that correspond precisely to the disjoint connected components of the structure at level $n$. Here, connectedness is defined in terms of edge-path connectedness: two vertices $x, y \in V_{n,i}$ are connected if there exists a sequence of hyperedges $(e_1, \dots, e_k) \in E_{n,i}^k$ such that $x \in e_1$, $y \in e_k$, and $e_j \cap e_{j+1} \neq \emptyset$ for all $j \in \{1, \dots, k-1\}$. Each $H_{n,i}$ captures one such maximal connected component.

### (ii) Top-Level Singularity and Connectedness
At the maximum depth level $N$, the set $\mathcal{H}_N$ consists of a single hypergraph $H_{N,1} = (V, E_{N,1})$ whose vertex set spans the entire space $V$, and this top-level hypergraph is edge-path connected.

### (iii) Base Level Structure
The base level $\mathcal{H}_0 = \{H_{0,i}\}_{i \in I_0}$ partitions the base elements into trivial singleton hypergraphs. Specifically, for each element $v \in V$, there exists a corresponding base hypergraph $H_{0,i} = (\{v\}, \{\{v\}\})$, such that the union of all base vertex sets recovers $V$:
```math
\bigcup_{i \in I_0} V_{0,i} = V

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
