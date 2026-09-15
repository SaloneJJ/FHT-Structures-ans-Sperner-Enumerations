# FHT-Structures-and-Sperner-Enumerations
Python scripts for the exact enumeration, structural generation, symmetry analysis and metric properties of Fractal Hyper-Tree (FHT) structures.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22163671.svg)](https://doi.org/10.5281/zenodo.22163671)

## Overview
This repository provides computational frameworks to enumerate, classify, and analyze Finite Fractal Hyper-Trees (FHTs). 

<ins>**Finite Fractal Hyper-Tree Structure**</ins></br>
Let $`\mathcal{F}`$ denote a Fractal Hyper-Tree (FHT) defined over a finite set of vertices $`V = \{v_1, v_2, \dots, v_{|V|}\}`$. It is structured as a graded family of hypergraph sets $`(\mathcal{H}_n)_{n=0}^N`$, where $`N \in \mathbb{N}`$ is a fixed maximum depth. For each level $`n`$, the set of hypergraphs is written as $`\mathcal{H}_n = \{H_{n,i}\}_{i \in I_n}`$, where each element is a pair $`H_{n,i} = (V_{n,i}, E_{n,i})`$ with $`V_{n,i} \subseteq V`$, and $I_n$ is a finite index set.

The family satisfies the following axiomatic conditions:</br>
**Partition by Connected Components:** </br>
For each level $`n \in \{0, \dots, N\}`$, the set $`\mathcal{H}_n = \{H_{n,i}\}_{i \in I_n}`$ consists of hypergraphs $`H_{n,i}`$ that correspond precisely to the disjoint connected components of the structure at level $`n`$. Here, connectedness is defined in terms of edge-path connectedness: two vertices $`x, y \in V_{n,i}`$ are connected if there exists a sequence of hyperedges $`(e_1, \dots, e_k) \in E_{n,i}^k`$ such that $`x \in e_1$, $y \in e_k`$, and $`e_j \cap e_{j+1} \neq \emptyset`$ for all $`j \in \{1, \dots, k-1\}`$. Each $`H_{n,i}`$ captures one such maximal connected component.</br>
**Top-Level Singularity and Connectedness:** </br>
At the maximum depth level $`N`$, the set $`\mathcal{H}_N`$ consists of a single hypergraph $`H_{N,1} = (V, E_{N,1})`$ whose vertex set spans the entire space $`V`$, and this top-level hypergraph is edge-path connected.</br>
**Base Level Structure:** </br>
The base level $`\mathcal{H}_0 = \{H_{0,i}\}_{i \in I_0}`$ partitions the base elements into trivial singleton hypergraphs. Specifically, for each element $`v \in V`$, there exists a corresponding base hypergraph $`H_{0,i} = (\{v\}, \{\{v\}\})`$, such that the union of all base vertex sets recovers $`V`$:

$$\bigcup_{i \in I_0} V_{0,i} = V$$

**Fractal Inclusion:** </br>
This axiom defines the specific structural class of FHTs by prohibiting partial intersections: for each level $`n \in \{1, \dots, N\}`$, for every $`H \in \mathcal{H}_n`$, for every $`e \in E_H`$, and for every $`H' \in \mathcal{H}_k`$ with $`0 \le k < n`$, if the hyperedge $`e`$ intersects the vertex set $`V_{H'}`$ ($`e \cap V_{H'} \neq \emptyset`$), then $`e`$ must strictly contain the entire vertex set $`V_{H'}`$:

$$e \cap V_{H'} \neq \emptyset \implies V_{H'} \subsetneq e$$

**Hierarchical Filiation:** </br>
For any level $`n \in \{1, \dots, N\}`$ and each $`H \in \mathcal{H}_n`$, there exists a lower-level hypergraph $`H' \in \mathcal{H}_{n-1}`$ and at least one hyperedge $`e \in E_H`$ such that the entire vertex set of $`H'`$ is strictly contained within $`e`$:

$$\exists H' \in \mathcal{H}_{n-1}, \quad \exists e \in E_H \quad \text{such that} \quad V_{H'} \subsetneq e$$

**Hyperedge Non-Triviality:** </br>
To prevent degenerate unary branching, for all $`n \in \{1, \dots, N\}`$, for every $`H \in \mathcal{H}_n`$, and for every hyperedge $`e \in E_H`$, the set $`C`$ of maximal lower-level hypergraphs strictly contained in $`e`$—defined as:

$$C = \left\{ H' \in \bigcup_{k=0}^{n-1} \mathcal{H}_k \;\middle|\; (V_{H'} \subsetneq e) \text{ and } ( \nexists H'' \in \bigcup_{k=0}^{n-1} \mathcal{H}_k \bigm| V_{H'} \subsetneq V_{H''} \subseteq e ) \right\}$$

satisfies the non-triviality condition of branching into at least two components: $`|C| \ge 2`$


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
