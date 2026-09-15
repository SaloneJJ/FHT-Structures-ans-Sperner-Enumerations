# FHT-Structures-and-Sperner-Enumerations
Python scripts for the exact enumeration, structural generation, symmetry analysis and metric properties of Fractal Hyper-Tree (FHT) structures.
This release introduces complete enumeration capabilities for Fractal Hyper-Trees (FHTs) alongside rigorous symmetry analysis, including automorphism group sizes, invariant (fixed) points, and vertex orbits.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22163671.svg)](https://doi.org/10.5281/zenodo.22163671)

## Overview
This repository provides computational frameworks to enumerate, classify, and analyze Finite Fractal Hyper-Trees (FHTs). 

Definition 1 (Finite Fractal Hyper-Tree Structure).  Let  denote a Fractal Hyper-Tree (FHT) defined over a finite set of vertices . It is structured as a graded family of hypergraph sets , where  is a fixed maximum depth. For each level , the set of hypergraphs is written as , where each element is a pair  with , and  is a finite index set.
The family satisfies the following axiomatic conditions:
Partition by Connected Components: For each level , the set  consists of hypergraphs  that correspond precisely to the disjoint connected components of the structure at level . Here, connectedness is defined in terms of edge-path connectedness: two vertices  are connected if there exists a sequence of hyperedges  such that , , and  for all . Each  captures one such maximal connected component.
Top-Level Singularity and Connectedness: At the maximum depth level , the set  consists of a single hypergraph  whose vertex set spans the entire space , and this top-level hypergraph is edge-path connected.
Base Level Structure: The base level  partitions the base elements into trivial singleton hypergraphs. Specifically, for each element , there exists a corresponding base hypergraph , such that the union of all base vertex sets recovers : 
Fractal Inclusion (Strict Hierarchy): This axiom defines the specific structural class of FHTs by prohibiting partial intersections: for each level , for every , for every , and for every  with , if the hyperedge  intersects the vertex set  (), then  must strictly contain the entire vertex set : 
Hierarchical Filiation: For any level  and each , there exists a lower-level hypergraph  and at least one hyperedge  such that the entire vertex set of  is strictly contained within : 
Hyperedge Non-Triviality (Maximal Children Coverage): To prevent degenerate unary branching, for all , for every , and for every hyperedge , the set  of maximal lower-level hypergraphs strictly contained in —defined as: $$C = \left\{ H' \in \bigcup_{k=0}^{n-1} \mathcal{H}_k \;\middle|\; (V_{H'} \subsetneq e) \text{ and } ( \nexists H'' \in \bigcup_{k=0}^{n-1} \mathcal{H}_k \bigm| V_{H'} \subsetneq V_{H''} \subseteq e ) \right\}$$ satisfies the non-triviality condition of branching into at least two components: 

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
