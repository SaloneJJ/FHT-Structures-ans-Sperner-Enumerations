# FHT-Structures-and-Sperner-Enumerations
Python scripts for the exact enumeration, structural generation, and symmetry analysis of Fractal Hyper-Tree (FHT) structures.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22163671.svg)](https://doi.org/10.5281/zenodo.22163671)

## Overview
This repository provides computational frameworks to enumerate, classify, and analyze Finite Fractal Hyper-Trees (FHTs). It covers both **flat Sperner hypergraph topologies ($S(v)$)** and **hierarchical imbricated structures ($I(v)$)** through formal integer partition compositions.

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
