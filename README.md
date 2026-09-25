# Fractal Hyper-Trees (FHT) — enumeration, symmetries, renormalized ultrametric

Code accompanying the preprint *Fractal Hyper-Trees: ...* (J.-J. Salone, 2026).

[![DOI](https://zenodo.org/badge/1350712518.svg)](https://doi.org/10.5281/zenodo.22163670)

Version 4.0.0: doi:10.5281/zenodo.22948863

## Installation
    python3 -m pip install -r requirements.txt

## Usage
    python3 fht.py                                   # interactive mode
    python3 fht.py --check                           # checks of the host lemmas
    python3 fht.py --vmax 9 --methods recurrence     # counts f(v), v <= 9
    python3 fht.py --vmax 4 --methods brute --list --sym --dist

## Correspondence with the paper
- brute force: Proposition "Recursive decomposition", Lemma "Grading", Tables of symmetry data and metric classes
- recurrence: Corollary "Recurrence, explicit form", Lemmas "single top point" and "two top points"

## Results
f(v), v = 1..9: 1, 1, 4, 21, 199, 16296, 489998312, 1392195548400053714,
789204635842035014359735420510795311 (see results/).

## Note on earlier versions
Versions <= 3.x used a different definition and over-counted I(v) for v >= 6.
