#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FHT v4.0 -- Fractal Hyper-Trees: enumeration, symmetries, renormalized ultrametric.
SALONE Jean-Jacques 2026 -- jean-jacques.salone@univ-antilles.fr

Interactive mode : python3 fht.py
Command line     : python3 fht.py --vmax 9 --methods recurrence
                   python3 fht.py --vmax 4 --methods brute,recurrence --list --sym --dist
                   python3 fht.py --check
"""

import argparse
import json
import sys
import time
from collections import Counter
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, product
from math import factorial, prod

# ======================================================================
# Guard-rails (hard limits, also enforced inside the functions)
# ======================================================================
LIMIT_BRUTE = 6        # brute force: v <= 6 (v = 6: long run, several GB of RAM)
LIMIT_RECURRENCE = 9   # recurrence: v <= 9 (S(v) is known up to 9 in A261006)
LIMIT_DETAILS = 6      # listing / symmetries / distances need the brute force
WARN_DETAILS = 5       # above this value the detailed output is very long
GENERIC_MAX_M = 7      # generic Burnside host counting: at most 7 host points
LEMMAS_FROM_M = 7      # hosts on >= 7 points: Lemmas lem:singletop / lem:twotop
TWO_TOP_MAX_U = 5      # Lemma lem:twotop: at most 5 low points
FLAT_INTERNAL_MAX = 6  # S(v) computed internally up to 6, then read from A261006

# OEIS A261006 (connected simplicial complexes = connected clutters = flat FHTs)
A261006 = {1: 1, 2: 1, 3: 3, 4: 14, 5: 157, 6: 15942, 7: 489980450,
           8: 1392195547909966848,
           9: 789204635842035012967539870068113408}
# OEIS A006602 (antichain covers of an unlabeled n-set; differences of A003182)
A006602 = {1: 1, 2: 2, 3: 5, 4: 20, 5: 180, 6: 16143, 7: 489996795,
           8: 1392195548399980210}


class GuardError(Exception):
    """Raised when a computation would exceed a guard-rail."""


def log(*args):
    print(*args, flush=True)


# ======================================================================
# Combinatorial helpers
# ======================================================================
def bits(mask):
    i = 0
    while mask:
        if mask & 1:
            yield i
        mask >>= 1
        i += 1


def popcount(x):
    return bin(x).count("1")


def partitions(n, maxpart=None):
    """Integer partitions of n, parts in non-increasing order."""
    maxpart = n if maxpart is None else maxpart
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


def part_counts(lam):
    r = {}
    for c in lam:
        r[c] = r.get(c, 0) + 1
    return r


def z(lam):
    """Order of the centralizer of a permutation of cycle type lam."""
    return prod(c ** k * factorial(k) for c, k in part_counts(lam).items())


def mult_partitions(F, k):
    """(lam, w): ways to choose k blocks among F classes with multiplicities lam."""
    for lam in partitions(k):
        l = len(lam)
        if l > F:
            continue
        w = prod(F - i for i in range(l))
        w //= prod(factorial(x) for x in part_counts(lam).values())
        yield lam, w


def perm_from_cycles(cycles, start, perm):
    p = start
    for c in cycles:
        for i in range(c):
            perm[p + i] = p + (i + 1) % c
        p += c
    return p


def image(mask, perm):
    out, i = 0, 0
    while mask:
        if mask & 1:
            out |= 1 << perm[i]
        mask >>= 1
        i += 1
    return out


def cycles_of(perm, mask=None):
    """Cycles (as lists) of a permutation given as a tuple, inside mask if given."""
    seen, out = set(), []
    for i in range(len(perm)):
        if (mask is not None and not (mask >> i) & 1) or i in seen:
            continue
        cyc, j = [], i
        while j not in seen:
            seen.add(j)
            cyc.append(j)
            j = perm[j]
        out.append(cyc)
    return out


def relabel(perm, W):
    """Restriction of perm to the invariant set W, relabelled on range(|W|)."""
    elems = list(bits(W))
    idx = {e: i for i, e in enumerate(elems)}
    return tuple(idx[perm[e]] for e in elems)


# ======================================================================
# Host counting 1: generic Burnside (fixed admissible hosts)
# ======================================================================
def connected(edges, full):
    reach, rest = edges[0], list(edges[1:])
    changed = True
    while changed and rest:
        changed, keep = False, []
        for e in rest:
            if e & reach:
                reach |= e
                changed = True
            else:
                keep.append(e)
        rest = keep
    return reach == full


_fixg_cache = {}


def fix_generic(m, top_cycles, low_cycles):
    """Admissible hosts on m points (top points 0..t-1) fixed by a permutation."""
    if m > GENERIC_MAX_M:
        raise GuardError(f"generic host counting is limited to {GENERIC_MAX_M} "
                         f"points (requested: {m})")
    key = (m, top_cycles, low_cycles)
    if key in _fixg_cache:
        return _fixg_cache[key]
    perm = {}
    t = perm_from_cycles(top_cycles, 0, perm)
    perm_from_cycles(low_cycles, t, perm)
    full, topmask = (1 << m) - 1, (1 << t) - 1
    seen, orbits = set(), []
    for s in range(1, full + 1):
        if s in seen or popcount(s) < 2 or not (s & topmask):
            continue
        orb, x = [], s
        while x not in seen:
            seen.add(x)
            orb.append(x)
            x = image(x, perm)
        orbits.append(orb)
    n = len(orbits)
    union = [0] * n
    for i, o in enumerate(orbits):
        for a in o:
            union[i] |= a
    compat = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if not any((a & b) in (a, b) for a in orbits[i] for b in orbits[j]):
                compat[i] |= 1 << j

    def dfs(cand, cover, edges):
        if cover != full:
            rest, c = 0, cand
            while c:
                lb = c & -c
                rest |= union[lb.bit_length() - 1]
                c ^= lb
            if cover | rest != full:
                return 0
        total = 1 if (cover == full and connected(edges, full)) else 0
        while cand:
            lb = cand & -cand
            i = lb.bit_length() - 1
            cand ^= lb
            total += dfs(cand & compat[i], cover | union[i], edges + orbits[i])
        return total

    res = dfs((1 << n) - 1, 0, []) if n else 0
    _fixg_cache[key] = res
    return res


# ======================================================================
# Host counting 2: Lemma lem:twotop (two top points), via up-sets
# ======================================================================
_UPSETS = {}


def upsets(k):
    """All up-sets of 2^{0..k-1}, as bitmasks over the 2^k subsets."""
    if k > TWO_TOP_MAX_U:
        raise GuardError(f"up-sets are enumerated on at most {TWO_TOP_MAX_U} points")
    if k not in _UPSETS:
        if k == 0:
            _UPSETS[k] = [0, 1]
        else:
            prev, half = upsets(k - 1), 1 << (k - 1)
            _UPSETS[k] = [f0 | (f1 << half) for f0 in prev for f1 in prev
                          if f0 & ~f1 == 0]
    return _UPSETS[k]


_inv_cache = {}


def inv_upsets(pi):
    """Up-sets of 2^{range(k)} invariant under the permutation pi."""
    if pi not in _inv_cache:
        k = len(pi)
        sub_map = [image(A, pi) for A in range(1 << k)]

        def fam_image(F):
            out = 0
            for A in bits(F):
                out |= 1 << sub_map[A]
            return out

        _inv_cache[pi] = [F for F in upsets(k) if fam_image(F) == F]
    return _inv_cache[pi]


def below_counts(Rs, Ss):
    """For each R in Rs, number of S in Ss with S contained in R."""
    try:
        import numpy as np
    except ImportError:
        return [sum(1 for S in Ss if S & ~R == 0) for R in Rs]
    S = np.array(Ss, dtype=np.uint64)
    out = []
    for i in range(0, len(Rs), 512):
        R = np.array(Rs[i:i + 512], dtype=np.uint64)
        mask = (S[None, :] & ~R[:, None]) == 0
        out.extend(int(x) for x in mask.sum(axis=1))
    return out


_M_cache, _N2_cache, _Ns_cache, _c_cache = {}, {}, {}, {}


def M_inv(pi):
    if pi not in _M_cache:
        _M_cache[pi] = len(inv_upsets(pi))
    return _M_cache[pi]


def sum_N2(pi):
    """Sum over invariant up-sets R of N_pi(R)^2."""
    if pi not in _N2_cache:
        R = inv_upsets(pi)
        _N2_cache[pi] = sum(c * c for c in below_counts(R, R))
    return _N2_cache[pi]


def sum_Nswap(pi):
    """Sum over pi-invariant up-sets R of N_{pi^2}(R)."""
    if pi not in _Ns_cache:
        pi2 = tuple(pi[pi[i]] for i in range(len(pi)))
        _Ns_cache[pi] = sum(below_counts(inv_upsets(pi), inv_upsets(pi2)))
    return _Ns_cache[pi]


def c_cover(pi):
    """Invariant antichains with union range(k); convention c = 1 if k = 0."""
    if pi in _c_cache:
        return _c_cache[pi]
    if len(pi) == 0:
        return 1
    cyc = [sum(1 << x for x in c) for c in cycles_of(pi)]
    total = 0
    for sel in range(1 << len(cyc)):
        W = 0
        for j in range(len(cyc)):
            if (sel >> j) & 1:
                W |= cyc[j]
        total += (-1) ** (len(cyc) - popcount(sel)) * M_inv(relabel(pi, W))
    _c_cache[pi] = total
    return total


_fixt_cache = {}


def fix_twotop(u, top_cycles, low_cycles):
    """Fix(g) for two top points, by Lemma lem:twotop."""
    if u > TWO_TOP_MAX_U:
        raise GuardError(f"Lemma lem:twotop is used for at most {TWO_TOP_MAX_U} "
                         f"low points (requested: {u})")
    key = (u, top_cycles, low_cycles)
    if key in _fixt_cache:
        return _fixt_cache[key]
    perm = {}
    perm_from_cycles(low_cycles, 0, perm)
    h = tuple(perm[i] for i in range(u))
    cyc_lists = cycles_of(h)
    cyc = [sum(1 << x for x in c) for c in cyc_lists]
    swap = (top_cycles == (2,))
    full = (1 << u) - 1

    F = 0
    for sel in range(1 << len(cyc)):
        W = sum(cyc[j] for j in range(len(cyc)) if (sel >> j) & 1)
        sign = (-1) ** (len(cyc) - popcount(sel))
        pi = relabel(h, W)
        a = (-1 + sum_Nswap(pi)) if swap else (1 - 2 * M_inv(pi) + sum_N2(pi))
        F += sign * a

    D = 0
    if not swap:
        for sel in range(1 << len(cyc)):
            U1 = sum(cyc[j] for j in range(len(cyc)) if (sel >> j) & 1)
            D += c_cover(relabel(h, U1)) * c_cover(relabel(h, full ^ U1))
    elif all(len(c) % 2 == 0 for c in cyc_lists):
        h2 = tuple(h[h[i]] for i in range(u))
        for sel in range(1 << len(cyc_lists)):
            U1 = 0
            for j, c in enumerate(cyc_lists):
                start = (sel >> j) & 1
                for x in c[start::2]:
                    U1 |= 1 << x
            D += c_cover(relabel(h2, U1))

    _fixt_cache[key] = F - D
    return F - D


# ======================================================================
# Host orbits H(alpha; beta)
# ======================================================================
def fix_count(m, top_cycles, low_cycles, use_lemmas=True):
    if use_lemmas and m >= LEMMAS_FROM_M and sum(top_cycles) == 2:
        return fix_twotop(m - 2, top_cycles, low_cycles)
    return fix_generic(m, top_cycles, low_cycles)


_H_cache = {}


def host_orbits(top_parts, low_parts, use_lemmas=True):
    key = (top_parts, low_parts, use_lemmas)
    if key in _H_cache:
        return _H_cache[key]
    m = sum(top_parts) + sum(low_parts)
    if (use_lemmas and m >= LEMMAS_FROM_M and top_parts == (1,)
            and len(low_parts) == 1 and low_parts[0] in A006602):
        res = A006602[low_parts[0]]                       # Lemma lem:singletop
    else:
        groups = [(k, True) for k in top_parts] + [(k, False) for k in low_parts]
        order = prod(factorial(k) for k, _ in groups)
        per_group = [[(lam, factorial(k) // z(lam)) for lam in partitions(k)]
                     for k, _ in groups]
        total = 0
        for combo in product(*per_group):
            cnt = prod(c for _, c in combo)
            topc = tuple(sorted((c for (lam, _), (_, top) in zip(combo, groups)
                                 if top for c in lam), reverse=True))
            lowc = tuple(sorted((c for (lam, _), (_, top) in zip(combo, groups)
                                 if not top for c in lam), reverse=True))
            total += cnt * fix_count(m, topc, lowc, use_lemmas)
        assert total % order == 0, "Burnside: non-integer orbit count"
        res = total // order
    _H_cache[key] = res
    return res


# ======================================================================
# Method 1: recurrence (Corollary cor:explicit)
# ======================================================================
def run_recurrence(vmax, use_lemmas=True):
    if vmax > LIMIT_RECURRENCE:
        raise GuardError(f"the recurrence is limited to v <= {LIMIT_RECURRENCE}")
    F = {1: {0: 1}}
    res = {1: {"S": 1, "I": 0, "f": 1, "by_depth": {0: 1}}}
    log("[recurrence] v=1: S=1  I=0  f=1")
    for v in range(2, vmax + 1):
        t0 = time.time()
        if v <= FLAT_INTERNAL_MAX:
            S = host_orbits((v,), (), use_lemmas)
            if S != A261006[v]:
                log(f"  WARNING: computed S({v})={S} differs from A261006 "
                    f"({A261006[v]})")
            src = "Burnside"
        else:
            S, src = A261006[v], "A261006"
        Fv = {1: S}
        types = [(s, d) for s in range(1, v) for d in sorted(F[s])]

        def process(chosen):
            if sum(k for _, _, k in chosen) < 2:
                return
            L = max(d for _, d, _ in chosen)
            if L == 0:
                return                                   # flat case, counted in S
            lists = [list(mult_partitions(F[s][d], k)) for s, d, k in chosen]
            for combo in product(*lists):
                w = prod(x[1] for x in combo)
                top = tuple(sorted((p for (lam, _), (_, d, _) in zip(combo, chosen)
                                    if d == L for p in lam), reverse=True))
                low = tuple(sorted((p for (lam, _), (_, d, _) in zip(combo, chosen)
                                    if d < L for p in lam), reverse=True))
                Fv[L + 1] = Fv.get(L + 1, 0) + w * host_orbits(top, low, use_lemmas)

        def rec(i, remaining, chosen):
            if remaining == 0:
                process(chosen)
                return
            if i == len(types):
                return
            s, d = types[i]
            for k in range(remaining // s, -1, -1):
                rec(i + 1, remaining - k * s, chosen + [(s, d, k)] if k else chosen)

        rec(0, v, [])
        F[v] = Fv
        f = sum(Fv.values())
        res[v] = {"S": S, "I": f - S, "f": f, "by_depth": dict(sorted(Fv.items()))}
        log(f"[recurrence] v={v}: S={S} ({src})  I={f - S}  f={f}  "
            f"depth={dict(sorted(Fv.items()))}  ({time.time() - t0:.1f}s)")
    return res


# ======================================================================
# Method 2: brute force (Proposition prop:recursive + pynauty)
# ======================================================================
def nauty_graph(v, edges):
    from pynauty import Graph
    n = v + len(edges)
    adj = {i: [] for i in range(n)}
    for j, e in enumerate(edges):
        for x in bits(e):
            adj[x].append(v + j)
            adj[v + j].append(x)
    return Graph(n, directed=False, adjacency_dict=adj,
                 vertex_coloring=[set(range(v)), set(range(v, n))])


def host_clutters(m):
    """All labelled connected clutters on m points, covering (no pruning)."""
    if m > LIMIT_BRUTE:
        raise GuardError(f"labelled clutters are enumerated on at most "
                         f"{LIMIT_BRUTE} points")
    full = (1 << m) - 1
    cands = [s for s in range(1, full + 1) if popcount(s) >= 2]
    out = []

    def bt(i, chosen, cover):
        if chosen and cover == full and connected(chosen, full):
            out.append(tuple(chosen))
        for j in range(i, len(cands)):
            s = cands[j]
            if all((s & e) != s and (s & e) != e for e in chosen):
                chosen.append(s)
                bt(j + 1, chosen, cover | s)
                chosen.pop()

    bt(0, [], 0)
    return out


def edge_index_map(v):
    """Index of each edge in the order (size, lexicographic): 1, 2, ..., 12, 13, ..."""
    order = [sum(1 << x for x in comb)
             for r in range(1, v + 1) for comb in combinations(range(v), r)]
    return {mask: i for i, mask in enumerate(order)}


def run_brute(vmax):
    if vmax > LIMIT_BRUTE:
        raise GuardError(f"the brute force is limited to v <= {LIMIT_BRUTE}")
    try:
        from pynauty import certificate
    except ImportError:
        raise GuardError("the brute force needs pynauty "
                         "(python3 -m pip install pynauty)")
    fht = {1: [{"edges": (1,), "depth": 0, "key": (0,)}]}
    res = {1: {"S": 1, "I": 0, "f": 1, "by_depth": {0: 1}}}
    log("[brute] v=1: S=1  I=0  f=1")
    hosts = {}
    for v in range(2, vmax + 1):
        t0 = time.time()
        idx = edge_index_map(v)
        seen, classes = {}, []
        plist = [p for p in partitions(v) if len(p) >= 2]
        for pi, p in enumerate(plist):
            m = len(p)
            if m not in hosts:
                hosts[m] = host_clutters(m)
            flat = (m == v)
            sizes = sorted(set(p), reverse=True)
            choice_lists = [list(combinations_with_replacement(
                range(len(fht[s])), p.count(s))) for s in sizes]
            for choice in product(*choice_lists):
                blocks, offset = [], 0                   # (mask, depth, edges)
                for s, idxs in zip(sizes, choice):
                    for i in idxs:
                        c = fht[s][i]
                        mask = ((1 << s) - 1) << offset
                        edges = [] if s == 1 else [e << offset for e in c["edges"]]
                        blocks.append((mask, c["depth"], edges))
                        offset += s
                L = max(b[1] for b in blocks)
                topmask = sum(1 << j for j, b in enumerate(blocks) if b[1] == L)
                base = [e for b in blocks for e in b[2]]
                for host in hosts[m]:
                    if any(not (h & topmask) for h in host):          # axiom (A5)
                        continue
                    new = [sum(blocks[j][0] for j in bits(h)) for h in host]
                    E = tuple(sorted(base + new))
                    key = (tuple(sorted(idx[e] for e in E)) if flat
                           else (pi, choice, len(host), host))
                    cert = certificate(nauty_graph(v, E))
                    if cert in seen:
                        cl = classes[seen[cert]]
                        assert cl["decomp"] == (p, choice), \
                            f"non-unique decomposition: {edges_str(v, E)}"
                        if key < cl["key"]:
                            cl["key"], cl["edges"] = key, E
                        continue
                    seen[cert] = len(classes)
                    classes.append({"edges": E, "depth": L + 1, "key": key,
                                    "decomp": (p, choice)})
        classes.sort(key=lambda c: (c["depth"], c["key"]))
        fht[v] = classes
        by_depth = dict(sorted(Counter(c["depth"] for c in classes).items()))
        S = by_depth.get(1, 0)
        res[v] = {"S": S, "I": len(classes) - S, "f": len(classes),
                  "by_depth": by_depth}
        log(f"[brute] v={v}: S={S}  I={len(classes) - S}  f={len(classes)}  "
            f"depth={by_depth}  ({time.time() - t0:.1f}s)")
    return fht, res


# ======================================================================
# Symmetries
# ======================================================================
def compose(p, q):
    return tuple(p[q[i]] for i in range(len(q)))


def group_closure(gens, n):
    ident = tuple(range(n))
    elements, frontier = {ident}, [ident]
    while frontier:
        new = []
        for a in frontier:
            for g in gens:
                b = compose(g, a)
                if b not in elements:
                    elements.add(b)
                    new.append(b)
        frontier = new
    return elements


def element_order(p):
    ident, q, k = tuple(range(len(p))), p, 1
    while q != ident:
        q, k = compose(p, q), k + 1
    return k


def group_name(elements, gens):
    """Name from the order and element-order statistics (reliable for degree <= 6)."""
    n = len(elements)
    o = Counter(element_order(p) for p in elements)
    ab = all(compose(a, b) == compose(b, a) for a in gens for b in gens)
    names = {1: "{e}", 2: "C2", 3: "C3", 5: "C5"}
    if n in names:
        return names[n]
    if n == 4:
        return "C4" if o[4] else "V4"
    if n == 6:
        return "C6" if ab else "S3"
    if n == 8:
        if ab:
            return "C8" if o[8] else ("C4 x C2" if o[4] else "C2^3")
        return "D4" if o[2] == 5 else "Q8"
    if n == 10:
        return "C10" if ab else "D5"
    if n == 12:
        if ab:
            return "C12" if o[12] else "C6 x C2"
        if o[2] == 3 and o[3] == 8:
            return "A4"
        return "D6" if o[2] == 7 else "Dic3"
    if n == 16 and not ab and o[2] == 11:
        return "D4 x C2"
    if n == 24 and o[2] == 9 and o[3] == 8 and o[4] == 6:
        return "S4"
    if n == 48 and not ab and o[2] == 19:
        return "S4 x C2"
    if n == 72 and not ab:
        return "S3 wr C2"
    if n == 120 and o[5] == 24:
        return "S5"
    if n == 720:
        return "S6"
    return f"group of order {n}"


def cycle_notation(p):
    cyc = [c for c in cycles_of(p) if len(c) > 1]
    if not cyc:
        return "id"
    return "".join("(" + " ".join(str(x + 1) for x in c) + ")" for c in cyc)


def symmetry_data(v, E):
    from pynauty import autgrp
    gens, g1, g2, orbits, _ = autgrp(nauty_graph(v, E))
    gens_v = [tuple(g[:v]) for g in gens]
    elements = group_closure(gens_v, v)
    size = int(round(g1 * 10 ** g2))
    assert size == len(elements), "group order mismatch"
    classes = {}
    for x in range(v):
        classes.setdefault(orbits[x], []).append(x)
    orb = sorted(classes.values())
    inv = [c[0] for c in orb if len(c) == 1]
    return {"order": len(elements), "name": group_name(elements, gens_v),
            "gens": [cycle_notation(g) for g in gens_v if g != tuple(range(v))],
            "orbits": orb, "k": len(orb), "inv": inv}


# ======================================================================
# Renormalized ultrametric (Definition def:renorm)
# ======================================================================
def component_tree(v, E):
    """Components of level >= 1, children, local host diameters, valuation."""
    edges = sorted(set(E), key=lambda e: (popcount(e), e))
    rank = {}
    for e in edges:                                   # Lemma lem:grading
        below = [rank[f] for f in rank if (f & e) == f and f != e]
        rank[e] = 1 + max(below, default=0)
    comps = []
    for n in sorted(set(rank.values())):
        groups = []
        for e in (e for e in edges if rank[e] == n):
            mask, es, keep = e, [e], []
            for gm, ge in groups:
                if gm & e:
                    mask, es = mask | gm, es + ge
                else:
                    keep.append((gm, ge))
            groups = keep + [(mask, es)]
        comps += [{"mask": gm, "level": n, "edges": ge, "children": []}
                  for gm, ge in groups]
    by_mask = {c["mask"]: c for c in comps}

    def parent(mask):
        cands = [c for c in comps if c["mask"] != mask and (c["mask"] & mask) == mask]
        return min(cands, key=lambda c: popcount(c["mask"]))["mask"] if cands else None

    for x in range(v):
        by_mask[parent(1 << x)]["children"].append(1 << x)
    for c in comps:
        c["parent"] = parent(c["mask"])
        if c["parent"] is not None:
            by_mask[c["parent"]]["children"].append(c["mask"])
    for c in comps:
        ch, k = c["children"], len(c["children"])
        adj = [set() for _ in range(k)]
        for e in c["edges"]:
            inside = [i for i, b in enumerate(ch) if (e & b) == b]
            for i in inside:
                adj[i].update(j for j in inside if j != i)
        diam = 0
        for i in range(k):
            dist, queue = {i: 0}, [i]
            for a in queue:
                for b in adj[a]:
                    if b not in dist:
                        dist[b] = dist[a] + 1
                        queue.append(b)
            diam = max(diam, max(dist.values()))
        c["diam"], c["s"] = diam, Fraction(1, diam + 1)
    for c in sorted(comps, key=lambda c: -popcount(c["mask"])):
        c["phi"] = c["s"] * (by_mask[c["parent"]]["phi"] if c["parent"] else 1)
    return {"v": v, "comps": comps, "by_mask": by_mask,
            "root": max(comps, key=lambda c: popcount(c["mask"]))["mask"]}


def renorm_matrix(tree):
    v, comps = tree["v"], tree["comps"]
    M = [[Fraction(0)] * v for _ in range(v)]
    for x in range(v):
        for y in range(x + 1, v):
            pair = (1 << x) | (1 << y)
            a = min((c for c in comps if (c["mask"] & pair) == pair),
                    key=lambda c: popcount(c["mask"]))
            M[x][y] = M[y][x] = a["phi"]
    return M


def dendro_key(tree, mask):
    """Canonical form of the weighted dendrogram (Proposition prop:metricclasses)."""
    if popcount(mask) == 1:
        return "x"
    c = tree["by_mask"][mask]
    return f"({c['diam']}:" + ",".join(sorted(dendro_key(tree, k)
                                               for k in c["children"])) + ")"


def metric_data(v, E):
    if v == 1:
        return {"key": "x", "matrix": [[Fraction(0)]], "skeleton": "empty",
                "D": "--"}
    tree = component_tree(v, E)
    sk = sorted(tree["comps"], key=lambda c: (popcount(c["mask"]), c["mask"]))
    return {"key": dendro_key(tree, tree["root"]), "matrix": renorm_matrix(tree),
            "skeleton": ", ".join(set_str(c["mask"]) for c in sk),
            "D": "(" + ",".join(str(c["diam"]) for c in sk) + ")"}


# ======================================================================
# Output
# ======================================================================
def set_str(mask):
    return "".join(str(x + 1) for x in bits(mask))


def edges_str(v, E):
    if v == 1:
        return "{1}"
    es = sorted(E, key=lambda e: (popcount(e), tuple(bits(e))))
    return "{" + ", ".join(set_str(e) for e in es) + "}"


def matrix_str(M, indent="      "):
    cells = [[str(x) for x in row] for row in M]
    w = max(len(c) for row in cells for c in row)
    return "\n".join(indent + "[ " + "  ".join(c.rjust(w) for c in row) + " ]"
                     for row in cells)


def run_details(fht, vdet, do_list, do_sym, do_dist):
    catalogue, gid = [], 0
    for v in range(1, vdet + 1):
        for c in fht[v]:
            gid += 1
            c["id"], c["v"] = gid, v
            catalogue.append(c)
    if do_sym:
        for c in catalogue:
            c["sym"] = symmetry_data(c["v"], c["edges"])
    classes = {}
    if do_dist:
        for c in catalogue:
            c["metric"] = metric_data(c["v"], c["edges"])
            key = (c["v"], c["metric"]["key"])
            if key not in classes:
                classes[key] = {"id": len(classes), "members": [], "rep": c}
            classes[key]["members"].append(c["id"])
            c["mclass"] = classes[key]["id"]

    if do_list or do_sym:
        log("\n" + "=" * 70 + "\n STRUCTURES\n" + "=" * 70)
        for c in catalogue:
            line = f"F_{c['id']:<5} v={c['v']}  N={c['depth']}  E={edges_str(c['v'], c['edges'])}"
            if do_dist:
                line += f"   [metric class M_{c['mclass']}]"
            log(line)
            if do_sym:
                s = c["sym"]
                orb = ", ".join("{" + ",".join(str(x + 1) for x in o) + "}"
                                for o in s["orbits"])
                inv = "{" + ",".join(str(x + 1) for x in s["inv"]) + "}"
                log(f"      Aut = {s['name']} (order {s['order']})  "
                    f"generators: {', '.join(s['gens']) or 'id'}")
                log(f"      orbits: {orb}   k = {s['k']}   Inv = {inv}")
        if do_sym:
            log("\n Symmetry classes by v:")
            for v in range(1, vdet + 1):
                cnt = Counter(c["sym"]["name"] for c in catalogue if c["v"] == v)
                log(f"   v={v}: " + ", ".join(f"{n}: {k}" for n, k in sorted(cnt.items())))

    if do_dist:
        log("\n" + "=" * 70 + "\n METRIC CLASSES (renormalized ultrametric)\n" + "=" * 70)
        for (v, _), cl in classes.items():
            r = cl["rep"]["metric"]
            log(f"\nM_{cl['id']}  v={v}  skeleton: {r['skeleton']}  D={r['D']}")
            log(f"   FHTs: {', '.join('F_' + str(i) for i in cl['members'])}")
            log(f"   distance matrix (labels of F_{cl['rep']['id']}):")
            log(matrix_str(r["matrix"]))
        log("\n Number of metric classes by v: " + ", ".join(
            f"v={v}: {sum(1 for (w, _) in classes if w == v)}"
            for v in range(1, vdet + 1)))


def run_checks():
    ok = True
    log("\n[check] Lemma lem:singletop against Burnside, u <= 6 "
        "(u = 6 may take a few minutes)")
    for u in range(1, 7):
        g = host_orbits((1,), (u,), use_lemmas=False)
        good = (g == A006602[u])
        ok &= good
        log(f"   u={u}: Burnside={g}  A006602={A006602[u]}  {'OK' if good else 'FAIL'}")
    log("[check] Lemma lem:twotop against Burnside, u <= 4")
    for u in range(0, 5):
        for topc in ((1, 1), (2,)):
            for lowc in partitions(u):
                a, b = fix_twotop(u, topc, lowc), fix_generic(u + 2, topc, lowc)
                ok &= (a == b)
                if a != b:
                    log(f"   FAIL u={u} top={topc} low={lowc}: lemma={a} Burnside={b}")
        log(f"   u={u}: done")
    log(f"[check] {'all checks passed' if ok else 'SOME CHECKS FAILED'}")
    return ok


# ======================================================================
# User interface
# ======================================================================
def ask_int(prompt, lo, hi, default):
    while True:
        s = input(f"{prompt} [{lo}-{hi}, default {default}]: ").strip()
        if not s:
            return default
        try:
            x = int(s)
        except ValueError:
            print("  Please enter an integer.")
            continue
        if lo <= x <= hi:
            return x
        print(f"  Refused: the value must lie between {lo} and {hi}.")


def ask_yes(prompt, default=False):
    d = "Y/n" if default else "y/N"
    while True:
        s = input(f"{prompt} [{d}]: ").strip().lower()
        if not s:
            return default
        if s in ("y", "yes", "o", "oui"):
            return True
        if s in ("n", "no", "non"):
            return False


def interactive():
    log("Limits: brute force v <= %d, recurrence v <= %d, details v <= %d."
        % (LIMIT_BRUTE, LIMIT_RECURRENCE, LIMIT_DETAILS))
    log("Methods: 1 = brute force, 2 = recurrence, 3 = both")
    choice = ask_int("Method", 1, 3, 2)
    methods = {1: {"brute"}, 2: {"recurrence"}, 3: {"brute", "recurrence"}}[choice]
    hi = max(LIMIT_BRUTE if "brute" in methods else 0,
             LIMIT_RECURRENCE if "recurrence" in methods else 0)
    vmax = ask_int("Maximal number of vertices v", 1, hi, min(5, hi))
    if vmax > LIMIT_DETAILS:
        log(f"  Note: listing, symmetries and distances are limited to v <= {LIMIT_DETAILS}.")
    return {"vmax": vmax, "methods": methods,
            "list": ask_yes("List the FHT structures?"),
            "sym": ask_yes("Compute and display the symmetries?"),
            "dist": ask_yes("Compute the renormalized distance matrices?"),
            "check": ask_yes("Check the lemmas against Burnside?"),
            "out": "fht_results.json"}


def parse_args():
    ap = argparse.ArgumentParser(description="FHT v4.0: enumeration, symmetries, metrics")
    ap.add_argument("--vmax", type=int, default=0, help="maximal number of vertices")
    ap.add_argument("--methods", default="recurrence",
                    help="brute, recurrence, or brute,recurrence")
    ap.add_argument("--list", action="store_true", help="list the structures")
    ap.add_argument("--sym", action="store_true", help="symmetries")
    ap.add_argument("--dist", action="store_true", help="renormalized distance matrices")
    ap.add_argument("--check", action="store_true", help="check the lemmas")
    ap.add_argument("--out", default="fht_results.json", help="JSON file for the counts")
    a = ap.parse_args()
    methods = {m.strip() for m in a.methods.split(",") if m.strip()}
    if not methods <= {"brute", "recurrence"}:
        sys.exit("Unknown method. Use brute, recurrence or brute,recurrence.")
    return {"vmax": a.vmax, "methods": methods, "list": a.list, "sym": a.sym,
            "dist": a.dist, "check": a.check, "out": a.out}


def main():
    log("=" * 70 + "\n FRACTAL HYPER-TREES v4.0\n" + "=" * 70)
    cfg = interactive() if len(sys.argv) == 1 else parse_args()
    details = cfg["list"] or cfg["sym"] or cfg["dist"]

    if cfg["check"]:
        run_checks()
    if cfg["vmax"] < 1:
        if not cfg["check"]:
            log("Nothing to do: give --vmax (and possibly --check).")
        return

    limits = {"brute": LIMIT_BRUTE, "recurrence": LIMIT_RECURRENCE}
    if cfg["vmax"] > max(limits[m] for m in cfg["methods"]):
        sys.exit(f"Refused: v = {cfg['vmax']} exceeds the limit of every chosen "
                 f"method ({', '.join(f'{m} <= {limits[m]}' for m in cfg['methods'])}).")
    v_brute = 0
    if "brute" in cfg["methods"]:
        v_brute = min(cfg["vmax"], LIMIT_BRUTE)
        if cfg["vmax"] > LIMIT_BRUTE:
            log(f"Warning: the brute force is limited to v <= {LIMIT_BRUTE}; "
                f"it will run up to {LIMIT_BRUTE}.")
    v_det = min(cfg["vmax"], LIMIT_DETAILS) if details else 0
    if details:
        if cfg["vmax"] > LIMIT_DETAILS:
            log(f"Warning: details are limited to v <= {LIMIT_DETAILS}.")
        if v_det > v_brute:
            log(f"Note: the details need the structures; the brute force will run "
                f"up to v = {v_det}.")
            v_brute = v_det
        if v_det > WARN_DETAILS:
            log("Warning: v = 6 gives 16296 structures; redirect the output to a file.")
    if v_brute == LIMIT_BRUTE:
        log("Warning: the brute force for v = 6 is long and needs several GB of RAM.")

    results = {}
    try:
        if "recurrence" in cfg["methods"]:
            results["recurrence"] = run_recurrence(cfg["vmax"])
        if v_brute:
            fht, results["brute"] = run_brute(v_brute)
            if details:
                run_details(fht, v_det, cfg["list"], cfg["sym"], cfg["dist"])
    except GuardError as err:
        sys.exit(f"Stopped by a guard-rail: {err}")

    if "recurrence" in results and "brute" in results:
        log("\n[comparison] brute force versus recurrence")
        for v in results["brute"]:
            a, b = results["brute"][v], results["recurrence"][v]
            same = a["f"] == b["f"] and a["by_depth"] == b["by_depth"]
            log(f"   v={v}: f={a['f']}  {'OK' if same else 'DIFFERENT'}")

    with open(cfg["out"], "w") as fh:
        json.dump({m: {str(v): {"S": str(r["S"]), "I": str(r["I"]), "f": str(r["f"]),
                                "by_depth": {str(d): str(n)
                                             for d, n in r["by_depth"].items()}}
                       for v, r in res.items()} for m, res in results.items()},
                  fh, indent=1)
    log(f"\nCounts written to {cfg['out']}.")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit("\nInterrupted.")
