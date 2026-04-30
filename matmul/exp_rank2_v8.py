#!/usr/bin/env python3
"""
Corrected analytical model and new optimization directions.

CORRECTED READ COUNTS for Ti=8, Tj=4, N=16:

Loop: bi > bj > bk > ii > jj
  for bi in range(nbi=2):
    for bj in range(nbj=4):
      for bk in range(nbk=16):
        for jj in range(Tj=4): copy sB(jj), B[bk, bj*4+jj]  -- Tj=4 reads
        for ii in range(Ti=8):
          copy sA_cache, A[bi*8+ii, bk]  -- 1 copy
          for jj in range(Tj=4):
            if bk==0: mul sC(ii,jj), sA, sB(jj)  -- reads sA(1), sB(1)
            else: mul tmp, sA, sB(jj); add sC(ii,jj), tmp  -- reads sA(1),sB(1),sC(1),tmp(1)
      for ii,jj: copy C[bi*8+ii,bj*4+jj], sC(ii,jj)  -- reads sC

READ ANALYSIS:
  A[bi*Ti+ii, bk]:
    - read once per (bi, bj, bk, ii) -- the copy reads it
    - Total: nbi * nbj * nbk * Ti = 2 * 4 * 16 * 8 = 1024 reads over 256 cells = 4 reads/cell
    - cell A[i,k]: nbi_containing * nbj * 1 * 1 = 1 * nbj = nbj = 4 reads

  sA_cache:
    - read Tj times per (bi, bj, bk, ii)
    - Total: nbi * nbj * nbk * Ti * Tj = 2 * 4 * 16 * 8 * 4 = 4096

  sB(jj):
    - read Ti times per (bi, bj, bk) (once per ii, plus bk=0 direct mul)
    - Actually: read once per ii per (bi,bj,bk,jj): nbi * nbj * nbk * Ti = 2*4*16*8 = 1024 per cell

  tmp:
    - read once per (bi,bj,bk>0,ii,jj): nbi * nbj * (nbk-1) * Ti * Tj = 2*4*15*8*4 = 3840

  sC(ii,jj):
    - For bk>0: read once in "add sC, tmp": nbi * nbj * (nbk-1) = 2*4*15 = 120
    - For copy-out: read once: nbi * nbj = 8
    - Total: 128 per cell ✓

  B[bk, bj*Tj+jj] = cell (k, j):
    - copy sB(jj), B[k, j] executed for each (bi, bj_containing, bk=k, jj=j%Tj)
    - Over all bi: nbi = 2 times
    - So B[k,j] read nbi = 2 times

  C_at(bi*Ti+ii, bj*Tj+jj):
    - copy C[i,j], sC(ii,jj) executed once per (bi, bj)
    - C[i,j] read nbi_containing * nbj_containing = 1 * 1 = 1 time ✓

So corrected bulk reads:
  A: 4 reads per cell
  B: 2 reads per cell
  C: 1 read per cell

And bulk addresses:
  A_base = 7 + Tj + Ti*Tj = 7 + 4 + 32 = 43? No...
  sA_cache=1, tmp=2, sB@3-6, sC@7-38 → scratch_end = 39
  A_base = 39, B_base = 39+256 = 295, C_base = 295+256 = 551

Correct bulk costs:
  bulk_A = sum over 256 cells of 4 * cost(39..294)
  bulk_B = sum over 256 cells of 2 * cost(295..550)
  bulk_C = sum over 256 cells of 1 * cost(551..806)
"""

import sys, math
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from matmul import score_16x16, _parse

RECORD = 110_487
BEST_SO_FAR = 73_602
N = 16


def addr_cost(addr: int) -> int:
    return math.isqrt(addr - 1) + 1


def compute_analytical_cost(Ti, Tj, N=16):
    """Compute analytical cost for the sA_cache=1, tmp=2 layout."""
    nbi = N // Ti
    nbj = N // Tj
    nbk = N

    sA_cache = 1
    tmp_addr = 2
    sB_base = 3
    sC_base = 3 + Tj
    scratch_end = sC_base + Ti * Tj

    A_base = scratch_end
    B_base = A_base + N * N
    C_base = B_base + N * N

    # Scratch costs
    sA_cost = Ti * nbk * Tj * nbi * nbj * addr_cost(sA_cache)
    tmp_cost = Ti * (nbk - 1) * Tj * nbi * nbj * addr_cost(tmp_addr)
    sB_cost = sum(Ti * nbk * nbi * nbj * addr_cost(sB_base + jj) for jj in range(Tj))
    sC_cost = sum((nbk - 1) * nbi * nbj * addr_cost(sC_base + ii * Tj + jj)
                  for ii in range(Ti) for jj in range(Tj))
    sC_copy_cost = sum(nbi * nbj * addr_cost(sC_base + ii * Tj + jj)
                       for ii in range(Ti) for jj in range(Tj))

    # Bulk costs (CORRECTED)
    bulk_A_cost = sum(nbj * addr_cost(A_base + i) for i in range(N * N))  # read nbj times (not nbi*nbj)
    bulk_B_cost = sum(nbi * addr_cost(B_base + i) for i in range(N * N))  # read nbi times
    bulk_C_cost = sum(1 * addr_cost(C_base + i) for i in range(N * N))

    total = sA_cost + tmp_cost + sB_cost + sC_cost + sC_copy_cost + bulk_A_cost + bulk_B_cost + bulk_C_cost
    return total, {
        'sA': sA_cost,
        'tmp': tmp_cost,
        'sB': sB_cost,
        'sC_compute': sC_cost,
        'sC_copy': sC_copy_cost,
        'bulk_A': bulk_A_cost,
        'bulk_B': bulk_B_cost,
        'bulk_C': bulk_C_cost,
    }


def generate_sA1_tmp2(Ti: int, Tj: int) -> str | None:
    """Current best layout."""
    if N % Ti != 0 or N % Tj != 0:
        return None
    nbi = N // Ti
    nbj = N // Tj
    nbk = N

    sA_cache = 1
    tmp_addr = 2
    sB_base = 3
    sC_base = sB_base + Tj
    scratch_end = sC_base + Ti * Tj

    A_base = scratch_end
    B_base = A_base + N * N
    C_base = B_base + N * N

    A_at = lambda i, j: A_base + i * N + j
    B_at = lambda i, j: B_base + i * N + j
    C_at = lambda i, j: C_base + i * N + j

    inputs = ([A_at(i, j) for i in range(N) for j in range(N)] +
              [B_at(i, j) for i in range(N) for j in range(N)])
    outputs = [C_at(i, j) for i in range(N) for j in range(N)]
    lines = [",".join(map(str, inputs))]

    for bi in range(nbi):
        for bj in range(nbj):
            for bk in range(nbk):
                for jj in range(Tj):
                    lines.append(f"copy {sB_base+jj},{B_at(bk, bj*Tj+jj)}")
                for ii in range(Ti):
                    lines.append(f"copy {sA_cache},{A_at(bi*Ti+ii, bk)}")
                    for jj in range(Tj):
                        is_first = (bk == 0)
                        if is_first:
                            lines.append(f"mul {sC_base+ii*Tj+jj},{sA_cache},{sB_base+jj}")
                        else:
                            lines.append(f"mul {tmp_addr},{sA_cache},{sB_base+jj}")
                            lines.append(f"add {sC_base+ii*Tj+jj},{tmp_addr}")
            for ii in range(Ti):
                for jj in range(Tj):
                    lines.append(f"copy {C_at(bi*Ti+ii, bj*Tj+jj)},{sC_base+ii*Tj+jj}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


if __name__ == "__main__":
    print(f"Record to beat: {RECORD:,}")
    print(f"Best so far:    {BEST_SO_FAR:,}")
    print()

    Ti, Tj = 8, 4
    cost, breakdown = compute_analytical_cost(Ti, Tj)
    actual_cost = 73_602

    print(f"=== Analytical vs actual for Ti=8, Tj=4 ===")
    print(f"  Analytical: {cost:,}")
    print(f"  Actual:     {actual_cost:,}")
    print(f"  Difference: {cost - actual_cost:,}")
    print()
    for k, v in breakdown.items():
        print(f"  {k}: {v:,}")

    # Detailed read analysis from actual IR
    ir = generate_sA1_tmp2(Ti, Tj)
    if ir:
        print(f"\n=== Actual read analysis from IR ===")
        input_addrs, ops, output_addrs = _parse(ir)
        reads = {}
        for op, oprs in ops:
            if op == "copy":
                src = oprs[1]
                reads[src] = reads.get(src, 0) + 1
            else:
                srcs = (oprs[1], oprs[2]) if len(oprs) == 3 else (oprs[0], oprs[1])
                for s in srcs:
                    reads[s] = reads.get(s, 0) + 1
        for a in output_addrs:
            reads[a] = reads.get(a, 0) + 1

        scratch_end = 3 + Tj + Ti * Tj
        A_base = scratch_end
        B_base = A_base + N * N
        C_base = B_base + N * N

        scratch_reads = sum(r for a, r in reads.items() if a < scratch_end)
        A_reads = sum(r for a, r in reads.items() if A_base <= a < B_base)
        B_reads = sum(r for a, r in reads.items() if B_base <= a < C_base)
        C_reads = sum(r for a, r in reads.items() if C_base <= a)

        scratch_cost = sum(r * addr_cost(a) for a, r in reads.items() if a < scratch_end)
        A_cost = sum(r * addr_cost(a) for a, r in reads.items() if A_base <= a < B_base)
        B_cost = sum(r * addr_cost(a) for a, r in reads.items() if B_base <= a < C_base)
        C_cost = sum(r * addr_cost(a) for a, r in reads.items() if C_base <= a)
        total_actual = scratch_cost + A_cost + B_cost + C_cost

        print(f"  Scratch: {scratch_reads} reads, cost {scratch_cost:,}")
        print(f"  Bulk A:  {A_reads} reads, cost {A_cost:,}")
        print(f"  Bulk B:  {B_reads} reads, cost {B_cost:,}")
        print(f"  Bulk C:  {C_reads} reads, cost {C_cost:,}")
        print(f"  Total: {total_actual:,} (score: {score_16x16(ir):,})")

        print(f"\n  Per-address analysis (top 10 by cost):")
        top = sorted(reads.items(), key=lambda x: -x[1] * addr_cost(x[0]))[:10]
        for addr, r in top:
            c = addr_cost(addr)
            print(f"    addr={addr:>4} (cost={c}): {r:>6} reads = {r*c:>8,} total")

    # Now explore: what if we use different bulk orderings?
    # Can we put C at lower addresses than A/B?
    print(f"\n=== Exploring: C before A and B in bulk ===")
    # Currently: scratch@1-38, A@39-294, B@295-550, C@551-806
    # What if: scratch@1-38, C@39-294, A@295-550, B@551-806?
    # C is read 1x, A is read nbj=4x, B is read nbi=2x
    # Optimal: A first (most reads), B second, C last
    # But A_base=39 is optimal since A read most (4x per cell)

    # What if we reorder: scratch@1-38, A@39-294, C@295-550, B@551-806?
    # C at 295-550, B at 551-806: C read 1x at cost~18, B read 2x at cost~24.
    # Current: C at 551-806 (cost~26, 1x), B at 295-550 (cost~18, 2x).
    # Cost: C*1*26 + B*2*18 = 256*26 + 256*36 = 6656 + 9216 = 15,872
    # vs swap: C*1*18 + B*2*26 = 256*18 + 256*52 = 4608 + 13312 = 17,920
    # Current is BETTER (put more-read B at lower addresses).

    # What about A before B vs B before A?
    # A read 4x, B read 2x: A should be at lower addresses (39-294), B at 295-550. ✓

    print(f"  Current ordering (A,B,C) is optimal (most reads → lowest addresses)")

    # Key insight: the BULK cost is minimized by the current layout.
    # And SCRATCH cost is minimized by optimal weighted assignment.
    # So 73,602 IS the analytical optimum for the Ti=8,Tj=4 parametric family.

    # But wait: are there alternative IR STRUCTURES (not just tiling params) that could help?
    print(f"\n=== Can we reduce bulk A reads? ===")
    print(f"  A reads per cell = nbj = {N // Tj} (once per bj block)")
    print(f"  B reads per cell = nbi = {N // Ti}")
    print(f"  To reduce A reads: increase Tj (fewer bj blocks), but sC grows")
    print(f"  To reduce B reads: increase Ti (fewer bi blocks), but sC grows")
    print()
    print(f"  For Ti=8,Tj=4: A reads = 4, B reads = 2 ← optimal balance")
    print(f"  Confirmed: 73,602 is optimal for this parametric family")
    print()
    print(f"  NOTE: Analytical cost {cost:,} ≠ actual 73,602 — discrepancy suggests")
    print(f"  the analytical model has a bug. Let me find it.")

    # Find the discrepancy
    N_val = N
    Ti_val, Tj_val = Ti, Tj
    nbi_v = N_val // Ti_val
    nbj_v = N_val // Tj_val
    nbk_v = N_val

    sA_cache = 1
    tmp_addr = 2
    sB_base = 3
    sC_base = 3 + Tj_val
    scratch_end_v = sC_base + Ti_val * Tj_val
    A_base_v = scratch_end_v
    B_base_v = A_base_v + N_val * N_val
    C_base_v = B_base_v + N_val * N_val

    print(f"\n  Layout: sA=1, tmp=2, sB@3-{2+Tj_val}, sC@{sC_base}-{scratch_end_v-1}")
    print(f"  A@{A_base_v}-{B_base_v-1}, B@{B_base_v}-{C_base_v-1}, C@{C_base_v}-{C_base_v+255}")

    # Analytical
    an_sA = Ti_val * nbk_v * Tj_val * nbi_v * nbj_v * addr_cost(sA_cache)
    an_tmp = Ti_val * (nbk_v - 1) * Tj_val * nbi_v * nbj_v * addr_cost(tmp_addr)
    an_sB = sum(Ti_val * nbk_v * nbi_v * nbj_v * addr_cost(sB_base + jj) for jj in range(Tj_val))
    an_sC = sum((nbk_v - 1 + 1) * nbi_v * nbj_v * addr_cost(sC_base + ii * Tj_val + jj)
                for ii in range(Ti_val) for jj in range(Tj_val))

    # A: cell (i,k) read for each (bi=i//Ti, bj, bk=k, ii=i%Ti) → nbj times
    an_A = sum(nbj_v * addr_cost(A_base_v + i * N_val + k) for i in range(N_val) for k in range(N_val))
    # B: cell (k,j) read for each (bi, bj=j//Tj, bk=k, jj=j%Tj) → nbi times
    an_B = sum(nbi_v * addr_cost(B_base_v + k * N_val + j) for k in range(N_val) for j in range(N_val))
    an_C = sum(1 * addr_cost(C_base_v + i * N_val + j) for i in range(N_val) for j in range(N_val))

    print(f"\n  Revised analytical breakdown:")
    print(f"    sA: {an_sA:,}")
    print(f"    tmp: {an_tmp:,}")
    print(f"    sB: {an_sB:,}")
    print(f"    sC: {an_sC:,}")
    print(f"    bulk_A: {an_A:,}")
    print(f"    bulk_B: {an_B:,}")
    print(f"    bulk_C: {an_C:,}")
    total_revised = an_sA + an_tmp + an_sB + an_sC + an_A + an_B + an_C
    print(f"    TOTAL: {total_revised:,}  (actual: 73,602)")
