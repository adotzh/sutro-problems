#!/usr/bin/env python3
"""
New optimization ideas from 73,602.

KEY INSIGHT: Can we reduce sA_cache reads by moving ii OUTSIDE bj?

Current loop: bi > bj > bk > ii > jj
  - sA_cache reads: Ti * nbk * Tj * nbi * nbj = 8*16*4*2*4 = 4096 (at cost 1)
  - Each A[bi*Ti+ii, bk] is loaded nbi*nbj*Tj = 2*4*4 = 32 times? No wait:
    For each (bi,bj): for each (bk,ii): copy sA_cache once per ii per bk
    Then read Tj times in inner jj loop.
    So loads: nbi*nbj*Ti*nbk = 2*4*8*16 = 1024 copies
    reads from sA_cache: nbi*nbj*Ti*nbk*Tj = 1024*4 = 4096

NEW LOOP: bi > ii > bj > bk > jj (ii moved outermost after bi)
  - For each (bi, ii): process all (bj, bk, jj) combos
  - sA_cache[ii, bk] stays same for ALL bj iterations:
    Load sA_cache once per (bi, ii, bk):  nbi * Ti * nbk = 2*8*16 = 256 copies
    Then read Tj * nbj times from sA_cache:  nbi * Ti * nbk * Tj * nbj = 256*4*4 = 4096 reads
  - BUT: sC can only hold one ROW of the output (Tj cells for fixed ii)
    sC_row: Tj = 4 cells at some address
    Need to save/restore sC_row for each (bj) transition — actually NO:
    For each (bi, ii, bj), we accumulate over bk, then write out:
      For each bk: copy sA_cache, A_at(bi*Ti+ii, bk); for jj: mul/add into sC_row
      After all bk: copy C[bi*Ti+ii, bj*Tj+jj] from sC_row
    This is exactly what we want! sA_cache is loaded nbi*Ti*nbk = 256 times total.
    But wait: for each (bi, ii, bj, bk): copy sA_cache still!
    Because for bj=0 bk=0: load A[bi*Ti+ii, 0] into sA_cache
    For bj=0 bk=1: load A[bi*Ti+ii, 1] into sA_cache
    ...
    For bj=1 bk=0: load A[bi*Ti+ii, 0] into sA_cache again!
    So loads = nbi * Ti * nbj * nbk = 2*8*4*16 = 1024 (same as current!)

  The key question: by moving ii outside bj, do A cells get loaded fewer times?
  A[bi*Ti+ii, bk] is loaded once per (bi, bj, bk, ii) regardless of loop order.
  So sA_cache LOADS = nbi * Ti * nbj * nbk = 2*8*4*16 = 1024 always.
  And sA_cache READS = loads * Tj = 4096 always (for Tj=4).

  NO BENEFIT from reordering! The sA cache reads are invariant.

ANOTHER IDEA: What if we combine sA_cache=1 with also caching sB element at addr 2?

Currently:
  - sA_cache at addr 1: each (bk,ii) produces 4 reads of addr 1
  - sB at addr 3-6: each sB(jj) produces 1 read per (ii,bk)
  - For each (bk,ii): read sA_cache 4 times (once per jj), read sB(jj) once each

If we also do "sB caching" per jj:
  - Move inner loop order to: bk > jj > ii
  - sB_cache at addr 2: for each jj: copy sB_cache, B_at(bk, bj*Tj+jj)
  - sA at addr 3-10 (Ti=8 cells): for each ii: read sA_cell(ii)
  - For each (bk,jj,ii): mul tmp, sA_cell(ii), sB_cache; add sC(ii,jj), tmp

  Analysis: sB_cache reads Ti * nbk * nbi * nbj = 8*16*2*4 = 1024 per outer jj loop?
  No: for each (bi,bj,bk,jj): Ti reads from sB_cache
  = nbi * nbj * nbk * Tj * Ti = 2*4*16*4*8 = 4096 (same!)

  sA scratchpad reads: Tj reads per (bi,bj,bk,ii) = nbi*nbj*nbk*Ti*Tj = 4096 (same!)

  The key savings: sA is now at addr 3-10 (cost 2-4) instead of... well sA_cache was at addr 1 (cost 1).
  sB_cache is now at addr 2 (cost 2) instead of addr 3-6 (cost 2-3).
  We're ADDING cost for sA and SLIGHTLY reducing cost for sB.
  WORSE overall.

REAL INSIGHT: Can we reduce tmp reads?

tmp reads = Ti * (nbk-1) * Tj * nbi * nbj = 3840 (bk=0 skips tmp)
What if we restructure to eliminate tmp entirely?

For bk=0: mul sC(ii,jj), sA, sB(jj)  -- direct init
For bk>0: mul tmp, sA, sB(jj); add sC(ii,jj), tmp

The add instruction: "add dest, src" = dest := dest + src (reads both dest and src)
Is there another way to accumulate?

What if we use "mul sC_new, sA, sB; add sC, sC_new" but sC_new is a separate cell?
That's the same as tmp.

What if we chain accumulation differently?
bk=0: mul sC, sA0, sB0
bk=1: need sC += sA1*sB0, i.e., tmp = sA1*sB0, sC += tmp. Unavoidable.

BUT WAIT: What if we use the sC cell as both accumulator AND temp for the LAST bk?
For the LAST bk (bk=nbk-1): we don't need the previous sC value after the add.
So for bk=nbk-1, we do: mul tmp, sA, sB; add sC, tmp (reads tmp once).
For an earlier bk: mul tmp, sA, sB; add sC, tmp.

No way to avoid tmp reads for bk>0.

IDEA: Use bk=0 special case more aggressively.
Currently: bk=0 does `mul sC(ii,jj), sA, sB(jj)` — 0 tmp reads.
This saves Ti*Tj*nbi*nbj = 32*8 = 256 tmp reads.

Could we do bk=0 AND bk=1 without tmp?
bk=0: mul sC, sA0, sB0  (sC = a0*b0)
bk=1: we need sC += a1*b1. Without tmp, we'd have to read sC to compute the new value.
  Option: add sC, sB1; mul sC, sA1 -- NO, this computes sA1*(sC + sB1) = sA1*(a0*b0 + b1), wrong!
  There's no IR sequence for sC += a1*b1 without a temp.

IDEA: Can we reduce nbk (number of bk blocks) by increasing Tk?
Tk=1: nbk=16, tmp reads = Ti*(nbk-1)*Tj*nbi*nbj = 8*15*4*2*4 = 3840
Tk=2: nbk=8, but sA_cache has 2 cells at addr 1-2, tmp at addr 3
  tmp reads = Ti*(nbk-1)*Tk*Tj*nbi*nbj = 8*7*2*4*2*4 = 3584
  But tmp is now at addr 3 (cost 2) vs addr 2 (cost 2): same.
  However sA_cache[0] at addr 1 (cost 1), sA_cache[1] at addr 2 (cost 2).
  Net: 3584*2 vs 3840*2. Better by 256 reads * cost 2 = 512 savings on tmp.
  But: sA_cache now has 2 cells. Cell[0] reads = Ti*nbk*(Tk*Tj)*nbi*nbj = 8*8*8*2*4 = 4096? No...
  For Tk=2: inner kk loop: for kk in range(Tk): read sA_cache[kk];
  Each sA_cache cell read Tj times per (bi,bj,bk,ii):
    nbi*nbj*Ti*nbk*Tj = 2*4*8*8*4 = 2048 reads from each sA_cache cell.
  sA_cache[0] at addr 1 (cost 1): 2048 reads = 2048
  sA_cache[1] at addr 2 (cost 2): 2048 reads = 4096
  vs current Tk=1: sA_cache[0] at addr 1: 4096 reads = 4096

  For Tk=2:
  sA_cache cost = 2048*1 + 2048*2 = 6144 vs 4096 for Tk=1. WORSE.
  Even though tmp_reads are fewer, sA_cache costs more.

OK let me try some truly new ideas:

1. DIRECT BULK ACCUMULATION: What if sC is at bulk addresses directly?
   - sC at addr 551+ (C output location): would read C at cost 24+. Much worse.

2. ZERO-COPY OUTPUT: Currently we copy sC → bulk_C at the end.
   Can we avoid this? No -- the output must be at specific addresses.

3. SHARED sA_cache AND sB_cache:
   For a 2x1 tile: Ti=2, Tj=1, nbi=8, nbj=16, nbk=16
   sA_cache=1 (1 cell), tmp=2, sB=3 (1 cell), sC@4-5 (2 cells at cost 2 each)
   For each (bi,bj,bk,ii):
     copy sA_cache, A_at(bi*2+ii, bk): 8*16*16*2 = 4096 loads
     copy sB, B_at(bk, bj+0): 8*16*16 = 2048 loads
   sA_cache reads: 2048*1*nbi*nbj = 8*16*2048? No...
   Ti=2, Tj=1, so Tj reads per (ii,bk): 1 read from sA_cache per ii per bk.
   sA_cache reads = Ti*nbk*Tj*nbi*nbj = 2*16*1*8*16 = 4096 at cost 1 = 4096
   tmp reads = Ti*(nbk-1)*Tj*nbi*nbj = 2*15*1*8*16 = 3840 at cost 2 = 7680
   sB reads = Ti*nbk*nbi*nbj = 2*16*8*16 = 4096 at cost 2 = 8192
   sC: 2 cells at addr 4-5 (cost 2 each), reads = nbk*nbi*nbj = 16*8*16 = 2048 each
   sC cost = 2*2048*2 = 8192
   Bulk A: 256 cells, each read nbi*nbj = 128 times at cost addr_39+
   ...
   This seems similar to what we had with Ti=8,Tj=4 in terms of the high-cost sC.

Let me try something completely different:
WHAT IF we have TWO sA_cache cells at addr 1 and 2, for TWO ii values simultaneously?
For Ti=8, Tj=4, process ii=0 and ii=1 together:
  copy sA_cache_0, A[bi*8+0, bk]  (addr 1)
  copy sA_cache_1, A[bi*8+1, bk]  (addr 2)
  for jj=0..3:
    mul tmp_0, sA_cache_0, sB(jj); add sC(0,jj), tmp_0
    mul tmp_1, sA_cache_1, sB(jj); add sC(1,jj), tmp_1

  We need TWO tmp cells! And they'd be at higher addresses (3,4)?
  tmp reads: same number, but at addr 3-4 (cost 2 each) vs current single tmp at addr 2 (cost 2).
  So same cost for tmp. But sA reads halved? No:
  sA_cache_0 reads = Tj per pair = 4, sA_cache_1 reads = Tj per pair = 4
  Total = 8 = Ti reads per bk (same as before).
  No benefit in terms of read counts.

TRULY DIFFERENT IDEA: Can we OMIT the sC scratchpad for the last (bi,bj) block?
No -- the IR must be fixed (not adaptive).

IDEA: What if we use a smaller sC and accumulate row by row?
For Ti=8, Tj=4: sC row = 4 cells (one ii row).
Layout: sA_cache=1, tmp=2, sB@3-6, sC_row@7-10 (4 cells, cost 3 each)
For each (bi,bj):
  Zero sC_full (in bulk? No...)
  For ii=0..7:
    Initialize sC_row to 0? But IR has no "zero" instruction.
    Actually for bk=0: mul sC_row(jj), sA_cache, sB(jj) -- initializes sC_row
    For bk>0: mul tmp, sA_cache, sB(jj); add sC_row(jj), tmp
    After all bk: copy C[bi*8+ii, bj*4+jj] from sC_row(jj)

  But we're initializing sC_row for EACH ii! Total sC_row cells = Tj = 4.
  sC_row reads per ii per (bi,bj): nbk reads (bk=0: mul into sC; bk>0: 1 read each; copy-out: 1)
  = nbk reads per ii per (bi,bj)
  Total: Ti * nbk * nbi * nbj = 8*16*2*4 = 1024 reads per sC_row cell
  vs current: 128 reads per sC cell (32 cells)
  Total sC reads: 4 cells * 1024 = 4096 vs 32 * 128 = 4096 (SAME!)

  But the sC_row cells are at addr 7-10 (cost 3), vs current sC at addr 7-38 (cost 3-7).
  sC_row cost = 4096 * 3 = 12,288
  vs current sC cost = 4096 * avg_cost where avg_cost > 3.

  Current sC cost: sum over 32 cells of (128 reads * cost(addr))
  = 128 * sum_{ii,jj} cost(7 + ii*4 + jj)
  = 128 * sum_{k=0}^{31} cost(7+k)
  = 128 * (cost(7)+cost(8)+...+cost(38))

  Let me compute: isqrt(6)+1=3, isqrt(7)+1=3, isqrt(8)+1=3, isqrt(9)+1=3, isqrt(10)+1=4
  cost(7)=3, cost(8)=3, cost(9)=3, cost(10)=4, cost(11)=4, ..., cost(16)=4, cost(17)=5, ...
  cost range 7-38: 3,3,3,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,7,7
  sum: 3*3 + 7*4 + 9*5 + 11*6 + 2*7 = 9+28+45+66+14 = 162
  sC cost current: 128 * 162 = 20,736

  sC_row cost: 4096 * cost(7-10) = 4096 * (3+3+3+4)/4 ≈ but unequal reads!
  Actually each sC_row cell gets ALL 4096 total reads equally: 1024 reads each.
  sC_row cost = 1024*(cost(7)+cost(8)+cost(9)+cost(10)) = 1024*(3+3+3+4) = 1024*13 = 13,312

  vs current sC: 20,736
  SAVINGS: 20,736 - 13,312 = 7,424!

  But we eliminated 32-4=28 sC cells, and the rest of the layout shifts.
  Scratch space now: sA_cache=1, tmp=2, sB@3-6, sC_row@7-10 → only 10 cells for scratch
  A_base = 11, A_base + 256 = 267, ... (much lower! bulk A and B at lower addresses)

  This is the ROW-CHUNKED idea! Let me implement it.
"""

import sys, math, itertools
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from matmul import score_16x16, _parse

RECORD = 110_487
BEST_SO_FAR = 73_602
N = 16


def addr_cost(addr: int) -> int:
    return math.isqrt(addr - 1) + 1


def cost_breakdown(ir: str) -> dict:
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
    tiers = {}
    for addr, count in reads.items():
        c = addr_cost(addr)
        t = tiers.setdefault(c, {"reads": 0, "cost": 0, "addrs": []})
        t["reads"] += count
        t["cost"] += count * c
        t["addrs"].append(addr)
    return tiers


def generate_row_chunked(Ti: int, Tj: int) -> str | None:
    """
    Row-chunked sC: only Tj sC cells at a time (one ii row), avoiding expensive high-address sC.

    Loop: bi > bj > ii > bk > jj
    Layout: sA_cache=1, tmp=2, sB@3..(Tj+2), sC_row@(Tj+3)..(Tj+2+Tj)

    Key: sC_row has only Tj cells (reused for each ii), at much lower addresses.
    Trade-off: A must be loaded nbi*nbj*Ti*nbk times (same), but sC addresses are cheaper.

    Wait - problem: if ii is INNER (bk > ii), we'd reload sA for each (bk, ii).
    If ii is OUTER (ii > bk), we can cache sA across bk... but sA changes with bk!

    For loop bi > bj > ii > bk > jj:
    - For each (ii, bk): load sA_cache = A[bi*Ti+ii, bk]
    - For each jj: mul/add into sC_row(jj)
    - After all bk: copy out sC_row → C row
    - Then move to ii+1 (sC_row is re-initialized at bk=0)

    sA_cache loads: nbi * nbj * Ti * nbk (same as before)
    sA_cache reads: nbi * nbj * Ti * nbk * Tj (same as before = 4096)

    sB loads: nbi * nbj * nbk * Tj (same as before: once per bk per jj)
    But wait: if bk is INNER to ii, do we reload sB for each ii?
    YES! For each (bi, bj, ii, bk, jj): copy sB(jj), B[bk, bj*Tj+jj]
    sB loads = nbi * nbj * Ti * nbk * Tj = 2*4*8*16*4 = 4096 (vs current 512!)
    MUCH WORSE.

    The issue: we can't have ii outermost AND avoid reloading sB.
    To avoid reloading sB: bk must be outermost relative to ii.

    SOLUTION: loop bi > bj > bk > ii > jj (current order), but use sC_row!
    Problem: with current bk > ii > jj loop, different ii share the same sC_row.
    We need to RESTORE sC_row when we move from ii=0 to ii=1 under the same bk.
    But sC_row holds the PARTIAL SUM for the current ii; when we switch to ii=1,
    sC_row must hold different values.

    So: we need SEPARATE sC cells for each ii, OR we need to save/restore sC_row.
    Saving/restoring sC_row to bulk would cost O(Tj * nbi * nbj * Ti * nbk) bulk reads. Terrible.

    CONCLUSION: sC_row approach requires ii OUTER to bk, which forces sB reload.
    Only useful if sB savings > sC savings, but sB savings = 0 (we reload more!).

    Let me TRY it anyway to get the actual score with the loop bi > bj > ii > bk > jj.
    """
    if N % Ti != 0 or N % Tj != 0:
        return None
    nbi = N // Ti
    nbj = N // Tj
    nbk = N  # Tk=1 always

    sA_cache = 1
    tmp_addr = 2
    sB_base = 3
    sC_row_base = sB_base + Tj  # = 3 + Tj

    # Layout: 1, 2, sB@3..(2+Tj), sC_row@(3+Tj)..(2+2*Tj)
    scratch_end = sC_row_base + Tj  # only Tj cells for sC_row

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

    # Loop: bi > bj > ii > bk > jj
    # sB is loaded for EVERY (ii, bk) -- expensive!
    for bi in range(nbi):
        for bj in range(nbj):
            for ii in range(Ti):
                for bk in range(nbk):
                    # Load sA_cache
                    lines.append(f"copy {sA_cache},{A_at(bi*Ti+ii, bk)}")
                    # Load sB (EXPENSIVE: done for every ii!)
                    for jj in range(Tj):
                        lines.append(f"copy {sB_base+jj},{B_at(bk, bj*Tj+jj)}")
                    # Accumulate
                    for jj in range(Tj):
                        is_first = (bk == 0)
                        if is_first:
                            lines.append(f"mul {sC_row_base+jj},{sA_cache},{sB_base+jj}")
                        else:
                            lines.append(f"mul {tmp_addr},{sA_cache},{sB_base+jj}")
                            lines.append(f"add {sC_row_base+jj},{tmp_addr}")
                # Write out sC_row for this ii
                for jj in range(Tj):
                    lines.append(f"copy {C_at(bi*Ti+ii, bj*Tj+jj)},{sC_row_base+jj}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


def generate_sA1_tmp2_v2(Ti: int, Tj: int) -> str | None:
    """
    Current best layout for reference: sA_cache=1, tmp=2, sB@3..(Tj+2), sC@(Tj+3)..
    Loop: bi > bj > bk > ii > jj
    """
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


def generate_sB_cached_per_jj(Ti: int, Tj: int) -> str | None:
    """
    Experimental: cache sB at addr 1, sA at addr 2-9 (Ti cells), sC at addr 10+
    Loop: bi > bj > bk > jj > ii

    For each (bk, jj): copy sB_cache = B[bk, bj*Tj+jj]; for ii: accumulate sC(ii,jj)
    sB_cache reads: Ti * nbk * Tj * nbi * nbj = 8*16*4*2*4 = 4096 at cost 1
    sA reads: Tj * nbk * nbi * nbj * Ti = same 4096 (but from addr 2-9, cost 2-4)

    This is worse since sA is now at higher addresses than sA_cache was.
    """
    if N % Ti != 0 or N % Tj != 0:
        return None
    nbi = N // Ti
    nbj = N // Tj
    nbk = N

    sB_cache = 1
    tmp_addr = 2
    sA_base = 3  # Ti cells: addr 3..(Ti+2)
    sC_base = 3 + Ti  # = Ti + 3

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
                # Load sA tile (one column of A tile)
                for ii in range(Ti):
                    lines.append(f"copy {sA_base+ii},{A_at(bi*Ti+ii, bk)}")
                # For each jj: cache sB, then accumulate over ii
                for jj in range(Tj):
                    lines.append(f"copy {sB_cache},{B_at(bk, bj*Tj+jj)}")
                    for ii in range(Ti):
                        is_first = (bk == 0)
                        if is_first:
                            lines.append(f"mul {sC_base+ii*Tj+jj},{sA_base+ii},{sB_cache}")
                        else:
                            lines.append(f"mul {tmp_addr},{sA_base+ii},{sB_cache}")
                            lines.append(f"add {sC_base+ii*Tj+jj},{tmp_addr}")
            for ii in range(Ti):
                for jj in range(Tj):
                    lines.append(f"copy {C_at(bi*Ti+ii, bj*Tj+jj)},{sC_base+ii*Tj+jj}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


def generate_dual_cache(Ti: int, Tj: int) -> str | None:
    """
    Both sA and sB cached: sA_cache=1, sB_cache=2, tmp=3, sC@4..
    Loop: bi > bj > bk > ii > jj (as usual)

    For each bk: load sB_row into individual cells (but we only have 1 sB_cache cell!)
    This only works if Tj=1 (one sB cell).

    For Tj=1: sB_cache=2, sA_cache=1, tmp=3, sC@4..(Ti+3)
    Loop: bi > bj > bk > ii:
      copy sA_cache, A[bi*Ti+ii, bk]  -- but this goes INSIDE bk, overwriting sA each time!
      Wait: for Tj=1, there's only ONE B element per (bk, bj).
      copy sB_cache, B[bk, bj*1+0]   -- load B
      for ii: copy sA_cache, A; mul sC(ii), sA_cache, sB_cache

    Wait: for Tj=1, mul sC(ii), sA_cache, sB_cache for bk=0,
          mul tmp, sA_cache, sB_cache; add sC(ii), tmp for bk>0.

    sA_cache at addr 1, sB_cache at addr 2, tmp at addr 3.
    sA reads: nbi*nbj*Ti*nbk*Tj = same 4096 but Tj=1, so = nbi*nbj*Ti*nbk = Ti*nbk*nbi*nbj
    For Ti=16, Tj=1: 16*16*1*16 = 4096 at cost 1.
    sB reads: nbi*nbj*Ti*nbk = 16*16*16*1 = 4096 at cost 2. (Tj=1, one read per ii)
    vs current (sB at addr 3-6, Tj=4): 1024 reads each = 4096 total, cost avg (2+2+3+3)/4*4096 = 2.5*4096
    So sB cost for dual_cache Tj=1: 4096*2 = 8192 vs 4096*2.5 = 10240. Saves 2048.

    But tmp at addr 3 vs addr 2: 3840 * (cost 2) vs (Ti=16,Tj=1: tmp reads = Ti*(nbk-1)*Tj*nbi*nbj = 16*15*1*1*16 = 3840)
    tmp at addr 3 (cost 2): 3840 * 2 = 7680. Same!

    sC: Ti cells at addr 4..Ti+3. For Ti=16: addr 4..19 (cost 2-5).
    vs current Ti=16,Tj=1 layout: sC at addr 4..19 (same!).

    So for Tj=1: sB_cache at addr 2 saves (2.5-2)*4096 = 2048 vs sB at addr 3.
    But tmp moves from 2 to 3 (same cost 2). Net: +2048 savings!

    But Ti=16, Tj=1 current score: 78,291 (much worse than Ti=8, Tj=4 at 73,602).
    Even with +2048 savings: 78,291 - 2048 = 76,243, still worse than 73,602.

    Can dual_cache work for Tj=2?
    sA_cache=1 (1 cell), sB_cache@2-3 (2 cells), tmp=4, sC@5..
    sB_cache reads: Ti*nbk*nbi*nbj = 8*16*2*8 = 2048 per cell, at cost 2-2: 2048*2+2048*2 = 8192
    vs current sB@3-4 (Tj=2): 2048 reads each at cost 2: 4096. Hmm wait, for Tj=2:
    current: sB@3-4, 1024 reads each at cost 2 total = 2048
    dual: sB_cache@2-3, 1024 reads each... same positions, same cost!
    tmp moves from addr 2 to addr 4 (cost 2): same cost.
    Net change: 0. No benefit.

    So dual_cache only helps for Tj=1 (slightly), but Tj=1 is already non-optimal.
    """
    # For completeness, let's try Tj=1 with sA_cache=1, sB_cache=2
    if Tj != 1 or N % Ti != 0:
        return None

    nbi = N // Ti
    nbj = N  # Tj=1
    nbk = N

    sA_cache = 1
    sB_cache = 2
    tmp_addr = 3
    sC_base = 4  # Ti cells

    scratch_end = sC_base + Ti
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
                lines.append(f"copy {sB_cache},{B_at(bk, bj)}")
                for ii in range(Ti):
                    lines.append(f"copy {sA_cache},{A_at(bi*Ti+ii, bk)}")
                    is_first = (bk == 0)
                    if is_first:
                        lines.append(f"mul {sC_base+ii},{sA_cache},{sB_cache}")
                    else:
                        lines.append(f"mul {tmp_addr},{sA_cache},{sB_cache}")
                        lines.append(f"add {sC_base+ii},{tmp_addr}")
            for ii in range(Ti):
                lines.append(f"copy {C_at(bi*Ti+ii, bj)},{sC_base+ii}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


def generate_sA_sC_reuse(Ti: int, Tj: int) -> str | None:
    """
    Idea: What if we DON'T use a separate sC scratchpad,
    and instead accumulate directly into bulk C?

    For the FIRST bk (bk=0): mul C_at(i,j), sA_cache, sB(jj) -- direct into bulk
    For bk>0: mul tmp, sA_cache, sB(jj); add C_at(i,j), tmp

    This reads C[i,j] at expensive bulk addresses (cost 24+) for all but bk=0.
    Total C reads: Ti*Tj*(nbk-1)*nbi*nbj = 32*15*8 = 3840 at cost ~26 avg = 99,840
    MUCH WORSE (current sC reads = 20,736).
    """
    return None  # Not worth implementing


if __name__ == "__main__":
    print(f"Record to beat: {RECORD:,}")
    print(f"Best so far:    {BEST_SO_FAR:,}")
    print()

    # Compare reference vs alternatives
    configs = [
        ("sA1_tmp2 Ti=8,Tj=4 (BEST)", lambda: generate_sA1_tmp2_v2(8, 4)),
        ("row_chunked Ti=8,Tj=4", lambda: generate_row_chunked(8, 4)),
        ("row_chunked Ti=8,Tj=8", lambda: generate_row_chunked(8, 8)),
        ("row_chunked Ti=16,Tj=4", lambda: generate_row_chunked(16, 4)),
        ("row_chunked Ti=4,Tj=4", lambda: generate_row_chunked(4, 4)),
        ("sB_cached_per_jj Ti=8,Tj=4", lambda: generate_sB_cached_per_jj(8, 4)),
        ("sB_cached_per_jj Ti=4,Tj=8", lambda: generate_sB_cached_per_jj(4, 8)),
        ("sB_cached_per_jj Ti=16,Tj=4", lambda: generate_sB_cached_per_jj(16, 4)),
        ("dual_cache Ti=16,Tj=1", lambda: generate_dual_cache(16, 1)),
        ("dual_cache Ti=8,Tj=1", lambda: generate_dual_cache(8, 1)),
    ]

    results = []
    for name, gen_fn in configs:
        ir = gen_fn()
        if ir is None:
            print(f"  {name}: SKIP")
            continue
        try:
            cost = score_16x16(ir)
            delta = RECORD - cost
            results.append((cost, name))
            marker = "  ***" if cost < BEST_SO_FAR else "  "
            print(f"{marker}{name}: cost={cost:,}  delta={delta:+,}")
        except ValueError as e:
            print(f"  {name}: ERROR - {e}")

    results.sort()
    print(f"\n=== Summary ===")
    for cost, name in results:
        print(f"  {cost:,}  {name}")

    # Save any new records
    for cost, name in results:
        if cost < BEST_SO_FAR:
            ir = None
            for n, gen_fn in configs:
                if n == name:
                    ir = gen_fn()
                    break
            if ir:
                out_path = Path(__file__).parent / "ir" / f"new_record_{cost}.ir"
                out_path.parent.mkdir(exist_ok=True)
                out_path.write_text(ir + "\n")
                print(f"\n*** SAVED: {out_path} ***")
