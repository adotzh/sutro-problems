#!/usr/bin/env python3
"""
Advanced optimizations to push below 82,477.

Key insight: In the Tk=1 approach with Ti=8, Tj=4:
  - sA reads: 8 cells * 256 reads = 2048 reads at cost 3-4 -> ~7,000 cost
  - sB reads: 4 cells * 256 reads = 1024 reads at cost 2-3 -> ~2,500 cost
  - sC reads: 32 cells * 120 reads = 3840 reads at cost 4-7 -> ~21,000 cost
  - tmp reads: 3840 reads at cost 1 -> 3,840 cost
  - Bulk reads: many at high cost -> ~47,000 cost

The bulk reads dominate! Each A cell read nbj=4 times, each B cell read nbi=2 times.

Bulk A cost: sum_{i,j} cost(A[i,j]) * nbj = 4 * sum over all A addrs
Bulk B cost: sum_{i,j} cost(B[i,j]) * nbi = 2 * sum over all B addrs

For Ti=8, Tj=4:
  A_base = 2+4+8+32 = 46, A_end = 46+256 = 301
  B_base = 302, B_end = 557
  C_base = 558, C_end = 813

The bulk costs depend on where we place A, B, C in memory.
Currently (for Ti=8, Tj=4, sB@2-5, sA@6-13, sC@14-45):
  scratch_end = 46, A@46-301, B@302-557, C@558-813

Can we rearrange to put B first (since it's read 2x with nbi=2)?

Wait: current approach puts A then B. But B is read more total times (nbi=2 across different bi).
Actually: A[i,j] is read nbj times (for j-direction blocking)
        B[i,j] is read nbi times (for i-direction blocking)

For Ti=8, Tj=4: nbj=4, nbi=2. So A is read MORE (4x) than B (2x)!
Put A at lower addresses for cheaper reads.

Current order: A before B — this is already optimal!

Can we reduce nbj? That requires larger Tj.
Ti=8, Tj=8: nbi=2, nbj=2. Both read 2x. But sC gets bigger (64 cells, addr 18-81, cost 5-10).
We computed this already: cost=87,358 -- worse.

Ti=8, Tj=16: nbi=2, nbj=1. B is read once! A is read once!
But sC = 128 cells — very expensive.
And nbi*nbj=2, so there are fewer (bi,bj) pairs to accumulate over.

Ti=16, Tj=16: nbi=1, nbj=1. Both read once!
sC = 256 cells — way too expensive.

Alternative: Can we precompute sums to reduce bulk reads?

Actually, let me think about something different:
The bulk B cost for Ti=8, Tj=4 is:
  nbi=2 reads of each B cell. B starts at addr ~302.
  sum of addr_cost(302..557) * 2 = ~22,000

What if we put B in LOWER addresses than A?
Current: scratch@2-45, A@46-301, B@302-557
New:     scratch@2-45, B@46-301, A@302-557
A is read nbj=4 times, B is read nbi=2 times.
By putting B first (lower addr): B reads get cheaper (read 2x), A reads get expensive (read 4x).
But A is read MORE times -- so keeping A first is better.

Let me try explicitly placing A, B, C in different orders.

Also: for the Tk=1 case, the inputs to the IR are:
  A at bulk_A locations, B at bulk_B locations
  But inputs can be declared at any locations!
  The input line says where A and B actually live.
  We can put B at lower addresses than A.

Wait, can we put B BEFORE A in the scratchpad?
Currently: inputs = [A_at(0,0), A_at(0,1), ..., A_at(N-1,N-1), B_at(0,0), ..., B_at(N-1,N-1)]
And the A_base, B_base, C_base determine where they go.

The cost depends on WHERE we read FROM (A_base, B_base, C_base).
Currently: A@46-301, B@302-557, C@558-813

For Ti=8, Tj=4, nbj=4, nbi=2:
- Each A cell is read nbj=4 times
- Each B cell is read nbi=2 times

Optimization: Put B at lower addresses than A!
New layout: scratch@2-45, B@46-301, A@302-557, C@558-813
A cells at 302-557 (cost ~18-24), read 4x each -> very expensive
B cells at 46-301 (cost 7-17), read 2x each -> cheaper

For Ti=4, Tj=4, nbj=4, nbi=4:
- Equal reads, so ordering A,B doesn't matter (symmetric).

For Ti=8, Tj=4: nbj > nbi, so A should be at lower cost.
Current order (A before B) is already correct.

Hmm. Let me think about other angles...

What if we use a completely non-contiguous memory layout?
Put the most frequently read A rows at the LOWEST available addresses?
The most frequently read A row is...all rows equally (each row is read 4 times).

Actually: A[bi*8..(bi+1)*8, j] for all bk (j=bk from 0 to 15), for all bj.
Each A cell is read 4 times (once per bj). All A cells equally.

What about REORDERING the bulk memory?
Instead of A in row-major order, store the column being loaded FIRST.
For Tk=1, each bk iteration loads A[:, bk] (column bk of A).
The order we access the columns is bk=0,1,...,15.

Currently A is stored row-major: A[0,0], A[0,1], ..., A[0,15], A[1,0], ...
Column bk accesses: A[0,bk], A[1,bk], ..., A[15,bk] -- strided reads.

If we store A in COLUMN-MAJOR order: A[0,0], A[1,0], ..., A[15,0], A[0,1], ...
Then column bk accesses A[bk*16..(bk+1)*16] -- contiguous!

But does reordering the bulk layout change costs?
The cost is ceil(sqrt(addr)), so it depends on WHERE the data lives, not when it's accessed.
In column-major: A[i,j] is at addr A_base + j*N + i
The same set of addresses is used, just accessed in different order.

What DOES matter: can we put A at lower addresses overall?
The answer is: we already start A right after the scratchpad.
The scratchpad ends at addr 46 for Ti=8, Tj=4.

Can we reduce the scratchpad size to push A to lower addresses?
Minimum scratchpad for Tk=1:
  - tmp: 1 cell
  - sA: Ti cells
  - sB: Tj cells
  - sC: Ti*Tj cells
  Total: 1 + Ti + Tj + Ti*Tj cells

For Ti=8, Tj=4: 1 + 8 + 4 + 32 = 45 cells, so A starts at addr 46. Can't reduce.

What if we eliminate sC and accumulate directly into A bulk?
No, we need to accumulate C values.

What if we use partial products without scratchpad?
The fundamental limit is: to avoid reading bulk repeatedly, we need scratchpad.

Let me try another approach: what if sA and sB are not in scratchpad but bulk memory
is laid out so that the FIRST COLUMN is at the lowest address?

Actually the whole point of scratchpad is that bulk reads are expensive.
The key trade-off is scratchpad size (determines address costs) vs number of bulk reads.

Let me try an asymmetric approach: cache B (read nbi=2 times) differently from A (read nbj=4 times).

NEW IDEA: For Tk=1 with Ti=8, Tj=4:
Instead of loading B fresh for each (bi,bj,bk) iteration,
could we cache a full "B column" across bi iterations?

B column at bk, bj: B[bk, bj*Tj:(bj+1)*Tj] -- 4 cells, loaded nbi=2 times.
If we cache B columns in a persistent scratchpad, each B column is loaded ONCE.

But to do this, we'd need to store all Tj*nbk = 4*16 = 64 B cells persistently.
Plus sA (8 cells) and sC (32 cells) -> total scratch 1+8+64+32=105 cells.
B would be at addr 2+8+32=42 to 2+8+32+64-1=105.
Cost of cached_B[bk,jj]: addr 42+bk*Tj+jj, cost ~7-11.

Instead of reading bulk B at cost 18-24, cached B at cost 7-11.
Savings: for each B read, save ~(18-10) = 8 per read.
Total B reads = Tj * nbi * nbj * nbk = 4 * 2 * 4 * 16 = 512 reads
Old cost: ~512 * 18 avg = 9,216
New cost: ~512 * 9 avg = 4,608

But loading cost: load all B once = Tj*nbj*nbk = 4*4*16 = 256 bulk B reads
plus storing from bulk to cached scratchpad = 256 reads from high addresses.
Net change: +256 * avg_bulk_B_cost - 512 * (avg_bulk_B_cost - avg_cached_B_cost)
            = 256 * 18 - 512 * 9 = 4,608 - 4,608 = 0?

Actually this isn't quite right. Let me think more carefully.

The key question: does caching the B-side help when it's read nbi=2 times?

Currently B cells are at addr ~302-557, cost ~18-24.
Each B cell read nbi=2 times: total 256*2 = 512 reads at cost ~18-24.

If we cache B in scratchpad at addr ~42-105 (cost 7-11):
- Load cost: 256 reads at cost 18-24 (loading from bulk B to cache B) = 4,608-6,144
- Usage: 512 reads at cost 7-11 instead of 18-24, saving ~7-17 per read = 3,584-8,704

Rough: savings ≈ 512*9 - 256*21 = 4,608 - 5,376 = -768 (WORSE!)

Because loading costs as much as the repeated access savings.

Let me try the analysis with exact numbers.

Actually: for B_base = 302 to 557:
sum_cost = sum(addr_cost(302+i) for i in range(256)) * nbi
= sum(addr_cost(302+i) for i in range(256)) * 2

sum(addr_cost(302+i) for i in range(256)) = ?

Also: if we make the scratchpad smaller (no sC) and put A/B at lower addr?
Without sC: no accumulation possible without reading bulk C.

OK I think we're at a local optimum with Ti=8, Tj=4, Tk=1.

Let me try a radically different approach: use NON-POWER-OF-2 splits.
N=16 only has divisors {1,2,4,8,16}, so no new options there.

But what about splitting the computation DIFFERENTLY?
"""

import sys, math, itertools
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from matmul import score_16x16, _parse

RECORD = 110_487
BEST_SO_FAR = 82_477
N = 16


def addr_cost(addr: int) -> int:
    return math.isqrt(addr - 1) + 1


def generate_tk1_B_cached(Ti: int, Tj: int) -> str | None:
    """
    Cache all of B in scratchpad (loaded once per bi block),
    then process all bj, bk with A from bulk.

    For Ti=8: nbi=2 bi blocks. For each bi block:
    1. Cache all of B (256 cells) into sB_full scratch
    2. For each bj, bk: load A column from bulk, outer product into sC tile

    This eliminates all B reads from bulk (except initial load).
    But sB_full is 256 cells = expensive!

    Instead: cache one BLOCK of B columns for each bj (not all of B).
    For a given bj: B[0..15, bj*Tj..(bj+1)*Tj] = N*Tj = 16*4 = 64 cells.
    Cache this into scratchpad, then process all bi, bk with this cached B.

    Loop: bj > bi > bk
    For each bj: cache B[:, bj*Tj:(bj+1)*Tj] -- 64 cells
    For each bi, bk: load A column, outer product into sC

    sB_full scratch: N*Tj = 64 cells
    sA: Ti cells
    sC: Ti*Tj cells

    Total scratch: 1 + 64 + Ti + Ti*Tj = 1 + 64 + 8 + 32 = 105 cells
    A_base = 106

    sB_full at addr 2..65 (cost 2-9 avg ~5)
    sA at addr 66..73 (cost 9)
    sC at addr 74..105 (cost 9-11)

    B full cache reads: 64 * nbi * nbk = 64 * 2 * 16 = 2048 reads at avg cost 5
    vs. original: 4 * nbi * nbj * nbk = 4 * 2 * 4 * 16 = 512 reads at avg cost 21

    Hmm, caching saves 512 - 2048 + (load savings)...
    Actually: original B bulk reads = 512 * avg_bulk_B = 512 * 21 = 10,752
    Cached approach:
    - Load B once per bj: nbj * N * Tj = 4 * 64 = 256 bulk reads (load cost)
    - Use cached B: nbi * nbj * nbk * Tj = 2 * 4 * 16 * 4 = 512 reads from cache at cost 5
    - But also load B from bulk for caching: 256 * avg_bulk_B_cost

    Wait, we still need to LOAD B from bulk into cache. That's unavoidable.
    The savings come from USING the cached version (cost 5) instead of bulk (cost 21).

    Reads FROM cache: nbi * nbk * N * Tj = ? No...

    Let me reconsider. Loop bj > bi > bk:
    - For each bj: load sB_full = B[0..15, bj*Tj..(bj+1)*Tj] into scratchpad
      This is 64 reads from bulk B. Done ONCE per bj (not per bi).
    - For each bi, bk: load A[:, bk] into sA (Ti reads from bulk A per (bi,bk))
    - For each ii, jj: sC[ii,jj] += sA[ii] * sB_full[bk, jj]

    The key: sB_full[bk, jj] is read at scratchpad cost (much less than bulk).
    Original: each B cell read nbi=2 times from bulk. With bj>bi>bk:
    Each B[bk, bj*Tj+jj] cell is read nbi=2 times from bulk WITH CACHING.
    But after caching: it's READ from the cache as sB_full[bk, jj].
    Total reads of sB_full cell[bk, jj] = nbi = 2 times (once per bi in the inner bi loop).

    Wait, that's the same number of reads, just at lower cost.
    But we ALSO have to READ it once from bulk to load into cache (per bj).
    So total reads = nbi (cached) + 1 (load) = 3 vs old 2 (bulk).

    Hmm: more reads overall, but cheaper. Net effect?
    Old: 2 reads from bulk B at avg cost 21 = 42 per cell
    New: 1 read from bulk B at avg cost 21 (load) + 2 reads from cache at avg cost 5 = 31 per cell

    Save: 11 per cell * 256 cells = 2,816 savings!

    But there's overhead:
    - sB_full is 64 cells, costs scratchpad space (pushes sA and sC to higher addresses)
    - The extra space makes sA more expensive and sC more expensive

    Let me implement and measure.
    """
    if N % Ti != 0 or N % Tj != 0:
        return None

    nbi = N // Ti
    nbj = N // Tj
    nbk = N  # Tk=1

    tmp = 1

    # Scratchpad layout for bj > bi > bk order with cached B:
    # sB_full: N*Tj cells (all rows for current bj column)
    # sA: Ti cells
    # sC: Ti*Tj cells
    sB_full_base = 2          # N*Tj = 64 cells
    sA_base = sB_full_base + N * Tj   # 66
    sC_base = sA_base + Ti              # 74 (for Ti=8)
    scratch_end = sC_base + Ti * Tj    # 106

    sB_full = lambda k, j: sB_full_base + k * Tj + j   # B[k, bj*Tj+j] cached
    sA_f = lambda ii: sA_base + ii
    sC_f = lambda ii, jj: sC_base + ii * Tj + jj

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

    for bj in range(nbj):
        # Cache all of B[:, bj*Tj:(bj+1)*Tj] -> sB_full
        for k in range(N):
            for j in range(Tj):
                lines.append(f"copy {sB_full(k,j)},{B_at(k, bj*Tj+j)}")

        for bi in range(nbi):
            for bk in range(nbk):
                # Load A column: A[bi*Ti..(bi+1)*Ti, bk] -> sA
                for ii in range(Ti):
                    lines.append(f"copy {sA_f(ii)},{A_at(bi*Ti+ii, bk)}")

                # Outer product: sC[ii,jj] += sA[ii] * sB_full[bk, jj]
                for ii in range(Ti):
                    for jj in range(Tj):
                        is_first = (bk == 0)
                        if is_first:
                            lines.append(f"mul {sC_f(ii,jj)},{sA_f(ii)},{sB_full(bk,jj)}")
                        else:
                            lines.append(f"mul {tmp},{sA_f(ii)},{sB_full(bk,jj)}")
                            lines.append(f"add {sC_f(ii,jj)},{tmp}")

            # Write sC to bulk C
            for ii in range(Ti):
                for jj in range(Tj):
                    lines.append(f"copy {C_at(bi*Ti+ii, bj*Tj+jj)},{sC_f(ii,jj)}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


def generate_tk1_A_cached(Ti: int, Tj: int) -> str | None:
    """
    Cache all of A[:, bk] for the current bi block in scratchpad.
    Loop order: bi > bk > bj

    For each bi:
      For each bk: load A column for ALL bi (Ti*N cells... too many)

    Actually: For each (bi, bk): load A column (Ti cells).
    This is already done in the standard approach!
    The issue is A column is loaded AGAIN for each bj.

    In bi > bj > bk order: A[bi*Ti:.., bk] is loaded nbj times per (bi,bk).
    In bi > bk > bj order: A[bi*Ti:.., bk] is loaded ONCE per (bi,bk)!
    But then sC[bi,bj] has partial results -- we need to store them in bulk.

    For bi > bk > bj order:
    - For each bk: compute sC_part[bi,bj][ii,jj] = sA[ii]*sB[bk,jj] once per bi
    - But sC needs to be Ti*Tj cells, shared across bk

    The issue: with bk as MIDDLE loop (bi > bk > bj):
    - For (bi, bk): load sA once
    - For each bj: load sB, compute partial outer product
    - sC for each (bi, bj) accumulates over bk

    Problem: when we're in the inner bj loop, sC for (bi, bj=0) is filled after bk=0 step,
    but we need to CONTINUE accumulating into it for bk=1, 2, ...
    So between bj iterations (for the same bk), sC must change.

    With Ti*Tj = 32 cells for sC (just one tile), we can only keep ONE (bi,bj) tile at a time.
    After processing bj=0 with bk=0, we write sC to somewhere, then process bj=1.
    Then for bk=1: we need to RELOAD sC for each (bi,bj) from where we stored it.

    If we store partial sC in bulk C: read back at high cost.
    This is Direction B from the problem statement.

    Let's implement and measure:
    - Loop: bi > bk > bj
    - For bk=0, bj=0: sC = outer product of first A column and first B row
    - After each bj: write sC to partial_C in bulk (separate from output C)
    - For bk>0: read partial_C[bi,bj] from bulk, add new outer product
    """
    if N % Ti != 0 or N % Tj != 0:
        return None

    nbi = N // Ti
    nbj = N // Tj
    nbk = N

    tmp = 1

    # Minimal scratchpad: sB first, sA next, sC last
    sB_base = 2
    sA_base = sB_base + Tj
    sC_base = sA_base + Ti
    scratch_end = sC_base + Ti * Tj

    sA_f = lambda ii: sA_base + ii
    sB_f = lambda jj: sB_base + jj
    sC_f = lambda ii, jj: sC_base + ii * Tj + jj

    A_base = scratch_end
    B_base = A_base + N * N
    C_base = B_base + N * N
    # We need a separate partial C region
    # Put partial C right after bulk C
    partC_base = C_base + N * N

    A_at = lambda i, j: A_base + i * N + j
    B_at = lambda i, j: B_base + i * N + j
    C_at = lambda i, j: C_base + i * N + j
    partC_at = lambda bi, bj, ii, jj: (partC_base +
                                        bi * nbj * Ti * Tj +
                                        bj * Ti * Tj +
                                        ii * Tj + jj)

    inputs = ([A_at(i, j) for i in range(N) for j in range(N)] +
              [B_at(i, j) for i in range(N) for j in range(N)])
    outputs = [C_at(i, j) for i in range(N) for j in range(N)]
    lines = [",".join(map(str, inputs))]

    for bi in range(nbi):
        for bk in range(nbk):
            # Load A column once (shared across all bj)
            for ii in range(Ti):
                lines.append(f"copy {sA_f(ii)},{A_at(bi*Ti+ii, bk)}")

            for bj in range(nbj):
                # Load B row for this (bk, bj)
                for jj in range(Tj):
                    lines.append(f"copy {sB_f(jj)},{B_at(bk, bj*Tj+jj)}")

                if bk == 0:
                    # Initialize sC with outer product
                    for ii in range(Ti):
                        for jj in range(Tj):
                            lines.append(f"mul {sC_f(ii,jj)},{sA_f(ii)},{sB_f(jj)}")
                else:
                    # Load partial C from "partial storage" and add
                    for ii in range(Ti):
                        for jj in range(Tj):
                            lines.append(f"copy {sC_f(ii,jj)},{partC_at(bi,bj,ii,jj)}")
                    for ii in range(Ti):
                        for jj in range(Tj):
                            lines.append(f"mul {tmp},{sA_f(ii)},{sB_f(jj)}")
                            lines.append(f"add {sC_f(ii,jj)},{tmp}")

                if bk < nbk - 1:
                    # Save partial sC to partial storage (will reload next bk)
                    for ii in range(Ti):
                        for jj in range(Tj):
                            lines.append(f"copy {partC_at(bi,bj,ii,jj)},{sC_f(ii,jj)}")
                else:
                    # Last bk: write to final C
                    for ii in range(Ti):
                        for jj in range(Tj):
                            lines.append(f"copy {C_at(bi*Ti+ii, bj*Tj+jj)},{sC_f(ii,jj)}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


def cost_breakdown(ir: str) -> dict:
    input_addrs, ops, output_addrs = _parse(ir)
    reads: dict[int, int] = {}
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
    tiers: dict[int, dict] = {}
    for addr, count in reads.items():
        c = addr_cost(addr)
        t = tiers.setdefault(c, {"reads": 0, "cost": 0, "addrs": []})
        t["reads"] += count
        t["cost"] += count * c
        t["addrs"].append(addr)
    return tiers


if __name__ == "__main__":
    print(f"Record to beat: {RECORD:,}")
    print(f"Best so far:    {BEST_SO_FAR:,}")
    print()

    print("=== B-cached approach (all B columns in scratchpad) ===")
    for Ti, Tj in [(8, 4), (4, 4), (4, 8), (8, 8)]:
        ir = generate_tk1_B_cached(Ti, Tj)
        if ir:
            try:
                cost = score_16x16(ir)
                delta = RECORD - cost
                marker = " *** BEATS BEST!" if cost < BEST_SO_FAR else ""
                print(f"  Ti={Ti:>2} Tj={Tj:>2}  cost={cost:>10,}  delta={delta:>+8,}{marker}")
            except ValueError as e:
                print(f"  Ti={Ti:>2} Tj={Tj:>2}  ERROR: {e}")

    print()
    print("=== A-cached (bi>bk>bj with partial C in bulk) ===")
    for Ti, Tj in [(8, 4), (4, 4), (4, 8), (16, 4), (8, 8)]:
        ir = generate_tk1_A_cached(Ti, Tj)
        if ir:
            try:
                cost = score_16x16(ir)
                delta = RECORD - cost
                marker = " *** BEATS BEST!" if cost < BEST_SO_FAR else ""
                print(f"  Ti={Ti:>2} Tj={Tj:>2}  cost={cost:>10,}  delta={delta:>+8,}{marker}")
            except ValueError as e:
                print(f"  Ti={Ti:>2} Tj={Tj:>2}  ERROR: {e}")

    print()
    print("=== B-cached breakdown for Ti=8, Tj=4 ===")
    ir = generate_tk1_B_cached(8, 4)
    if ir:
        try:
            cost = score_16x16(ir)
            print(f"Cost: {cost:,}  delta={RECORD-cost:+,}")
            tiers = cost_breakdown(ir)
            total = sum(t["cost"] for t in tiers.values())
            for c in sorted(tiers):
                t = tiers[c]
                addrs = sorted(t["addrs"])
                pct = t["cost"] / total * 100
                print(f"  cost={c:>2}  addrs={addrs[0]:>3}-{addrs[-1]:>3}  "
                      f"reads={t['reads']:>7,}  cost={t['cost']:>8,}  {pct:>5.1f}%")
        except ValueError as e:
            print(f"ERROR: {e}")
