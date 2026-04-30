#!/usr/bin/env python3
"""
Explore: Full A-row caching to reduce bulk A reads from 4x to 1x.

IDEA: For each bi block (and each bk), cache ALL A[bi*Ti+ii, bk] for ii=0..Ti-1
in scratchpad. Then loop over ALL bj blocks reading from scratchpad.

Layout:
  sA_row: Ti*N cells (full A tile row) at addr 1..Ti*N
  But Ti=8, N=16: 128 cells at addr 1-128 (cost 1-12)
  sB: Tj=4 cells at addr 129-132 (cost 12)
  tmp: addr 133 (cost 12)
  sC: Ti*Tj=32 cells at addr 134-165 (cost 12-13)

  Then for each (bi, bj):
    For bk=0..N-1: load sB row; for ii: mul/add into sC using sA_row[ii*N+bk]
    copy out sC

A is loaded ONCE per bi block (not per bj). So A reads = nbi (2 for Ti=8).
But the sA_row addresses are high: 1..128, cost 1-12. Avg cost much higher.
And sB/tmp/sC at even higher addresses.

ANALYTICAL COMPARISON:
Current (Ti=8, Tj=4, sA1_tmp2):
  sA cost: 4096 * 1 = 4,096
  bulk_A cost: 256 * 4 * avg_cost(39..294) = 1024 * avg(~13) ≈ 13,328

Full A-row cache:
  sA_row reads: 4096 total (same), but at addr 1-128 (cost 1-12, avg 7)
  sA_row cost ≈ 4096 * 7 = 28,672 (MUCH WORSE)
  bulk_A reads: nbi * N^2 = 2 * 256 = 512 (vs current 1024)
  bulk_A cost ≈ 512 * 13 = 6,656 (savings: ~6,672)
  Net change: +28,672 - 13,328 - 6,672 = +8,672. WORSE.

The fundamental problem: sA_row cells are at medium-cost addresses (1-128 avg 7)
while current sA_cache is at cost 1. The extra scratchpad reads hurt more than
the savings on bulk A.

ALTERNATIVE: Cache only ONE ROW of A (Ti cells for one bk at a time).
That's what current sA_cache does (single cell at addr 1).

DIFFERENT APPROACH: Reduce tmp cost by eliminating tmp.

For the 2-operand 'add', what addresses can we use?
add sC(ii,jj), sA_cache means sC(ii,jj) = sC(ii,jj) + sA_cache.
This is NOT what we want (we want sC += sA*sB).

WAIT — NEW INSIGHT: What if we use 3-operand add?
add sC(ii,jj), sC(ii,jj), [precomputed product]
But we don't have a 3-operand add that reads 3 values... actually we DO:
'add dest, src1, src2' reads src1 AND src2, stores src1+src2.
So 'add sC, sC, product_addr' would read sC AND product_addr.
But product must be stored somewhere first.

There's no way around computing the product and storing it before adding to sC.

NEW IDEA: What if sC uses SUBTRACTION instead of addition to avoid expensive reads?
sC starts at some large value and we subtract from it? No, this doesn't make sense
for arbitrary matrix values.

NEW IDEA: Can we REORDER BULK DATA to place frequently-accessed elements at
lower addresses?

For A[i,k]: element accessed with frequency nbj = 4 (read 4 times per cell).
All 256 A cells have the same read frequency (4). Current placement: A@39-294.
This is already the lowest available range.

For B[k,j]: element accessed with frequency nbi = 2.
All 256 B cells have frequency 2. Current placement: B@295-550.
To put B at lower addresses than A, we'd need to swap.
But A is read more often (4x vs 2x), so A should be at lower addresses. ✓

WHAT IF we use a DIFFERENT Ti and Tj that reduces the address range?
E.g., with very large Ti*Tj, scratch grows but sC cells are more spread out.
With very small Ti*Tj, scratch is tiny but more bj iterations → higher bulk costs.

Let me do a FULL SWEEP with the analytical model to find the global optimum:
"""

import sys, math
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from matmul import score_16x16

RECORD = 110_487
BEST_SO_FAR = 73_602
N = 16


def addr_cost(addr: int) -> int:
    return math.isqrt(addr - 1) + 1


def generate_sA1_tmp2(Ti: int, Tj: int) -> str | None:
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


def generate_sA1_sB_first(Ti: int, Tj: int) -> str | None:
    """
    Alternative: put sB BEFORE sC in scratch. Already explored (this IS the current layout).
    Current: sA=1, tmp=2, sB@3-6, sC@7-38. Exactly this.
    """
    return generate_sA1_tmp2(Ti, Tj)  # same


def generate_transpose_inner(Ti: int, Tj: int) -> str | None:
    """
    Swap inner jj and ii: loop order bk > jj > ii.
    sB(jj) is loaded first, then for each ii: sA_cache, mul/add.
    This changes which items share sA_cache (now sB is "outer cache").
    But as analyzed: sA reads are same regardless of ii/jj order.
    """
    if N % Ti != 0 or N % Tj != 0:
        return None
    nbi = N // Ti
    nbj = N // Tj
    nbk = N

    sA_cache = 1
    tmp_addr = 2
    sB_base = 3  # Only 1 sB cell needed (cache jj)? No: we need Tj cells for efficiency.
    # Actually: if jj is inner and loads sB each time, we'd load sB more.
    # Keep Tj sB cells loaded for each bk. Then for each ii, read them.
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
                # Load sB first
                for jj in range(Tj):
                    lines.append(f"copy {sB_base+jj},{B_at(bk, bj*Tj+jj)}")
                # Inner loop: ii first, then jj (SAME AS CURRENT — no change)
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


def generate_no_sB_scratchpad(Ti: int, Tj: int) -> str | None:
    """
    What if we DON'T cache sB in scratchpad, but load B directly from bulk?
    For each (bk, ii, jj): load B from bulk, multiply by sA_cache, add to sC.

    This saves sB scratchpad (no Tj cells at addr 3-6), so sC starts at addr 3.
    But now B is read Ti times per (bi,bj,bk) instead of 1 time.
    B reads per cell: nbi * Ti (instead of nbi = 2). For Ti=8: 16 reads per B cell.
    bulk_B cost: 256 * 16 * avg_cost(295...) ≈ MUCH WORSE.

    Save: sB cost = 10,240.
    Extra: 14 * B_reads_cost ≈ 14 * 10,738 / 2 = 75,166. WORSE by 64,926.
    """
    return None  # Much worse, don't implement


def generate_full_loop_A_cache(Ti: int, Tj: int) -> str | None:
    """
    Cache ENTIRE A tile (Ti*N cells) to reduce A bulk reads to 1x.
    Layout: sA_tile@1..(Ti*N), tmp@(Ti*N+1), sB@(Ti*N+2)..(Ti*N+1+Tj), sC@...

    For each bi: load all A[bi*Ti+ii, k] for ii=0..Ti-1, k=0..N-1.
    Then for each bj, bk, ii, jj: multiply from cached sA.

    A read: nbi times per cell (vs current 4 times).
    sA_tile reads: Ti * nbk * Tj * nbi * nbj = same 4096 (same total sA reads).
    But at higher addresses: Ti*N = 128 cells at addr 1-128 (avg cost 7 vs 1).
    sA_tile cost ≈ 4096 * 7 = 28,672 (vs current 4,096). WORSE by 24,576.
    Bulk A savings: (4-1) * 1024 * avg_cost(A) ≈ 3 * 13,328/4 * 3 ≈ 9,996.
    Net: much worse.
    """
    return None


def generate_col_major_sC(Ti: int, Tj: int) -> str | None:
    """
    Store sC in column-major order: sC[ii,jj] at sC_base + jj*Ti + ii.
    This doesn't change addresses used, just mapping. Same cost.
    """
    return None  # Same cost, different mapping


def generate_sA_col_cache(Ti: int, Tj: int) -> str | None:
    """
    NEW IDEA: Cache a COLUMN of A (Ti elements for fixed bk) across bj blocks.
    For each (bi, bk): load sA_col[ii] = A[bi*Ti+ii, bk] for all ii.
    Then for each bj: load sB, for each ii,jj: multiply.

    Layout: sA_col@1..Ti, tmp@Ti+1, sB@Ti+2..(Ti+1+Tj), sC@...
    For Ti=8: sA_col@1-8, tmp@9, sB@10-13, sC@14-45.

    A reads: A[bi*Ti+ii, bk] loaded nbi * nbk = 2*16 = 32 times per cell (not nbj*nbi*nbk = 128 or 4).
    Wait: loaded once per (bi, bk, ii). Over all: nbi * nbk * Ti = 2*16*8 = 256 cells × 1 load each.
    But A has 256 cells. Each loaded ONCE! (vs current 4 times). HUGE SAVINGS!

    But sA_col reads: for each bj: for each ii,jj: read sA_col[ii].
    = nbi * nbk * nbj * Ti * Tj = 2*16*4*8*4 = 4096 reads.
    At addr 1-8 (cost 1-3). Cost ≈ 4096 * avg_cost(1-8) = 4096 * 2.25 ≈ 9,216.
    vs current sA_cache: 4096 * 1 = 4,096. WORSE by 5,120.

    Bulk A savings: (4-1) * bulk_A_cost = 3 * 13,328/4 * 4 = 3 * 13,328... wait:
    Current: A read 4 times per cell = 1024 total reads, cost 13,328.
    New: A read 1 time per cell = 256 total reads, cost 13,328/4 ≈ 3,332.
    Savings: ~9,996.

    But sA_col extra cost vs sA_cache: 9,216 - 4,096 = 5,120 extra.
    tmp moves from addr 2 to addr 9 (Ti+1=9): cost 3 vs 2.
    tmp extra cost: 3840 * (3-2) = 3,840.
    sB moves from addr 3-6 to addr 10-13: cost 4 each vs 2-3.
    sB extra cost: 1024 * (4+4+4+4 - 2-2-3-3) = 1024 * 6 = 6,144.
    sC moves from addr 7-38 to addr 14-45: slightly higher addresses.
    sC extra cost: let me compute.

    Current sC@7-38: cost 3-7.
    New sC@14-45: cost 4-7.
    sC extra cost: 128 * sum_{k=14}^{45} cost(k) - 128 * sum_{k=7}^{38} cost(k)
    sum(7-38): 3*3+7*4+9*5+11*6+2*7 = 9+28+45+66+14 = 162
    sum(14-45): cost(14)=4,...,cost(16)=4(3 cells), cost(17)=5,...,cost(25)=5(9),
               cost(26)=6,...,cost(36)=6(11), cost(37)=7,...,cost(45)=7(9)
    = 3*4+9*5+11*6+9*7 = 12+45+66+63 = 186
    sC extra cost: 128*(186-162) = 128*24 = 3,072.

    TOTAL EXTRA: 5,120 (sA_col) + 3,840 (tmp) + 6,144 (sB) + 3,072 (sC) = 18,176.
    SAVINGS: 9,996 (bulk_A).
    NET: +18,176 - 9,996 = +8,180. WORSE.

    So the column cache approach is WORSE. The per-element sA_cache at addr 1 is better
    because it keeps the most-read item at the cheapest address.

    HOWEVER: what if Ti is LARGER? For Ti=16: sA_col@1-16 (cost 1-4 avg ~3), more spread.
    Even worse.

    CONCLUSION: The single sA_cache cell at addr 1 is the best possible for A caching.
    """
    return None


if __name__ == "__main__":
    print(f"Record to beat: {RECORD:,}")
    print(f"Best so far:    {BEST_SO_FAR:,}")
    print()

    # Full analytical sweep to find any better (Ti, Tj) combinations
    print("=== Analytical sweep: all Ti, Tj divisors of 16 ===")
    results = []
    divisors = [d for d in range(1, N+1) if N % d == 0]
    for Ti in divisors:
        for Tj in divisors:
            nbi = N // Ti
            nbj = N // Tj
            nbk = N

            sB_base = 3
            sC_base = 3 + Tj
            scratch_end = sC_base + Ti * Tj
            A_base = scratch_end
            B_base = A_base + N * N
            C_base = B_base + N * N

            # Scratchpad cost
            sA_cost = Ti * nbk * Tj * nbi * nbj * addr_cost(1)   # 4096 * 1 always
            tmp_cost = Ti * (nbk - 1) * Tj * nbi * nbj * addr_cost(2)  # 3840 * 2 always
            sB_cost = sum(Ti * nbk * nbi * nbj * addr_cost(sB_base + jj) for jj in range(Tj))
            sC_cost = sum(nbk * nbi * nbj * addr_cost(sC_base + ii * Tj + jj)
                          for ii in range(Ti) for jj in range(Tj))

            # Bulk costs (CORRECTED: A read nbj times, B read nbi times)
            bulk_A_cost = nbj * sum(addr_cost(A_base + i) for i in range(N * N))
            bulk_B_cost = nbi * sum(addr_cost(B_base + i) for i in range(N * N))
            bulk_C_cost = sum(addr_cost(C_base + i) for i in range(N * N))

            total = sA_cost + tmp_cost + sB_cost + sC_cost + bulk_A_cost + bulk_B_cost + bulk_C_cost
            results.append((total, Ti, Tj))

    results.sort()
    print("Top 10:")
    for cost, Ti, Tj in results[:10]:
        print(f"  Ti={Ti:>2} Tj={Tj:>2}  analytical={cost:,}  delta={RECORD-cost:+,}")

    print()

    # Verify by actual scoring
    print("=== Actual scores for top configs ===")
    for cost, Ti, Tj in results[:5]:
        ir = generate_sA1_tmp2(Ti, Tj)
        if ir:
            actual = score_16x16(ir)
            print(f"  Ti={Ti:>2} Tj={Tj:>2}: analytical={cost:,}  actual={actual:,}  "
                  f"(diff={actual-cost:+})")

    # Breakdown for Ti=8, Tj=4
    print()
    print("=== Breakdown for Ti=8, Tj=4 ===")
    Ti, Tj = 8, 4
    nbi = N // Ti
    nbj = N // Tj
    nbk = N

    sB_base = 3
    sC_base = 3 + Tj
    scratch_end = sC_base + Ti * Tj
    A_base = scratch_end
    B_base = A_base + N * N
    C_base = B_base + N * N

    sA_cost = 4096
    tmp_cost = 3840 * 2
    sB_cost = sum(Ti * nbk * nbi * nbj * addr_cost(sB_base + jj) for jj in range(Tj))
    sC_cost = sum(nbk * nbi * nbj * addr_cost(sC_base + ii * Tj + jj)
                  for ii in range(Ti) for jj in range(Tj))
    bulk_A_cost = nbj * sum(addr_cost(A_base + i) for i in range(N * N))
    bulk_B_cost = nbi * sum(addr_cost(B_base + i) for i in range(N * N))
    bulk_C_cost = sum(addr_cost(C_base + i) for i in range(N * N))

    total = sA_cost + tmp_cost + sB_cost + sC_cost + bulk_A_cost + bulk_B_cost + bulk_C_cost
    print(f"  sA_cost:    {sA_cost:>8,}  ({sA_cost/total*100:.1f}%)")
    print(f"  tmp_cost:   {tmp_cost:>8,}  ({tmp_cost/total*100:.1f}%)")
    print(f"  sB_cost:    {sB_cost:>8,}  ({sB_cost/total*100:.1f}%)")
    print(f"  sC_cost:    {sC_cost:>8,}  ({sC_cost/total*100:.1f}%)")
    print(f"  bulk_A:     {bulk_A_cost:>8,}  ({bulk_A_cost/total*100:.1f}%)")
    print(f"  bulk_B:     {bulk_B_cost:>8,}  ({bulk_B_cost/total*100:.1f}%)")
    print(f"  bulk_C:     {bulk_C_cost:>8,}  ({bulk_C_cost/total*100:.1f}%)")
    print(f"  TOTAL:      {total:>8,}")
