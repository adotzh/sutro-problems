#!/usr/bin/env python3
"""
Direction B: A-tile caching across bj (bi>bk>bj loop order).

In the standard bi>bj>bk loop, A[bi,bk] is reloaded for every bj (4x redundant
for the square T=4 case). With Tk=1, A[bi,bk] is just a single column of A,
loaded for every (bi, bj, bk) — that's nbj*nbk = Tj_blocks * N times.

Fix: restructure to bi>bk>bj so A is loaded once per (bi, bk).

For the Tk=1 case (Ti=8, Tj=4):
  Standard:  nbi*nbj*nbk loads of sA = 2*4*16 = 128 loads of 8 cells = 1024 bulk reads
  Cached:    nbi*nbk loads of sA = 2*16 = 32 loads of 8 cells = 256 bulk reads
  Savings:   768 bulk A reads

Cost: for bk>0, when we move to a new bj with the same bk, we need to reload sC
from the partial result stored in bulk C. This costs Ti*Tj * addr_cost(C_addr).

Let's compute if this is a win.

Also implement Direction B for the original T=4 square tiling, and for Ti=8, Tj=4, Tk=1.
"""

import sys, math
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from matmul import score_16x16, _parse

RECORD = 110_487
BEST_SO_FAR = 82_477
N = 16


def addr_cost(addr: int) -> int:
    return math.isqrt(addr - 1) + 1


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


def generate_tk1_acache(Ti: int, Tj: int, slot_order=(1, 0, 2)) -> str | None:
    """
    Tk=1, bi>bk>bj loop order (cache A-tile across bj iterations).

    Layout: sB first (smaller), then sA, then sC.
    For sC: needs to read partial sums from bulk C when bj > 0.

    Key change:
    - for bi in range(nbi):
        for bk in range(nbk):  # bk is middle loop
          load sA once
          for bj in range(nbj):  # bj is inner loop
            load sB
            for ii in range(Ti):
              for jj in range(Tj):
                sC[ii,jj] += sA[ii] * sB[jj]  (first bk: initialize)
            write back sC to bulk C (always, since next bj needs a fresh sC)

    Wait, the accumulation over bk means we can't write sC after each bk.
    We need to accumulate over ALL bk.

    Two strategies:
    A) bi>bk>bj: After inner bj loop, sC has partial sum for (bi, bj=0..nbj-1, bk).
       But sC only has one tile's worth, so for each (bi,bk), we do:
         for bj: load sB, compute outer product, ACCUMULATE INTO BULK C directly
       This requires reading bulk C for bj > 0 (to add to existing partial).

    B) bi>bk>bj with partial C in scratchpad (requires Ti*Tj cells per (bi,bj)):
       Not feasible unless we have multiple sC tiles.

    Let's implement Strategy A:
    - for bi in range(nbi):
        for bk in range(nbk):
          load sA (once per bk)
          for bj in range(nbj):
            load sB
            for ii, jj:
              if bk == 0:
                C_at[bi*Ti+ii, bj*Tj+jj] += sA[ii] * sB[jj]  (direct to bulk C)
              else:
                read C_at[bi*Ti+ii, bj*Tj+jj], add sA[ii]*sB[jj], write back

    This avoids the scratchpad sC entirely but requires reading bulk C N times.
    Cost of reading bulk C: for each (bi,bj,bk) with bk>0: Ti*Tj reads at bulk C cost.
    """
    if N % Ti != 0 or N % Tj != 0:
        return None

    nbi = N // Ti
    nbj = N // Tj
    nbk = N  # Tk=1

    tmp = 1

    # Scratchpad: just sA and sB (no sC)
    # slot_order for (sA=0, sB=1): which comes first
    # Using sB first (cheaper for small Tj):
    sB_base = 2
    sA_base = sB_base + Tj
    scratch_end = sA_base + Ti

    sA = lambda ii: sA_base + ii
    sB = lambda jj: sB_base + jj

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
        for bk in range(nbk):
            # Load column of A once (shared across all bj)
            for ii in range(Ti):
                lines.append(f"copy {sA(ii)},{A_at(bi*Ti+ii, bk)}")

            for bj in range(nbj):
                # Load row of B
                for jj in range(Tj):
                    lines.append(f"copy {sB(jj)},{B_at(bk, bj*Tj+jj)}")

                # Accumulate outer product directly into bulk C
                for ii in range(Ti):
                    for jj in range(Tj):
                        if bk == 0:
                            # Initialize bulk C directly
                            lines.append(f"mul {C_at(bi*Ti+ii, bj*Tj+jj)},{sA(ii)},{sB(jj)}")
                        else:
                            # mul tmp = sA*sB, then add to bulk C
                            lines.append(f"mul {tmp},{sA(ii)},{sB(jj)}")
                            lines.append(f"add {C_at(bi*Ti+ii, bj*Tj+jj)},{tmp}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


def generate_tk1_acache_with_sC(Ti: int, Tj: int) -> str | None:
    """
    A-caching with scratchpad sC.

    bi>bk>bj loop, but we keep sC in scratchpad.
    Problem: sC accumulates over bk, but for different bj, we need different sC.
    Solution: for each (bi, bj), we have a separate sC.

    But with Ti=8, Tj=4: each sC is 32 cells. For nbj=4 tiles: 128 cells — too expensive.

    Alternative: Only cache A for fixed bj (bi>bj>bk but pre-cache A):
    - for bi, bj:
        pre-load: A[bi*Ti.., :] into a larger scratchpad (Ti*N cells) -- too big
    Not feasible.

    Better: Only cache within a single bj at a time.
    bi > bj > bk with bA cached once:
    - for bi in range(nbi):
        for bj in range(nbj):
          # Cache A-column for all bk
          # A[bi*Ti+ii, bk] for all ii, bk: Ti*nbk cells = Ti*N cells too big
          for bk in range(nbk):
            ... standard

    We can't fit all A into scratch without it being too expensive.

    Let's instead try: cache A for a GROUP of bk iterations.
    If we group bk into chunks of size Gk, within each group A stays cached.

    For Ti=8, Tj=4, Tk=1, Gk=2:
    - sA: Ti * Gk = 16 cells (can hold 2 consecutive A columns)
    - sB: Tj * Gk = 8 cells
    - sC: Ti * Tj = 32 cells

    Actually the simplest fix: in the standard bi>bj>bk loop, preload the
    entire A-row (Ti*N cells) into a bigger scratchpad before the bj loop.

    But Ti*N = 128 cells — addresses up to 2+128+128+32 = ~290.
    Each address would cost ~17. Too expensive for reads.

    Let's just measure the direct-to-bulk-C approach:
    """
    return generate_tk1_acache(Ti, Tj)


def generate_tk1_hybrid(Ti: int, Tj: int, bk_group: int = 4) -> str | None:
    """
    Hybrid approach: group bk iterations together to amortize A reload.

    For bk groups of size Gk:
    - Load A strip: Ti * Gk cells into scratchpad (wider A tile)
    - Load B strip: Gk * Tj cells into scratchpad (taller B tile)
    - Accumulate Ti * Tj sC

    This is equivalent to using Tk=Gk in the original non-square tiling.
    Let's try Tk=2, Tk=4, etc.
    """
    Tk = bk_group
    if N % Ti != 0 or N % Tj != 0 or N % Tk != 0:
        return None

    nbi = N // Ti
    nbj = N // Tj
    nbk = N // Tk

    tmp = 1

    # sB first (smaller if Tj <= Ti)
    if Tj <= Ti:
        sB_base = 2
        sA_base = sB_base + Tk * Tj
        sC_base = sA_base + Ti * Tk
    else:
        sA_base = 2
        sB_base = sA_base + Ti * Tk
        sC_base = sB_base + Tk * Tj

    scratch_end = sC_base + Ti * Tj

    sA = lambda ii, kk: sA_base + ii * Tk + kk
    sB = lambda kk, jj: sB_base + kk * Tj + jj
    sC = lambda ii, jj: sC_base + ii * Tj + jj

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
                for ii in range(Ti):
                    for kk in range(Tk):
                        lines.append(f"copy {sA(ii,kk)},{A_at(bi*Ti+ii, bk*Tk+kk)}")
                for kk in range(Tk):
                    for jj in range(Tj):
                        lines.append(f"copy {sB(kk,jj)},{B_at(bk*Tk+kk, bj*Tj+jj)}")
                for ii in range(Ti):
                    for jj in range(Tj):
                        for kk in range(Tk):
                            is_first = (bk == 0) and (kk == 0)
                            if is_first:
                                lines.append(f"mul {sC(ii,jj)},{sA(ii,kk)},{sB(kk,jj)}")
                            else:
                                lines.append(f"mul {tmp},{sA(ii,kk)},{sB(kk,jj)}")
                                lines.append(f"add {sC(ii,jj)},{tmp}")
            for ii in range(Ti):
                for jj in range(Tj):
                    lines.append(f"copy {C_at(bi*Ti+ii, bj*Tj+jj)},{sC(ii,jj)}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


if __name__ == "__main__":
    print(f"Record to beat: {RECORD:,}")
    print(f"Best so far:    {BEST_SO_FAR:,}")
    print()

    # Strategy A: direct to bulk C (no sC scratchpad)
    print("=== Strategy A: Direct accumulation into bulk C (bi>bk>bj) ===")
    for Ti, Tj in [(8, 4), (4, 4), (16, 4), (8, 2), (4, 8)]:
        if N % Ti != 0 or N % Tj != 0:
            continue
        ir = generate_tk1_acache(Ti, Tj)
        if ir is None:
            continue
        try:
            cost = score_16x16(ir)
            delta = RECORD - cost
            best_marker = " *** BEATS BEST!" if cost < BEST_SO_FAR else ""
            print(f"  Ti={Ti:>2} Tj={Tj:>2}  cost={cost:>10,}  delta={delta:>+8,}{best_marker}")
        except ValueError as e:
            print(f"  Ti={Ti:>2} Tj={Tj:>2}  ERROR: {e}")

    print()

    # Strategy B: hybrid with larger Tk (group bk)
    print("=== Strategy B: Larger Tk groups (Ti=8, Tj=4 base) ===")
    for Ti in [8, 4, 16]:
        for Tj in [4, 2, 8]:
            for Tk in [1, 2, 4, 8]:
                if N % Ti != 0 or N % Tj != 0 or N % Tk != 0:
                    continue
                scratch = Ti * Tk + Tk * Tj + Ti * Tj
                ir = generate_tk1_hybrid(Ti, Tj, bk_group=Tk)
                if ir is None:
                    continue
                try:
                    cost = score_16x16(ir)
                    delta = RECORD - cost
                    best_marker = " *** BEATS BEST!" if cost < BEST_SO_FAR else ""
                    print(f"  Ti={Ti:>2} Tj={Tj:>2} Tk={Tk:>2}  scratch={scratch:>3}  "
                          f"cost={cost:>10,}  delta={delta:>+8,}{best_marker}")
                except ValueError as e:
                    print(f"  Ti={Ti:>2} Tj={Tj:>2} Tk={Tk:>2}  ERROR: {e}")

    print()

    # Cost breakdown for the best direct-acache approach
    print("=== Cost breakdown: Ti=8 Tj=4 direct-to-bulk-C ===")
    ir = generate_tk1_acache(8, 4)
    if ir:
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

    # Also compare vs best (Ti=8, Tj=4, Tk=1 scratchpad)
    print()
    print("=== Cost breakdown: Ti=8 Tj=4 Tk=1 (scratchpad, BEST) ===")
    from exp_nonsquare_v2 import generate_tk1_optimized
    ir = generate_tk1_optimized(8, 4, slot_order=(1, 0, 2))
    if ir:
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
