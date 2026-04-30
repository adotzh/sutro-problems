#!/usr/bin/env python3
"""
Further optimization from 73,602.

Current best: Ti=8, Tj=4, sA_cache=1, tmp=2, cost=73,602
Layout: sA_cache=1, tmp=2, sB@3-6, sC@7-38, A@39-294, B@295-550, C@551-806

Now try:
1. Different sB base (sB at 3 is already optimal since 1,2 are used)
2. What if we put some sB cells interleaved with sC to lower average cost?
   No - sB and sC must be contiguous for the indexing.
3. What if sB and sC ordering is swapped? sC@3..Ti*Tj+2, sB@Ti*Tj+3..?
   For Ti=8,Tj=4: sC@3-34 (cost 2-6), sB@35-38 (cost 6). sB becomes expensive.
   sB reads 1024 each at cost 6 = 6,144 total.
   vs current sB@3-6 at cost 2-3 avg 2.75: 1024*2.75 = 2,816.
   Penalty: 3,328. Much worse.

4. What if we use sA_cache for EACH ROW of sC accumulation instead of per-cell?
   sA_cache stores A[ii, bk] for one ii at a time.
   sB has 4 cells for the current bk row of B.
   sC has 32 cells for the accumulated tile.

   Actually: What if we do OUTER LOOP DIFFERENTLY?
   Current: bi > bj > bk > ii > jj loop.

   Alternative: bi > bj > ii > bk > jj
   For each (bi, bj, ii): accumulate sC_row[jj] = sum_{bk} A[bi*Ti+ii, bk] * B[bk, bj*Tj+jj]
   sC_row: only Tj = 4 cells at addr 3..6 (cost 2-3)!
   sA_cache: 1 cell at addr 1 (cost 1)
   sB: could be eliminated (read directly from bulk B)
   tmp: 2 (cost 2)

   For each (ii, bk): copy sA_cache, A_at(i, bk); then for each jj: do the accumulation.
   If sB is eliminated (read B_bulk directly in mul):
     - sB reads: Ti * nbi * nbj * nbk * cost_B_bulk = 8 * 2 * 4 * 16 * avg_B = 8192 * 20 ≈ 163,840. TOO EXPENSIVE.

   Keep sB: sB at addr 3..6 (4 cells), load B[bk, bj*Tj:..] once per bk.
   sC_row at addr 7..10 (cost 4 each) -- 4 cells, re-initialized for each ii.
   tmp at addr 2.

   Loop: bi > bj > ii > bk > jj
   For each (bi, bj, ii):
     for bk:
       copy sA_cache(=1), A_at(bi*Ti+ii, bk)
       for jj: copy sB(jj), B_at(bk, bj*Tj+jj)  -- Load B EACH TIME for each ii too!
       OR load B once per bk at the start and use for all ii... but then bk is outer.

   The issue: if ii is between bj and bk, we need to load B for each (bj, ii, bk).
   B[bk, bj*Tj+jj] doesn't depend on ii, so it could be loaded once per (bj, bk).
   But then we need bk to be OUTERMOST within bj, which conflicts with ii being inner.

   Loop: bi > bj > bk > ii > jj (CURRENT)
   - B loaded once per (bj, bk): ✓ (Tj loads per bk block)
   - A loaded Ti times per (bj, bk): Ti * Tj * nbk loads total per (bi,bj)

   Loop: bi > bj > ii > bk > jj
   - B loaded Tj times per (ii, bk): nbi * nbj * Ti * nbk * Tj loads = massive!
     = 2 * 4 * 8 * 16 * 4 = 4096 copies of B vs current 2 * 4 * 16 * 4 = 512 copies!
   - A loaded once per (ii, bk): nbi * nbj * Ti * nbk = 2*4*8*16 = 1024 (same as current).
   MUCH WORSE.

So the current bi > bj > bk > ii > jj loop order is optimal.

Let me try something very different: MULTI-LEVEL caching.

Idea: for each (bi, bj) outer block, precompute the ENTIRE bk sum:
  sum_k A[bi*Ti+ii, k] * B[k, bj*Tj+jj] for all (ii, jj)
This IS the C tile computation. No further optimization here.

Let me look at the actual cost breakdown to understand what to optimize:
  addr 2 (tmp): 3840 reads, cost 7,680 (10.4%)
  addr 1 (sA_cache): 4096 reads, cost 4,096 (5.6%)
  addr 3-6 (sB): 1024 reads each, cost 2,560 each = 10,240 total (13.9%)
  addr 7-38 (sC): 128 reads each, cost varies, total 20,736 (28.1%)
  Bulk A: 256 cells * 4 reads = 1024 reads at cost 7-18, total ≈ 13,100 (17.7%)
  Bulk B: 256 cells * 2 reads = 512 reads at cost 18-24, total ≈ 10,180 (13.8%)
  Bulk C: 256 reads at cost 24-29, total ≈ 6,750 (9.1%)

The biggest costs:
1. sC (28.1%): 20,736 - can't reduce
2. tmp/sA_cache (16%): 7,680 + 4,096 = 11,776 - swapping helps
3. sB (13.9%): 10,240 - can we reduce sB reads?

sB reads: Ti * nbi * nbj * nbk = 8 * 2 * 4 * 16 = 1024 per cell.
Total sB reads = Tj * 1024 = 4096.
sB is at addr 3-6 (cost 2-3 avg 2.5): total = 4096 * 2.5 = 10,240.

Can we reduce sB reads? No - each sB cell must be read Ti=8 times per bk.

Can we reduce sB address cost? sB is already at addr 3-6 (cost 2-3).
Addr 3 = cost 2, addr 4 = cost 2, addr 5 = cost 3, addr 6 = cost 3.
If sB was at addr 1-4 instead: but addr 1 is sA_cache and addr 2 is tmp.
We'd need to move sA_cache and tmp to higher addresses.

What if: sA_cache=1, tmp=2, sB=3 (1 cell for Tj=1), sC@4... -- this is Tj=1!
Already showed Ti=16,Tj=1 case: 97,338. Much worse (too many nbj iterations).

Wait, let me check: can we put ALL FOUR sB cells at addr 2-5?
If tmp=2 AND sB[0]=2: conflict! When we write to tmp (mul 2,...), we'd overwrite sB[0].
But sB[0] is already loaded before the inner ii loop, and we don't update it within ii.
Actually: mul tmp, sA_cache, sB(0) uses dest=tmp=2, src1=sA_cache=1, src2=sB(0).
If sB(0)=2 (same as tmp), then src2=dest=2. This means we're doing "mul 2, 1, 2":
  dest=2, src1=1, src2=2.
  Read src1=1 (sA_cache), read src2=2 (current sB[0]).
  Compute 1*sB[0] = sA_cache * sB[0].
  WRITE to dest=2, OVERWRITING sB[0]!
  Then "add sC(ii,jj), 2" reads the product (correct).
  But next time we need sB[0], it's been overwritten by the product!

So we can't alias sB[0] with tmp without careful control.

UNLESS: after each (ii, jj=0) pair, we reload sB[0]? But that doubles B loads. Worse.

Actually wait: can we make sB[0] at addr 1 (sA_cache) and use 2 for both tmp and sB[0]?
No, conflicts everywhere.

INSIGHT: The only way to reduce sB cost is to have sB cells at CHEAPER addresses.
The cheapest available addresses after (sA_cache=1, tmp=2) are 3, 4, 5, 6.
These are already the cheapest possible for sB.

IDEA: What if sB_base = 2 and we use a DIFFERENT addr for tmp?
  sA_cache=1, sB@2-5, tmp=? where tmp ≠ 1,2,3,4,5 and ≠ sC addresses.
  tmp could be a sC address? No -- when we write tmp in mul, we'd overwrite sC.
  tmp must be a dedicated, non-sC, non-sA_cache, non-sB address.

  So: sA_cache=1, sB@2-5, sC@6-37, tmp=38 (cost 7).
  tmp reads: 3840 * 7 = 26,880 vs current 3840 * 2 = 7,680.
  Penalty: 19,200. MUCH WORSE.

  Alternative: sA_cache=1, sB@2-5, tmp=7, sC@8-39 (shifting sC by 1)?
  But then tmp=7 is in the "sC region" if we're not careful.
  We can't use a sC address as tmp because mul overwrites dest.

  What about: sA_cache=1, tmp=7, sB@2-5, sC@8-39?
  tmp at cost 3 (addr 7): 3840 * 3 = 11,520 vs current 7,680.
  sB at addr 2-5 vs 3-6: saves (2+2+3+3)-(2+3+3+3) = 10-11 = -1. Not much.
  Overall: worse by 3,840.

  So tmp must stay at addr 2 (cheapest available after sA_cache=1).

ANOTHER IDEA: What if we SHARE tmp with one of the sC cells?
After we write to tmp (bk>0 mul), we then read from tmp (in add sC, tmp).
These happen in sequence: mul THEN add.
Could we use sC(ii=0, jj=0) as tmp?

  for bk=0: mul sC(0,0), sA_cache, sB(0) -- initialized sC(0,0)
  for bk=1: mul sC(0,0), sA_cache, sB(0) -- OVERWRITES sC(0,0) with product (WRONG!)

  We can't alias sC cells with tmp because mul overwrites the destination.

THE CURRENT LAYOUT APPEARS OPTIMAL. Let me try some variations of Ti and Tj with the new (sA=1, tmp=2) assignment.

Also: can we reduce the OUTER COPY overhead?
At end of each (bi,bj): copy sC to bulk C: Ti*Tj = 32 copies, each reading sC.
Total: 32 * nbi * nbj = 32 * 2 * 4 = 256 reads from sC.
These are already counted in the sC reads (part of the 128 per cell).
Can we eliminate these? No -- we must write C outputs somewhere.

What if C output is at LOWER addresses (to reduce output read cost)?
Output reads at C bulk (addr 551-806, cost 24-29).
Can we put C at addr 39-294 (A's location)? Then A would need to move.
But A is read more times (4x) vs C (1x output), so A should be at lower addresses.

FINAL CHECK: Can we use non-divisors for Ti and Tj?
N=16, divisors are {1,2,4,8,16}. For non-divisors, the tiling would be non-uniform.
We'd need special-casing for boundary tiles. Complex to implement and likely not beneficial.

Let me just run an exhaustive search with the new (sA_cache=1, tmp=2) layout
over all valid Ti,Tj pairs and also try Tk>1 variants.
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


def generate_sA1_tmp2(Ti: int, Tj: int, Tk: int = 1) -> str | None:
    """
    Optimal: sA_cache=1, tmp=2, sB@3..(Tk*Tj+2), sC@(Tk*Tj+3)..(Tk*Tj+Ti*Tj+2)
    For Tk>1: cache Tk consecutive sA entries.
    """
    if N % Ti != 0 or N % Tj != 0 or N % Tk != 0:
        return None
    nbi = N // Ti
    nbj = N // Tj
    nbk = N // Tk

    sA_cache_base = 1  # Tk cells at addr 1..Tk
    tmp = Tk + 1        # addr Tk+1

    # Wait: for Tk=1: sA_cache=1 (1 cell), tmp=2, sB@3..
    # For Tk=2: sA_cache=1,2 (2 cells), tmp=3, sB@4..
    # But if sA_cache occupies addr 1..Tk and tmp is at Tk+1:
    # Tk=1: sA_cache=1, tmp=2
    # Tk=2: sA_cache=1..2, tmp=3 -- but then tmp reads cost 2 (addr 3)
    #        vs Tk=1 where tmp=2 (cost 2)
    # So Tk>1 doesn't help: sA_cache has Tk cells at addr 1..Tk (avg cost increases),
    # tmp stays at Tk+1 (cost increases), sB starts at Tk+2 (cost increases), etc.
    # Tk=1 is optimal for this layout family.

    if Tk == 1:
        sA_cache = 1
        tmp_addr = 2
        sB_base = 3
        sC_base = sB_base + Tj  # = 3 + Tj
    else:
        # For Tk>1, sA_cache occupies multiple addrs
        sA_cache_base = 1  # Tk addrs
        tmp_addr = Tk + 1
        sB_base = Tk + 2
        sC_base = sB_base + Tk * Tj

    scratch_end = sC_base + Ti * Tj

    if Tk == 1:
        sA_fn = lambda kk=0: sA_cache  # always addr 1
        sB_fn = lambda kk, jj: sB_base + jj  # kk always 0
        sC_fn = lambda ii, jj: sC_base + ii * Tj + jj
    else:
        sA_fn = lambda kk: sA_cache_base + kk
        sB_fn = lambda kk, jj: sB_base + kk * Tj + jj
        sC_fn = lambda ii, jj: sC_base + ii * Tj + jj

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
                # Load sB tile
                for kk in range(Tk):
                    for jj in range(Tj):
                        lines.append(f"copy {sB_fn(kk,jj)},{B_at(bk*Tk+kk, bj*Tj+jj)}")

                for ii in range(Ti):
                    if Tk == 1:
                        lines.append(f"copy {sA_fn()},{A_at(bi*Ti+ii, bk)}")
                        for jj in range(Tj):
                            is_first = (bk == 0)
                            if is_first:
                                lines.append(f"mul {sC_fn(ii,jj)},{sA_fn()},{sB_fn(0,jj)}")
                            else:
                                lines.append(f"mul {tmp_addr},{sA_fn()},{sB_fn(0,jj)}")
                                lines.append(f"add {sC_fn(ii,jj)},{tmp_addr}")
                    else:
                        for kk in range(Tk):
                            lines.append(f"copy {sA_fn(kk)},{A_at(bi*Ti+ii, bk*Tk+kk)}")
                        for jj in range(Tj):
                            for kk in range(Tk):
                                is_first = (bk == 0) and (kk == 0)
                                if is_first:
                                    lines.append(f"mul {sC_fn(ii,jj)},{sA_fn(kk)},{sB_fn(kk,jj)}")
                                else:
                                    lines.append(f"mul {tmp_addr},{sA_fn(kk)},{sB_fn(kk,jj)}")
                                    lines.append(f"add {sC_fn(ii,jj)},{tmp_addr}")

            for ii in range(Ti):
                for jj in range(Tj):
                    lines.append(f"copy {C_at(bi*Ti+ii, bj*Tj+jj)},{sC_fn(ii,jj)}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


if __name__ == "__main__":
    print(f"Record to beat: {RECORD:,}")
    print(f"Best so far:    {BEST_SO_FAR:,}")
    print()

    print("=== Full sweep: sA_cache=1, tmp=2, all Ti,Tj,Tk ===")
    best_cost = BEST_SO_FAR
    best_ir = None
    best_config = None
    results = []

    divisors = [d for d in range(1, N+1) if N % d == 0]

    for Ti in divisors:
        for Tj in divisors:
            for Tk in [1]:  # Only Tk=1 for now (Tk>1 is handled differently)
                ir = generate_sA1_tmp2(Ti, Tj, Tk)
                if ir is None:
                    continue
                try:
                    cost = score_16x16(ir)
                    results.append((cost, Ti, Tj, Tk))
                    if cost < best_cost:
                        best_cost = cost
                        best_ir = ir
                        best_config = (Ti, Tj, Tk)
                        print(f"  *** NEW BEST! Ti={Ti} Tj={Tj} Tk={Tk}  "
                              f"cost={cost:,}  delta={RECORD-cost:+,}")
                except ValueError:
                    pass

    results.sort()
    print(f"\nTop 10:")
    for cost, Ti, Tj, Tk in results[:10]:
        print(f"  Ti={Ti:>2} Tj={Tj:>2} Tk={Tk:>2}  cost={cost:>10,}  delta={RECORD-cost:>+8,}")

    print("\n=== Cost breakdown for best (Ti=8, Tj=4) ===")
    ir = generate_sA1_tmp2(8, 4)
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

    if best_ir and best_cost < BEST_SO_FAR:
        out_path = Path(__file__).parent / "ir" / f"new_record_{best_cost}.ir"
        out_path.parent.mkdir(exist_ok=True)
        out_path.write_text(best_ir + "\n")
        print(f"\nSaved: {out_path}")

    print(f"\nFinal best: {best_cost:,}  (record: {RECORD:,}  delta={RECORD-best_cost:+,})")
