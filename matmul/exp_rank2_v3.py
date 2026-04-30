#!/usr/bin/env python3
"""
Direction: Further optimize the sA-cache approach.

Current best: Ti=8, Tj=4, sA-cache at addr 2, cost=73,858
Layout: tmp=1, sA_cache=2, sB@3-6, sC@7-38, A@39-294, B@295-550, C@551-806

The cost breakdown showed:
  cost=2: addr 2-4, reads=6,144, cost=12,288 (16.6%)  <- sA_cache and first 3 sB cells
  cost=3: addr 5-9, reads=2,432, cost=7,296  (9.9%)   <- rest of sB and first sC cells
  cost=4: addr 10-16, reads=896, cost=3,584  (4.9%)
  cost=5: addr 17-25, reads=1152, cost=5,760 (7.8%)
  cost=6: addr 26-36, reads=1408, cost=8,448 (11.4%)

Idea 1: Can we also cache sB in a two-level loop?
  Current inner loop (per bk):
    for jj in Tj: copy sB(jj) from bulk B
    for ii in Ti: copy sA_cache from bulk A; for jj: mul/add

  The sB load (Tj copies from bulk B) is Tj=4 reads per bk block.
  sB is read Ti=8 times per jj value per block.

  If we cache sB[jj] at addr 3 (overwriting each time for each jj):
    for jj: copy addr3, B_at(bk, bj*Tj+jj)
    for ii: copy sA_cache, A_at(...)
    mul addr4, sA_cache, addr3  <- no separate sB array needed
    ...

  But then we load sB once per jj, and sA_cache once per ii: total Ti+Tj loads per bk.
  Currently: sA_cache once per ii + sB(jj) once per jj = Ti+Tj loads also.
  Same number. But now we read from addr 3 (cost 2) instead of sB[jj] (cost 2-3).

  Wait: in current approach:
  - sB is loaded ONCE per bk iteration: Tj copies from B_bulk
  - sB[jj] is READ Ti times per bk (from addr 3..6, cost 2-3)

  With "inline sB" (no sB scratchpad):
  - sB[jj] would be read from B_bulk EACH TIME: Ti reads per jj per bk (from B_bulk)
  - B_bulk at addr 295-550 (cost 18-24)

  That's MUCH worse. The whole point is to avoid repeated B_bulk reads.

  The current scheme is optimal: load sB once, then read Tj-cell scratchpad Ti times.

Idea 2: What if we process multiple ii rows at once (SIMD-like)?
  E.g., copy sA_cache0, A[0,bk] and sA_cache1, A[1,bk] before the jj loop.
  Then for each jj: two muls into sC[0,jj] and sC[1,jj].

  This doubles the number of sA_cache slots needed (2 addrs instead of 1).
  But saves: 2 muls sharing the SAME sB[jj] read.

  Wait: in our current code for inner jj loop:
    "mul tmp, sA_cache, sB(jj)" reads sA_cache AND sB(jj).
    If we process two ii at once: we'd need TWO muls:
      "mul tmp1, sA_cache0, sB(jj)"  reads sA_cache0 and sB(jj)
      "mul tmp2, sA_cache1, sB(jj)"  reads sA_cache1 and sB(jj)

    sB(jj) is read TWICE for two muls. We can't "share" the read.
    There's no instruction that does "mul BOTH, sA0_and_sA1, sB_once".

  So processing multiple rows doesn't save sB reads. Same cost.

Idea 3: Reorganize the inner loop to reduce sC reads.
  Current: for bk=0: mul sC directly; for bk>0: mul tmp; add sC, tmp.
  The "add sC, tmp" reads BOTH sC and tmp.

  Alternative: accumulate tmp across all jj iterations for a fixed ii and bk?
  No: we need separate accumulation for each (ii,jj) target.

Idea 4: What if sC is stored MORE COMPACTLY?
  For each (ii,jj), we need a persistent accumulation cell.
  The 32 cells are required. No way to reduce.

Idea 5: Use a DIFFERENT accumulation strategy.
  Current: first bk -> direct write to sC; others -> mul tmp; add sC, tmp.
  Alternative: use sC itself as the "tmp" destination, with in-place add.

  For bk=0: mul sC(ii,jj), sA_cache, sB(jj)  [no change]
  For bk>0:
    mul sC(ii,jj), sC(ii,jj), sB(jj)  -- WRONG! This would compute sC*sB, not sA*sB+sC.

  Can't do this.

Idea 6: Use MULTIPLE sA caches to reduce repeated loads.
  For each bk iteration, each sA[ii] value is the same for all jj.
  We load it once (per bk) into addr 2 before the jj loop.
  This is already optimal for Ti=1 per-ii copy.

Idea 7: What if we preload sA values into MULTIPLE cheap addresses?
  If we had Ti cheap addresses (2..Ti+1), we could load ALL sA values there.
  This is exactly what the original sA scratchpad does!
  But at Ti=8: sA at addr 2-9 (cost 2-3 avg 2.5).
  vs sA_cache at addr 2 (cost 2) each time.

  Actually the current approach IS better because we use addr 2 for ALL Ti reads,
  whereas the old scratchpad had sA at addr 2-9 (avg cost 2.5, not 2).

  The KEY INSIGHT: the cache approach uses ONE cheap address for ALL sA values,
  at the cost of loading it once per ii per bk.

Idea 8: What if we can reduce the NUMBER of sA cache loads?
  Current: Ti * nbi * nbj * nbk = 8 * 2 * 4 * 16 = 1024 cache load operations.
  Each reads from bulk A (at cost 7-18) -- but this cost is unavoidable.
  The cache load is: "copy addr2, A_at(bi*Ti+ii, bk)"
  Cost = addr_cost(A_at(bi*Ti+ii, bk)).
  This is the SAME as the bulk A read cost we'd pay anyway.

  Can we reduce the NUMBER of loads? Only if we cache more than one sA value.
  If we have Ti=8 cheap addresses (2..9), and preload ALL 8 sA values for the current bk:
  - Load: 8 reads from bulk A (same cost)
  - Then for each jj: read sA[ii] from addr ii+2 (cost 2-3) vs addr 2 (cost 2).
  This is the original sA scratchpad approach! Already tried, costs more.

Idea 9: Can we use addr 1 for sB instead of tmp?
  If sB[0] is at addr 1 (cheapest!) and tmp is at addr 2:
  - tmp reads in add: cost 2 per add (vs 1 per add currently)
  - sB[0] reads: cost 1 per read (vs cost 2 currently for addr 3)

  For Ti=8, Tj=4:
  sB[0] reads = Ti * nbi * nbj * nbk = 8 * 128 = 1024 at cost 2 = 2048 currently.
  If at addr 1: 1024 * 1 = 1024. Savings: 1024.

  tmp reads = (nbk-1)*Ti*Tj*nbi*nbj = 3840 at cost 1 = 3840 currently.
  If tmp at addr 2: 3840 * 2 = 7680. Penalty: 3840.

  Net: WORSE by 2816.

Idea 10: TRY ALL POSSIBLE addr assignments for (tmp, sA_cache, sB cells, sC cells).
  The address space has the constraint: tmp ≠ sA_cache ≠ sB[j] ≠ sC[ii,jj].
  What's the OPTIMAL assignment?

  The read patterns:
  - tmp: (nbk-1)*Ti*Tj*nbi*nbj = 3840 reads
  - sA_cache: Tj * Ti * nbi * nbj * nbk = Tj*256 = 1024 reads
  - sB[jj]: Ti * nbi * nbj * nbk = 256 per cell, Tj=4 cells total = 1024 reads
  - sC[ii,jj]: nbk * nbi * nbj = 128 per cell, Ti*Tj=32 cells total = 4096 reads
  - Bulk A: nbj per cell = 4 per cell (256 cells)
  - Bulk B: nbi per cell = 2 per cell (256 cells)
  - Bulk C: 1 per cell (256 cells, at output)

  To minimize total cost: assign cheapest addresses to most-read things.
  Sorted by reads (descending):
    1. sC cells (4096 total): should be cheapest → addr 1..32
    2. tmp (3840 reads): addr 33 (next)
    3. sA_cache (1024 reads): addr 34
    4. sB cells (1024 total): addr 35..38
    Then bulk...

  But wait: tmp must be different from sC addresses (since we write to sC AND read from it).
  And if tmp is in addr 1..32 (sC region), we'd overwrite sC values!

  Actually: we WRITE to tmp (it's the destination of mul), so it can't be a sC address.
  Similarly: sA_cache gets written by "copy addr2, A_bulk".

  Can sC and tmp share addresses? No -- when we do "mul tmp, sA, sB; add sC, tmp",
  tmp and sC(ii,jj) must be DIFFERENT addresses.

  The OPTIMAL assignment assuming no conflicts:
  - tmp at addr 1 (cheapest non-sC address)
  - sA_cache at addr 2
  - sB at addr 3..Tj+2
  - sC at addr Tj+3..Tj+Ti*Tj+2

  But this puts sC AFTER sB, with sC starting at addr Tj+3.
  For Ti=8, Tj=4: sC @ 7..38. This is current layout!

  Alternative: put sC BEFORE sB after tmp and sA_cache:
  - tmp=1, sA_cache=2, sC@3..Ti*Tj+2, sB@Ti*Tj+3..Ti*Tj+Tj+2
  - For Ti=8,Tj=4: sC@3-34, sB@35-38 (cost 6)

  sC reads: 4096 reads at avg cost(3..34) = 5.0 avg = 20,480
  vs current: 4096 reads at avg cost(7..38) = 5.0625 avg = 20,736
  Savings: 256

  sB reads: 1024 reads at avg cost(35..38) = 6 avg = 6,144
  vs current: 1024 reads at avg cost(3..6) = 2.5 avg = 2,560
  Penalty: 3,584

  Net: WORSE by 3,328.

  So the current layout IS optimal for this category of approaches!

Summary: Within the current sA-cache paradigm, Ti=8, Tj=4 at 73,858 appears optimal.

Let me try to think outside this box entirely.

RADICAL IDEA: What if we precompute DIFFERENCES of A values?
For Strassen-like approach, instead of computing C += a[i]*b[j] for each (i,j,k),
compute something like: (a[0]+a[1])*(b[0]+b[1]) = a[0]*b[0] + a[0]*b[1] + a[1]*b[0] + a[1]*b[1]
in one operation? No -- this is still 4 products in the cost model (arithmetic is free, reads cost).

ACTUALLY: "(a[0]+a[1])*(b[0]+b[1])" computes a SUM PRODUCT.
In the Dally model: add addr_sum_A, sA(0), sA(1) costs addr_cost(sA(0)) + addr_cost(sA(1)).
Then mul addr_M1, addr_sum_A, addr_sum_B costs addr_cost(addr_sum_A) + addr_cost(addr_sum_B).

The sum addr_sum_A is a derived quantity, stored at some address.
If addr_sum_A is at addr 2 (cost 2 per read): we SAVE on the reads of the individual components.

For Strassen 2x2: M1 = (A[0,0]+A[1,1]) * (B[0,0]+B[1,1])
  Standard: 4 products reading A[0,0], A[1,1], B[0,0], B[1,1] each twice
  Strassen: compute A[0,0]+A[1,1] and B[0,0]+B[1,1], then multiply.
    Cost of computing sums: addr_cost(sA[0,0]) + addr_cost(sA[1,1]) + addr_cost(sB[0,0]) + addr_cost(sB[1,1])
    Cost of mul: addr_cost(sum_A) + addr_cost(sum_B)

  For sA at addresses 2-9 (cost 2-3) and sB at 10-13 (cost 4):
    Standard 4 products reading sA(0,0)+sA(1,1) = cost(2)+cost(9) = 2+3=5
    Reading sB(0,0)+sB(1,1) = cost(10)+cost(13) = 4+4=8
    Total: 5+8 = 13 per 2 components.
    Strassen: adds sA costs (2+3=5) + sB costs (4+4=8) + then reads sums at cost 2+2=4
    Total: 5+8+4 = 17 -- MORE reads!

This is why Strassen is worse in the Dally model.

Let me try one more approach: use a DIFFERENT approach to handle the accumulation.
Instead of mul+add separately, what if we use the fact that:
  a * b + c = a * b + c
can be computed as:
  t = mul(a, b)  [reads a, b]
  c += t         [reads c, t]

OR:
  a * (b + c/a) -- only works if a != 0 and c/a is precomputed

Not helpful in general.

OPTIMIZATION IDEA: Reduce sC reads by doing DIRECT accumulation with in-place mul-add.

In the current code for bk>0:
  "mul tmp, sA_cache, sB(jj)"   -- reads sA_cache and sB(jj)
  "add sC(ii,jj), tmp"           -- reads sC(ii,jj) and tmp

Total reads per iteration (non-first): addr_cost(2) + addr_cost(sB(jj)) + addr_cost(sC(ii,jj)) + addr_cost(1)

What if we could do: "add sC(ii,jj), sA_cache * sB(jj)"?
The ISA supports: mul dest, src1, src2 and add dest, src1, src2.
But NOT fused multiply-add.

HOWEVER: We could do "mul sC(ii,jj), sA_cache, sB(jj)" -- this writes to sC(ii,jj) but doesn't ADD.
It OVERWRITES! Not what we want.

UNLESS: we only call this for the first iteration, and use add for subsequent ones.
That's exactly what we do with direct_first=True! Already implemented.

NEW INSIGHT: What if we use ADDITION FIRST then multiply?
"add tmp, sA_cache, 0" -- can't do this with 0.
We'd need a zero address.

Actually: There IS a trick. For bk=0 first product: "mul sC(ii,jj), sA_cache, sB(jj)".
For bk>0: "mul tmp, sA_cache, sB(jj); add sC(ii,jj), tmp".

What if for bk>0, instead of reading sC(ii,jj) in the add, we use in-place:
"add sC(ii,jj), tmp" is "add sC(ii,jj), sC(ii,jj), tmp" which reads sC TWICE? No.

Actually "add dest, src" is shorthand for "add dest, dest, src" which reads both src1=dest and src2.
So it reads sC(ii,jj) (as src1) and tmp (as src2).

There's no way to avoid reading sC in the add. It's the accumulation step.

Let me try a smarter loop unrolling to reduce sC reads:
Current: for each bk, sC(ii,jj) is read once in the add.
What if we compute all bk contributions to sC(ii,jj) at once?
  For fixed (ii, jj): sC(ii,jj) = sum_{bk} sA[ii] * sB[jj] -- but sA and sB change per bk!

We can't precompute this more efficiently.

OK, let me try yet another angle: the INIT overhead.
For the FIRST bk=0: we use "mul sC(ii,jj), sA_cache, sB(jj)" -- no tmp read.
For ALL other bk: mul tmp + add sC: 2 reads per iteration.

What if we initialize sC differently?
  Before the bk loop: copy sC(ii,jj), const0 -- set to 0.
  But the IR doesn't have constants! We can't write 0.

Alternative: use the first bk value to initialize with direct mul (already done).

FINAL IDEA: Can we save on the COPY INSTRUCTIONS (loading A and B into scratch)?

Currently: copy sA_cache, A_bulk [costs addr_cost(A_bulk)]
           copy sB(jj), B_bulk [costs addr_cost(B_bulk)]

These are UNAVOIDABLE -- we must read from bulk.

The total cost of these copies:
  A loads: Ti * nbi * nbj * nbk * avg(A_bulk_cost) = 8*2*4*16*10 = 10,240
  B loads: Tj * nbi * nbj * nbk * avg(B_bulk_cost) = 4*2*4*16*20 = 10,240

These match what we see in the cost breakdown.

WHAT IF: we use a "smarter" A storage layout where frequently accessed A cells are closer?
We access A[i, bk] for each (bi, bj, bk, ii). Each A cell is read nbj times.
All A cells are read equally, so no cell is "more frequent" than others.

The bulk A cost = nbj * sum_{i,j} addr_cost(A_base + i*N + j)
This sum is minimized when A_base is as small as possible.
Current A_base = 39 (scratch_end for Ti=8,Tj=4).
Can't reduce further without making scratch more expensive.

WHAT IF: we try to store A row-major but with FEWER cells (e.g., only the A cells we actually use)?
No -- we use all 256 cells of A.

I'm converging on the conclusion that 73,858 may be near-optimal for this approach.

Let me try a variation: What if we use Tk=2 with sA_cached?
For Tk=2: sA needs Ti*Tk = 16 cells (instead of Ti=8 cells).
But with sA_cached approach: we cache 2 A elements at a time (for kk=0 and kk=1)?
Actually Tk=2 means: for each (bk, ii): load A[bi*Ti+ii, bk*2] and A[bi*Ti+ii, bk*2+1].
Then for each jj: compute sC(ii,jj) += sA[0]*sB[0,jj] + sA[1]*sB[1,jj].

We'd need 2 cache addresses (addr 2 and 3) for the two sA elements.
sB would be 2*Tj = 8 cells at addr 4-11.
sC at addr 12-43.
A_base = 44.

vs Tk=1 Ti=8,Tj=4: A_base=39.
A_base is 5 higher for Tk=2 -> slightly more expensive A reads.
Also sC is higher: 12-43 vs 7-38.
Probably worse.
"""

import sys, math
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from matmul import score_16x16, _parse

RECORD = 110_487
BEST_SO_FAR = 73_858
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


def generate_sA_cached_tk_any(Ti: int, Tj: int, Tk: int) -> str | None:
    """
    sA_cached approach extended to Tk>1.
    Cache Tk consecutive sA[ii,0..Tk-1] elements at addr 2..Tk+1.
    sB: Tk*Tj cells at addr Tk+2..
    sC: Ti*Tj cells at addr ...
    """
    if N % Ti != 0 or N % Tj != 0 or N % Tk != 0:
        return None
    nbi = N // Ti
    nbj = N // Tj
    nbk = N // Tk

    tmp = 1
    # For Tk>1: need Tk sA_cache cells (one per kk value)
    sA_cache_base = 2  # Tk cells at addr 2..Tk+1
    sB_base = sA_cache_base + Tk  # Tk*Tj cells
    sC_base = sB_base + Tk * Tj  # Ti*Tj cells
    scratch_end = sC_base + Ti * Tj

    sA_cache = lambda kk: sA_cache_base + kk
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
                # Load sB tile: B[bk*Tk:(bk+1)*Tk, bj*Tj:(bj+1)*Tj]
                for kk in range(Tk):
                    for jj in range(Tj):
                        lines.append(f"copy {sB(kk,jj)},{B_at(bk*Tk+kk, bj*Tj+jj)}")

                # For each ii: cache Tk A elements, then compute Ti*Tj partial products
                for ii in range(Ti):
                    for kk in range(Tk):
                        lines.append(f"copy {sA_cache(kk)},{A_at(bi*Ti+ii, bk*Tk+kk)}")

                    for jj in range(Tj):
                        for kk in range(Tk):
                            is_first = (bk == 0) and (kk == 0)
                            if is_first:
                                lines.append(f"mul {sC(ii,jj)},{sA_cache(kk)},{sB(kk,jj)}")
                            else:
                                lines.append(f"mul {tmp},{sA_cache(kk)},{sB(kk,jj)}")
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

    print("=== sA-cached Tk>1 ===")
    best_cost = BEST_SO_FAR
    best_ir = None
    best_config = None

    divisors = [d for d in range(1, N+1) if N % d == 0]
    results = []

    for Ti in divisors:
        for Tj in divisors:
            for Tk in divisors:
                scratch = Tk + Tk*Tj + Ti*Tj
                if scratch > 100:
                    continue
                ir = generate_sA_cached_tk_any(Ti, Tj, Tk)
                if ir is None:
                    continue
                try:
                    cost = score_16x16(ir)
                    results.append((cost, Ti, Tj, Tk))
                    if cost < best_cost:
                        best_cost = cost
                        best_ir = ir
                        best_config = (Ti, Tj, Tk)
                        print(f"  *** NEW BEST! Ti={Ti} Tj={Tj} Tk={Tk}  cost={cost:,}  delta={RECORD-cost:+,}")
                except ValueError:
                    pass

    results.sort()
    print(f"\nTop 10:")
    for cost, Ti, Tj, Tk in results[:10]:
        scratch = Tk + Tk*Tj + Ti*Tj
        print(f"  Ti={Ti:>2} Tj={Tj:>2} Tk={Tk:>2}  scratch={scratch:>3}  cost={cost:>10,}  delta={RECORD-cost:>+8,}")

    if best_ir and best_cost < BEST_SO_FAR:
        Ti, Tj, Tk = best_config
        out_path = Path(__file__).parent / "ir" / f"new_record_{best_cost}.ir"
        out_path.parent.mkdir(exist_ok=True)
        out_path.write_text(best_ir + "\n")
        print(f"\nSaved: {out_path}")

    print(f"\nFinal best: {best_cost:,}  (record: {RECORD:,}  delta={RECORD-best_cost:+,})")
