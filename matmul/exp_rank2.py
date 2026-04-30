#!/usr/bin/env python3
"""
Rank-2 update idea: instead of rank-1 updates (Tk=1), use rank-2 updates (Tk=2)
but with the SMARTER placement (sB first then sA).

And further: try to reduce the number of tmp reads by batching operations.

Also: consider that for Tk=2, we can use TWO outer products simultaneously,
potentially reusing one sA or sB load across both products.

For rank-2 update: C += A[:,0]*B[0,:] + A[:,1]*B[1,:]
  = two outer products summed

Standard Tk=1 (two separate iterations):
  iter 1: load sA (Ti cells) + load sB (Tj cells), compute Ti*Tj products
  iter 2: load sA (Ti cells) + load sB (Tj cells), compute Ti*Tj products
  Total loads: 2*Ti + 2*Tj = 2*(8+4) = 24 cells
  Total inner products: 2*Ti*Tj = 64

For rank-2 (Tk=2, keep BOTH columns of A and BOTH rows of B in scratchpad):
  load sA (2*Ti = 16 cells) + load sB (2*Tj = 8 cells), compute Ti*Tj products (2 per output)
  Total loads: 24 cells (same as two Tk=1 iterations)
  Total inner products: Ti*Tj*2 = 64 (same)
  But: sA is at addr 2+8+16=26 (2*Tj+2*Ti+2) -- much higher than Tk=1!

So Tk=2 doesn't reduce loads, just changes scratchpad addresses. We already showed Tk=2 is worse.

DIFFERENT IDEA: For the inner (ii,jj) loop, sA[ii] is read Tj times.
What if we cache sA[ii] in a dedicated tmp register before the jj loop?
  "copy tmp_ii, sA[ii]" then use tmp_ii Tj times.
  Cost: 1 read from sA[ii] + Tj reads from tmp_ii
  = addr_cost(sA[ii]) + Tj * addr_cost(tmp_ii)
  vs. Tj reads from sA[ii]: Tj * addr_cost(sA[ii])

For sA[ii] at addr 6-13 (cost 3-4):
  Current: Tj * addr_cost(sA[ii]) = 4 * 3.5 = 14 avg per ii
  With tmp cache at addr 1: 3.5 + 4*1 = 7.5 avg per ii

SAVINGS: 14 - 7.5 = 6.5 per ii per (bi,bj,bk)!
Total: 6.5 * Ti * nbi * nbj * nbk = 6.5 * 8 * 2 * 4 * 16 = 6,656
Minus: extra copy reads: Ti * nbi * nbj * nbk = 8 * 128 = 1,024 reads at cost 3.5 avg = 3,584

Wait: the SAVINGS come from reading sA[ii] only ONCE and then reading tmp Tj-1 MORE times.
Old: Tj reads of sA[ii] at cost 3.5
New: 1 read of sA[ii] at cost 3.5 + Tj reads of tmp at cost 1 + 1 copy op?

Actually: "copy tmp, sA[ii]" costs addr_cost(sA[ii]) (the SRC is what costs).
Then "mul dest, tmp, sB[jj]" costs addr_cost(tmp) + addr_cost(sB[jj]).

Old: "mul dest, sA[ii], sB[jj]" costs addr_cost(sA[ii]) + addr_cost(sB[jj]).

So by caching sA[ii] in tmp:
  Instead of addr_cost(sA[ii]) per mul, we pay:
  (1/Tj) * addr_cost(sA[ii]) + addr_cost(tmp) per mul (amortized over Tj jj values)

For sA[ii] at cost 3.5 avg, Tj=4, tmp at cost 1:
Old: 3.5 per mul
New: 3.5/4 + 1 = 1.875 per mul
Savings per mul: 1.625
Total muls = Ti * Tj * nbi * nbj * nbk = 8*4*2*4*16 = 4096
Total savings: 1.625 * 4096 = 6,656

But we also need to DO the copy operations:
Each "copy tmp, sA[ii]" is an extra instruction reading sA[ii].
Extra reads: Ti * nbi * nbj * nbk = 8 * 128 = 1024 reads at cost 3.5 avg = 3,584

Net savings: 6,656 - 3,584 = 3,072! That's a real improvement!

BUT: there's a problem. We have only ONE tmp address (addr 1).
To cache sA[0], sA[1], ..., sA[Ti-1] for the inner jj loop:
  We need Ti separate tmp addresses (one per ii value).
  But if we process jj INNERMOST (for a given ii), we only need ONE tmp:
    for ii in range(Ti):
      copy tmp, sA[ii]  # cache sA[ii]
      for jj in range(Tj):
        mul something, tmp, sB[jj]
        ...

This works! We only need tmp to hold ONE sA[ii] at a time (since ii is outer).
The inner loop is: for ii: {copy tmp, sA[ii]} then {for jj: use tmp and sB[jj]}

Let me verify this gives the savings we expect.

Actually wait: in the current implementation, the inner loop is:
  for ii in range(Ti):
    for jj in range(Tj):
      if bk==0: mul sC(ii,jj), sA(ii), sB(jj)
      else:     mul tmp, sA(ii), sB(jj); add sC(ii,jj), tmp

When bk>0: mul reads tmp (addr 1), sA(ii) (addr ~6-13), sB(jj) (addr 2-5).
The "tmp" in "mul tmp, ..." is the DESTINATION, not a source read!
Only src1 and src2 are read.

So "mul {tmp},{sA(ii)},{sB(jj)}" reads sA(ii) and sB(jj), writes tmp.
Then "add {sC(ii,jj)},{tmp}" reads sC(ii,jj) and tmp.

The tmp READ happens in the ADD instruction, not the MUL instruction.

To cache sA[ii]: we'd need to COPY sA[ii] somewhere cheaper before the jj loop:
  copy tmp2, sA[ii]    # reads sA[ii], writes tmp2 (some cheap addr)
  for jj:
    mul tmp, tmp2, sB[jj]   # reads tmp2 (cost 1) and sB[jj] (cost 2-3)
    add sC(ii,jj), tmp       # reads sC(ii,jj) and tmp

But we need a SECOND tmp address (tmp2) for the cached sA[ii].
And this tmp2 would be read Tj times per ii per (bi,bj,bk).

If tmp2 is at addr 1 and the original tmp is somewhere else:
  tmp2 = 1 (cheapest)
  original_tmp = ??? (needs to be at some address)

We could use tmp=1 for "sA cache" and put the accumulation tmp at some address after sC.
  addr 1: sA_cache (used to cache sA[ii] for Tj inner iterations)
  No addr for accumulation tmp -- but we can use sC_base + Ti*Tj as accumulation tmp!

Actually: if the dest of mul is at addr 1 (same as current tmp), and we cache sA[ii]
at addr 1 too... there's a conflict.

Let me think: to cache sA[ii] at addr 1:
  Step 1: "copy 1, sA(ii)" -- reads sA(ii), writes to addr 1
  Step 2: "mul X, 1, sB(jj)" for each jj -- reads addr 1 (cost 1) and sB(jj)
  Step 3: add to sC somehow

But what is X (the destination of mul)? We need it somewhere to then add to sC.
If X = sC(ii,jj) (direct write on first bk):
  "mul sC(ii,jj), 1, sB(jj)" -- reads 1 (cost 1) and sB(jj) (cost 2-3), writes sC.
  This is the first bk iteration! And it works!

For subsequent bk iterations:
  "copy 1, sA(ii)" -- reads sA(ii), writes to addr 1
  "mul TEMP, 1, sB(jj)" -- reads addr 1 and sB(jj), writes TEMP
  "add sC(ii,jj), TEMP" -- reads sC(ii,jj) and TEMP

Where is TEMP? We need a second tmp location.
Let's use addr 2 as TEMP if sB starts at addr 3:
  tmp (sA cache) = 1
  TEMP = 2 (also used by sB[0] if we're not careful)

Alternatively: dedicate addr N+Ti+Tj+Ti*Tj+2 (right after scratchpad) as TEMP.
But that would be at a HIGH address (cost 8 for current layout).

Actually: can we eliminate TEMP entirely?
After "copy 1, sA(ii)", we can do "mul sC(ii,jj), 1, sB(jj)" for jj=0 (first product).
Then for jj=1..Tj-1 and for bk>0:
  "mul TEMP, 1, sB(jj)"
  "add sC(ii,jj), TEMP"

We still need TEMP. Unless: we can use addr 1 (tmp) for TEMP too.
  After "copy 1, sA(ii)":
  for jj in range(Tj):
    if first accumulation for (ii,jj):
      mul sC(ii,jj), 1, sB(jj)   # 1=sA[ii] cached, writes directly to sC
    else:
      mul 1, 1, sB(jj)             # reads addr 1 (=sA[ii]), writes to addr 1!
      add sC(ii,jj), 1             # reads sC and addr 1 (now = sA[ii]*sB[jj])

CONFLICT: "mul 1, 1, sB(jj)" reads addr 1 (as src1=sA[ii]) and WRITES to addr 1.
The read happens before the write, so the result is: addr1 = addr1 * sB[jj] = sA[ii]*sB[jj].
Then "add sC, 1" correctly reads this product.

But now addr 1 = sA[ii]*sB[jj], not sA[ii] anymore!
For the NEXT jj, "mul 1, 1, sB(jj+1)" would compute sA[ii]*sB[jj] * sB[jj+1] which is WRONG.

So we can't reuse addr 1 as both sA_cache and tmp.

SOLUTION: Use TWO different cheap addresses:
  addr 1: sA cache (copy sA[ii] here once)
  addr 2: TEMP for products (mul 2, 1, sB[jj])

But addr 2 is currently the first sB address!
With sB at addr 2-5 (Tj=4 cells): we'd lose one sB slot if addr 2 is TEMP.

Alternative layout:
  addr 1: TEMP (for products)
  addr 2: sA cache
  sB: addr 3-6 (Tj cells)
  sA: NO LONGER NEEDED -- we copy into addr 2 each time
  sC: addr 7-38

Wait: we don't need a SEPARATE sA scratchpad if we just copy sA[ii] into addr 2 each time!
But we still need to load sA[ii] from bulk A into addr 2. We save by only loading it ONCE
per ii (instead of Tj times).

Let me redesign:
  addr 1: product TEMP
  addr 2: current sA[ii] cache
  sB: addr 3..Tj+2 (Tj cells)
  sC: addr Tj+3..Tj+Ti*Tj+2 (Ti*Tj cells)
  A_base: after sC

For Ti=8, Tj=4: scratch = 2 + 4 + 32 = 38 cells starting at addr 2+1 = 2+38 = 40?
Wait: tmp=1 (1 cell), sA_cache=2 (1 cell), sB=3-6 (4 cells), sC=7-38 (32 cells)
scratch_end = 39, A_base = 39.

Inner loop per (bi,bj,bk,ii):
  copy 2, A_at(bi*Ti+ii, bk)     # cache sA[ii] into addr 2 (reads bulk A)
  [we NO LONGER need an explicit sA scratchpad]
  for jj in range(Tj):
    sB(jj) = addr 3+jj
    is_first = (bk==0)
    if is_first:
      mul sC(ii,jj), 2, sB(jj)   # reads addr 2 (cost 2) and sB(jj) (cost 2-3)
    else:
      mul 1, 2, sB(jj)            # reads addr 2 (cost 2) and sB(jj), writes to tmp at addr 1
      add sC(ii,jj), 1            # reads sC(ii,jj) and addr 1 (cost 1)

Cost comparison per (bi,bj,bk,ii):
  Old: copy (addr Ti+2+jj, A_bulk) * Tj times (during sA load) -- one copy per COLUMN bk, per ii
       mul sC(ii,jj), sA(ii), sB(jj) or mul tmp, sA(ii), sB(jj):
       reads sA(ii) (cost 3-4) and sB(jj) (cost 2-3) per mul

  New: copy (addr 2, A_bulk) * 1 time per ii
       mul sC(ii,jj), 2, sB(jj) or mul 1, 2, sB(jj):
       reads addr 2 (cost 2) and sB(jj) (cost 2-3) per mul

Wait: in the old code, the sA LOAD (copy sA(ii), A_bulk) is a COPY of the BULK A cell.
The sA(ii) address is the DESTINATION (free). The A_bulk address is the SOURCE (costs).
So OLD load cost = addr_cost(A_bulk) -- same as NEW load cost.

The SAVINGS come from:
  OLD: mul reads sA(ii) at cost 3-4
  NEW: mul reads addr 2 at cost 2

Each mul saves (3.5 - 2) = 1.5 reads cost.
But each jj iteration now: for the non-first bk case:
  OLD: mul tmp, sA(ii), sB(jj) [reads sA(ii)+sB(jj)]; add sC(ii,jj), tmp [reads sC+tmp]
  NEW: mul 1, 2, sB(jj) [reads addr2+sB(jj)]; add sC(ii,jj), 1 [reads sC+addr1]

OLD tmp reads: cost 1 per add (addr 1)
NEW addr1 reads: cost 1 per add (addr 1) -- SAME

OLD sA reads per (bi,bj,bk,ii): addr_cost(sA(ii)) * Tj (one per jj mul)
NEW addr2 reads per (bi,bj,bk,ii): addr_cost(2) * Tj = 2 * Tj

For sA(ii) at cost 3-4 avg 3.5: OLD = 3.5 * 4 = 14; NEW = 2 * 4 = 8. SAVES 6!

But NEW also has EXTRA copy: "copy addr2, A_bulk" once per (bi,bj,bk,ii).
OLD has: "copy sA(ii), A_bulk" once per (bi,bj,bk,ii) -- SAME NUMBER OF COPIES!
The difference: destination.
  OLD dest: sA(ii) at addr Ti+2+jj -- FREE (writes are free)
  NEW dest: addr 2 -- FREE (writes are free)

Both are just copies reading from A_bulk. Same cost!

So: savings come PURELY from reading addr 2 instead of sA(ii) in the mul instruction.
Savings per mul = cost(sA(ii)) - cost(2) = 3.5 - 2 = 1.5 per mul
Total muls = Ti * Tj * nbi * nbj * nbk = 8*4*2*4*16 = 4096
Total savings = 1.5 * 4096 = 6,144

But there's one issue: for the FIRST bk iteration:
  OLD: mul sC(ii,jj), sA(ii), sB(jj) [reads sA(ii)+sB(jj)]
  NEW: mul sC(ii,jj), 2, sB(jj) [reads addr2+sB(jj)] -- same savings

And sC is at addr 7-38 (vs old addr 14-45) since we removed the sA scratchpad!
sC cost savings: 32 cells, reads 16*(nbi*nbj) = 128 times each.
Old sC addrs: 14-45 (cost 4-7), avg 5.5
New sC addrs: 7-38 (cost 3-7), avg: let me compute.

For Ti=8,Tj=4: new sC_base=7, cells at 7-38
  cost sequence: addr_cost(7)=3, addr_cost(8)=3, ..., addr_cost(9)=3, addr_cost(10)=4, ...
  7,8,9 -> cost 3 (3 cells)
  10-16 -> cost 4 (7 cells)
  17-25 -> cost 5 (9 cells)
  26-36 -> cost 6 (11 cells)
  37-38 -> cost 7 (2 cells)
  Total 3+7+9+11+2 = 32 cells ✓
  Sum = 3*3 + 7*4 + 9*5 + 11*6 + 2*7 = 9+28+45+66+14 = 162
  Avg = 162/32 = 5.0625

vs Old sC addrs 14-45 avg:
  14-16 -> cost 4 (3 cells)
  17-25 -> cost 5 (9 cells)
  26-36 -> cost 6 (11 cells)
  37-45 -> cost 7 (9 cells)
  Sum = 3*4 + 9*5 + 11*6 + 9*7 = 12+45+66+63 = 186
  Avg = 186/32 = 5.8125

sC savings: (5.8125 - 5.0625) * 4096 = 0.75 * 4096 = 3,072

Total estimated savings: 6,144 (mul) + 3,072 (sC) = 9,216!

But WAIT: sC has 32 cells (addr 7-38), which is LARGER than old layout (addr 14-45 also 32 cells).
The sC addresses are LOWER, which is what we want.

BUT: by eliminating sA scratchpad (Ti=8 cells), A_base moves from 46 to 38+1 = 39!
Bulk A savings: A starts at 39 instead of 46 -- 7 addresses lower.
sum(addr_cost(39+i) for i in 0..255) vs sum(addr_cost(46+i) for i in 0..255)
= 7 fewer high addresses at start.
Each of 256 A cells read nbj=4 times.
Savings per cell: addr_cost(39+i) vs addr_cost(46+i) -- at i=0: cost(39)=7 vs cost(46)=7.
Hmm, cost(39) = ceil(sqrt(39)) = ceil(6.24) = 7. cost(46) = ceil(sqrt(46)) = ceil(6.78) = 7.
So the first several cells have the SAME cost! The savings from lower A_base are minimal here.

Let me just implement and measure.
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


def generate_cached_sA(Ti: int, Tj: int) -> str | None:
    """
    Eliminate sA scratchpad. Cache current sA[ii] in addr 2 for the jj loop.
    tmp (for products) is at addr 1.
    sB at addr 3..Tj+2.
    sC at addr Tj+3..Tj+Ti*Tj+2.
    A_base at Tj+Ti*Tj+3.
    """
    if N % Ti != 0 or N % Tj != 0:
        return None
    nbi = N // Ti
    nbj = N // Tj
    nbk = N

    tmp = 1       # product tmp
    sA_cache = 2  # cache current A[ii, bk] here

    sB_base = 3
    sC_base = sB_base + Tj
    scratch_end = sC_base + Ti * Tj

    sB = lambda jj: sB_base + jj
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
                # Load sB row: B[bk, bj*Tj:(bj+1)*Tj]
                for jj in range(Tj):
                    lines.append(f"copy {sB(jj)},{B_at(bk, bj*Tj+jj)}")

                # For each ii: cache sA[ii] in addr sA_cache, then process all jj
                for ii in range(Ti):
                    # Cache A[bi*Ti+ii, bk] in sA_cache
                    lines.append(f"copy {sA_cache},{A_at(bi*Ti+ii, bk)}")

                    for jj in range(Tj):
                        is_first = (bk == 0)
                        if is_first:
                            # Direct mul into sC
                            lines.append(f"mul {sC(ii,jj)},{sA_cache},{sB(jj)}")
                        else:
                            # Use tmp for product, then accumulate
                            lines.append(f"mul {tmp},{sA_cache},{sB(jj)}")
                            lines.append(f"add {sC(ii,jj)},{tmp}")

            # Write sC to bulk C
            for ii in range(Ti):
                for jj in range(Tj):
                    lines.append(f"copy {C_at(bi*Ti+ii, bj*Tj+jj)},{sC(ii,jj)}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


def generate_cached_sA_v2(Ti: int, Tj: int) -> str | None:
    """
    Same idea but with sB CACHED at low addresses for ALL bk values,
    eliminating the need to reload sB Tj times per bk.

    Wait: sB changes for each bk (different row of B).
    We can't cache ALL sB rows -- too many cells.
    """
    return generate_cached_sA(Ti, Tj)


def generate_double_cache(Ti: int, Tj: int) -> str | None:
    """
    Cache BOTH sA and sB differently:
    - sA_cache at addr 2 (1 cell, for current ii)
    - sB_cache at addr 3 (1 cell, for current jj)  -- NO BENEFIT since sB loaded from bulk anyway
    Actually: if we cache sB[jj] in addr 3 and read from there:
      addr 3 is cheaper than sB_base+jj when Tj > 1?
      If sB_base = 3 (after sA_cache=2), then sB[0] = addr 3 anyway.
      No benefit.
    """
    return generate_cached_sA(Ti, Tj)


if __name__ == "__main__":
    print(f"Record to beat: {RECORD:,}")
    print(f"Best so far:    {BEST_SO_FAR:,}")
    print()

    print("=== sA-cached approach (no sA scratchpad) ===")
    best_cost = BEST_SO_FAR
    best_ir = None
    best_config = None

    for Ti in [1, 2, 4, 8, 16]:
        for Tj in [1, 2, 4, 8, 16]:
            if N % Ti != 0 or N % Tj != 0:
                continue
            ir = generate_cached_sA(Ti, Tj)
            if ir is None:
                continue
            try:
                cost = score_16x16(ir)
                delta = RECORD - cost
                marker = " *** BEATS BEST!" if cost < best_cost else ""
                sC_base = 3 + Tj
                sC_end = sC_base + Ti * Tj - 1
                A_base = sC_end + 1
                print(f"  Ti={Ti:>2} Tj={Tj:>2}  sC@{sC_base:>3}-{sC_end:>3}  "
                      f"A_base={A_base:>3}  cost={cost:>10,}  delta={delta:>+8,}{marker}")
                if cost < best_cost:
                    best_cost = cost
                    best_ir = ir
                    best_config = (Ti, Tj)
            except ValueError as e:
                print(f"  Ti={Ti:>2} Tj={Tj:>2}  ERROR: {e}")

    if best_ir and best_config:
        Ti, Tj = best_config
        print(f"\n=== Breakdown for best config: Ti={Ti}, Tj={Tj} ===")
        print(f"Cost: {best_cost:,}  delta={RECORD-best_cost:+,}")
        tiers = cost_breakdown(best_ir)
        total = sum(t["cost"] for t in tiers.values())
        for c in sorted(tiers):
            t = tiers[c]
            addrs = sorted(t["addrs"])
            pct = t["cost"] / total * 100
            print(f"  cost={c:>2}  addrs={addrs[0]:>3}-{addrs[-1]:>3}  "
                  f"reads={t['reads']:>7,}  cost={t['cost']:>8,}  {pct:>5.1f}%")

        if best_cost < RECORD:
            out_path = Path(__file__).parent / "ir" / f"new_record_{best_cost}.ir"
            out_path.parent.mkdir(exist_ok=True)
            out_path.write_text(best_ir + "\n")
            print(f"\nSaved: {out_path}")

    print(f"\nFinal best: {best_cost:,}  (record: {RECORD:,}  delta={RECORD-best_cost:+,})")
