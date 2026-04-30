#!/usr/bin/env python3
"""
Direction: sA-cache v2. Explore variations on the 73,858 result.

Best: Ti=8, Tj=4, cached sA at addr 2, cost=73,858
Layout: tmp=1, sA_cache=2, sB@3-6, sC@7-38, A@39-294, B@295-550, C@551-806

Key insight: by caching sA[ii] at addr 2 (cost 2 per read), we read addr 2 instead of
old sA addresses (3-4 cost), saving 1.5-2 per inner product.

Now explore:
1. Different Ti, Tj combinations (already done in v1, Ti=8,Tj=4 wins)
2. Different layout: can addr 2 and addr 3 be even cheaper? No, addr 2 is already cost 2.
3. Cache BOTH sA and sB?
4. Different inner loop order?
5. Asymmetric access: cache sA at addr 2, put sB at addr 3 (same as now)
   What if Tj=2? sC would start at addr 5 instead of 7.
6. What about Ti=8, Tj=2 (sC@5-20, A@21-276, B@277-532)?
   Already shows 78,550 -- worse.
   But let me see if Ti=8,Tj=2 beats Ti=8,Tj=4 with the NEW layout optimized.

Wait: the results show:
  Ti=8, Tj=2: cost=78,550
  Ti=8, Tj=4: cost=73,858
  Ti=16,Tj=2: cost=78,547
  Ti=16,Tj=4: cost=75,843

All worse than Ti=8,Tj=4. Let me try other options:
- Can we also cache sB?
- Or use a different sC layout?

Also: what if we use the SAME trick for sB (cache sB[jj] at addr 3 for outer loop)?
If ii is innermost (and jj is outermost), then sB[jj] is read Ti times.
We could cache sB[jj] at addr 2 for the ii loop:
  for jj: copy addr2, sB[jj]; for ii: mul ..., sA[ii], addr2

But then sA[ii] (cost 3-4) is read Ti times instead of jj times.
For Ti=8: sA read 8 times per (jj,bk) -- more reads from expensive sA.
For Tj=4: sB (old) read 4 times per (ii,bk) -- same as current sA reads.

Let's see which is better:
Current (sA cached at 2, sB at 3-6):
  sA_cache reads: Tj * nbi * nbj * nbk per (ii) = Tj * 128 per ii
  = 4 * (nbi*nbj*nbk=128) reads at cost 2 per (ii) per block
  Total: Ti * Tj * 128 * 2 = 8*4*128*2 = 8,192

Alternative (sB cached at 2, sA at 3-10):
  sB_cache reads: Ti * nbi * nbj * nbk per (jj) = Ti * 128 per jj
  = 8 * 128 reads at cost 2 per jj per block
  Total: Tj * Ti * 128 * 2 = 4*8*128*2 = 8,192

SAME COST. But old sA (addr 3-10) vs old sB (addr 3-6):
  sA has Ti=8 cells at addr 3-10 (cost 2-3-3-3-3-3-4-4 = all between 2-4)
  sB has Tj=4 cells at addr 3-6 (cost 2-3-3-3)

For sB-cached version:
  sA reads per mul: addr_cost(sA[ii]) where sA is at 3-10 (cost 2-4, avg 3)
  For each (jj, bk, bi, bj): Ti muls reading sA[ii]
  Total sA reads = Tj * Ti * nbi * nbj * nbk * cost_avg_sA
  = 4 * 8 * 128 * 3 = 12,288

vs current sA_cache reads = 8,192 -- current is BETTER.

So the current approach (cache sA, sB in fixed positions) is optimal.

Actually I should also try:
- Cache sA at addr 2, sB at addr 3: current approach
- Cache sA at addr 2, sB at addr 3..Tj+2: current approach ✓
- Put sC before sB after sA_cache: sA_cache@2, sC@3.., sB at higher?

Wait: if sC comes before sB:
  sA_cache=2, sC@3..(Ti*Tj+2), sB@(Ti*Tj+3)..(Ti*Tj+Tj+2)
  For Ti=8, Tj=4: sC@3-34, sB@35-38
  sB reads: addr 35-38 (cost 6) * 256 per cell = 6*256*4 = 6,144 (bad!)
  sC cost: lower starting addr but... let me compute.

Alternative ordering:
  sA_cache=2, sB@3-6, sC@7-38: CURRENT (cost 73,858)
  sA_cache=2, sC@3-34, sB@35-38: worse (sB at high addr)

What about:
  sA_cache=2, sB@3-6, sC in DIFFERENT layout?
  sC column-major vs row-major? Addresses are the same 7-38.
  Just access order changes, not cost.

NEW IDEA: Use addr 2 as sA_cache AND eliminate the explicit sA copy
by computing mul directly.

Actually: what if we use addr 2 NOT for sA_cache but to hold a RUNNING SUM
across the jj loop? That doesn't help...

ANOTHER IDEA: Can we also avoid the `copy sA_cache, A_at(...)` instruction
by computing it as part of the first mul?

"mul sC(0,0), A_at(...), sB(0)" -- this reads from BULK A directly!
But bulk A is at high addresses (39+), cost 7+. This is what we already do
for the FIRST bk iteration (direct mul into sC).

Hmm, for the first bk=0:
  OLD: copy sA_cache, A_at(bi*Ti+ii, 0); mul sC(ii,0), sA_cache, sB(0); mul sC(ii,1), sA_cache, sB(1); ...
  EVEN SIMPLER: for bk=0: mul sC(ii,jj), A_at(bi*Ti+ii, bk), sB(jj) -- reads bulk A directly!
    No sA_cache needed for first bk (no accumulation needed for first bk).

This saves 1 copy per ii per bk=0 block: Ti * nbi * nbj = 8 * 2 * 4 = 64 copy ops
And instead reads from bulk A in mul: 64 muls reading bulk A (cost 7+).
But we ALREADY read bulk A for the copy! So net:
  OLD bk=0, ii=0: copy sA_cache, A_at(0,0) [reads A_at]; mul sC(0,0), sA_cache, sB(0) [reads addr2, sB(0)]
  NEW bk=0, ii=0: mul sC(0,0), A_at(0,0), sB(0) [reads A_at, sB(0)]
  Saves: reading addr2 (cost 2). Costs: reading A_at instead of addr2 in mul.
  Net: A_at cost - addr2 cost = 7-2 = 5 per bk=0 mul.
  For all bk=0 muls: Ti * Tj * nbi * nbj = 8*4*2*4 = 256 muls: 5*256 = 1,280 WORSE!

So this optimization doesn't help (even though it eliminates the copy, it reads bulk A in mul).

Let me now try a different sA_cache location:
Can we put sA_cache at an address LOWER than 2?
Addr 1 is tmp (needed for product accumulation).
Can sA_cache = 1 and tmp = somewhere_else?

If tmp = sC_base (addr 7 for Ti=8,Tj=4):
  sA_cache=1 (cost 1), sB@2-5 (cost 2-3), sC@6-37 (cost 3-7), tmp=38 (cost 7)

But wait: tmp is used as the destination of "mul tmp, sA_cache, sB(jj)" then as source of "add sC, tmp".
If tmp=38 (cost 7): each add reads tmp at cost 7.
Old tmp at 1: each add reads at cost 1.
tmp reads = (nbk-1)*Ti*Tj*nbi*nbj = 15*8*4*2*4 = 3840 reads.
Cost penalty: 3840 * 6 = 23,040. Much worse.

What if we arrange so tmp=1 is the sA cache?
  tmp (product) = some_other_addr (e.g., 3)
  sA_cache = 1 (cost 1)
  sB @ 2 and 4-6 (use addr 2 for one sB, skip 3=tmp)... awkward.

Actually: what if sA_cache=1 and we DON'T use a separate tmp?
In the inner loop for bk>0:
  mul sA_cache, sA_cache, sB(jj)   -- reads sA_cache(=1) twice, overwrites addr 1
  add sC(ii,jj), sA_cache           -- reads sC and addr 1 (now=A*B_jj)

But then sA_cache is DESTROYED after the first jj! We can't use the same sA_cache for jj+1.

Unless: we restore it after each mul!
  copy sA_cache_copy, sA_cache  -- reads sA_cache, saves copy
  mul sA_cache, sA_cache, sB(jj)  -- overwrites sA_cache with A*B_jj
  add sC(ii,jj), sA_cache
  copy sA_cache, sA_cache_copy  -- restore sA_cache from copy

This requires an ADDITIONAL address (sA_cache_copy) and 2 extra copies per jj.
Not worth it.

CONCLUSION: The current approach (sA_cache at addr 2, tmp at addr 1) is optimal
for this family. The cost is 73,858.

Let me now think about whether we can improve further by combining the sA_cache idea
with the slot_order optimization.

For the current best (Ti=8,Tj=4):
  sA_cache @ 2 (reads cost 2)
  sB @ 3-6 (reads cost 2-3 avg 2.75)

What if we put sB FIRST at addr 2-5 (cost 2-3) and sA_cache at addr 6 (cost 3)?
Then: sB reads cost 2-3 (same), sA_cache reads cost 3 (vs 2).
Net: sA_cache reads: 4096 reads at cost 3 vs 2 = 4096 more cost. WORSE.

What about: NO sA_cache, NO sB cache, just use the sA_single approach from exp_nonsquare_v4?
The row-chunked approach was Ti=4,Tj=4 at cost 143,331 -- much worse.

Can we COMBINE: sA_cached at addr 2 AND also cache sB at addr 3?
For the current Ti=8,Tj=4: sB is at addr 3-6 (4 cells).
  sB[0] is at addr 3 (cost 2), sB[1] at addr 4 (cost 2), etc.
  These are ALREADY cheap! No benefit from caching sB.

What if Tj is larger (e.g., Tj=8): sB at addr 3-10 (cost 2-4).
  Ti=8,Tj=8: sC@11-74, A@75-330, B@331-586 -- A_base=75, B_base=331
  nbj=2, nbi=2, nbk=16.
  B read 2x (nbi=2). B at addr 331-586 (cost 19-25, avg 22). 256 cells * 2 * 22 = 11,264.
  A read 2x (nbj=2). A at addr 75-330 (cost 9-19, avg 13.5). 256 * 2 * 13.5 = 6,912.
  sC: 64 cells at 11-74 (cost 4-9, avg ~6.5). 64 * 128 reads. Hmm.

vs Ti=8,Tj=4:
  B read 2x at 295-550 (avg 20). 256*2*20 = 10,240
  A read 4x at 39-294 (avg 10). 256*4*10 = 10,240

Ti=8,Tj=8 might be interesting. Already measured as 77,243 -- worse.

Let me think about Ti=16,Tj=4 (cost 75,843):
  nbi=1, nbj=4, nbk=16.
  sC: 64 cells at 7-70 (cost 3-9, avg 6). 64 * 16 * 1 * 4 reads = 64 * 64 = 4096. Cost = 4096*6 = 24,576.
  A_base = 71 (cost 9). A read 4x (nbj=4). 256*4*cost_A. cost_A at 71-326 = avg 10. = 10,240.
  B_base = 327. B read 1x (nbi=1). 256*1*avg18.5 = 4,736.
  Total est: 24,576 + 10,240 + 4,736 + ... = ~39,552 + other = ~75,843 ✓

For Ti=8,Tj=4:
  sC: 32 cells at 7-38. 32*128 reads. avg cost 5. Cost = 32*128*5 = 20,480.
  A_base=39, B_base=295. A read 4x (nbj=4). 256*4*10 = 10,240. B read 2x. 256*2*20 = 10,240.
  Total: 20,480 + 10,240 + 10,240 + ... ≈ 40,960 + other ≈ 73,858 ✓

Can we find Ti,Tj where all three costs (sC, A, B) are simultaneously lower?
- Small sC: need Ti*Tj small => small Ti*Tj => more (bi,bj) pairs => more sC reads
  But total sC reads = 4096 * avg_cost(sC) always.
  Wait: for sA-cached: sC reads per cell = nbk * nbi * nbj = 16 * (N/Ti) * (N/Tj)
  Total sC reads = Ti*Tj * 16 * (N/Ti) * (N/Tj) = 16 * N^2 = 4096 (always!)
  So total sC reads = 4096. The COST depends on sC address.
  To minimize sC cost: put sC at lowest addresses.

For sA-cached: sC starts at addr 3 + Tj.
  To minimize sC_base: minimize Tj. Smallest Tj=1.
  Ti=16, Tj=1: sC_base=4, sC@4-19 (cost 2-5, avg 3.5). 4096 * 3.5 = 14,336.
  A_base=20. A read 16x (nbj=16). 256*16*avg9 = 36,864.
  B_base=276. B read 1x. 256*avg17 = 4,352.
  sA_cache reads: Tj*nbi*nbj*nbk = 1*1*16*16=256 reads per (ii per bi-block).
    Actually: Ti * (nbi*nbj*nbk*Tj) = 16 * (1*16*16*1) = 16*256 = 4096 reads at cost 2.
    Total sA_cache cost = 4096 * 2 = 8,192.
  sB reads: sB has 1 cell at addr 3 (cost 2). Reads = nbi*nbj*nbk*Ti = 1*16*16*16 = 4096 at cost 2.
    Total sB cost = 4096 * 2 = 8,192.
  tmp reads: (nbk-1)*Ti*Tj*nbi*nbj = 15*16*1*1*16 = 3840 at cost 1 = 3840.

  Total = 14,336 + 36,864 + 4,352 + 8,192 + 8,192 + 3840 + bulk_C
  = 75,776 + bulk_C
  bulk_C = sum(addr_cost(276+256+i) for i in range(256)) = sum(addr_cost(532+i) for i in range(256))
  ≈ 6,000
  Total ≈ 81,776

vs Ti=8,Tj=4: 73,858. Ti=16,Tj=1 is worse. ✓ (We measured Ti=16,Tj=1 as 97,338.)

So Ti=8,Tj=4 seems to be a good sweet spot.

Let me try other Ti,Tj combinations more systematically with the sA_cached layout.
And also try Ti=8,Tj=4 with REVERSED order: sC before sB.
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


def generate_sA_cached_reversed_sC(Ti: int, Tj: int) -> str | None:
    """
    sA_cache at addr 2, sC BEFORE sB (to get sC at lower addresses).
    sC@3..(Ti*Tj+2), sB@(Ti*Tj+3)..(Ti*Tj+Tj+2).

    But sB reads: Ti * nbi * nbj * nbk per cell, at addr Ti*Tj+3 to Ti*Tj+Tj+2.
    For Ti=8, Tj=4: sC@3-34, sB@35-38 (cost 6 each).
    sB: 4 cells * 256 reads * cost 6 = 6,144.
    vs current sB: 4 cells * 256 reads at cost 2.5 = 2,560.
    Penalty: 3,584.

    But sC starts at 3 instead of 7: 4 fewer address units.
    sC savings: 4096 reads * avg_cost_savings.
    avg cost at 3-34 vs 7-38: slight improvement for first few cells.
    Probably not worth it.
    """
    if N % Ti != 0 or N % Tj != 0:
        return None
    nbi = N // Ti
    nbj = N // Tj
    nbk = N

    tmp = 1
    sA_cache = 2
    sC_base = 3           # sC first
    sB_base = sC_base + Ti * Tj  # sB after sC
    scratch_end = sB_base + Tj

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
                for jj in range(Tj):
                    lines.append(f"copy {sB(jj)},{B_at(bk, bj*Tj+jj)}")
                for ii in range(Ti):
                    lines.append(f"copy {sA_cache},{A_at(bi*Ti+ii, bk)}")
                    for jj in range(Tj):
                        is_first = (bk == 0)
                        if is_first:
                            lines.append(f"mul {sC(ii,jj)},{sA_cache},{sB(jj)}")
                        else:
                            lines.append(f"mul {tmp},{sA_cache},{sB(jj)}")
                            lines.append(f"add {sC(ii,jj)},{tmp}")
            for ii in range(Ti):
                for jj in range(Tj):
                    lines.append(f"copy {C_at(bi*Ti+ii, bj*Tj+jj)},{sC(ii,jj)}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


def generate_sA_cached_standard(Ti: int, Tj: int) -> str | None:
    """Standard: sA_cache at 2, sB@3..(Tj+2), sC@(Tj+3).."""
    if N % Ti != 0 or N % Tj != 0:
        return None
    nbi = N // Ti
    nbj = N // Tj
    nbk = N

    tmp = 1
    sA_cache = 2
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
                for jj in range(Tj):
                    lines.append(f"copy {sB(jj)},{B_at(bk, bj*Tj+jj)}")
                for ii in range(Ti):
                    lines.append(f"copy {sA_cache},{A_at(bi*Ti+ii, bk)}")
                    for jj in range(Tj):
                        is_first = (bk == 0)
                        if is_first:
                            lines.append(f"mul {sC(ii,jj)},{sA_cache},{sB(jj)}")
                        else:
                            lines.append(f"mul {tmp},{sA_cache},{sB(jj)}")
                            lines.append(f"add {sC(ii,jj)},{tmp}")
            for ii in range(Ti):
                for jj in range(Tj):
                    lines.append(f"copy {C_at(bi*Ti+ii, bj*Tj+jj)},{sC(ii,jj)}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


def generate_sB_cached_standard(Ti: int, Tj: int) -> str | None:
    """Variation: cache sB[jj] at addr 2, sA at 3..(Ti+2), sC@(Ti+3).."""
    if N % Ti != 0 or N % Tj != 0:
        return None
    nbi = N // Ti
    nbj = N // Tj
    nbk = N

    tmp = 1
    sB_cache = 2
    sA_base = 3
    sC_base = sA_base + Ti
    scratch_end = sC_base + Ti * Tj

    sA = lambda ii: sA_base + ii
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
                # Load sA row
                for ii in range(Ti):
                    lines.append(f"copy {sA(ii)},{A_at(bi*Ti+ii, bk)}")
                # For each jj: cache sB[jj] at addr 2, then process all ii
                for jj in range(Tj):
                    lines.append(f"copy {sB_cache},{B_at(bk, bj*Tj+jj)}")
                    for ii in range(Ti):
                        is_first = (bk == 0)
                        if is_first:
                            lines.append(f"mul {sC(ii,jj)},{sA(ii)},{sB_cache}")
                        else:
                            lines.append(f"mul {tmp},{sA(ii)},{sB_cache}")
                            lines.append(f"add {sC(ii,jj)},{tmp}")
            for ii in range(Ti):
                for jj in range(Tj):
                    lines.append(f"copy {C_at(bi*Ti+ii, bj*Tj+jj)},{sC(ii,jj)}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


def generate_both_cached(Ti: int, Tj: int) -> str | None:
    """
    Cache BOTH: sA at addr 2, sB at addr 3.
    No sA or sB scratchpad at all!
    Inner loops restructured: for each (ii, jj): copy sA into 2, copy sB into 3, then mul.

    But we'd need to re-copy BOTH for each (ii, jj) combination.
    For ii: copy sA(ii) at addr 2
    For jj: copy sB(jj) at addr 3
    mul ..., 2, 3

    This requires nbi*nbj*nbk * Ti * Tj copies of sA (one per (ii,jj) combo)
    Plus nbi*nbj*nbk * Ti * Tj copies of sB.
    vs current: nbi*nbj*nbk * (Ti + Tj) copies.

    Ti*Tj vs Ti+Tj: for Ti=8, Tj=4: 32 vs 12. Much MORE copies!
    Bad idea.
    """
    return None


if __name__ == "__main__":
    print(f"Record to beat: {RECORD:,}")
    print(f"Best so far:    {BEST_SO_FAR:,}")
    print()

    print("=== sA-cached standard vs reversed sC ===")
    best_cost = BEST_SO_FAR
    best_ir = None
    best_config = None

    for Ti in [1, 2, 4, 8, 16]:
        for Tj in [1, 2, 4, 8, 16]:
            if N % Ti != 0 or N % Tj != 0:
                continue
            for name, gen in [("sA-cache", generate_sA_cached_standard),
                               ("sB-cache", generate_sB_cached_standard),
                               ("sC-first", generate_sA_cached_reversed_sC)]:
                ir = gen(Ti, Tj)
                if ir is None:
                    continue
                try:
                    cost = score_16x16(ir)
                    delta = RECORD - cost
                    if cost < best_cost:
                        best_cost = cost
                        best_ir = ir
                        best_config = (Ti, Tj, name)
                        print(f"  *** NEW BEST! Ti={Ti:>2} Tj={Tj:>2} {name}  "
                              f"cost={cost:>10,}  delta={delta:>+8,}")
                    elif cost < BEST_SO_FAR + 2000:
                        print(f"  Ti={Ti:>2} Tj={Tj:>2} {name}  cost={cost:>10,}  delta={delta:>+8,}")
                except ValueError as e:
                    print(f"  Ti={Ti:>2} Tj={Tj:>2} {name}  ERROR: {e}")

    print(f"\nBest: {best_cost:,}  delta={RECORD-best_cost:+,}")

    if best_ir and best_cost < BEST_SO_FAR:
        Ti, Tj, name = best_config
        out_path = Path(__file__).parent / "ir" / f"new_record_{best_cost}.ir"
        out_path.parent.mkdir(exist_ok=True)
        out_path.write_text(best_ir + "\n")
        print(f"Saved: {out_path}")

    # Also analyze all results to understand the landscape
    print("\n=== Full sA-cache sweep ===")
    results = []
    for Ti in [1, 2, 4, 8, 16]:
        for Tj in [1, 2, 4, 8, 16]:
            if N % Ti != 0 or N % Tj != 0:
                continue
            ir = generate_sA_cached_standard(Ti, Tj)
            if ir:
                try:
                    cost = score_16x16(ir)
                    results.append((cost, Ti, Tj))
                    delta = RECORD - cost
                    print(f"  Ti={Ti:>2} Tj={Tj:>2}  cost={cost:>10,}  delta={delta:>+8,}")
                except ValueError:
                    pass

    results.sort()
    print(f"\nTop 5:")
    for cost, Ti, Tj in results[:5]:
        print(f"  Ti={Ti:>2} Tj={Tj:>2}  cost={cost:>10,}  delta={RECORD-cost:>+8,}")

    print(f"\nFinal best: {best_cost:,}  (record: {RECORD:,}  delta={RECORD-best_cost:+,})")
