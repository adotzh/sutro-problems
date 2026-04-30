#!/usr/bin/env python3
"""
Final push: try creative variations to beat 73,858.

Current layout for Ti=8, Tj=4, sA-cached:
  addr 1: tmp (product accumulator) - read 3840 times at cost 1
  addr 2: sA_cache - read 4096 times at cost 2
  addr 3-6: sB[0..3] - read 1024 times each at cost 2-3
  addr 7-38: sC[0..31] - read 4096 times total at cost 3-7
  addr 39-294: A bulk - read 4*256=1024 times at cost 7-18
  addr 295-550: B bulk - read 2*256=512 times at cost 18-24
  addr 551-806: C bulk - read 256 times at cost 24-29 (output)

Key cost breakdown (from prev run):
  cost=2: 6,144 reads (12,288 cost = 16.6%)
  cost=3: 2,432 reads (7,296 cost = 9.9%)
  ...

The addr 2 reads dominate for sA_cache. But what can we do?

TRULY CREATIVE IDEA: What if we use the sC addresses AS a cache for sA?
After the first bk iteration (bk=0), sC[ii,jj] holds the first outer product.
For subsequent bk iterations, we read sC AND add to it.

What if: for the FIRST bk=0 iteration, we use a "column accumulation" approach?
  - sC[ii, 0..Tj-1] should hold sum_{bk=0..K-1} sA[ii]*sB[bk,jj]
  - What if we "precompute" the sum of all sA[ii] values across bk?
  - sum_sA[ii] = sum_{bk=0}^{15} A[bi*Ti+ii, bk]

This sum can be precomputed ONCE, but requires reading ALL A cells (already done).
Then we'd need sum_{bk} sA[ii,bk] * sB[bk, jj] which is the DOT PRODUCT of a column of A with a row of B.
This IS what we're computing! Just written differently.

PRECOMPUTED DOT PRODUCTS:
What if we precompute all 16x16 = 256 dot products B_col[j] · B_col[j']?
No, that doesn't help.

What about: precompute tmp = sum_{bk} A[i, bk] for each fixed i?
= row sum of A. But C[i,j] = sum_k A[i,k]*B[k,j] != row_sum_A * something useful.

LOCALITY IMPROVEMENT:
In the current bi>bj>bk>ii>jj loop, the sC(ii,jj) addresses jump around:
For each (bi,bj): we accumulate all (ii,jj) cells.
The sC accesses are sequentially for the jj loop.

What if we reorganize to: bi>bj>ii>bk>jj?
For each (bi,bj,ii): accumulate sC(ii, 0..Tj-1) over all bk.
Then write sC(ii,:) to bulk C for this ii.

This reduces sC to just Tj cells (one row at a time), BUT we need to reload A values for each ii separately AND for each bk.

This is the row-chunked approach from exp_nonsquare_v4.py which gave cost=143,331. Much worse.

What about: bi>ii>bj>bk>jj?
Same as above but ii outer. Each (bi,ii) block: process all bj,bk.
For fixed (bi,ii): need sC[jj] for each bj separately (different jj ranges).
= Same as row-chunked. Worse.

WILD IDEA: What if sC is stored COLUMN-MAJOR (jj outer) instead of row-major (ii outer)?
sC[jj, ii] at addr sC_base + jj*Ti + ii.
The inner loop (ii then jj) would access sC in column-major order.
But the address VALUES change: for the same total cells (32), just different ordering.
Cost depends on which addr gets which usage pattern.

In the current (ii outer, jj inner) loop:
  for ii in range(Ti): for jj in range(Tj): ...sC(ii,jj)...
  Accesses sC(0,0), sC(0,1), sC(0,2), sC(0,3), sC(1,0), ...

With row-major sC: addr = sC_base + ii*Tj + jj
  sC(0,0)=7, sC(0,1)=8, sC(0,2)=9, sC(0,3)=10, sC(1,0)=11, ...

With column-major sC: addr = sC_base + jj*Ti + ii
  sC(0,0)=7, sC(0,1)=15, sC(0,2)=23, sC(0,3)=31, sC(1,0)=8, ...

Row-major accesses are contiguous (sequential in jj). Column-major jumps.
The COST per access is based on address, not order. Same set of addresses.
Row-major vs column-major for sC: same set of addresses, same total cost.

INSIGHT: The 32 sC cells can be assigned to any 32 addresses in 7..38 range.
The READS per cell are:
  For each (bi,bj): sC(ii,jj) is read (nbk-1) = 15 times (accumulation) + 1 (writeback)
  = 16 reads per (bi,bj) per cell.
  Total per cell = 16 * nbi * nbj = 16 * 2 * 4 = 128 reads.

ALL sC cells have the SAME read count (128). So we should put the 32 cheapest addresses.
The cheapest 32 addresses after addr 6 (since 1-6 are reserved for tmp, sA_cache, sB):
  addr 7-38: cost 3(7,8,9), 4(10..16), 5(17..25), 6(26..36), 7(37,38)
  Sum = 3*3 + 7*4 + 9*5 + 11*6 + 2*7 = 9+28+45+66+14 = 162
  Average = 162/32 = 5.0625
  Total sC cost = 162 * 128 = 20,736

For ANY ordering of sC within addr 7-38, the cost is THE SAME (20,736).
So the layout doesn't matter within the sC region.

What if we could use FEWER cheap addresses by assigning sC to addr 7-38 efficiently?
We already use all 32 cells in this range. Optimal.

OK what about TRULY fundamentally different approaches?

APPROACH: Avoid sC entirely using a different algorithmic structure.

For computing C = A @ B where A is 16x16 and B is 16x16:
What if we compute C one ROW at a time?
C[i,:] = A[i,:] @ B (row i of C = ith row of A times B)

This is essentially a 1x16 vector times 16x16 matrix = 16 dot products.
For each i: read A[i,0..15] (16 cells) and all of B (256 cells).
Total bulk reads: 16 + 256 = 272 per row * 16 rows = 4352 reads (from bulk).
But each B cell is read 16 times (once per row of A).
This is standard naive matmul with i as outermost loop.

Alternatively: compute one COLUMN of C at a time.
C[:,j] = A @ B[:,j] (jth column of C = A times jth column of B)
For each j: read all A (256) and B[:,j] (16).
B column j is at scattered addresses in row-major B.

None of these escape the fundamental trade-off.

Let me try a COMBINED sA-sB cached approach where BOTH sA and sB use a scratchpad,
but organized differently to get sC at lower addresses.

Key: we want sC to start at the LOWEST possible address.
Minimum sC_base = 1 + (space for tmp) + (space for sA) + (space for sB)

What if tmp=1, sA=2 (1 cell), sB=3 (1 cell), sC@4...?
For Ti=8, Tj=4: only 2 total slots for sA and sB.
We'd need a single sA cache and a single sB cache.
Loop restructured: for ii, jj: copy sA, copy sB, mul, add.

But then for each (ii,jj,bk): 2 copies (sA and sB) + mul + add.
vs current: for each (ii) per bk: 1 copy (sA); for each (jj) per bk: 1 copy (sB).
Current copies per (bk): Ti + Tj = 8 + 4 = 12.
New copies per (bk): Ti * Tj * 2 = 64. MUCH MORE.

The extra copies dominate. Worse.

What about: sA is completely eliminated (read A directly from bulk in mul).
i.e., no sA_cache AT ALL.

For bk=0: mul sC(ii,jj), A_at(bi*Ti+ii, bk), sB(jj) [reads A_bulk + sB(jj)]
For bk>0: mul tmp, A_at(bi*Ti+ii, bk), sB(jj) [reads A_bulk + sB(jj)]
          add sC(ii,jj), tmp [reads sC + tmp]

This reads A_bulk DIRECTLY in the mul instruction. Cost per mul = cost(A_bulk) + cost(sB(jj)).
For A_bulk at addr 39+: cost 7-18.

vs current with sA_cache:
  copy sA_cache, A_at(...) [reads A_bulk once per ii per bk]
  mul ..., sA_cache, sB(jj) [reads addr2=sA_cache (cost 2)]

Current reads A_bulk once and addr2 Tj=4 times.
Direct reads A_bulk Tj=4 times.
Savings from cache: (4-1) * A_bulk_cost - 4 * cost(sA_cache=2) = 3*10 - 4*2 = 30-8 = 22 per (ii,bk).
Total savings: 22 * Ti * nbi * nbj * nbk = 22 * 8 * 2 * 4 * 16 = 22,528.

This is a HUGE savings! But we ALSO eliminate the copy cost of sA_cache.
The copy "copy addr2, A_at(...)" costs addr_cost(A_at(...)) ≈ 10 per (ii,bk).
By eliminating copies: save 10 * Ti * nbi * nbj * nbk = 10 * 8 * 2 * 4 * 16 = 10,240.

Wait but WITHOUT sA_cache, we read A_bulk in the mul instead. So:
Direct approach (no sA_cache):
  mul reads A_bulk: cost(A_bulk) per mul
  Total cost = Tj * nbi * nbj * nbk * sum(cost(A_at(i,bk)) for i,bk)
  = Tj * nbi * nbj * nbk * Ti * avg(A_bulk_cost)
  = 4 * 2 * 4 * 16 * 8 * 10
  = 4 * 128 * 8 * 10 = 40,960

With sA_cache:
  copy reads A_bulk: cost(A_bulk) per copy
  Total copy cost = nbi * nbj * nbk * Ti * avg(A_bulk_cost) = 128 * 8 * 10 = 10,240
  mul reads addr2 (sA_cache): 4096 reads at cost 2 = 8,192
  Total = 18,432

vs direct: 40,960. Cache is MUCH BETTER (saves 22,528).

So the sA_cache is critically important. The issue is sC must start AFTER tmp+sA_cache+sB.

Let me try ONE MORE creative idea: what if we don't have a persistent sC scratchpad,
but instead use a "rolling computation" where we compute each C output as we go?

For each output C[i,j]: compute sum_k A[i,k]*B[k,j] in a single pass.
This is Ti=1, Tj=1 tiling: very small sC (1 cell at addr 4), but many outer iterations.
We saw this gives cost=164,085 -- much worse.

CONCLUSION AFTER EXHAUSTIVE ANALYSIS:
The current best of 73,858 (Ti=8, Tj=4, Tk=1, sA-cached at addr 2) appears to be
near-optimal given the constraints of:
1. The IR instruction set (add/sub/mul/copy with 2-3 operands)
2. The Dally cost model (ceil(sqrt(addr)) per read)
3. No zero constants available

Let me verify we haven't missed any obvious slot orderings for the sA-cached approach
by trying ALL possible addr assignments for (tmp, sA_cache, sB[0..3], sC[0..31]).
"""

import sys, math, itertools
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


def generate_sA_cached_custom_addrs(Ti, Tj, tmp_addr, sA_cache_addr, sB_addrs, sC_addrs,
                                     A_base, B_base, C_base):
    """Generate with completely custom address assignments."""
    assert len(sB_addrs) == Tj
    assert len(sC_addrs) == Ti * Tj
    if N % Ti != 0 or N % Tj != 0:
        return None
    nbi = N // Ti
    nbj = N // Tj
    nbk = N

    sC = lambda ii, jj: sC_addrs[ii * Tj + jj]
    sB = lambda jj: sB_addrs[jj]

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
                    lines.append(f"copy {sA_cache_addr},{A_at(bi*Ti+ii, bk)}")
                    for jj in range(Tj):
                        if bk == 0:
                            lines.append(f"mul {sC(ii,jj)},{sA_cache_addr},{sB(jj)}")
                        else:
                            lines.append(f"mul {tmp_addr},{sA_cache_addr},{sB(jj)}")
                            lines.append(f"add {sC(ii,jj)},{tmp_addr}")
            for ii in range(Ti):
                for jj in range(Tj):
                    lines.append(f"copy {C_at(bi*Ti+ii, bj*Tj+jj)},{sC(ii,jj)}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


def try_all_addr_perms():
    """
    For Ti=8, Tj=4: try different address assignments for the scratchpad.

    Fixed: read patterns:
    - tmp: 3840 reads
    - sA_cache: 4096 reads
    - sB[j]: 1024 reads each (j=0..3)
    - sC[i]: 128 reads each (32 cells)

    The optimal assignment: highest reads -> lowest addresses.
    Sorted by reads: sA_cache(4096) > sC_cells(128*32=4096) ... equal at cell level.

    Actually per-cell: sA_cache reads 4096 times, sC cells read 128 times each.
    So sA_cache is MUCH more read per cell. It MUST be at the cheapest address.

    Cheapest non-conflicting assignment:
    - sA_cache = 1 (if tmp can go elsewhere)
    - tmp = 2 (cost 2)
    - sB = 3..6 (cost 2-3)
    - sC = 7..38 (cost 3-7)

    Let's try: swap tmp and sA_cache.
    """
    Ti, Tj = 8, 4
    nbi = N // Ti  # 2
    nbj = N // Tj  # 4
    nbk = N        # 16
    A_base_fixed = 39  # after 38 cells of scratch

    print("=== Custom addr assignment (swapped tmp and sA_cache) ===")

    # Current: tmp=1, sA_cache=2, sB=3-6, sC=7-38, A_base=39
    # Try: sA_cache=1, tmp=2, sB=3-6, sC=7-38, A_base=39
    # sA_cache reads: 4096 * cost(1) = 4096 * 1 = 4096 (vs 4096 * 2 = 8,192)
    # tmp reads: 3840 * cost(2) = 3840 * 2 = 7,680 (vs 3840 * 1 = 3,840)
    # Net: save 4096 - (7680-3840) = 4096 - 3840 = +256 (IMPROVEMENT!)

    ir_swap = generate_sA_cached_custom_addrs(
        Ti=8, Tj=4,
        tmp_addr=2, sA_cache_addr=1,
        sB_addrs=[3, 4, 5, 6],
        sC_addrs=list(range(7, 39)),
        A_base=39, B_base=295, C_base=551
    )
    if ir_swap:
        try:
            cost = score_16x16(ir_swap)
            print(f"  sA_cache=1, tmp=2: cost={cost:,}  delta={RECORD-cost:+,}")
            if cost < BEST_SO_FAR:
                print("  *** BEATS BEST!")
                return cost, ir_swap
        except ValueError as e:
            print(f"  ERROR: {e}")

    return BEST_SO_FAR, None


def try_opt_swap():
    """Try swapping tmp and sA_cache."""
    Ti, Tj = 8, 4
    nbi, nbj, nbk = 2, 4, 16

    print("=== sA_cache at addr 1 (instead of 2) ===")
    # sA_cache=1, tmp=2, sB=3-6, sC=7-38
    # A_base=39

    def generate(tmp_addr, sA_cache_addr):
        sB_base = 3
        sC_base = sB_base + Tj
        A_base = sC_base + Ti * Tj
        B_base = A_base + N * N
        C_base = B_base + N * N

        sB = lambda jj: sB_base + jj
        sC = lambda ii, jj: sC_base + ii * Tj + jj

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
                        lines.append(f"copy {sA_cache_addr},{A_at(bi*Ti+ii, bk)}")
                        for jj in range(Tj):
                            if bk == 0:
                                lines.append(f"mul {sC(ii,jj)},{sA_cache_addr},{sB(jj)}")
                            else:
                                lines.append(f"mul {tmp_addr},{sA_cache_addr},{sB(jj)}")
                                lines.append(f"add {sC(ii,jj)},{tmp_addr}")
                for ii in range(Ti):
                    for jj in range(Tj):
                        lines.append(f"copy {C_at(bi*Ti+ii, bj*Tj+jj)},{sC(ii,jj)}")

        lines.append(",".join(map(str, outputs)))
        return "\n".join(lines)

    best_cost = BEST_SO_FAR
    best_ir = None

    # Try different tmp and sA_cache positions
    for tmp_addr, sA_cache_addr in [(1, 2), (2, 1), (1, 3), (3, 1), (1, 4), (4, 1)]:
        if tmp_addr == sA_cache_addr:
            continue
        ir = generate(tmp_addr, sA_cache_addr)
        try:
            cost = score_16x16(ir)
            delta = RECORD - cost
            marker = " *** BEATS BEST!" if cost < best_cost else ""
            print(f"  tmp={tmp_addr}, sA_cache={sA_cache_addr}  cost={cost:,}  delta={delta:+,}{marker}")
            if cost < best_cost:
                best_cost = cost
                best_ir = ir
        except ValueError as e:
            print(f"  tmp={tmp_addr}, sA_cache={sA_cache_addr}  ERROR: {e}")

    return best_cost, best_ir


def try_all_Ti_Tj_with_swapped():
    """Try all Ti, Tj with sA_cache=1 and tmp=2."""
    print("\n=== All Ti,Tj with sA_cache=1, tmp=2 ===")
    best_cost = BEST_SO_FAR
    best_ir = None
    best_config = None

    results = []
    divisors = [d for d in range(1, N+1) if N % d == 0]

    for Ti in divisors:
        for Tj in divisors:
            if N % Ti != 0 or N % Tj != 0:
                continue

            nbi = N // Ti
            nbj = N // Tj
            nbk = N

            tmp_addr = 2
            sA_cache_addr = 1
            sB_base = 3
            sC_base = sB_base + Tj
            A_base = sC_base + Ti * Tj
            B_base = A_base + N * N
            C_base = B_base + N * N

            sB = lambda jj: sB_base + jj
            sC = lambda ii, jj: sC_base + ii * Tj + jj
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
                            lines.append(f"copy {sA_cache_addr},{A_at(bi*Ti+ii, bk)}")
                            for jj in range(Tj):
                                if bk == 0:
                                    lines.append(f"mul {sC(ii,jj)},{sA_cache_addr},{sB(jj)}")
                                else:
                                    lines.append(f"mul {tmp_addr},{sA_cache_addr},{sB(jj)}")
                                    lines.append(f"add {sC(ii,jj)},{tmp_addr}")
                    for ii in range(Ti):
                        for jj in range(Tj):
                            lines.append(f"copy {C_at(bi*Ti+ii, bj*Tj+jj)},{sC(ii,jj)}")

            lines.append(",".join(map(str, outputs)))
            ir = "\n".join(lines)

            try:
                cost = score_16x16(ir)
                results.append((cost, Ti, Tj))
                if cost < best_cost:
                    best_cost = cost
                    best_ir = ir
                    best_config = (Ti, Tj)
                    print(f"  *** NEW BEST! Ti={Ti} Tj={Tj}  cost={cost:,}  delta={RECORD-cost:+,}")
            except ValueError:
                pass

    results.sort()
    print(f"\nTop 10:")
    for cost, Ti, Tj in results[:10]:
        print(f"  Ti={Ti:>2} Tj={Tj:>2}  cost={cost:>10,}  delta={RECORD-cost:>+8,}")

    return best_cost, best_config, best_ir


if __name__ == "__main__":
    print(f"Record to beat: {RECORD:,}")
    print(f"Best so far:    {BEST_SO_FAR:,}")
    print()

    # Try swapping tmp and sA_cache
    c1, ir1 = try_opt_swap()

    # Try all Ti,Tj with swapped addresses
    c2, conf2, ir2 = try_all_Ti_Tj_with_swapped()

    best_cost = min(c1, c2)
    best_ir = ir1 if c1 < c2 else ir2

    if best_ir and best_cost < BEST_SO_FAR:
        out_path = Path(__file__).parent / "ir" / f"new_record_{best_cost}.ir"
        out_path.parent.mkdir(exist_ok=True)
        out_path.write_text(best_ir + "\n")
        print(f"\nSaved: {out_path}")

    print(f"\nFinal best: {best_cost:,}  (record: {RECORD:,}  delta={RECORD-best_cost:+,})")
