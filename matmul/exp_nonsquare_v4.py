#!/usr/bin/env python3
"""
Direction A v4: Advanced optimizations for non-square Tk=1 tiling.

Best so far: Ti=8, Tj=4, Tk=1, cost=82,477
Layout: sB@2-5, sA@6-13, sC@14-45

The cost breakdown shows:
  - sC reads dominate: addr 14-45, cost 4-7, ~37% of total
  - sA reads: addr 6-13, cost 3-4, ~23%
  - sB reads: addr 2-5, cost 2-3, ~18%
  - tmp reads (addr 1): ~5%
  - Bulk reads: ~17%

Ideas to reduce sC cost:
1. Reduce the number of sC accumulation reads
   - Currently: (nbk - 1) = 15 accumulations per (ii,jj) per (bi,bj) block
   - Can't reduce this without changing algorithm

2. Move sC to lower addresses
   - sC needs Ti*Tj = 32 cells. Currently at 14-45.
   - Could we place just sC starting at addr 2?
   - sC@2-33, sA@34-41, sB@42-45
   - sC costs: 2-6 (vs current 4-7)
   - sA costs: 6-7 (vs current 3-4) -- WORSE for sA
   - sB costs: 7 (vs current 2-3) -- MUCH WORSE for sB

3. Better: use a layout where sC gets addresses with lower cost
   - sB=4 cells first (addr 2-5, cost 2-3)
   - sC=32 cells next (addr 6-37, cost 3-7) -- 2 cells at cost 3!
   - sA=8 cells next (addr 38-45, cost 7)
   - vs current: sB@2-5, sA@6-13, sC@14-45

   Let's compare:
   Current:   sA@6-13 (cost 3-4), sB@2-5 (cost 2-3), sC@14-45 (cost 4-7)
   Alt slot:  sB@2-5 (cost 2-3), sC@6-37 (cost 3-7), sA@38-45 (cost 7)

4. Try Ti=16, Tj=4, Tk=1:
   - sB=4 cells, sA=16 cells, sC=64 cells
   - scratch = 84 cells, starting at addr 2+84=86
   - sC@22-85 (cost 5-10) -- gets expensive
   - But there's only 1 bi block (N/Ti=1), so we save on nbi*nbj repeats

5. Key insight: The cost of sC reads depends on:
   a) Number of accumulations per sC cell = nbk - 1 = 15 (for Tk=1, N=16)
   b) Address cost of each sC cell
   c) Number of (bi,bj) pairs = nbi * nbj

   Total sC accum cost = sum_{ii,jj} addr_cost(sC(ii,jj)) * (nbk - 1) * nbi * nbj

   For Ti=8, Tj=4: nbi=2, nbj=4, nbk=16
   sC cells: 32, addr 14-45, avg cost ~5.5
   sC accum reads = 32 * 15 * 2 * 4 = 3840
   sC cost = ~5.5 * 3840 ≈ 21,120

6. What if we use a bigger Tj to get fewer nbj repeats?
   Ti=8, Tj=8: nbi=2, nbj=2, nbk=16
   sC cells: 64, addr needs 2+8+8+64=82 cells starting at addr 2+16=18
   sC@18-81 (cost 5-10) -- expensive
   sC accum reads = 64 * 15 * 2 * 2 = 3840 (same number)
   But higher cost per read.

7. What if we do multi-level? Ti=8, Tj=4, but do chunks of bk (Tk=2)?
   This was tried in exp_acache.py -- Ti=8, Tj=4, Tk=2 gave cost=93,545 (worse).
   The issue is sA+sB get bigger (Ti*Tk + Tk*Tj = 16+8=24) and push sC higher.

Let me try a completely different approach:
- Split sC into multiple small pieces using a "register file" model
- Or: use partial sC accumulation with writebacks

Actually, let me look at the cost breakdown more carefully:
- addr 14-16 (3 cells): cost 4, reads ~2432
- addr 17-25 (9 cells): cost 5, reads ~1152
- addr 26-36 (11 cells): cost 6, reads ~1408
- addr 37-45 (9 cells): cost 7, reads ~1168

The high-cost sC reads (cost 6-7) come from addr 26-45. Can we reduce the number of cells there?

What if Ti=4, Tj=2, Tk=1?
- sC = 8 cells at addr 2+2+4=8 to addr 15 (cost 3-4)
- But nbi=4, nbj=8, nbk=16: many more blocks
- sC accum reads = 8 * 15 * 4 * 8 = 3840 (still same!)
- But lower cost per cell: avg ~3.5 vs 5.5

Actually! The number of sC reads depends on:
total_reads = sC_cells * (nbk-1) * nbi * nbj
= (Ti*Tj) * (N/Tk - 1) * (N/Ti) * (N/Tj)
= Ti*Tj * (16-1) * (16/Ti) * (16/Tj)
= Ti*Tj * 15 * 256 / (Ti*Tj)
= 15 * 256 = 3840

So the number of sC accumulation reads is ALWAYS 3840, regardless of Ti,Tj!
The only thing that changes is the ADDRESS COST of those reads.

So we want: minimize sum_{ii,jj} addr_cost(sC(ii,jj)) * (nbk-1) * nbi*nbj
= 3840 * avg_cost(sC cells)

To minimize avg cost of sC cells:
- Put sC as early as possible (lowest addresses)
- Minimize the number of other things before sC

The minimum possible sC start = 2 + 0 (if sC is first)
Then sA and sB go after sC, which makes them more expensive.

Trade-off: sA cells read more times than sC cells?
sA reads per cell = nbj * nbk * nbi * Tj (where Tj = Tj reads per (ii,kk) per block)
Wait: sA[ii] is read Tj times per (bi,bj,bk) block, total nbi*nbj*nbk blocks
= Tj * nbi * nbj * nbk = Tj * (N/Ti) * (N/Tj) * N = N^2 = 256 per cell

sB reads per cell = Ti * nbi * nbj * nbk = N^2 = 256 per cell

sC reads per cell = (nbk-1) * nbi * nbj = (N-1) * (N/Ti) * (N/Tj)
For Ti=8, Tj=4: (15) * 2 * 4 = 120 per cell

So sA and sB cells are read 256 times each vs sC cells 120 times.
This means we should prioritize low addresses for sA and sB over sC!

Current best: sB@2-5 (4 cells * 256 reads), sA@6-13 (8 cells * 256 reads), sC@14-45
The ordering sB,sA,sC makes sB (smallest, most reads per cell) cheapest.

But wait: sC cells are read 120 times vs sA/sB 256 times each.
Also sC has 32 cells vs sA 8 and sB 4. The total reads:
- sA: 8 cells * 256 = 2048 reads
- sB: 4 cells * 256 = 1024 reads
- sC: 32 cells * 120 = 3840 reads

So sC has MORE total reads! Even though each cell is read less often, there are more cells.

The key is to minimize: sum_addr (addr_cost * reads_for_that_addr)

For sA: 256 reads each; for sB: 256 reads each; for sC: 120 reads each.
Since 256 > 120, we should put sA and sB at lower addresses than sC.
This is already what we're doing (sB@2-5, sA@6-13, sC@14-45).

But what if we had fewer sC cells? That would push sC addresses lower.
sC cells = Ti * Tj. Minimize Ti * Tj while keeping performance.

For smaller Ti*Tj: smaller tiles, more tile iterations, same total sC reads (3840 always).
But sA,sB would be smaller too, fitting in lower addresses.

Let's compute the optimal Ti,Tj:
Want to minimize: sum_addr cost(addr) * reads(addr)
= sum_{sA addrs} cost * 256 + sum_{sB addrs} cost * 256 + sum_{sC addrs} cost * 120

For given Ti, Tj (with sB first, then sA, then sC):
sB: Tj cells at addr 2,3,...,Tj+1
sA: Ti cells at addr Tj+2,...,Tj+Ti+1
sC: Ti*Tj cells at addr Tj+Ti+2,...,Ti+Tj+Ti*Tj+1
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


def predict_cost_tk1(Ti: int, Tj: int, slot_order: tuple = (1, 0, 2)) -> dict:
    """
    Analytically predict cost for Tk=1 config.
    slot_order: position of (sA=0, sB=1, sC=2) in scratchpad.
    """
    nbi = N // Ti
    nbj = N // Tj
    nbk = N  # Tk=1

    sizes = [Ti, Tj, Ti * Tj]  # sA, sB, sC sizes
    order = sorted(range(3), key=lambda s: slot_order[s])
    bases = [0, 0, 0]
    cur = 2
    for s in order:
        bases[s] = cur
        cur += sizes[s]

    sA_base, sB_base, sC_base = bases

    # sA reads per cell: Tj times per (bi,bj,bk) block, total nbi*nbj*nbk blocks
    reads_per_sA = Tj * nbi * nbj * nbk

    # sB reads per cell: Ti times per block
    reads_per_sB = Ti * nbi * nbj * nbk

    # sC reads per cell: (nbk-1) times per (bi,bj) (accum) + 1 (copy to bulk)
    reads_per_sC_accum = (nbk - 1) * nbi * nbj
    reads_per_sC_copy = nbi * nbj
    reads_per_sC = reads_per_sC_accum + reads_per_sC_copy

    # tmp reads: (nbk-1) times per (ii,jj) per (bi,bj) = (nbk-1)*Ti*Tj*nbi*nbj total
    reads_tmp = (nbk - 1) * Ti * Tj * nbi * nbj
    cost_tmp = addr_cost(1) * reads_tmp

    # Bulk A reads: Ti cells per (bi,bk), loaded nbi*nbj*nbk times
    # But each A column (bi,bk) is loaded nbj times (wasteful in bi>bj>bk order)
    reads_bulk_A = Ti * nbi * nbj * nbk  # same formula as sA load copies
    # Bulk B reads: Tj per (bk,bj), loaded nbi*nbj*nbk times
    reads_bulk_B = Tj * nbi * nbj * nbk

    # Bulk C reads: Ti*Tj per (bi,bj) (output reads at end)
    reads_bulk_C = Ti * Tj * nbi * nbj

    A_base = cur
    B_base = A_base + N * N
    C_base = B_base + N * N

    cost_sA = sum(addr_cost(sA_base + i) * reads_per_sA for i in range(Ti))
    cost_sB = sum(addr_cost(sB_base + i) * reads_per_sB for i in range(Tj))
    cost_sC = sum(addr_cost(sC_base + i) * reads_per_sC for i in range(Ti * Tj))

    # Bulk costs (approximate: assume uniform distribution within each array)
    cost_bulk_A = sum(addr_cost(A_base + i * N + j) for i in range(N) for j in range(N)) * nbi * nbj
    # Wait: each bulk A cell is read only when that specific (bi,bk) block is active
    # A_at(bi*Ti+ii, bk) is read nbi*nbj times (once per bj) for each (bi,bk)
    # Total reads for A_at(i,j) = nbj (one per bj iteration)
    cost_bulk_A = sum(addr_cost(A_base + i * N + j) * nbj for i in range(N) for j in range(N))

    # B_at(bk, bj*Tj+jj) is read nbi times (once per bi) for each (bk,bj)
    cost_bulk_B = sum(addr_cost(B_base + i * N + j) * nbi for i in range(N) for j in range(N))

    # C_at is read once at output
    cost_bulk_C = sum(addr_cost(C_base + i * N + j) for i in range(N) for j in range(N))

    total = cost_sA + cost_sB + cost_sC + cost_tmp + cost_bulk_A + cost_bulk_B + cost_bulk_C

    return {
        'cost_sA': cost_sA, 'cost_sB': cost_sB, 'cost_sC': cost_sC,
        'cost_tmp': cost_tmp, 'cost_bulk_A': cost_bulk_A,
        'cost_bulk_B': cost_bulk_B, 'cost_bulk_C': cost_bulk_C,
        'total': total,
        'sA_base': sA_base, 'sB_base': sB_base, 'sC_base': sC_base,
        'A_base': A_base
    }


def generate_tk1(Ti: int, Tj: int, slot_order=(1, 0, 2)) -> str | None:
    """Generate Tk=1 tiled IR."""
    if N % Ti != 0 or N % Tj != 0:
        return None
    nbi = N // Ti
    nbj = N // Tj
    nbk = N

    tmp = 1
    sizes = [Ti, Tj, Ti * Tj]
    order = sorted(range(3), key=lambda s: slot_order[s])
    bases = [0, 0, 0]
    cur = 2
    for s in order:
        bases[s] = cur
        cur += sizes[s]

    sA_base, sB_base, sC_base = bases
    scratch_end = cur

    sA = lambda ii: sA_base + ii
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
                for ii in range(Ti):
                    lines.append(f"copy {sA(ii)},{A_at(bi*Ti+ii, bk)}")
                for jj in range(Tj):
                    lines.append(f"copy {sB(jj)},{B_at(bk, bj*Tj+jj)}")
                for ii in range(Ti):
                    for jj in range(Tj):
                        if bk == 0:
                            lines.append(f"mul {sC(ii,jj)},{sA(ii)},{sB(jj)}")
                        else:
                            lines.append(f"mul {tmp},{sA(ii)},{sB(jj)}")
                            lines.append(f"add {sC(ii,jj)},{tmp}")
            for ii in range(Ti):
                for jj in range(Tj):
                    lines.append(f"copy {C_at(bi*Ti+ii, bj*Tj+jj)},{sC(ii,jj)}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


def analytical_sweep():
    """Use analytical model to find best Ti, Tj, slot_order."""
    print("=== Analytical cost prediction for all Tk=1 configs ===")
    results = []
    slot_perms = list(itertools.permutations(range(3)))

    for Ti in [1, 2, 4, 8, 16]:
        for Tj in [1, 2, 4, 8, 16]:
            if N % Ti != 0 or N % Tj != 0:
                continue
            for slot in slot_perms:
                pred = predict_cost_tk1(Ti, Tj, slot)
                results.append((pred['total'], Ti, Tj, slot, pred))

    results.sort()
    print(f"{'Ti':>4} {'Tj':>4} {'slot':>10} {'sA_base':>8} {'sB_base':>8} {'sC_base':>8} {'predicted':>12}")
    for total, Ti, Tj, slot, pred in results[:20]:
        print(f"  {Ti:>4} {Tj:>4} {str(slot):>10} "
              f"{pred['sA_base']:>8} {pred['sB_base']:>8} {pred['sC_base']:>8} "
              f"{total:>12,.0f}")
    return results


def verify_top_configs(results):
    """Score the top analytically predicted configs."""
    print("\n=== Verifying top predicted configs ===")
    best_cost = BEST_SO_FAR
    best_config = None
    best_ir = None

    for total, Ti, Tj, slot, pred in results[:30]:
        ir = generate_tk1(Ti, Tj, slot)
        if ir is None:
            continue
        try:
            cost = score_16x16(ir)
            delta = RECORD - cost
            marker = " *** BEATS BEST!" if cost < best_cost else ""
            print(f"  Ti={Ti:>2} Tj={Tj:>2} slot={slot}  "
                  f"predicted={total:>10,.0f}  actual={cost:>10,}  delta={delta:>+8,}{marker}")
            if cost < best_cost:
                best_cost = cost
                best_config = (Ti, Tj, slot)
                best_ir = ir
        except ValueError as e:
            print(f"  Ti={Ti:>2} Tj={Tj:>2} slot={slot}  ERROR: {e}")

    return best_cost, best_config, best_ir


def try_alternative_accumulation():
    """
    Alternative: accumulate into sC differently.
    In the current approach, sC accumulates over bk.
    What if we use a larger sC (storing multiple partial results)?

    Or: use a two-level tiling where the outer bk level is also tiled.
    """
    print("\n=== Two-level bk tiling for Tk=1 ===")
    # Ti=8, Tj=4, but bk is split into groups of Gk
    # Each group: accumulate into sC over Gk outer products
    # Between groups: sC is written to / read from bulk C

    def gen_two_level(Ti, Tj, Gk):
        if N % Ti != 0 or N % Tj != 0 or N % Gk != 0:
            return None
        nbi = N // Ti
        nbj = N // Tj
        n_bk_groups = N // Gk  # number of bk groups

        tmp = 1
        # sB first, then sA, then sC
        sB_base = 2
        sA_base = sB_base + Tj
        sC_base = sA_base + Ti
        scratch_end = sC_base + Ti * Tj

        sA = lambda ii: sA_base + ii
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
                for bk_group in range(n_bk_groups):
                    for bk_inner in range(Gk):
                        bk = bk_group * Gk + bk_inner
                        for ii in range(Ti):
                            lines.append(f"copy {sA(ii)},{A_at(bi*Ti+ii, bk)}")
                        for jj in range(Tj):
                            lines.append(f"copy {sB(jj)},{B_at(bk, bj*Tj+jj)}")
                        for ii in range(Ti):
                            for jj in range(Tj):
                                is_first = (bk == 0)
                                if is_first:
                                    lines.append(f"mul {sC(ii,jj)},{sA(ii)},{sB(jj)}")
                                else:
                                    lines.append(f"mul {tmp},{sA(ii)},{sB(jj)}")
                                    lines.append(f"add {sC(ii,jj)},{tmp}")

                    # After each Gk group, write sC to bulk C if not last group
                    # Actually this doesn't make sense -- sC needs to persist across groups
                    # So just continue (same as flat bk loop).
                    pass

                # End of all bk: write sC to bulk C
                for ii in range(Ti):
                    for jj in range(Tj):
                        lines.append(f"copy {C_at(bi*Ti+ii, bj*Tj+jj)},{sC(ii,jj)}")

        lines.append(",".join(map(str, outputs)))
        return "\n".join(lines)

    # This is the same as the standard Ti=8,Tj=4,Tk=1 -- the Gk doesn't change anything
    # (the bk_inner loop is equivalent to bk loop with different indexing)
    pass


def try_sC_chunked():
    """
    Split sC into rows: process one row of sC at a time.
    Instead of a full Ti×Tj sC tile, keep only Tj cells in sC at once.
    After processing all bk for a single ii row, write that row to bulk C.
    Then move to next ii.

    This reduces sC from Ti*Tj cells to Tj cells, at the cost of loading sA[ii]
    for each ii separately (instead of all Ti rows at once).

    Ti=8, Tj=4: sC reduces from 32 to 4 cells
    Layout: sB@2-5 (4 cells), sC@6-9 (4 cells), sA@10 (1 cell used at a time)
    sC addr: 6-9, cost 3 each (vs 4-7 in current)

    But: we now process one row ii at a time
    - For each (bi, bj, ii): loop over bk and accumulate sC[ii, jj] for all jj
    - After bk loop: write sC row to bulk C
    - Repeat for next ii

    Tradeoff: sA[ii] is now loaded N*nbj*Ti = 16*4*8 = 512 times (vs N*nbj times)
    Actually wait: for each (bi, bj, ii, bk): load 1 element of sA into sC_row_tmp
    Hmm, let's reconsider.

    Plan: Reorganize outer loops as bi > bj > ii > bk > jj (process one row of sC at a time)
    - sA: just 1 cell (current a[ii] value)
    - sB: Tj cells (a row of B)
    - sC_row: Tj cells (accumulated dot products for row ii of C tile)
    - After inner bk loop: write sC_row to bulk C

    This eliminates the need for Ti*Tj sC scratchpad.
    """
    def gen_row_chunked(Ti, Tj):
        if N % Ti != 0 or N % Tj != 0:
            return None
        nbi = N // Ti
        nbj = N // Tj
        nbk = N

        tmp = 1
        # sA_single: 1 cell (just the current a[ii, bk])
        # sB: Tj cells
        # sC_row: Tj cells
        # Layout: sB@2..Tj+1, sC_row@Tj+2..2Tj+1, sA_single@2Tj+2
        sB_base = 2
        sC_base = sB_base + Tj
        sA_single = sC_base + Tj  # just 1 cell for current a[ii] element

        scratch_end = sA_single + 1

        sB = lambda jj: sB_base + jj
        sC_r = lambda jj: sC_base + jj  # sC row for current ii

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
                for ii in range(Ti):
                    # Process one row: C[bi*Ti+ii, bj*Tj:..] += A[bi*Ti+ii, :] @ B[:, bj*Tj:..]
                    for bk in range(nbk):
                        # Load A[bi*Ti+ii, bk] into sA_single
                        lines.append(f"copy {sA_single},{A_at(bi*Ti+ii, bk)}")
                        # Load B row: B[bk, bj*Tj..bj*Tj+Tj] into sB
                        for jj in range(Tj):
                            lines.append(f"copy {sB(jj)},{B_at(bk, bj*Tj+jj)}")
                        # Outer product: sC_row[jj] += sA_single * sB[jj]
                        for jj in range(Tj):
                            if bk == 0:
                                lines.append(f"mul {sC_r(jj)},{sA_single},{sB(jj)}")
                            else:
                                lines.append(f"mul {tmp},{sA_single},{sB(jj)}")
                                lines.append(f"add {sC_r(jj)},{tmp}")
                    # Write sC row to bulk C
                    for jj in range(Tj):
                        lines.append(f"copy {C_at(bi*Ti+ii, bj*Tj+jj)},{sC_r(jj)}")

        lines.append(",".join(map(str, outputs)))
        return "\n".join(lines)

    print("=== Row-chunked sC (sC has Tj cells instead of Ti*Tj) ===")
    best_cost = BEST_SO_FAR
    best_ir = None
    best_config = None

    for Ti in [1, 2, 4, 8, 16]:
        for Tj in [1, 2, 4, 8, 16]:
            if N % Ti != 0 or N % Tj != 0:
                continue
            ir = gen_row_chunked(Ti, Tj)
            if ir is None:
                continue
            try:
                cost = score_16x16(ir)
                delta = RECORD - cost
                marker = " *** BEATS BEST!" if cost < best_cost else ""
                Tj2 = Tj  # sC and sB size
                sB_base = 2
                sC_base = sB_base + Tj2
                sA_single = sC_base + Tj2
                print(f"  Ti={Ti:>2} Tj={Tj:>2}  sB@{sB_base}-{sB_base+Tj-1}  "
                      f"sC@{sC_base}-{sC_base+Tj-1}  sA@{sA_single}  "
                      f"cost={cost:>10,}  delta={delta:>+8,}{marker}")
                if cost < best_cost:
                    best_cost = cost
                    best_ir = ir
                    best_config = (Ti, Tj)
            except ValueError as e:
                print(f"  Ti={Ti:>2} Tj={Tj:>2}  ERROR: {e}")

    return best_cost, best_config, best_ir


if __name__ == "__main__":
    print(f"Record to beat: {RECORD:,}")
    print(f"Best so far:    {BEST_SO_FAR:,}")
    print()

    # 1. Analytical sweep
    results = analytical_sweep()

    # 2. Verify top analytical predictions
    best_cost, best_config, best_ir = verify_top_configs(results)

    # 3. Try row-chunked approach
    rc_cost, rc_config, rc_ir = try_sC_chunked()
    if rc_cost < best_cost:
        best_cost, best_config, best_ir = rc_cost, rc_config, rc_ir

    if best_ir and best_cost < RECORD:
        out_path = Path(__file__).parent / "ir" / f"new_record_{best_cost}.ir"
        out_path.parent.mkdir(exist_ok=True)
        out_path.write_text(best_ir + "\n")
        print(f"\nSaved: {out_path}")

    print(f"\nFinal best: {best_cost:,}  (record: {RECORD:,}  delta={RECORD-best_cost:+,})")
