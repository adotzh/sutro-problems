#!/usr/bin/env python3
"""
Layout optimization v2: Minimize bulk read cost by optimal placement.

For Ti=8, Tj=4, Tk=1:
  - Bulk A: read 4x each (nbj=4 times)
  - Bulk B: read 2x each (nbi=2 times)
  - Bulk C: read 1x (output)
  - scratch_end = 46

Current layout: A@46-301, B@302-557, C@558-813

Since A is read 4x and B is read 2x, A should be at lower addresses (already the case).
But can we overlap A and B access patterns to reduce cost?

Insight: The bulk A and B arrays are at fixed relative positions in memory,
but their order can be swapped. The input declaration line specifies WHICH
addresses hold A and B values.

Key: We can declare A to be at ANY set of addresses!
Currently: A at linear range [46..301], B at [302..557].
But what if we interleave A and B, putting A's most-accessed elements first?

For Tk=1 loop: A[i, bk] for bk=0..15 and i in [bi*Ti..(bi+1)*Ti].
The access PATTERN for A is: for each (bi,bj,bk), we read A[bi*Ti:bi*Ti+Ti, bk].
So A column bk is accessed nbi*nbj = 2*4 = 8 times.
All A columns are accessed equally.

For B: B[bk, bj*Tj:bj*Tj+Tj] for each (bi,bj,bk).
Each B cell is accessed nbi = 2 times.

So within A: all cells read equally.
Within B: all cells read equally.

But A is read 4x and B is read 2x total, so we should interleave A,B,A,B at low addresses.

Actually the simplest: just put A,B,C in the SAME order (row-major) but with optimal spacing.

Let me try: put A and B interleaved.
For row i: A_row_i at [A_base+i*N .. A_base+i*N+N-1]
           B_row_i at [B_base+i*N .. B_base+i*N+N-1]

Interleaved: A_row_0, B_row_0, A_row_1, B_row_1, ...
This doesn't help directly since total addresses are the same.

Alternative: Put them contiguously but in different order:
1. A (read 4x) first: addr 46-301
2. B (read 2x) next: addr 302-557
3. C last: addr 558-813

This is already optimal!

Wait, but what if we PUT C in the SAME address range as the scratchpad results?
C is only READ once (at output). So C at high addresses is fine — we READ it once.
But C is also WRITTEN: bulk C addresses cost reads when we do the final copy FROM sC.

Actually: "copy dest, src" costs addr_cost(src). The dest is free (writes are free).
So when we "copy C_at(i,j), sC(ii,jj)", we PAY for reading sC_at(ii,jj), NOT C_at(i,j).
And the output "read" at the end costs addr_cost(output_addr).

So C addresses only matter for:
1. The output read at end (once per C cell)
2. Nothing else (writes to C are free)

So we want to put C at LOWER addresses to reduce the output read cost!
Currently C is at addr 558-813, output cost = sum(addr_cost(558+i) for i in range(256)).

What if we put C BEFORE A and B?
New layout: scratch@2-45, C@46-301, A@302-557, B@558-813

Then C addresses are cheap (7-18 range), but A and B are at 302-813 (18-29 range).
A is read 4x from high addresses -> very expensive.
B is read 2x from even higher addresses -> very expensive.

Net: definitely worse.

What if we put A, C, B?
A@46-301 (read 4x), C@302-557 (read 1x for output), B@558-813 (read 2x)
C output reads are at 302-557 (cost 18-24). vs current 558-813 (24-29).
B reads at 558-813 (cost 24-29) vs current 302-557 (18-24).
A stays at 46-301 (cost 7-18, read 4x) -- unchanged.

Change: C output saves (reduces cost by ~6 per cell * 256) = ~1,500
        B reads increase by ~6 per cell * 256 * 2 = ~3,000
Net: WORSE by ~1,500.

What about B, A, C? (B at lowest, then A, then C)
B@46-301 (read 2x), A@302-557 (read 4x), C@558-813 (read 1x)
B reads at 46-301 (cost 7-18, read 2x): total = 2 * sum(addr_cost(46+i) for i in range(256))
A reads at 302-557 (cost 18-24, read 4x): total = 4 * sum(addr_cost(302+i) for i in range(256))
C output at 558-813 (cost 24-29): total = sum(addr_cost(558+i) for i in range(256))

Current (A, B, C):
B reads at 302-557 (cost 18-24): 2 * sum(302-557) = smaller * 2
A reads at 46-301 (cost 7-18): 4 * sum(46-301) = smaller * 4

For B,A,C vs A,B,C:
A reads: from 302-557 instead of 46-301 -> MUCH MORE EXPENSIVE
B reads: from 46-301 instead of 302-557 -> MUCH CHEAPER
Net: 4*A_reads more vs 2*B_reads savings. With same addr ranges swapped, A loss > B gain.

Conclusion: Current A,B,C order is optimal for our access pattern (A read 4x, B 2x).

Let me try a completely different angle: What if we reduce nbj (B reads) by using larger Tj?
Ti=8, Tj=8: nbj=2. B read 2x, A read 2x.
sC = 64 cells at addr 18-81 (cost 5-10) -- more expensive.

What if we think about the problem differently:
The total cost = f(scratchpad_cost) + g(bulk_cost)
- scratchpad cost ∝ scratchpad_size and read_frequency
- bulk cost ∝ bulk_size and read_multiplier (nbj for A, nbi for B)

For the TOTAL COST to decrease, we need to reduce one of these without increasing the other too much.

One unexplored idea: use MULTIPLE sC tiles simultaneously!
Instead of one Ti×Tj sC tile, keep two: sC0 for bj=0 and sC1 for bj=1.
Then we can process bj=0 and bj=1 together within each bk iteration.

But this would double the sC scratchpad size, making sC more expensive.
The savings: we load A column ONCE per bk instead of twice (since both bj tiles share the same bk).
Savings on sA reads = reduced, but sC doubled in size and cost.

Let me quantify:
Current Ti=8, Tj=4: nbi=2, nbj=4, nbk=16
With 2 sC tiles (for bj=0 and bj=1): process 2 bj at once.
sC0 and sC1: 2 * 32 = 64 cells
sB for bj=0: 4 cells, sB for bj=1: 4 cells = 8 cells
sA: 8 cells

Total scratch: 1 + 8 + 8 + 64 = 81 cells -> A_base at 82

sA reads: reduced by 2 (load once for 2 bj values)
= sA reads = Tj_total * nbi * nbj * nbk / 2 (... actually same total work)

Hmm, the total mathematical work is fixed. Loading A once for 2 bj just means
we process 2 outer products per (bi,bk) iteration.

Cost comparison:
Current (single sC):
  sA reads: Ti * nbj * nbk * nbi = 8 * 4 * 16 * 2 = 1024 reads at cost 3-4 (avg 3.5) = 3,584
  sC reads: Ti*Tj * (nbk-1) * nbi * nbj = 32 * 15 * 2 * 4 = 3840 reads at cost 4-7 (avg 5.5) = 21,120

Dual-sC:
  sA reads: Ti * nbj/2 * nbk * nbi = 8 * 2 * 16 * 2 = 512 reads at cost ~9 (A_base=82) = 4,608
  sC reads: 2*Ti*Tj * (nbk-1) * nbi * nbj/2 = 64 * 15 * 2 * 2 = 3840 reads at cost 9-15 avg = 46,200!

MUCH WORSE: sC addresses double and costs 9-15 instead of 4-7.

This confirms: the best strategy is to keep sC small and cheap.

The fundamental constraint: sC_size = Ti*Tj must be at low addresses.
For Ti=8, Tj=4: 32 cells, starting at addr 14 (after 12 cells for sA, sB, tmp).
sC costs: addr 14-45, cost 4-7.

The ONLY way to reduce this is to make sC smaller (use smaller Ti*Tj).
But smaller Ti*Tj means more (bi,bj) pairs, and we saw from the analysis that
the total sC reads = 3840 regardless of Ti,Tj.

WAIT. Let me reconsider. The number 3840 = sum(addr_cost(sC(ii,jj)) * reads_per_cell).
For different Ti,Tj: the ADDRESSES change, hence the COST changes even though reads_count is same.

For Ti=4, Tj=4: sC at addr 10-25 (cost 4-5), avg cost = 4.375
  sC total cost = 16 * (nbk-1=15) * nbi*nbj = 16 * 15 * 4 * 4 = 3840 reads * avg_cost 4.375 = 16,800

For Ti=8, Tj=4: sC at addr 14-45 (cost 4-7), avg cost = 5.5
  sC total cost = 32 * 15 * 2 * 4 = 3840 reads * avg_cost 5.5 = 21,120

Hmm wait, for Ti=4, Tj=4: nbi=4, nbj=4, nbk=16
  sC reads per cell = (nbk-1) * nbi * nbj = 15 * 4 * 4 = 240
  Total sC reads = 16 cells * 240 = 3840 ✓

But the total cost depends on the SUM over all sC cells of cost(addr)*reads_per_cell:
For Ti=4,Tj=4: sum = sum_{i=0}^{15} addr_cost(10+i) * 240

Let me compute which Ti,Tj gives the CHEAPEST sC cost:
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


def compute_expected_cost(Ti, Tj, slot_order=(1, 0, 2)):
    """Compute the scratchpad cost analytically."""
    nbi = N // Ti
    nbj = N // Tj
    nbk = N

    sizes = [Ti, Tj, Ti * Tj]
    order = sorted(range(3), key=lambda s: slot_order[s])
    bases = [0, 0, 0]
    cur = 2
    for s in order:
        bases[s] = cur
        cur += sizes[s]

    sA_base, sB_base, sC_base = bases[0], bases[1], bases[2]
    A_base = cur
    B_base = A_base + N * N
    C_base = B_base + N * N

    # sA: each cell read Tj times per (bi,bj,bk), total nbi*nbj*nbk blocks
    sA_cost = sum(addr_cost(sA_base + i) * Tj * nbi * nbj * nbk for i in range(Ti))

    # sB: each cell read Ti times per (bi,bj,bk), total nbi*nbj*nbk blocks
    sB_cost = sum(addr_cost(sB_base + i) * Ti * nbi * nbj * nbk for i in range(Tj))

    # sC accum: each cell read (nbk-1) times per (bi,bj) plus 1 writeback
    sC_reads_per_cell = (nbk - 1 + 1) * nbi * nbj  # accum + writeback
    sC_cost = sum(addr_cost(sC_base + i) * sC_reads_per_cell for i in range(Ti * Tj))

    # tmp reads: (nbk-1)*Ti*Tj times per (bi,bj)
    tmp_cost = addr_cost(1) * (nbk - 1) * Ti * Tj * nbi * nbj

    # Bulk A: each cell read nbj times (in bi>bj>bk order, A is loaded for every bj)
    bulk_A_cost = sum(addr_cost(A_base + i) * nbj for i in range(N * N))

    # Bulk B: each cell read nbi times
    bulk_B_cost = sum(addr_cost(B_base + i) * nbi for i in range(N * N))

    # Bulk C: read once at output
    bulk_C_cost = sum(addr_cost(C_base + i) for i in range(N * N))

    total = sA_cost + sB_cost + sC_cost + tmp_cost + bulk_A_cost + bulk_B_cost + bulk_C_cost

    return {
        'sA': sA_cost, 'sB': sB_cost, 'sC': sC_cost, 'tmp': tmp_cost,
        'bulk_A': bulk_A_cost, 'bulk_B': bulk_B_cost, 'bulk_C': bulk_C_cost,
        'total': total,
        'A_base': A_base, 'sC_base': sC_base, 'scratch_end': cur
    }


def generate_tk1(Ti, Tj, slot_order=(1, 0, 2)):
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


def analyze_costs():
    """Break down costs for the best configuration."""
    print("=== Cost analysis for Ti=8, Tj=4, Tk=1 ===")
    pred = compute_expected_cost(8, 4)
    print(f"Predicted: {pred['total']:,}")
    print(f"  sA: {pred['sA']:,} ({pred['sA']/pred['total']*100:.1f}%)")
    print(f"  sB: {pred['sB']:,} ({pred['sB']/pred['total']*100:.1f}%)")
    print(f"  sC: {pred['sC']:,} ({pred['sC']/pred['total']*100:.1f}%)")
    print(f"  tmp: {pred['tmp']:,} ({pred['tmp']/pred['total']*100:.1f}%)")
    print(f"  bulk_A: {pred['bulk_A']:,} ({pred['bulk_A']/pred['total']*100:.1f}%)")
    print(f"  bulk_B: {pred['bulk_B']:,} ({pred['bulk_B']/pred['total']*100:.1f}%)")
    print(f"  bulk_C: {pred['bulk_C']:,} ({pred['bulk_C']/pred['total']*100:.1f}%)")

    print("\n=== Compare with Ti=4, Tj=4, Tk=1 ===")
    pred44 = compute_expected_cost(4, 4, (0, 1, 2))
    print(f"Predicted: {pred44['total']:,}")
    print(f"  sA: {pred44['sA']:,}, sB: {pred44['sB']:,}, sC: {pred44['sC']:,}")
    print(f"  tmp: {pred44['tmp']:,}")
    print(f"  bulk_A: {pred44['bulk_A']:,}, bulk_B: {pred44['bulk_B']:,}")

    print(f"\nWhy Ti=8,Tj=4 beats Ti=4,Tj=4:")
    print(f"  sA+sB: {pred['sA']+pred['sB']:,} vs {pred44['sA']+pred44['sB']:,}  diff={pred44['sA']+pred44['sB']-(pred['sA']+pred['sB']):+,}")
    print(f"  sC: {pred['sC']:,} vs {pred44['sC']:,}  diff={pred44['sC']-pred['sC']:+,}")
    print(f"  tmp: {pred['tmp']:,} vs {pred44['tmp']:,}  diff={pred44['tmp']-pred['tmp']:+,}")
    print(f"  bulk_A: {pred['bulk_A']:,} vs {pred44['bulk_A']:,}  diff={pred44['bulk_A']-pred['bulk_A']:+,}")
    print(f"  bulk_B: {pred['bulk_B']:,} vs {pred44['bulk_B']:,}  diff={pred44['bulk_B']-pred['bulk_B']:+,}")


def try_all_possible_layouts():
    """
    Try ALL possible memory orderings for A, B, C bulk arrays.
    For each Ti, Tj: try A before B, B before A, etc.
    Also try interleaving individual rows of A and B.
    """
    print("\n=== Different A/B bulk orderings ===")
    best_cost = BEST_SO_FAR
    best_ir = None

    def gen_bulk_order(Ti, Tj, A_first=True):
        """Generate with A before B (A_first=True) or B before A (False)."""
        if N % Ti != 0 or N % Tj != 0:
            return None
        nbi = N // Ti
        nbj = N // Tj
        nbk = N

        tmp = 1
        # Best slot order: sB first, then sA, then sC
        sB_base = 2
        sA_base = sB_base + Tj
        sC_base = sA_base + Ti
        scratch_end = sC_base + Ti * Tj

        sA = lambda ii: sA_base + ii
        sB = lambda jj: sB_base + jj
        sC = lambda ii, jj: sC_base + ii * Tj + jj

        if A_first:
            A_base = scratch_end
            B_base = A_base + N * N
        else:
            B_base = scratch_end
            A_base = B_base + N * N
        C_base = max(A_base, B_base) + N * N

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

    for Ti, Tj in [(8, 4), (4, 4), (8, 8)]:
        for A_first in [True, False]:
            ir = gen_bulk_order(Ti, Tj, A_first)
            if ir:
                try:
                    cost = score_16x16(ir)
                    delta = RECORD - cost
                    marker = " *** BEATS BEST!" if cost < best_cost else ""
                    nbi = N // Ti
                    nbj = N // Tj
                    print(f"  Ti={Ti:>2} Tj={Tj:>2} A_first={str(A_first):<5}  "
                          f"(A read {nbj}x, B read {nbi}x)  "
                          f"cost={cost:>10,}  delta={delta:>+8,}{marker}")
                    if cost < best_cost:
                        best_cost = cost
                        best_ir = ir
                except ValueError as e:
                    print(f"  ERROR: {e}")

    return best_cost, best_ir


def try_col_major_storage():
    """
    Store A in column-major order: A[0..N-1, bk] is contiguous for each bk.
    This puts frequently accessed COLUMNS of A at low addresses.
    For Tk=1: we access A[:, bk] for bk=0..15.
    In col-major: A[:,0] at addr A_base, A[:,1] at A_base+N, etc.
    """
    print("\n=== Column-major A storage ===")

    def gen_colmaj_A(Ti, Tj):
        if N % Ti != 0 or N % Tj != 0:
            return None
        nbi = N // Ti
        nbj = N // Tj
        nbk = N

        tmp = 1
        sB_base = 2
        sA_base = sB_base + Tj
        sC_base = sA_base + Ti
        scratch_end = sC_base + Ti * Tj

        sA = lambda ii: sA_base + ii
        sB = lambda jj: sB_base + jj
        sC = lambda ii, jj: sC_base + ii * Tj + jj

        A_base = scratch_end  # A in column-major
        B_base = A_base + N * N
        C_base = B_base + N * N

        # Column-major: A[i, j] at addr A_base + j*N + i
        A_at_col = lambda i, j: A_base + j * N + i
        # Row-major (standard): B[i, j] at addr B_base + i*N + j
        B_at = lambda i, j: B_base + i * N + j
        C_at = lambda i, j: C_base + i * N + j

        # Input declaration order must match: A first, B second
        # For column-major A: declare A in column-major order
        inputs = ([A_at_col(i, j) for i in range(N) for j in range(N)] +
                  [B_at(i, j) for i in range(N) for j in range(N)])
        outputs = [C_at(i, j) for i in range(N) for j in range(N)]
        lines = [",".join(map(str, inputs))]

        for bi in range(nbi):
            for bj in range(nbj):
                for bk in range(nbk):
                    # In col-major: A[bi*Ti+ii, bk] = A_at_col(bi*Ti+ii, bk) = A_base + bk*N + bi*Ti+ii
                    # These are CONTIGUOUS! A column bk starts at A_base + bk*N
                    for ii in range(Ti):
                        lines.append(f"copy {sA(ii)},{A_at_col(bi*Ti+ii, bk)}")
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

    best_cost = BEST_SO_FAR
    for Ti, Tj in [(8, 4), (4, 4), (4, 8), (8, 8)]:
        ir = gen_colmaj_A(Ti, Tj)
        if ir:
            try:
                cost = score_16x16(ir)
                delta = RECORD - cost
                marker = " *** BEATS BEST!" if cost < best_cost else ""
                print(f"  Ti={Ti:>2} Tj={Tj:>2} col-major A  cost={cost:>10,}  delta={delta:>+8,}{marker}")
                if cost < best_cost:
                    best_cost = cost
            except ValueError as e:
                print(f"  Ti={Ti:>2} Tj={Tj:>2}  ERROR: {e}")

    return best_cost


if __name__ == "__main__":
    print(f"Record to beat: {RECORD:,}")
    print(f"Best so far:    {BEST_SO_FAR:,}")
    print()

    analyze_costs()
    best_cost1, best_ir1 = try_all_possible_layouts()
    best_cost2 = try_col_major_storage()

    best_cost = min(best_cost1, best_cost2)
    if best_cost < BEST_SO_FAR:
        print(f"\n*** New best: {best_cost:,}  delta={RECORD-best_cost:+,}")

    print(f"\nFinal best: {best_cost:,}  (record: {RECORD:,}  delta={RECORD-best_cost:+,})")
