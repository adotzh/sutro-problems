#!/usr/bin/env python3
"""
Combined optimizations: try every plausible idea to beat 82,477.

Key analysis from exp_layout_v2.py for Ti=8, Tj=4, Tk=1:
  sA: 14,336  (17.4%)  -- read 256x each, addr 6-13 (cost 3-4)
  sB:  9,216  (11.2%)  -- read 256x each, addr 2-5  (cost 2-3)
  sC: 23,808  (28.9%)  -- read 120x each, addr 14-45 (cost 4-7)
  tmp: 3,840   (4.7%)  -- read 3840 total, addr 1 (cost 1)
  bulk_A: 13,636 (16.5%) -- each cell read 4x (nbj=4)
  bulk_B: 10,822 (13.1%) -- each cell read 2x (nbi=2)
  bulk_C:  6,819  (8.3%) -- each cell read 1x (output)

To beat 82,477, we need to save ~3,000.

Idea 1: REDUCE sC cost
  sC is at addr 14-45. If we can put sC at addr 10-41 (save ~4 per sC cell at worst):
  - Need to reduce pre-sC scratchpad from 12 cells to 8 cells
  - Currently: tmp(1) + sB(4) + sA(8) = 13 cells total, so sC starts at 14
  - To save: reduce sA to 4 cells? But Ti=8...
  - Or: eliminate one level of indirection

Idea 2: REDUCE bulk_A cost
  Bulk A is at addr 46-301 (scratch_end=46 for Ti=8,Tj=4,Tk=1).
  Each A cell read 4 times (nbj=4).
  What if we use a smaller sC to push A lower?
  - Ti=4, Tj=4: scratch_end=26, A@26-281, B@282-537 --> LOWER A costs!

  For Ti=4,Tj=4: A@26-281 vs Ti=8,Tj=4: A@46-301
  bulk_A difference = 4 * (sum(addr_cost(26+i) for i in range(256)) -
                             sum(addr_cost(46+i) for i in range(256)))
  Let me compute: lower A base = cheaper A reads.

Idea 3: Use COLUMN-MAJOR layout for B
  For Tk=1 outer loop, each bk loads B[bk, bj*Tj:...].
  These are scattered in row-major B.
  In col-major B: B[bk, j] = B_base + j*N + bk (contiguous for fixed j).
  For the access pattern: B[bk, bj*Tj+jj] = B_base + (bj*Tj+jj)*N + bk
  These are NOT contiguous for different jj values.
  So col-major B is worse for our access pattern.

Idea 4: INTERLEAVE A and B in address space
  Goal: reduce avg address of both A and B by interleaving.
  A[0,0], B[0,0], A[0,1], B[0,1], ..., A[N-1,N-1], B[N-1,N-1]
  A is still read 4x and B 2x. But now they start at lower addresses.
  Interleaved starts at scratch_end=46.
  A[i,j] at addr 46 + 2*(i*N+j) (even addrs)
  B[i,j] at addr 46 + 2*(i*N+j) + 1 (odd addrs)
  A has 256 cells, addresses 46,48,50,...,556 (even)
  B has 256 cells, addresses 47,49,51,...,557 (odd)
  Both arrays are in [46..557].
  OLD: A in [46..301], B in [302..557].
  NEW: A in [46..556], B in [47..557].
  NEW A addresses are HIGHER (avg ~302 vs avg ~174) -- WORSE for A (read 4x)!
  NEW B addresses are LOWER (avg ~301 vs avg ~430) -- BETTER for B (read 2x).
  Net: likely worse since A penalty > B gain.

Idea 5: Try Tk=1 with OPTIMAL tmp placement
  Currently tmp=1 (cost 1). That's already minimal.

Idea 6: Combine Ti=8,Tj=4 with a BIGGER sC size to reduce outer reads
  More sC cells means we can handle a larger (bi,bj) combined tile.
  But bigger sC = higher addresses = more expensive.

Idea 7: Use the scratchpad MORE CLEVERLY
  For the best config: when does tmp cost more than it saves?
  Each (mul tmp, sA, sB) + (add sC, tmp) reads tmp from addr 1 (cost 1).
  Alternatively: (mul sC, sA, sB) for first product -- already done!
  The 3840 tmp reads at cost 1 = 3840 total cost.
  Could we eliminate MORE tmp reads? Only by using direct mul into sC.
  Currently: first bk product goes direct to sC, rest use tmp.
  What if we used an in-place accumulate differently?

  Wait: "add sC, tmp" uses 2-operand form: reads sC (as destination) and tmp.
  That means sC is ALSO read during this operation!
  sC read count = (nbk-1) accumulations per (ii,jj) per (bi,bj) = 15*2*4 = 120x
  Plus: the sC write at the end of each (bi,bj) is also a read: 2*4 = 8x
  So total sC reads = 128 per cell. But our analysis shows 120+8 = 128? But we said 120+1=121?

  Let me recheck: "add {sC(ii,jj)},{tmp}" expands to "add dest, dest, src"
  - This reads sC(ii,jj) (as src1=dest) AND reads tmp (as src2)
  - Writes sC(ii,jj)
  So yes: each accumulation reads sC once.

  "copy {C_at(...)},{sC(ii,jj)}" reads sC once more.

  Total sC reads per (ii,jj) per (bi,bj): (nbk-1) + 1 = nbk = 16
  Total sC reads per cell: 16 * nbi * nbj = 16 * 2 * 4 = 128

  Wait, I had 120 earlier. Let me recheck:
  - For bk=0: mul direct into sC (no sC read for first product)
  - For bk=1..nbk-1: mul tmp; add sC,tmp (reads sC)
  - Total accumulation reads per (ii,jj): (nbk-1) = 15
  - Plus 1 for copy sC to bulk C
  Total = 16 reads per (ii,jj) per (bi,bj). With nbi=2,nbj=4: 16*8 = 128 per sC cell.

  My prediction formula had 120+1=121... let me fix that.
  Actually: reads_per_sC_accum = (nbk - 1) = 15 (accumulations)
  reads_per_sC_copy = 1 (writeback to bulk C)
  reads_per_sC_per_tile = 16 per (bi,bj)
  Total per cell = 16 * nbi * nbj = 16 * 2 * 4 = 128

  Hmm, but the cost breakdown shows sC total = 23,808.
  If 128 reads per cell for 32 cells: avg_cost = 23,808 / (128*32) = 23,808/4096 = 5.81.
  That matches the avg of addr_cost(14..45) ~ 5.4..5.7.

OK the formula was: (nbk-1) * nbi * nbj + nbi*nbj = nbk * nbi*nbj = 16*2*4 = 128.
That's correct.

Idea 8: Use REGISTER REUSE for sA
  In the inner loop (ii,jj), sA[ii] is read Tj=4 times (once per jj).
  If we loaded sA[ii] into tmp, then read tmp instead of sA[ii]:
  But tmp has only 1 address, and sA has multiple values.
  We can't cache all sA in a single tmp.

Idea 9: Try output from LOWER bulk addresses
  If C was at addr 46-301 instead of 558-813, output cost would be lower.
  But then A or B can't be there.
  What if: C at 46-301, then A,B at 302-813?
  Input layout: inputs = A_addrs + B_addrs; outputs = C_addrs
  We could put C at LOWER addresses if we don't need to use those for A or B.

  Since C is only read at output (not during computation), let's try:
  C at low addrs, A+B at higher addrs.
  scratch_end=46, C@46-301, A@302-557, B@558-813
  Then: A reads cost 18-24 (4x) = expensive, B reads cost 24-29 (2x) = very expensive.
  Output read from C@46-301 (cost 7-18, 1x) = saves vs current 558-813 (24-29, 1x).
  Total C output savings: ~6 per cell * 256 = 1,536
  A read penalty: ~6 per cell * 256 * 4 = 6,144
  B read penalty: ~6 per cell * 256 * 2 = 3,072
  Net: MUCH WORSE.

Idea 10: OPTIMAL ORDERING of the A bulk array
  Within A, some cells are accessed more "clustered" than others.
  The access pattern is: for each (bi,bj,bk): read A[bi*Ti:(bi+1)*Ti, bk]
  Each cell A[i,j] is read nbj=4 times (same for all cells).
  So there's no optimization within A.

Conclusion: The current Ti=8, Tj=4, Tk=1 appears near-optimal for this approach.

Let me try a few more creative ideas:

Idea 11: MULTIPLE tmp addresses
  Currently: 1 tmp at addr 1. What if we use 2 tmps?
  tmp1 = 1, tmp2 = addr of something else unused.
  The inner loop: for each (ii,jj,bk), we compute:
    mul tmp, sA[ii], sB[jj]
    add sC[ii,jj], tmp
  With 2 tmps, we could pipeline:
    mul tmp1, sA[0,0], sB[0,0]
    mul tmp2, sA[0,1], sB[0,1]  # overlapping?
  But this doesn't reduce reads -- each mul still reads sA and sB once.

Idea 12: Precompute sA[ii]*sA[ii]... No, that doesn't help.

Idea 13: Strassen specifically for outer products
  For rank-1 update: C += a * b^T
  Standard: Ti*Tj reads of sA + Ti*Tj reads of sB = 2*Ti*Tj reads per outer product.
  With Strassen: can't reduce for rank-1 updates (each output element needs both a[i] and b[j]).

Idea 14: BLOCK-STRASSEN at the tile level
  Instead of 4 outer blocks per (bi,bj) (4 uses of inner product with 4 A columns and 4 B rows),
  use Strassen's algorithm at the 2x2 tile-BLOCK level.
  E.g., treat the 16x16 matrix as four 8x8 sub-matrices.
  Apply Strassen to get 7 sub-matmuls instead of 8.

  The 8x8 sub-matmuls are themselves done with the Ti=8,Tj=4 tiling.
  Strassen saves 1 of the 8 sub-matmuls (8*8*8 inner products -> 7*8*8).
  Each 8x8 sub-matmul costs (with Ti=8,Tj=4 tiling): ~10,000 cost.
  Saving 1 sub-matmul = saving ~10,000.
  But Strassen at the 2x2 block level requires computing linear combinations of blocks
  (e.g., A[0,0]+A[1,1]) which need extra scratchpad.

  This could be complex to implement correctly. Let me think...

  Actually the Strassen approach at 16x16 level means:
  - Decompose 16x16 into four 8x8 matrices: A00,A01,A10,A11 and B00,B01,B10,B11
  - Compute 7 "products" M1..M7 (each is an 8x8 matmul)
  - Combine to get C00,C01,C10,C11

  Each M_i is an 8x8 matmul (perhaps with linear combination of A-blocks as input).
  The key: instead of 8 sub-matmuls, we do 7.
  But each M_i's input is a LINEAR COMBINATION of A-sub-blocks (sums/differences).
  These linear combinations cost reads to compute.

  For the Dally model: reads of the linear-combined blocks cost based on WHERE those temps are.
  If we store the combined blocks in scratchpad: depends on scratchpad addresses.

  This is complex. Let me implement a simplified version.
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


def generate_block_strassen_8x8():
    """
    Apply Strassen at the 2x2 block level (each block is 8x8).
    Then use Ti=8, Tj=4, Tk=1 tiling for each 8x8 sub-matmul.

    Strassen formulas for C = A @ B (2x2 blocks, each block is 8x8):
    Let H = 8 (half-size). N = 16 = 2*H.
    A = [[A00, A01], [A10, A11]] (each A_ij is H x H = 8x8)
    B = [[B00, B01], [B10, B11]] (each B_ij is H x H = 8x8)

    M1 = (A00 + A11) @ (B00 + B11)
    M2 = (A10 + A11) @ B00
    M3 = A00 @ (B01 - B11)
    M4 = A11 @ (B10 - B00)
    M5 = (A00 + A01) @ B11
    M6 = (A10 - A00) @ (B00 + B01)
    M7 = (A01 - A11) @ (B10 + B11)

    C00 = M1 + M4 - M5 + M7
    C01 = M3 + M5
    C10 = M2 + M4
    C11 = M1 - M2 + M3 + M6

    We need to:
    1. Compute the 8x8 block sums/differences (e.g., A00+A11, A10-A00)
    2. Store them in scratchpad
    3. Do 7 sub-matmuls of 8x8 matrices
    4. Combine into C

    The main scratchpad usage:
    - For each sub-matmul: need sA (8 cells), sB (4 cells), sC (32 cells) at low addr
    - But we also need temp 8x8 storage for the input combinations

    Temporary 8x8 storage: 64 cells each, multiple pairs needed.
    This could be very expensive. Let's estimate.

    Actually a key issue: in Strassen, the "operands" to each sub-matmul are
    the SUM/DIFFERENCE of two 8x8 matrices. These combined matrices need to
    be stored SOMEWHERE to serve as the "bulk A" and "bulk B" for each sub-matmul.

    If we place these combo matrices right after the scratchpad:
    Each combo matrix: 64 cells. We need up to 4 of them simultaneously.
    4 * 64 = 256 additional cells -> A_base much higher.

    This would make all sub-matmul bulk reads more expensive.

    Let me try to estimate if Strassen at the block level helps.

    Standard 7 sub-matmuls: each costs ~cost(Ti=8,Tj=4,Tk=1,N=8)
    But N=8 means: nbi=1, nbj=2, nbk=8, scratch same, A@46-109, B@110-173, C@174-237
    Much cheaper A,B,C reads!
    Cost of 8x8 sub-matmul with Ti=8,Tj=4,Tk=1:
    """

    # First, let's compute the cost of an 8x8 sub-matmul
    def cost_8x8_submatmul(A_base, B_base, C_base, scratch_end=46):
        """Estimate cost of 8x8 matmul using Ti=8,Tj=4,Tk=1 with given bulk positions."""
        n = 8
        Ti, Tj = 8, 4
        nbi = n // Ti  # 1
        nbj = n // Tj  # 2
        nbk = n        # 8

        # Scratchpad at fixed positions
        sB_base = 2
        sA_base = sB_base + Tj   # 6
        sC_base = sA_base + Ti   # 14
        # (scratch_end = 46 fixed)

        sA = lambda ii: sA_base + ii
        sB = lambda jj: sB_base + jj
        sC = lambda ii, jj: sC_base + ii * Tj + jj

        A_at = lambda i, j: A_base + i * n + j
        B_at = lambda i, j: B_base + i * n + j
        C_at = lambda i, j: C_base + i * n + j

        sA_cost = sum(addr_cost(sA(ii)) * Tj * nbi * nbj * nbk for ii in range(Ti))
        sB_cost = sum(addr_cost(sB(jj)) * Ti * nbi * nbj * nbk for jj in range(Tj))
        sC_cost = sum(addr_cost(sC(ii,jj)) * nbk * nbi * nbj for ii in range(Ti) for jj in range(Tj))
        tmp_cost = addr_cost(1) * (nbk - 1) * Ti * Tj * nbi * nbj
        bulk_A_cost = sum(addr_cost(A_at(i,j)) * nbj for i in range(n) for j in range(n))
        bulk_B_cost = sum(addr_cost(B_at(i,j)) * nbi for i in range(n) for j in range(n))
        bulk_C_cost = sum(addr_cost(C_at(i,j)) for i in range(n) for j in range(n))

        return sA_cost + sB_cost + sC_cost + tmp_cost + bulk_A_cost + bulk_B_cost + bulk_C_cost

    # Standard 16x16 with 8 sub-matmuls (no Strassen):
    # C00 = A00@B00 + A01@B10
    # C01 = A00@B01 + A01@B11
    # C10 = A10@B00 + A11@B10
    # C11 = A10@B01 + A11@B11
    # Each A_ij is at rows [0..7 or 8..15], cols [0..7 or 8..15]
    # With bulk A at 46-301 (row-major): A[i,j] at 46+i*16+j
    # A00 = A[0..7, 0..7]: address range within 46-301
    # A01 = A[0..7, 8..15]: same rows, different cols

    # For the 8x8 sub-matmul of A00 (rows 0..7, cols 0..7):
    # A00[i,j] = A_full[i,j] at addr 46+i*16+j for i in [0..7], j in [0..7]
    # But this is NOT contiguous (there are 8 columns, then 8 more in the row).

    # This is the key issue: the 8x8 sub-blocks of A,B are not contiguous in memory.
    # To use them as input to sub-matmuls, we'd need to either:
    # a) Store A,B in special layouts
    # b) Copy them to contiguous scratchpad regions first

    # If we store the 8x8 sub-blocks in scratchpad (64 cells each):
    # We need 4 A sub-blocks and 4 B sub-blocks = 8 * 64 = 512 cells
    # This is way too expensive.

    # Conclusion: Strassen at 2x2 block level requires too much scratchpad overhead.
    print("Block-Strassen (2x2 at 8x8 level) analysis:")
    print("  Requires storing 8x8 block sums in scratchpad (512+ cells)")
    print("  Would push all other scratchpad to expensive addresses")
    print("  Estimated to be much worse than current best.")


def generate_tk1_interleaved_load(Ti: int, Tj: int) -> str | None:
    """
    Try loading sA and sB interleaved within each bk iteration,
    then computing the partial outer product before loading next bk's A,B.
    The hope: better data reuse allows eliminating tmp reads somehow.

    Actually: the current code already does load sA, load sB, then compute.
    The order doesn't change reads.
    """
    return None  # No benefit


def generate_best_alternative() -> str:
    """
    Based on all analysis, try a few more configurations that might help:
    1. Ti=8, Tj=4, Tk=1 with different outer loop order that reduces bulk loads
    2. Loop: bi > bk > bj (A cached across bj, but no partial C needed since Tk=1)
       Wait: bk is nbk=16. For bi>bk>bj: load A once per (bi,bk), use for all bj.
       Total A loads: nbi * nbk = 2 * 16 = 32 (vs current 2 * 4 * 16 = 128).
       But we still need sC to accumulate over bk, and for different bj, we need different sC.
       Unless: we DON'T use sC, and write directly to bulk C for each (bi,bj).
       For each (bi,bj), C[bi*Ti:(bi+1)*Ti, bj*Tj:(bj+1)*Tj] = sum_{bk} sA[bk]*sB[bk]^T
       We can accumulate directly in bulk C without an sC scratchpad!
       For bj=0, bk=0: write product to C (no need to read C)
       For bj=0, bk=1: read C, add, write C
       etc.
    """
    # This is the bi>bk>bj with direct C update approach
    # We tried this in exp_acache.py and it was WORSE (147,365)
    # The reason: bulk C reads at addr 558-813 cost 24-29 each
    # Each C cell is read nbk-1 = 15 times to add partial products.
    # 256 * 15 * 29 = 111,360 just for C reads -- much worse!
    return None


def try_asymmetric_ab_storage():
    """
    Since A is read more than B, try putting A at lower addresses by making
    the B array smaller or finding other tricks.

    NOVEL IDEA: What if A occupies the same addresses as B, using overlapping?
    No -- addresses must be distinct.

    What if we use a DIFFERENT memory region for outputs?
    Currently inputs = A_addrs + B_addrs, outputs = C_addrs.
    We declare A,B,C in non-overlapping regions.

    What if we put C at ADDR 1 (and tmp elsewhere)?
    C needs 256 cells. If C starts at 1: C at 1-256, tmp at 257.
    Output cost: sum(addr_cost(1+i) for i in range(256)) = ~1,752 (very cheap!)
    vs current sum(addr_cost(558+i) for i in range(256)) = ~6,819
    Savings: 5,067

    But: tmp at 257 costs 17 per read. Total tmp reads = 3840. tmp cost = 65,280!
    vs current tmp at 1: tmp cost = 3,840.
    Lost: 61,440. Net: MUCH WORSE.

    What if we use tmp=2 (cost 2)?
    Then C starts at 3, A,B bulk starts at C_end=259.
    C at 3-258 (cost 2-17, avg ~9): output cost = 9 * 256 = 2,304
    tmp at 2 (cost 2): tmp reads = 3840 * 2 = 7,680
    A at 259-514, B at 515-770 (cost 17-28)
    This changes many costs. Let me compute.
    """
    print("\n=== C at low addresses experiment ===")

    def gen_low_C(Ti, Tj, C_start):
        """Put C at addr C_start, tmp at C_start-1."""
        if N % Ti != 0 or N % Tj != 0:
            return None
        nbi = N // Ti
        nbj = N // Tj
        nbk = N

        tmp_addr = C_start - 1

        # sB at C_start+256, then sA, then sC
        sB_base = C_start + N * N
        sA_base = sB_base + Tj
        sC_base = sA_base + Ti
        scratch_end = sC_base + Ti * Tj

        sA = lambda ii: sA_base + ii
        sB = lambda jj: sB_base + jj
        sC = lambda ii, jj: sC_base + ii * Tj + jj

        A_base = scratch_end
        B_base = A_base + N * N
        C_base = C_start  # C is at low addresses

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
                                lines.append(f"mul {tmp_addr},{sA(ii)},{sB(jj)}")
                                lines.append(f"add {sC(ii,jj)},{tmp_addr}")
                for ii in range(Ti):
                    for jj in range(Tj):
                        lines.append(f"copy {C_at(bi*Ti+ii, bj*Tj+jj)},{sC(ii,jj)}")

        lines.append(",".join(map(str, outputs)))
        return "\n".join(lines)

    best_cost = BEST_SO_FAR
    best_ir = None
    for Ti, Tj in [(8, 4), (4, 4)]:
        for C_start in [2, 3, 5, 10]:
            ir = gen_low_C(Ti, Tj, C_start)
            if ir:
                try:
                    cost = score_16x16(ir)
                    delta = RECORD - cost
                    marker = " *** BEATS BEST!" if cost < best_cost else ""
                    print(f"  Ti={Ti:>2} Tj={Tj:>2} C_start={C_start:>3}  "
                          f"cost={cost:>10,}  delta={delta:>+8,}{marker}")
                    if cost < best_cost:
                        best_cost = cost
                        best_ir = ir
                except ValueError:
                    pass

    return best_cost, best_ir


def try_mixed_precision_tiles():
    """
    Idea: use different Ti for different rows.
    Process rows 0-7 with Ti=8 and rows 8-15 with Ti=8.
    This is the same as Ti=8 overall.

    What about ASYMMETRIC tiling: first part with smaller sC?
    E.g., process 4 rows of A at a time but accumulate for ALL 16 cols of B.
    Ti=4, Tj=16, Tk=1: nbi=4, nbj=1, nbk=16
    sC: 64 cells... too big.
    """
    pass


def try_sC_column_optimized():
    """
    Instead of row-major sC layout, use column-major.
    sC[ii,jj] = sC_base + jj*Ti + ii (column-major)
    This doesn't change addresses but changes which cells are read together.
    In our loop (ii outer, jj inner): sC[ii,jj] is read in ii-first order.
    Column-major: addr = sC_base + jj*Ti + ii -- higher-address cells read first.
    Row-major: addr = sC_base + ii*Tj + jj -- sequential reads.
    The cost depends on the ADDRESS, not the order of reads.
    Column-major doesn't change total cost.
    """
    pass


if __name__ == "__main__":
    print(f"Record to beat: {RECORD:,}")
    print(f"Best so far:    {BEST_SO_FAR:,}")
    print()

    generate_block_strassen_8x8()

    best_cost, best_ir = try_asymmetric_ab_storage()

    print(f"\nFinal best: {best_cost:,}  (record: {RECORD:,}  delta={RECORD-best_cost:+,})")
