#!/usr/bin/env python3
"""
Direction C: Strassen's algorithm within tiles.

Strassen computes 2x2 matmul with 7 multiplications instead of 8.
Standard 2x2: C = A @ B requires 8 reads of sA and 8 reads of sB.
Strassen: 7 products M1..M7, each reads sA and sB twice => 14 reads total
  but produces 4 C values from 7 products instead of 8.

Since reads are NOT free, the question is:
  Do we save reads on the sA/sB side, vs the extra overhead of
  computing linear combinations (adds of sA, sB, and M values)?

Strassen's formulas for C = A @ B (2×2 blocks):
  A = [[a00, a01], [a10, a11]]
  B = [[b00, b01], [b10, b11]]

  M1 = (a00 + a11) * (b00 + b11)
  M2 = (a10 + a11) * b00
  M3 = a00 * (b01 - b11)
  M4 = a11 * (b10 - b00)
  M5 = (a00 + a01) * b11
  M6 = (a10 - a00) * (b00 + b01)
  M7 = (a01 - a11) * (b10 + b11)

  C[0,0] = M1 + M4 - M5 + M7
  C[0,1] = M3 + M5
  C[1,0] = M2 + M4
  C[1,1] = M1 - M2 + M3 + M6

Read count (per 2x2 block):
  Standard: 8 reads of A + 8 reads of B = 16 sA/sB reads
  Strassen: M1: 2A + 2B reads (but we precompute a00+a11, b00+b11)
            Actually: 7 product reads * 2 = 14 reads
            But we also read sums: if precomputed, add overhead

The issue: in the standard approach, each of the 4 C values is:
  C[i,j] = sum of 2 products = 4 reads of sA + 4 reads of sB
  For the full 4x4 inner loop: each sA[ii,kk] read Tj*nbk*nbj times

Let me think about this differently for the Tk=1 best case.

For Ti=8, Tj=4, Tk=1:
- Each bk iteration: compute 8*4=32 outer products (sA[ii] * sB[jj])
- That's 32 reads of sA + 32 reads of sB per bk block
- Total: 32 * nbk * nbi * nbj = 32 * 16 * 2 * 4 = 4096 reads each

Strassen doesn't directly apply here because we're doing a rank-1 update (Tk=1 outer product).
Strassen is more relevant for the square T×T case.

Let me implement Strassen for the T=4 square tiling:
- Tile is 4×4
- Within each tile, use 2×2 Strassen sub-blocks
- The 4×4 tile becomes four 2×2 blocks

But for T=4, Tk=4: the 4×4 matmul uses 64 sA+sB reads.
With 2×2 Strassen blocks: (4/2)^3 = 8 Strassen calls, each 7 muls vs 8 = saves 8 muls.
But in our model, mul is FREE. The savings come from READS, not muls.

Standard: each sA[ii,kk] is read Tj times, each sB[kk,jj] read Ti times.
With Strassen 2x2 sub-blocks within 4x4:
- We process 4 output 2x2 blocks using Strassen's 7 products per (bi_block, bj_block, bk_block)
- Each M_i reads 2 sA and 2 sB entries (after adding them together)
- But the additions themselves cost reads too!

Let me implement it concretely for the T=4 case with Tk=1 inner blocks.

Actually, Strassen is complex and the savings are unclear. Let me first try a simpler variant:
use 2x2 Strassen within each Ti=4, Tj=4 outer product block.

For Ti=4, Tj=4, Tk=1: each (bi,bj,bk) does 1 outer product: C += col * row^T.
- Standard: reads sA 4 times (once per ii) and sB 4 times (once per jj) = not quite.
  Actually for each (ii,jj): read sA[ii] and sB[jj] -> 32 reads per block (8*4+8*4... no)
  Wait: 4*4 = 16 (ii,jj) pairs, each reads sA[ii] and sB[jj] -> 16+16 = 32 reads
  But sA[ii] is read 4 times (once per jj), sB[jj] is read 4 times (once per ii).
  With Tj=4: sA has 4 cells, each read 4 times = 16 sA reads.
  With Ti=4... wait, for Ti=8, Tj=4: sA has 8 cells each read 4 times = 32 sA reads,
  sB has 4 cells each read 8 times = 32 sB reads per block.

Strassen on outer products doesn't naturally apply.

Let me try Strassen on the SQUARE T=4 (standard tiling) case:
The Ti=4, Tj=4, Tk=4 square tiling has an inner 4x4x4 matmul per tile.
Apply Strassen recursively to get 4x4 = 2x(2x2 matmul) with Strassen.

Actually let me try a completely different approach for Direction C:
Apply Strassen to the full 16x16 problem directly, breaking it into 8x8 sub-problems.
This is complex. Let me implement it properly.

For simplicity, focus on how Strassen reduces sC READS (since sC is at addr 10+ and costs 4-7).
In the standard approach:
  sC[ii,jj] is read (nbk*Tk - 1) times (for the add accumulations)
  + once for writeback to bulk C
  Total: nbk*Tk reads per sC cell.

With Strassen, if we can reduce the number of distinct (ii,jj) sC accumulation targets,
we'd reduce sC reads. But Strassen's intermediate values M1..M7 need to be stored
somewhere and then combined to form C — adding extra reads.

Let me just implement it and compare:
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


def generate_strassen_2x2_in_4x4(
    sA_locs,   # dict: (ii,kk) -> addr
    sB_locs,   # dict: (kk,jj) -> addr
    sC_locs,   # dict: (ii,jj) -> addr
    M_locs,    # list of 7 addr for M1..M7
    tmp_addr,  # address for tmp computations
    first_write,  # bool: first bk iteration (initialize sC)
) -> list:
    """
    Generate Strassen IR instructions for computing sC += sA @ sB (2x2 case).

    sA indexed as [ii,kk] for ii,kk in {0,1}
    sB indexed as [kk,jj] for kk,jj in {0,1}
    sC indexed as [ii,jj] for ii,jj in {0,1}

    Strassen formulas (using our sA/sB/sC naming):
      M1 = (sA[0,0] + sA[1,1]) * (sB[0,0] + sB[1,1])
      M2 = (sA[1,0] + sA[1,1]) * sB[0,0]
      M3 = sA[0,0] * (sB[0,1] - sB[1,1])
      M4 = sA[1,1] * (sB[1,0] - sB[0,0])
      M5 = (sA[0,0] + sA[0,1]) * sB[1,1]
      M6 = (sA[1,0] - sA[0,0]) * (sB[0,0] + sB[0,1])
      M7 = (sA[0,1] - sA[1,1]) * (sB[1,0] + sB[1,1])

      sC[0,0] += M1 + M4 - M5 + M7
      sC[0,1] += M3 + M5
      sC[1,0] += M2 + M4
      sC[1,1] += M1 - M2 + M3 + M6
    """
    lines = []
    a00 = sA_locs[(0, 0)]
    a01 = sA_locs[(0, 1)]
    a10 = sA_locs[(1, 0)]
    a11 = sA_locs[(1, 1)]
    b00 = sB_locs[(0, 0)]
    b01 = sB_locs[(0, 1)]
    b10 = sB_locs[(1, 0)]
    b11 = sB_locs[(1, 1)]
    M = M_locs  # M[0] = addr for M1, ..., M[6] = addr for M7

    # M1 = (sA[0,0] + sA[1,1]) * (sB[0,0] + sB[1,1])
    # Need a tmp for (a00+a11) and (b00+b11)
    # Use M[0] as temp for (a00+a11), another tmp for (b00+b11)
    # But we only have one tmp address for intermediate calculations
    # Strategy: use M addresses as temporaries during computation

    t = tmp_addr

    # M1: (a00+a11) * (b00+b11)
    lines.append(f"add {M[0]},{a00},{a11}")     # M[0] = a00 + a11  (reads a00,a11)
    lines.append(f"add {t},{b00},{b11}")         # t = b00 + b11     (reads b00,b11)
    lines.append(f"mul {M[0]},{M[0]},{t}")       # M1 = M[0] * t     (reads M[0],t)

    # M2: (a10+a11) * b00
    lines.append(f"add {M[1]},{a10},{a11}")      # M[1] = a10 + a11  (reads a10,a11)
    lines.append(f"mul {M[1]},{M[1]},{b00}")     # M2 = M[1] * b00   (reads M[1],b00)

    # M3: a00 * (b01 - b11)
    lines.append(f"sub {M[2]},{b01},{b11}")      # M[2] = b01 - b11  (reads b01,b11)
    lines.append(f"mul {M[2]},{a00},{M[2]}")     # M3 = a00 * M[2]   (reads a00,M[2])

    # M4: a11 * (b10 - b00)
    lines.append(f"sub {M[3]},{b10},{b00}")      # M[3] = b10 - b00  (reads b10,b00)
    lines.append(f"mul {M[3]},{a11},{M[3]}")     # M4 = a11 * M[3]   (reads a11,M[3])

    # M5: (a00+a01) * b11
    lines.append(f"add {M[4]},{a00},{a01}")      # M[4] = a00 + a01  (reads a00,a01)
    lines.append(f"mul {M[4]},{M[4]},{b11}")     # M5 = M[4] * b11   (reads M[4],b11)

    # M6: (a10-a00) * (b00+b01)
    lines.append(f"sub {M[5]},{a10},{a00}")      # M[5] = a10 - a00  (reads a10,a00)
    lines.append(f"add {t},{b00},{b01}")         # t = b00 + b01     (reads b00,b01)
    lines.append(f"mul {M[5]},{M[5]},{t}")       # M6 = M[5] * t     (reads M[5],t)

    # M7: (a01-a11) * (b10+b11)
    lines.append(f"sub {M[6]},{a01},{a11}")      # M[6] = a01 - a11  (reads a01,a11)
    lines.append(f"add {t},{b10},{b11}")         # t = b10 + b11     (reads b10,b11)
    lines.append(f"mul {M[6]},{M[6]},{t}")       # M7 = M[6] * t     (reads M[6],t)

    # Now combine: sC[ii,jj] += Mxx
    c00 = sC_locs[(0, 0)]
    c01 = sC_locs[(0, 1)]
    c10 = sC_locs[(1, 0)]
    c11 = sC_locs[(1, 1)]

    if first_write:
        # sC[0,0] = M1 + M4 - M5 + M7
        lines.append(f"add {c00},{M[0]},{M[3]}")  # c00 = M1 + M4     (reads M[0],M[3])
        lines.append(f"sub {c00},{M[4]}")          # c00 -= M5         (reads M[4])
        lines.append(f"add {c00},{M[6]}")          # c00 += M7         (reads M[6])

        # sC[0,1] = M3 + M5
        lines.append(f"add {c01},{M[2]},{M[4]}")  # c01 = M3 + M5     (reads M[2],M[4])

        # sC[1,0] = M2 + M4
        lines.append(f"add {c10},{M[1]},{M[3]}")  # c10 = M2 + M4     (reads M[1],M[3])

        # sC[1,1] = M1 - M2 + M3 + M6
        lines.append(f"sub {c11},{M[0]},{M[1]}")  # c11 = M1 - M2     (reads M[0],M[1])
        lines.append(f"add {c11},{M[2]}")          # c11 += M3         (reads M[2])
        lines.append(f"add {c11},{M[5]}")          # c11 += M6         (reads M[5])
    else:
        # Accumulate: sC[ii,jj] += ...
        lines.append(f"add {t},{M[0]},{M[3]}")     # t = M1 + M4       (reads M[0],M[3])
        lines.append(f"sub {t},{M[4]}")            # t -= M5           (reads M[4])
        lines.append(f"add {t},{M[6]}")            # t += M7           (reads M[6])
        lines.append(f"add {c00},{t}")             # c00 += t          (reads t)

        lines.append(f"add {t},{M[2]},{M[4]}")    # t = M3 + M5       (reads M[2],M[4])
        lines.append(f"add {c01},{t}")             # c01 += t          (reads t)

        lines.append(f"add {t},{M[1]},{M[3]}")    # t = M2 + M4       (reads M[1],M[3])
        lines.append(f"add {c10},{t}")             # c10 += t          (reads t)

        lines.append(f"sub {t},{M[0]},{M[1]}")    # t = M1 - M2       (reads M[0],M[1])
        lines.append(f"add {t},{M[2]}")            # t += M3           (reads M[2])
        lines.append(f"add {t},{M[5]}")            # t += M6           (reads M[5])
        lines.append(f"add {c11},{t}")             # c11 += t          (reads t)

    return lines


def generate_strassen_tiled_4x4() -> str:
    """
    Apply Strassen's algorithm within each 4×4 tile.

    The 4×4 tile matmul is decomposed into four 2×2 sub-matmuls.
    Each 2×2 sub-matmul uses Strassen's 7-product algorithm.

    Layout:
      tmp    at 1
      sA     at 2-17   (4×4 = 16 cells, cost 2-5)
      sB     at 18-33  (4×4 = 16 cells, cost 5-6)
      sC     at 34-49  (4×4 = 16 cells, cost 6-7)
      M1..M7 at 50-56  (7 Strassen intermediates, cost 8)
      Bulk A at 57..
    """
    n, T = 16, 4
    nb = n // T

    tmp = 1
    sA_base = 2
    sB_base = sA_base + T * T   # 18
    sC_base = sB_base + T * T   # 34
    M_base  = sC_base + T * T   # 50
    scratch_end = M_base + 7    # 57

    A_base = scratch_end
    B_base = A_base + n * n
    C_base = B_base + n * n

    sA = lambda ii, kk: sA_base + ii * T + kk
    sB = lambda kk, jj: sB_base + kk * T + jj
    sC = lambda ii, jj: sC_base + ii * T + jj
    M = [M_base + m for m in range(7)]

    A_at = lambda i, j: A_base + i * n + j
    B_at = lambda i, j: B_base + i * n + j
    C_at = lambda i, j: C_base + i * n + j

    inputs = ([A_at(i, j) for i in range(n) for j in range(n)] +
              [B_at(i, j) for i in range(n) for j in range(n)])
    outputs = [C_at(i, j) for i in range(n) for j in range(n)]

    lines = [",".join(map(str, inputs))]

    for bi in range(nb):
        for bj in range(nb):
            for bk in range(nb):
                # Load A-tile
                for ii in range(T):
                    for kk in range(T):
                        lines.append(f"copy {sA(ii,kk)},{A_at(bi*T+ii, bk*T+kk)}")
                # Load B-tile
                for kk in range(T):
                    for jj in range(T):
                        lines.append(f"copy {sB(kk,jj)},{B_at(bk*T+kk, bj*T+jj)}")

                # Apply Strassen in four 2x2 sub-blocks
                # The 4x4 tile is split into 4 quadrants of size 2x2:
                # A = [[A00, A01], [A10, A11]] where A_ij is 2x2
                # B = [[B00, B01], [B10, B11]] where B_ij is 2x2
                # C = A @ B
                # Each 2x2 sub-block index (sri, srk) for A, (srk, srj) for B

                # For the outer block structure of the TILE:
                # ii in [0,1] = tile-row offset (2 rows each), jj in [0,1] = tile-col offset
                # kk in [0,1] = tile-inner offset
                # We have 2x2 outer blocks, each is a 2x2 sub-matmul

                for sri in range(2):  # tile row block (0 or 1, each 2 rows)
                    for srj in range(2):  # tile col block
                        # C[sri*2:sri*2+2, srj*2:srj*2+2] = sum_srk A[sri*2:..,srk*2:..] @ B[srk*2:..,srj*2:..]
                        # This is a 2x2 matmul accumulated over srk=0,1
                        for srk in range(2):
                            # Build the sA_locs and sB_locs for this 2x2 sub-block
                            sA_sub = {
                                (r, c): sA(sri * 2 + r, srk * 2 + c)
                                for r in range(2) for c in range(2)
                            }
                            sB_sub = {
                                (r, c): sB(srk * 2 + r, srj * 2 + c)
                                for r in range(2) for c in range(2)
                            }
                            sC_sub = {
                                (r, c): sC(sri * 2 + r, srj * 2 + c)
                                for r in range(2) for c in range(2)
                            }
                            first_write = (bk == 0 and srk == 0)
                            sub_lines = generate_strassen_2x2_in_4x4(
                                sA_sub, sB_sub, sC_sub, M, tmp, first_write)
                            lines.extend(sub_lines)

            # Write sC to bulk C
            for ii in range(T):
                for jj in range(T):
                    lines.append(f"copy {C_at(bi*T+ii, bj*T+jj)},{sC(ii,jj)}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


def generate_strassen_tiled_optimized() -> str:
    """
    Strassen with optimized M placement.

    Put M1..M7 at lowest available addresses.
    Since sA,sB,sC occupy 48 cells (addr 2-49), M goes to 50-56.

    Alternative: use a different layout to put M at lower cost addresses.

    Key observation: M values are read ~16 times each (4 C cells * 4 contributions * nbk*nbi*nbj).

    Wait, actually M values are computed and used within the same bk block.
    Each Mi is computed once and read multiple times to form the C sub-block outputs.

    Let's try putting sC AFTER M, so M gets lower addresses.
    Or use Tk=1 variant with Strassen's 2x2 inner product (outer product Strassen doesn't exist).

    Actually, for Tk=1 (outer product), each element C[ii,jj] = a[ii] * b[jj].
    There's no Strassen variant for rank-1 updates — you can't reduce reads.

    Strassen only helps when Tk > 1 (actual inner dimension to sum over).

    For the square T=4, Tk=4 tiling:
    The bottleneck is sC reads (at addr 34-49, cost 6-7).
    Strassen reduces the number of accumulation steps, potentially reducing sC reads.

    Standard 4x4: 64 inner products, 64-16 = 48 sC reads (one per accumulation) + 16 writes
    Strassen 2x2 sub-blocks: 4 outer blocks * 2 (srk=0,1) * 7 M computations = 56 sub-muls
    But Strassen requires more sC reads for combining the M values.

    Let's just measure.
    """
    return generate_strassen_tiled_4x4()


def generate_strassen_4x4_lowM() -> str:
    """
    Strassen with M1..M7 at low addresses.
    Move sC to higher addresses, put M before sC.

    Layout:
      tmp   at 1       (cost 1)
      M1..M7 at 2-8   (cost 2-3)  <- M values at lowest addresses
      sA    at 9-24   (4x4, cost 3-5)
      sB    at 25-40  (4x4, cost 5-7)
      sC    at 41-56  (4x4, cost 7)
    """
    n, T = 16, 4
    nb = n // T

    tmp = 1
    M_base  = 2
    sA_base = M_base + 7   # 9
    sB_base = sA_base + T * T  # 25
    sC_base = sB_base + T * T  # 41
    scratch_end = sC_base + T * T  # 57

    A_base = scratch_end
    B_base = A_base + n * n
    C_base = B_base + n * n

    sA = lambda ii, kk: sA_base + ii * T + kk
    sB = lambda kk, jj: sB_base + kk * T + jj
    sC = lambda ii, jj: sC_base + ii * T + jj
    M = [M_base + m for m in range(7)]

    A_at = lambda i, j: A_base + i * n + j
    B_at = lambda i, j: B_base + i * n + j
    C_at = lambda i, j: C_base + i * n + j

    inputs = ([A_at(i, j) for i in range(n) for j in range(n)] +
              [B_at(i, j) for i in range(n) for j in range(n)])
    outputs = [C_at(i, j) for i in range(n) for j in range(n)]

    lines = [",".join(map(str, inputs))]

    for bi in range(nb):
        for bj in range(nb):
            for bk in range(nb):
                for ii in range(T):
                    for kk in range(T):
                        lines.append(f"copy {sA(ii,kk)},{A_at(bi*T+ii, bk*T+kk)}")
                for kk in range(T):
                    for jj in range(T):
                        lines.append(f"copy {sB(kk,jj)},{B_at(bk*T+kk, bj*T+jj)}")

                for sri in range(2):
                    for srj in range(2):
                        for srk in range(2):
                            sA_sub = {
                                (r, c): sA(sri * 2 + r, srk * 2 + c)
                                for r in range(2) for c in range(2)
                            }
                            sB_sub = {
                                (r, c): sB(srk * 2 + r, srj * 2 + c)
                                for r in range(2) for c in range(2)
                            }
                            sC_sub = {
                                (r, c): sC(sri * 2 + r, srj * 2 + c)
                                for r in range(2) for c in range(2)
                            }
                            first_write = (bk == 0 and srk == 0)
                            sub_lines = generate_strassen_2x2_in_4x4(
                                sA_sub, sB_sub, sC_sub, M, tmp, first_write)
                            lines.extend(sub_lines)

            for ii in range(T):
                for jj in range(T):
                    lines.append(f"copy {C_at(bi*T+ii, bj*T+jj)},{sC(ii,jj)}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


if __name__ == "__main__":
    print(f"Record to beat: {RECORD:,}")
    print(f"Best so far:    {BEST_SO_FAR:,}")
    print()

    print("=== Strassen in 4×4 tiles ===")

    ir_std = generate_strassen_tiled_4x4()
    try:
        cost_std = score_16x16(ir_std)
        delta = RECORD - cost_std
        print(f"  Strassen (M after sC): cost={cost_std:,}  delta={delta:+,}")
        tiers = cost_breakdown(ir_std)
        total = sum(t["cost"] for t in tiers.values())
        for c in sorted(tiers):
            t = tiers[c]
            addrs = sorted(t["addrs"])
            pct = t["cost"] / total * 100
            print(f"    cost={c:>2}  addrs={addrs[0]:>3}-{addrs[-1]:>3}  "
                  f"reads={t['reads']:>7,}  cost={t['cost']:>8,}  {pct:>5.1f}%")
    except ValueError as e:
        print(f"  ERROR: {e}")

    print()
    ir_low = generate_strassen_4x4_lowM()
    try:
        cost_low = score_16x16(ir_low)
        delta = RECORD - cost_low
        print(f"  Strassen (M before sA): cost={cost_low:,}  delta={delta:+,}")
        tiers = cost_breakdown(ir_low)
        total = sum(t["cost"] for t in tiers.values())
        for c in sorted(tiers):
            t = tiers[c]
            addrs = sorted(t["addrs"])
            pct = t["cost"] / total * 100
            print(f"    cost={c:>2}  addrs={addrs[0]:>3}-{addrs[-1]:>3}  "
                  f"reads={t['reads']:>7,}  cost={t['cost']:>8,}  {pct:>5.1f}%")
    except ValueError as e:
        print(f"  ERROR: {e}")

    print()
    print("=== Comparison with standard T=4 ===")
    # Standard T=4 for reference
    from exp_layout_opt import generate_tiled_16x16_opt1
    ir_opt1 = generate_tiled_16x16_opt1()
    cost_opt1 = score_16x16(ir_opt1)
    print(f"  Standard T=4:            cost={cost_opt1:,}  delta={RECORD-cost_opt1:+,}")
    print(f"  Best non-square Tk=1:    cost={BEST_SO_FAR:,}  delta={RECORD-BEST_SO_FAR:+,}")

    print()
    print("=== Analysis: can we improve Strassen with better M placement? ===")
    # The key is M1..M7 are read multiple times to form sC:
    # In first_write: M read 2-3 times each (non-trivially)
    # In accumulate: M read more times
    # How many times total?
    print("\n  Counting M reads per Strassen call:")
    print("  M1 -> c00, c11 reads: 2 times (first), + ~2 accumulate")
    print("  M2 -> c10, c11 reads: 2 times")
    print("  M3 -> c01, c11 reads: 2 times")
    print("  M4 -> c00, c10 reads: 2 times")
    print("  M5 -> c00, c01 reads: 2 times")
    print("  M6 -> c11: 1 time")
    print("  M7 -> c00: 1 time")
    print("  (These are per 2x2 sub-block call)")
    print("  Total calls per (bi,bj): 4 outer * 2 srk = 8 calls")
    print("  Total calls overall: 8 * nb*nb = 8*16 = 128 calls")
    print()
    print("  If M is at cost-2 addresses (addr 2-8): M reads cost 2 each")
    print("  vs standard tmp at cost-1 (addr 1)")

    if best_cost_found := min([cost for cost in [cost_std if 'cost_std' in dir() else 999999,
                                                   cost_low if 'cost_low' in dir() else 999999]
                                if cost < BEST_SO_FAR], default=None):
        print(f"\n  *** Strassen improvement found: {best_cost_found:,}")
        if best_cost_found == cost_std:
            out_path = Path(__file__).parent / "ir" / f"new_record_{best_cost_found}.ir"
            out_path.write_text(ir_std + "\n")
        else:
            out_path = Path(__file__).parent / "ir" / f"new_record_{best_cost_found}.ir"
            out_path.write_text(ir_low + "\n")
        print(f"  Saved: {out_path}")
    else:
        print(f"\n  Strassen did not beat best ({BEST_SO_FAR:,}).")
