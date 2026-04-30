#!/usr/bin/env python3
"""
Targeted experiment: reduce sC reads.

For Ti=8, Tj=4, Tk=1 (best: 82,477):
  sC: 23,808 (28.9% of total)
  sC cells: 32, each read 128 times, at addr 14-45 (cost 4-7)

The ONLY way to reduce sC cost:
1. Reduce number of reads per sC cell = reduce nbk = use larger Tk
   But larger Tk makes sA+sB bigger -> pushes sC to higher addrs
2. Use fewer sC cells (smaller Ti*Tj) but same total work
   = use smaller tile and more tiles, same 3840 total reads
   But smaller tiles have smaller sA,sB, so sC starts at lower addr

The fundamental question: what's the OPTIMAL sC start address?
sC_start = tmp_addr_space + sA_size + sB_size + 1

For Ti=8, Tj=4: sC_start = 1 + 4 + 8 + 1 = 14 (plus tmp at 1)
For Ti=4, Tj=4: sC_start = 1 + 4 + 4 + 1 = 10

For Ti=4,Tj=4: sC_start=10, 16 cells at cost 4-5 avg 4.4, 240 reads/cell
  Total = 16 * 240 * 4.4 = 16,896 (approximately)

For Ti=8,Tj=4: sC_start=14, 32 cells at cost 4-7 avg 5.7, 128 reads/cell
  Total = 32 * 128 * 5.7 = 23,347 (approximately)

But wait: for Ti=4,Tj=4, sC reads/cell = nbk*nbi*nbj = 16*4*4 = 256
  vs Ti=8,Tj=4: 16*2*4 = 128

So Ti=4,Tj=4 has DOUBLE the reads/cell but at lower cost!
Ti=4,Tj=4: 16 * 256 * 4.4 = 18,022
Ti=8,Tj=4: 32 * 128 * 5.7 = 23,347

Ti=4,Tj=4 has lower sC cost! But Ti=8,Tj=4 beats it overall because:
- Bulk_B savings: B read 2x vs 4x -- 10,310 cheaper
- This outweighs sC penalty (23,347-18,022=5,325) and sA penalty.

Can we COMBINE these advantages?
- Use Ti=8, Tj=4 for the tiling (keeps B read 2x)
- BUT have a smaller sC by using a multi-pass approach?

MULTI-PASS IDEA:
Process the (bi,bj) tile in TWO PASSES:
Pass 1: Only process ii=0..3 (half the rows), accumulate into sC_half (16 cells)
Pass 2: Process ii=4..7 (other half), accumulate into sC_half (16 cells)

Layout:
  sB_base = 2 (4 cells)
  sA_half_base = 6 (4 cells, only half of sA)
  sC_half_base = 10 (16 cells at cost 4-5)
  scratch_end = 26

This is equivalent to Ti=4, Tj=4 but with a different OUTER loop structure
that processes bi in TWO halves. The bulk A start is at addr 26 instead of 46.

For BOTH halves combined:
sA_half (4 cells): read nbj*nbk per half = 4*16 per (bi_half, bj) = 64 * 4 (halves) * 4 (bj) = 1024 reads
Hmm, this is Ti=4,Tj=4 tiling with nbi=4. Same as Ti=4,Tj=4!

Wait but with Ti=8,Tj=4 the A_base=46. With Ti=4,Tj=4 the A_base=26.
So we could start A at 26 (lower cost) if we use Ti=4,Tj=4.
But then B is at 282-537 (read 4x) vs Ti=8,Tj=4 B at 302-557 (read 2x).
The B cost DOUBLES for Ti=4,Tj=4 despite B being at lower addresses.

This is the fundamental trade-off we've already analyzed.

NEW IDEA: What if we use TWO different tile strategies for different parts of the matrix?

For rows bi=0 (i=0..7): B rows need to be read for ALL bj (4 times).
For rows bi=1 (i=8..15): B rows need to be read again.

What if we process bi=0 and bi=1 in PARALLEL (interleaved)?
For each (bj, bk): load B once (shared between bi=0 and bi=1), load A for both bi=0 and bi=1.
This reduces B reads to 1x!

But: we need TWO sC tiles simultaneously (one for bi=0, one for bi=1).
sC_0 for bi=0: 32 cells
sC_1 for bi=1: 32 cells
Total sC: 64 cells

Layout:
  sB: 4 cells at addr 2-5
  sA0: 8 cells at addr 6-13 (A rows for bi=0)
  sA1: 8 cells at addr 14-21 (A rows for bi=1)
  sC0: 32 cells at addr 22-53
  sC1: 32 cells at addr 54-85
  scratch_end = 86

sB (4 cells): read Ti*2 * nbj * nbk = 16 * 4 * 16 = 1024 times at cost 2-3 avg 2.5 = 2,560
sA0 (8 cells): read Tj * nbj * nbk = 4 * 4 * 16 = 256 times at cost 3-4 avg 3.5 = 896
sA1 (8 cells): read 256 times at cost 4 = 1,024
sC0 (32 cells): read nbk*nbj = 16*4 = 64 times at cost 5-8 avg 6.5 = 13,312
sC1 (32 cells): read 64 times at cost 8-10 avg 9 = 18,432
tmp (1 cell): read (nbk-1)*Ti*Tj*nbj = 15*8*4*4 = 1920 times (for bi=0)
              + 1920 times (for bi=1) = 3840 times at cost 1 = 3,840

Hmm this is getting complex. The sC1 is at addr 54-85 (cost 8-10) -- WORSE.

What if we use a SMALLER Tj for each sA/sC pair?
sB: 2 cells (Tj=2)
sA0: 8, sA1: 8
sC0: 16 (Ti=8, Tj=2)
sC1: 16

But then nbj=8 (since Tj=2 now), which means MORE bj iterations.

This is really a parameterization choice and might not help.

Let me think about this more fundamentally.

For the CURRENT BEST (Ti=8, Tj=4, Tk=1):
- The TOTAL work is fixed: 16^3 = 4096 multiply-accumulate ops
- Each involves reading one sA cell and one sB cell (in scratchpad) plus updating one sC cell
- The addresses of these determines cost

The total read cost is:
  Sum over all ops of addr_cost(src1) + addr_cost(src2)

For inner loop (mul tmp/sC, sA, sB):
  - First bk (direct): reads sA(ii) + sB(jj)
  - Other bk: reads sA(ii) + sB(jj) + reads tmp + reads sC (for add)

Total reads:
- sA reads: sum over (bi,bj,bk,ii,jj) of addr_cost(sA(ii)) = N^3 * avg(sA_cost) / Ti
  Wait: for each (bi,bj,bk): Ti*Tj operations, each reads sA(ii) once.
  sA(ii) is read Tj times per (bi,bj,bk) block.
  Total sA reads = Tj * nbi * nbj * nbk = Tj * (N/Ti) * (N/Tj) * N = N^2 reads
  = 256 reads per sA cell.

- sB reads: similarly 256 reads per sB cell.

- sC reads: for each non-first bk op: add sC(ii,jj), tmp reads sC(ii,jj)
  Total per (ii,jj) per (bi,bj) = nbk-1 = 15 reads
  Plus 1 for writeback = 16 reads per (bi,bj)
  Total per cell = 16 * nbi * nbj = 16 * (N/Ti) * (N/Tj)

- tmp reads: for each non-first bk op: add sC reads tmp
  Total = (nbk-1) * Ti * Tj * nbi * nbj = 15 * Ti * Tj * N^2/(Ti*Tj) = 15 * N^2 = 61,440? No.
  = 15 * Ti * Tj * (N/Ti) * (N/Tj) = 15 * N^2 = too many reads.
  Wait: (nbk-1) = 15 per (ii,jj) per (bi,bj) = 15 * Ti*Tj * nbi*nbj
  For Ti=8,Tj=4,nbi=2,nbj=4: 15*32*8 = 3,840. Yes, 3840 tmp reads total.

The TOTAL READS are fixed at:
  sA reads: 256 total per sA cell * Ti cells = 256*Ti total
  sB reads: 256 per sB cell * Tj cells = 256*Tj total
  sC reads: 16*(N/Ti)*(N/Tj) per sC cell * Ti*Tj cells = 16*N^2 = 4096 total
  tmp reads: 15*N^2/(N/Tk) * Ti*Tj * (N/Ti)*(N/Tj) ... = 15 * Ti*Tj * (N/Ti)*(N/Tj) = 15*N^2? No.

  Let me just count for Ti=8,Tj=4,Tk=1:
  total sC reads = 32 * 128 = 4096 (as above)

  ALWAYS: sC total reads = Ti*Tj * 16*(N/Ti)*(N/Tj) = 16*N^2/N * N... = 16*N^2? No.
  = Ti*Tj * 16 * (N/Ti) * (N/Tj) = Ti*Tj * 16 * N^2 / (Ti*Tj) = 16*N^2 = 4096. YES!

  So total sC reads = 4096 ALWAYS, regardless of Ti, Tj!
  The COST changes because addresses change.

So the question is: what addresses minimize total read cost?

We want to minimize:
  256*Ti * avg_cost(sA) + 256*Tj * avg_cost(sB) + 4096 * avg_cost(sC)
  + 3840 * cost(tmp=1) + N^2 * nbj * avg_cost(A_bulk) + N^2 * nbi * avg_cost(B_bulk)
  + N^2 * avg_cost(C_bulk)

where:
  sB is at addr 2 to Tj+1 (smallest first)
  sA is at addr Tj+2 to Tj+Ti+1
  sC is at addr Tj+Ti+2 to Tj+Ti+Ti*Tj+1
  A_bulk starts at addr Tj+Ti+Ti*Tj+2

The sA+sB part:
  256*Ti * avg(sA_addrs) + 256*Tj * avg(sB_addrs)
  = 256 * (Ti * avg(sA) + Tj * avg(sB))

If we put sB first then sA:
  sB: addr 2..Tj+1, sA: addr Tj+2..Tj+Ti+1
  avg(sB) = mean of addr_cost(2..Tj+1)
  avg(sA) = mean of addr_cost(Tj+2..Tj+Ti+1)

If we put sA first then sB:
  sA: addr 2..Ti+1, sB: addr Ti+2..Ti+Tj+1

Which is better?
  Ti > Tj typically (8 > 4), so sA has more cells.
  Putting smaller array first keeps it at lower addresses AND leaves the larger array at HIGHER addrs.

Wait, for sA first (Ti=8):
  sA at 2-9 (cost 2-3-3-3-3-3-4-4): sum_cost = 2+2+2+3+3+3+4+4 = 23
  sB at 10-13 (cost 4-4-4-4): sum_cost = 16
  Total: 256*Ti*23/Ti + 256*Tj*16/Tj = 256*(23/8 * 8 + 16/4 * 4)...

OK let me just compute directly.
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


def compute_sa_sb_cost(Ti, Tj, sA_first=False):
    """Compute scratchpad component costs for sA before sB (sA_first=True) or sB first."""
    if sA_first:
        sA_base, sB_base = 2, 2 + Ti
    else:
        sB_base, sA_base = 2, 2 + Tj

    nbi = N // Ti
    nbj = N // Tj
    nbk = N
    sC_base = (2 + Ti + Tj) if not sA_first else (2 + Ti + Tj)
    # sC always starts after both sA and sB, regardless of order
    # Actually: sC starts at 2 + Ti + Tj (same total size)

    sA_cost = sum(addr_cost(sA_base + i) * 256 for i in range(Ti))
    sB_cost = sum(addr_cost(sB_base + i) * 256 for i in range(Tj))
    sC_cost = sum(addr_cost(sC_base + i) * (16 * nbi * nbj) for i in range(Ti * Tj))

    return sA_cost + sB_cost, sC_cost, sA_base, sB_base, sC_base


def generate_tk1_optimal(Ti, Tj, slot_order=(1, 0, 2)):
    """Generate Tk=1 with given slot order."""
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


def exhaustive_search_with_asymmetric_scratch():
    """
    Try ALL Ti, Tj (divisors of 16) with ALL slot orderings.
    Also try Tk != 1 for completeness.
    This is a complete search of the parametric space.
    """
    print("=== Exhaustive search: all Ti,Tj,Tk with all slot orders ===")
    results = []
    slot_perms = list(itertools.permutations(range(3)))
    tried = 0
    best_cost = BEST_SO_FAR

    divisors = [d for d in range(1, N+1) if N % d == 0]

    for Ti in divisors:
        for Tj in divisors:
            for Tk in divisors:
                # Compute scratch size
                scratch_size = Ti * Tk + Tk * Tj + Ti * Tj
                if scratch_size > 128:
                    continue  # Too large

                for slot in slot_perms:
                    # Generate the IR
                    ir = generate_nonsquare_full(Ti, Tj, Tk, slot)
                    if ir is None:
                        continue
                    tried += 1
                    try:
                        cost = score_16x16(ir)
                        results.append((cost, Ti, Tj, Tk, slot))
                        if cost < best_cost:
                            best_cost = cost
                            print(f"  *** NEW BEST! Ti={Ti} Tj={Tj} Tk={Tk} slot={slot}  "
                                  f"cost={cost:,}  delta={RECORD-cost:+,}")
                    except ValueError:
                        pass

    results.sort()
    print(f"\nTried {tried} configs. Best: {best_cost:,}  delta={RECORD-best_cost:+,}")
    print("\nTop 10:")
    for cost, Ti, Tj, Tk, slot in results[:10]:
        print(f"  Ti={Ti:>2} Tj={Tj:>2} Tk={Tk:>2} slot={slot}  cost={cost:,}  delta={RECORD-cost:+,}")

    return results, best_cost


def generate_nonsquare_full(Ti, Tj, Tk, slot_order=(0, 1, 2)):
    """Full non-square generator with arbitrary Tk and slot order."""
    if N % Ti != 0 or N % Tj != 0 or N % Tk != 0:
        return None

    nbi = N // Ti
    nbj = N // Tj
    nbk = N // Tk

    tmp = 1
    sizes = [Ti * Tk, Tk * Tj, Ti * Tj]  # sA, sB, sC
    order = sorted(range(3), key=lambda s: slot_order[s])
    bases = [0, 0, 0]
    cur = 2
    for s in order:
        bases[s] = cur
        cur += sizes[s]

    sA_base, sB_base, sC_base = bases
    scratch_end = cur

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

    # Run exhaustive search
    results, best_cost = exhaustive_search_with_asymmetric_scratch()

    if best_cost < BEST_SO_FAR and results:
        # Save the best IR
        cost_found = results[0][0]
        Ti, Tj, Tk, slot = results[0][1], results[0][2], results[0][3], results[0][4]
        ir = generate_nonsquare_full(Ti, Tj, Tk, slot)
        if ir:
            out_path = Path(__file__).parent / "ir" / f"new_record_{cost_found}.ir"
            out_path.parent.mkdir(exist_ok=True)
            out_path.write_text(ir + "\n")
            print(f"\nSaved: {out_path}")

    print(f"\nFinal: best={best_cost:,}  record={RECORD:,}  delta={RECORD-best_cost:+,}")
