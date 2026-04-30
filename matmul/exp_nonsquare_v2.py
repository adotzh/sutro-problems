#!/usr/bin/env python3
"""
Direction A v2: Deeper exploration of non-square tiles.

The v1 sweep found Ti=4, Tj=4, Tk=1 gives cost=84,627 (beats record by 25,860!).

Key insight: Tk=1 means the A-tile is just a column vector (Ti×1) and B-tile is
a row vector (1×Tj). The sA and sB scratchpads are tiny (4 cells each), so sC
at addr 10-25 is MUCH cheaper than with T=4 (where sC was at addr 34-49).

This v2 script:
1. Focuses on all Tk=1 configs and does detailed analysis
2. Tries alternative scratchpad layouts to put sC at even lower addrs
3. Tries different inner loop orderings for Tk=1 configs
4. Tries non-standard scratch ordering (sC first, then sA/sB)
"""

import sys, math
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from matmul import score_16x16, _parse

RECORD = 110_487
BEST_SO_FAR = 84_627  # from v1
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


def generate_nonsquare_custom(
    Ti: int, Tj: int, Tk: int,
    slot_order: tuple = (0, 1, 2),   # permutation: 0=sA, 1=sB, 2=sC
    inner_order: tuple = (0, 1, 2),  # permutation of (ii, jj, kk)
    outer_order: tuple = (0, 1, 2),  # permutation of (bi, bj, bk)
    direct_first: bool = True,
) -> str | None:
    """
    Generate tiled 16x16 IR with non-square tiles and configurable loop/layout order.
    slot_order[s] = which slot index (0, 1, or 2) to assign to sA, sB, sC.
    """
    if N % Ti != 0 or N % Tj != 0 or N % Tk != 0:
        return None

    nbi = N // Ti
    nbj = N // Tj
    nbk = N // Tk

    tmp = 1

    # Each slot has a different size: sA=Ti*Tk, sB=Tk*Tj, sC=Ti*Tj
    # slot_order[0]=position of sA, slot_order[1]=position of sB, slot_order[2]=position of sC
    sizes = [Ti * Tk, Tk * Tj, Ti * Tj]

    # Build slot base addresses
    # slot_order[s] = position (0,1,2) for array s (s=0:sA, s=1:sB, s=2:sC)
    # Sort arrays by position to get order
    order = sorted(range(3), key=lambda s: slot_order[s])
    bases = [0, 0, 0]
    cur = 2  # start after tmp
    for s in order:
        bases[s] = cur
        cur += sizes[s]

    sA_base, sB_base, sC_base = bases[0], bases[1], bases[2]
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

    # Map outer loop
    outer_ranges = [range(nbi), range(nbj), range(nbk)]  # 0=bi, 1=bj, 2=bk

    for b0 in outer_ranges[outer_order[0]]:
        for b1 in outer_ranges[outer_order[1]]:
            for b2 in outer_ranges[outer_order[2]]:
                bvals = [b0, b1, b2]
                # outer_order[0] is the variable in position 0 of bvals
                # We need: outer_order[i] tells us which variable (bi/bj/bk) runs at nesting level i
                # So if outer_order = (0,1,2): bi is outermost, bj middle, bk inner
                # bvals[i] = value of the variable at level i
                # outer_order[i] tells us which variable (0=bi,1=bj,2=bk) is at level i
                bi_val = bvals[list(outer_order).index(0)]
                bj_val = bvals[list(outer_order).index(1)]
                bk_val = bvals[list(outer_order).index(2)]

                # Load A-tile
                for ii in range(Ti):
                    for kk in range(Tk):
                        lines.append(
                            f"copy {sA(ii,kk)},{A_at(bi_val*Ti+ii, bk_val*Tk+kk)}")
                # Load B-tile
                for kk in range(Tk):
                    for jj in range(Tj):
                        lines.append(
                            f"copy {sB(kk,jj)},{B_at(bk_val*Tk+kk, bj_val*Tj+jj)}")

                # Inner contraction
                inner_ranges = [range(Ti), range(Tj), range(Tk)]  # 0=ii, 1=jj, 2=kk
                first_set = set()

                for i0 in inner_ranges[inner_order[0]]:
                    for i1 in inner_ranges[inner_order[1]]:
                        for i2 in inner_ranges[inner_order[2]]:
                            ivals = [i0, i1, i2]
                            ii_val = ivals[list(inner_order).index(0)]
                            jj_val = ivals[list(inner_order).index(1)]
                            kk_val = ivals[list(inner_order).index(2)]

                            is_first = (bk_val == 0) and ((ii_val, jj_val) not in first_set)
                            if is_first:
                                first_set.add((ii_val, jj_val))

                            if is_first and direct_first:
                                lines.append(
                                    f"mul {sC(ii_val,jj_val)},{sA(ii_val,kk_val)},{sB(kk_val,jj_val)}")
                            else:
                                lines.append(
                                    f"mul {tmp},{sA(ii_val,kk_val)},{sB(kk_val,jj_val)}")
                                if is_first:
                                    lines.append(f"copy {sC(ii_val,jj_val)},{tmp}")
                                else:
                                    lines.append(f"add {sC(ii_val,jj_val)},{tmp}")

            # Write sC to bulk C (after all bk for this bi,bj)
            # (outer_order tells us which loop variable corresponds to bi, bj, bk)
            # When we exit the b2 loop, we've finished all bk for current (b0,b1)
            # This write happens after exiting the innermost b2 loop
            for ii in range(Ti):
                for jj in range(Tj):
                    # bi_val and bj_val are from the outer two loops
                    lines.append(
                        f"copy {C_at(bi_val*Ti+ii, bj_val*Tj+jj)},{sC(ii,jj)}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


def generate_tk1_optimized(Ti: int, Tj: int, slot_order=(0,1,2)) -> str | None:
    """
    Specialized generator for Tk=1 case.
    With Tk=1:
    - sA is just Ti×1 = Ti cells (column of A)
    - sB is just 1×Tj = Tj cells (row of B)
    - sC is Ti×Tj cells (full C tile, 16 cells for Ti=Tj=4)
    - The kk loop has range(1) so it's trivial
    - Each (bi,bj,bk) block just does one outer product: sC += sA * sB^T
    """
    if N % Ti != 0 or N % Tj != 0:
        return None

    Tk = 1
    nbi = N // Ti
    nbj = N // Tj
    nbk = N  # since Tk=1, nbk=16

    tmp = 1

    # slot_order: position of (sA, sB, sC) in scratchpad
    sizes = [Ti, Tj, Ti * Tj]  # sA=Ti*1, sB=1*Tj, sC=Ti*Tj
    order = sorted(range(3), key=lambda s: slot_order[s])
    bases = [0, 0, 0]
    cur = 2
    for s in order:
        bases[s] = cur
        cur += sizes[s]

    sA_base, sB_base, sC_base = bases[0], bases[1], bases[2]
    scratch_end = cur

    sA = lambda ii: sA_base + ii          # sA[ii] = col of A block
    sB = lambda jj: sB_base + jj          # sB[jj] = row of B block
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
            for bk in range(nbk):  # bk goes from 0 to 15
                # Load single column of A: A[bi*Ti..(bi+1)*Ti, bk] -> sA
                for ii in range(Ti):
                    lines.append(f"copy {sA(ii)},{A_at(bi*Ti+ii, bk)}")
                # Load single row of B: B[bk, bj*Tj..(bj+1)*Tj] -> sB
                for jj in range(Tj):
                    lines.append(f"copy {sB(jj)},{B_at(bk, bj*Tj+jj)}")

                # Outer product: sC[ii,jj] += sA[ii] * sB[jj]
                for ii in range(Ti):
                    for jj in range(Tj):
                        is_first = (bk == 0)
                        if is_first:
                            # Direct mul into sC
                            lines.append(f"mul {sC(ii,jj)},{sA(ii)},{sB(jj)}")
                        else:
                            lines.append(f"mul {tmp},{sA(ii)},{sB(jj)}")
                            lines.append(f"add {sC(ii,jj)},{tmp}")

            # Write sC to bulk C
            for ii in range(Ti):
                for jj in range(Tj):
                    lines.append(f"copy {C_at(bi*Ti+ii, bj*Tj+jj)},{sC(ii,jj)}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


import itertools

def sweep_tk1_all_layouts():
    """Try all valid Tk=1 configs with all slot orderings."""
    print("=== Tk=1 configurations with all slot orderings ===")
    best_cost = BEST_SO_FAR
    best_config = None
    best_ir = None
    tried = 0

    slot_perms = list(itertools.permutations(range(3)))

    for Ti in [1, 2, 4, 8, 16]:
        for Tj in [1, 2, 4, 8, 16]:
            if N % Ti != 0 or N % Tj != 0:
                continue
            scratch_size = Ti + Tj + Ti * Tj
            if scratch_size > 64:
                continue

            for slot_order in slot_perms:
                ir = generate_tk1_optimized(Ti, Tj, slot_order=slot_order)
                if ir is None:
                    continue
                tried += 1
                try:
                    cost = score_16x16(ir)
                    delta = RECORD - cost
                    if cost < best_cost:
                        best_cost = cost
                        best_config = (Ti, Tj, slot_order)
                        best_ir = ir
                        print(f"  *** NEW BEST! Ti={Ti} Tj={Tj} slot={slot_order}  cost={cost:,}  delta={delta:+,}")
                    elif cost <= BEST_SO_FAR + 500:
                        print(f"  Ti={Ti} Tj={Tj} slot={slot_order}  cost={cost:,}  delta={delta:+,}")
                except ValueError as e:
                    print(f"  Ti={Ti} Tj={Tj} slot={slot_order}  ERROR: {e}")

    print(f"\nTried {tried} configs. Best: {best_cost:,}  delta={RECORD-best_cost:+,}")
    return best_cost, best_config, best_ir


def analyze_best():
    """Deep analysis of the Ti=4 Tj=4 Tk=1 winner."""
    print("\n=== Deep Analysis: Ti=4, Tj=4, Tk=1 ===")
    ir = generate_tk1_optimized(4, 4)
    cost = score_16x16(ir)
    print(f"Cost: {cost:,}  delta={RECORD-cost:+,}")

    tiers = cost_breakdown(ir)
    total = sum(t["cost"] for t in tiers.values())
    print(f"\nCost breakdown:")
    for c in sorted(tiers):
        t = tiers[c]
        addrs = sorted(t["addrs"])
        pct = t["cost"] / total * 100
        print(f"  cost={c:>2}  addrs={addrs[0]:>3}-{addrs[-1]:>3}  "
              f"reads={t['reads']:>7,}  cost={t['cost']:>8,}  {pct:>5.1f}%")
    print(f"  TOTAL: {total:,}")

    # Layout explanation
    Ti, Tj = 4, 4
    sA_base, sB_base, sC_base = 2, 6, 10
    print(f"\nLayout:")
    print(f"  tmp:  addr 1   (cost 1)")
    print(f"  sA:   addr {sA_base}-{sA_base+Ti-1}   (col of A, costs {[math.isqrt(sA_base+i-1)+1 for i in range(Ti)]})")
    print(f"  sB:   addr {sB_base}-{sB_base+Tj-1}  (row of B, costs {[math.isqrt(sB_base+i-1)+1 for i in range(Tj)]})")
    print(f"  sC:   addr {sC_base}-{sC_base+Ti*Tj-1}  (Ti*Tj acc, costs {[math.isqrt(sC_base+i-1)+1 for i in range(Ti*Tj)]})")


def try_bigger_tiles_tk1():
    """Try non-power-of-2 tile sizes (if N=16 is divisible)."""
    print("\n=== Extended Ti,Tj sweep for Tk=1 (relaxed scratch constraint) ===")
    best_cost = BEST_SO_FAR
    results = []

    # N=16 divisors include 1,2,4,8,16 only
    # But we could try asymmetric: Ti != Tj
    for Ti in [1, 2, 4, 8, 16]:
        for Tj in [1, 2, 4, 8, 16]:
            if N % Ti != 0 or N % Tj != 0:
                continue
            scratch = Ti + Tj + Ti * Tj
            ir = generate_tk1_optimized(Ti, Tj)
            if ir is None:
                continue
            try:
                cost = score_16x16(ir)
                results.append((cost, Ti, Tj))
                delta = RECORD - cost
                print(f"  Ti={Ti:>2} Tj={Tj:>2}  scratch={scratch:>3}  cost={cost:>10,}  delta={delta:>+8,}")
                if cost < best_cost:
                    best_cost = cost
            except ValueError:
                pass

    print(f"\nBest in sweep: {best_cost:,}")
    return results


if __name__ == "__main__":
    print(f"Record to beat: {RECORD:,}")
    print(f"Best so far:    {BEST_SO_FAR:,}")
    print()

    # 1. Analyze the best config
    analyze_best()

    # 2. Try all slot orderings for Tk=1
    print()
    best_cost, best_config, best_ir = sweep_tk1_all_layouts()

    # 3. Extended sweep (relaxed scratch constraint)
    try_bigger_tiles_tk1()

    if best_ir and best_cost < RECORD:
        Ti, Tj, slot = best_config
        out_path = Path(__file__).parent / "ir" / f"new_record_{best_cost}.ir"
        out_path.parent.mkdir(exist_ok=True)
        out_path.write_text(best_ir + "\n")
        print(f"\nSaved: {out_path}")

    print(f"\nFinal best: {best_cost:,}  (record: {RECORD:,}  delta={RECORD-best_cost:+,})")
