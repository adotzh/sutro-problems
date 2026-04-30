#!/usr/bin/env python3
"""
Direction A v3: Systematic exploration around the new best (Ti=8, Tj=4, Tk=1, slot=(1,0,2), cost=82,477).

Key insight from v2:
- Ti=8, Tj=4, Tk=1 with slot order sB,sA,sC:
  sB: addr 2-5   (Tj=4 cells, cost 2-3)
  sA: addr 6-13  (Ti=8 cells, cost 3-4)
  sC: addr 14-45 (Ti*Tj=32 cells, cost 4-7)

- vs Ti=4, Tj=4, Tk=1 with slot order sA,sB,sC:
  sA: addr 2-5   (4 cells)
  sB: addr 6-9   (4 cells)
  sC: addr 10-25 (16 cells)

Now:
1. Try ALL Ti,Tj (even non-divisors) by testing all combos with smarter scratchpad
2. Try ALL slot orderings for more Ti,Tj combos
3. Try reordering the outer loop (bi>bj>bk vs other orders)
4. Key question: can we put sC at lower addresses by using Tj<4?
   - If sC has Ti*Tj cells, smaller Tj means smaller sC and lower addresses
   - But we need to read sC more times (more bj iterations)
   - Trade-off: fewer sC cells but more sC reads
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


def generate_tk1_all_options(
    Ti: int, Tj: int,
    slot_order: tuple = (0, 1, 2),   # position of (sA, sB, sC) in scratchpad
    outer_order: tuple = (0, 1, 2),  # which of bi(0),bj(1),bk(2) is outermost..innermost
    direct_first: bool = True,
) -> str | None:
    """
    Tk=1 specialized generator with all layout options.
    slot_order[s] = which scratchpad position for sA(0), sB(1), sC(2).
    outer_order: permutation of (0,1,2) meaning the outer loop nest order (bi=0, bj=1, bk=2).

    Key: the sC writeback happens when bk loop ends (after nbk iterations).
    If bk is not the innermost loop, the logic changes.
    This version only supports outer_order where bk is innermost (outer_order[2]=2).
    """
    if N % Ti != 0 or N % Tj != 0:
        return None

    nbi = N // Ti
    nbj = N // Tj
    nbk = N  # Tk=1

    tmp = 1

    # Sizes of the three scratchpad arrays
    sizes = [Ti, Tj, Ti * Tj]  # sA, sB, sC

    # Place them in the order given by slot_order
    # slot_order[s] = position (0 means first in scratchpad, 2 means last)
    order = sorted(range(3), key=lambda s: slot_order[s])
    bases = [0, 0, 0]
    cur = 2
    for s in order:
        bases[s] = cur
        cur += sizes[s]

    sA_base, sB_base, sC_base = bases[0], bases[1], bases[2]
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

    # Only support bi>bj>bk (outer_order = (0,1,2)) and bj>bi>bk (outer_order=(1,0,2))
    # because sC is accumulated over bk and written back after all bk
    # For bi>bj>bk:
    if outer_order == (0, 1, 2):
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
                                if direct_first:
                                    lines.append(f"mul {sC(ii,jj)},{sA(ii)},{sB(jj)}")
                                else:
                                    lines.append(f"mul {tmp},{sA(ii)},{sB(jj)}")
                                    lines.append(f"copy {sC(ii,jj)},{tmp}")
                            else:
                                lines.append(f"mul {tmp},{sA(ii)},{sB(jj)}")
                                lines.append(f"add {sC(ii,jj)},{tmp}")
                for ii in range(Ti):
                    for jj in range(Tj):
                        lines.append(f"copy {C_at(bi*Ti+ii, bj*Tj+jj)},{sC(ii,jj)}")
    elif outer_order == (1, 0, 2):
        # bj > bi > bk
        for bj in range(nbj):
            for bi in range(nbi):
                for bk in range(nbk):
                    for ii in range(Ti):
                        lines.append(f"copy {sA(ii)},{A_at(bi*Ti+ii, bk)}")
                    for jj in range(Tj):
                        lines.append(f"copy {sB(jj)},{B_at(bk, bj*Tj+jj)}")
                    for ii in range(Ti):
                        for jj in range(Tj):
                            if bk == 0:
                                if direct_first:
                                    lines.append(f"mul {sC(ii,jj)},{sA(ii)},{sB(jj)}")
                                else:
                                    lines.append(f"mul {tmp},{sA(ii)},{sB(jj)}")
                                    lines.append(f"copy {sC(ii,jj)},{tmp}")
                            else:
                                lines.append(f"mul {tmp},{sA(ii)},{sB(jj)}")
                                lines.append(f"add {sC(ii,jj)},{tmp}")
                for ii in range(Ti):
                    for jj in range(Tj):
                        lines.append(f"copy {C_at(bi*Ti+ii, bj*Tj+jj)},{sC(ii,jj)}")
    else:
        return None  # unsupported

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


def sweep_comprehensive():
    """Full systematic sweep: all Ti, Tj, slot orders, outer orders."""
    print("=== Comprehensive Tk=1 sweep ===")
    print(f"Record to beat: {RECORD:,}  Best so far: {BEST_SO_FAR:,}")
    print()

    best_cost = BEST_SO_FAR
    best_config = None
    best_ir = None
    all_results = []
    slot_perms = list(itertools.permutations(range(3)))
    valid_outer = [(0, 1, 2), (1, 0, 2)]  # only bk-inner supported

    for Ti in [1, 2, 4, 8, 16]:
        for Tj in [1, 2, 4, 8, 16]:
            if N % Ti != 0 or N % Tj != 0:
                continue

            for slot_order in slot_perms:
                for outer_order in valid_outer:
                    for direct in [True, False]:
                        ir = generate_tk1_all_options(
                            Ti, Tj, slot_order=slot_order,
                            outer_order=outer_order, direct_first=direct)
                        if ir is None:
                            continue
                        try:
                            cost = score_16x16(ir)
                            all_results.append((cost, Ti, Tj, slot_order, outer_order, direct))
                            if cost < best_cost:
                                best_cost = cost
                                best_config = (Ti, Tj, slot_order, outer_order, direct)
                                best_ir = ir
                                delta = RECORD - cost
                                print(f"  *** NEW BEST! Ti={Ti} Tj={Tj} slot={slot_order} "
                                      f"outer={outer_order} direct={direct}  "
                                      f"cost={cost:,}  delta={delta:+,}")
                        except ValueError:
                            pass

    all_results.sort()
    print(f"\nTop 10 configurations:")
    for cost, Ti, Tj, slot, outer, direct in all_results[:10]:
        delta = RECORD - cost
        print(f"  Ti={Ti:>2} Tj={Tj:>2} slot={slot} outer={outer} direct={direct}  "
              f"cost={cost:>10,}  delta={delta:>+8,}")

    return best_cost, best_config, best_ir


def analyze_layout_impact():
    """Show how different slot orderings change address costs."""
    print("\n=== Layout impact for Ti=8, Tj=4, Tk=1 ===")
    Ti, Tj = 8, 4
    sizes = [Ti, Tj, Ti * Tj]
    names = ['sA', 'sB', 'sC']

    slot_perms = list(itertools.permutations(range(3)))
    results = []

    for slot_order in slot_perms:
        order = sorted(range(3), key=lambda s: slot_order[s])
        bases = [0, 0, 0]
        cur = 2
        for s in order:
            bases[s] = cur
            cur += sizes[s]

        sA_base, sB_base, sC_base = bases

        ir = generate_tk1_all_options(Ti, Tj, slot_order=slot_order)
        if ir is None:
            continue
        try:
            cost = score_16x16(ir)
            results.append((cost, slot_order, sA_base, sB_base, sC_base))
        except ValueError:
            pass

    results.sort()
    for cost, slot, sA_b, sB_b, sC_b in results[:6]:
        delta = RECORD - cost
        print(f"  slot={slot}  sA@{sA_b}-{sA_b+Ti-1}  sB@{sB_b}-{sB_b+Tj-1}  "
              f"sC@{sC_b}-{sC_b+Ti*Tj-1}  cost={cost:>10,}  delta={delta:>+8,}")


def try_mixed_inner_order():
    """Try different inner loop orderings (ii vs jj first)."""
    print("\n=== Inner loop order impact ===")

    def gen_with_inner(Ti, Tj, slot_order, jj_first=False):
        """Generate with jj-outer or ii-outer inner loop."""
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
        sA = lambda ii: sA_base + ii
        sB = lambda jj: sB_base + jj
        sC = lambda ii, jj: sC_base + ii * Tj + jj

        A_base = cur
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
                    if not jj_first:
                        pairs = [(ii, jj) for ii in range(Ti) for jj in range(Tj)]
                    else:
                        pairs = [(ii, jj) for jj in range(Tj) for ii in range(Ti)]
                    for ii, jj in pairs:
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

    best_slot = (1, 0, 2)
    for Ti, Tj in [(8, 4), (4, 4), (4, 8)]:
        for jj_first in [False, True]:
            ir = gen_with_inner(Ti, Tj, best_slot, jj_first=jj_first)
            if ir:
                try:
                    cost = score_16x16(ir)
                    print(f"  Ti={Ti} Tj={Tj} jj_first={jj_first}  cost={cost:,}  delta={RECORD-cost:+,}")
                except ValueError as e:
                    print(f"  Ti={Ti} Tj={Tj} jj_first={jj_first}  ERROR: {e}")


if __name__ == "__main__":
    analyze_layout_impact()
    try_mixed_inner_order()

    best_cost, best_config, best_ir = sweep_comprehensive()

    if best_ir and best_cost < RECORD:
        out_path = Path(__file__).parent / "ir" / f"new_record_{best_cost}.ir"
        out_path.parent.mkdir(exist_ok=True)
        out_path.write_text(best_ir + "\n")
        print(f"\nSaved: {out_path}")

    print(f"\nFinal best: {best_cost:,}  (record: {RECORD:,}  delta={RECORD-best_cost:+,})")
