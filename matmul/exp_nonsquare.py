#!/usr/bin/env python3
"""
Direction A: Non-square tile exploration for 16×16 matmul.

Instead of square T×T tiles, try independent tile dimensions:
  Ti = rows of A-tile / C-tile
  Tj = cols of B-tile / C-tile
  Tk = inner/reduction dimension

Layout:
  tmp   at addr 1
  sA    at addr 2 .. (1 + Ti*Tk)            (Ti×Tk cells)
  sB    at addr (2+Ti*Tk) .. (1+Ti*Tk+Tk*Tj) (Tk×Tj cells)
  sC    at addr (2+Ti*Tk+Tk*Tj) .. (1+Ti*Tk+Tk*Tj+Ti*Tj) (Ti×Tj cells)
  Bulk starts right after scratchpad.

Key insight: With non-square tiles we can choose Ti, Tj, Tk independently,
potentially placing sC at lower addresses than the T=4 square case,
or reducing the number of expensive scratchpad reads.
"""

import sys, math
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from matmul import score_16x16, _parse

RECORD = 110_487
N = 16


def addr_cost(addr: int) -> int:
    return math.isqrt(addr - 1) + 1


def cost_breakdown(ir: str) -> dict:
    """Return per-region read stats for an IR string."""
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


def generate_nonsquare(Ti: int, Tj: int, Tk: int) -> str | None:
    """
    Generate tiled 16x16 IR with non-square tile dimensions Ti×Tk for A,
    Tk×Tj for B, Ti×Tj for C.

    Requires N divisible by Ti, Tj, Tk.
    Ti*Tk + Tk*Tj + Ti*Tj <= 64 to keep scratchpad cheap.
    """
    if N % Ti != 0 or N % Tj != 0 or N % Tk != 0:
        return None

    nbi = N // Ti   # number of blocks in i direction
    nbj = N // Tj   # number of blocks in j direction
    nbk = N // Tk   # number of blocks in k direction

    tmp = 1
    # sA: Ti*Tk cells
    sA_base = 2
    # sB: Tk*Tj cells
    sB_base = sA_base + Ti * Tk
    # sC: Ti*Tj cells
    sC_base = sB_base + Tk * Tj

    scratch_end = sC_base + Ti * Tj  # first address after scratchpad

    A_base = scratch_end
    B_base = A_base + N * N
    C_base = B_base + N * N

    sA = lambda ii, kk: sA_base + ii * Tk + kk    # row-major in ii,kk
    sB = lambda kk, jj: sB_base + kk * Tj + jj    # row-major in kk,jj
    sC = lambda ii, jj: sC_base + ii * Tj + jj    # row-major in ii,jj

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
                # Load A-tile: A[bi*Ti..(bi+1)*Ti, bk*Tk..(bk+1)*Tk] -> sA
                for ii in range(Ti):
                    for kk in range(Tk):
                        lines.append(
                            f"copy {sA(ii,kk)},{A_at(bi*Ti+ii, bk*Tk+kk)}")
                # Load B-tile: B[bk*Tk..(bk+1)*Tk, bj*Tj..(bj+1)*Tj] -> sB
                for kk in range(Tk):
                    for jj in range(Tj):
                        lines.append(
                            f"copy {sB(kk,jj)},{B_at(bk*Tk+kk, bj*Tj+jj)}")

                # Inner contraction: sC[ii,jj] += sA[ii,kk] * sB[kk,jj]
                for ii in range(Ti):
                    for jj in range(Tj):
                        for kk in range(Tk):
                            is_first = (bk == 0) and (kk == 0)
                            if is_first:
                                # Direct mul into sC (no tmp read)
                                lines.append(
                                    f"mul {sC(ii,jj)},{sA(ii,kk)},{sB(kk,jj)}")
                            else:
                                lines.append(
                                    f"mul {tmp},{sA(ii,kk)},{sB(kk,jj)}")
                                lines.append(f"add {sC(ii,jj)},{tmp}")

            # Write sC tile to bulk C (after all bk for this (bi,bj))
            for ii in range(Ti):
                for jj in range(Tj):
                    lines.append(
                        f"copy {C_at(bi*Ti+ii, bj*Tj+jj)},{sC(ii,jj)}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


def analyze_config(Ti, Tj, Tk):
    """Compute expected cost analytically."""
    nbi = N // Ti
    nbj = N // Tj
    nbk = N // Tk

    sA_base = 2
    sB_base = sA_base + Ti * Tk
    sC_base = sB_base + Tk * Tj

    # Number of times each region is read
    # sA: read Ti*Tk times per (ii,jj,kk) block = Ti*Tk reads per inner loop
    # Each (bi,bk) block's sA tile is shared across all bj
    # sA load: nbi * nbj * nbk * Ti * Tk reads (one copy per element per block)
    # sA use: nbi * nbj * nbk * Ti * Tj * Tk reads (one per inner iter, 2 reads: sA and sB)
    # Actually sA is read once per (ii, kk, jj) = Ti*Tk*Tj times per block * nbi*nbj*nbk blocks

    # Analytical: sA addr sA(ii,kk) is read Tj times per inner jj loop, once per (bi,bk,bj)
    # total sA reads = nbi * nbj * nbk * Ti * Tk * Tj  ... but that's wrong
    # The inner loop is: for ii in Ti: for jj in Tj: for kk in Tk:
    # sA(ii,kk) is read once per (ii,jj,kk) => Ti*Tj*Tk times per (bi,bj,bk) block
    # Total sA reads = nbi * nbj * nbk * Ti * Tj * Tk (split per addr: Tj*Tk per addr ii)

    # For each sA addr (ii fixed, kk fixed): read Tj*Tk times per (bi,bj,bk) block * nbi*nbj*nbk
    # Actually: per (bi,bj,bk) block: for each ii, sA(ii,kk) is read Tj times (once per jj)
    # per addr sA(ii,kk): Tj reads per block * nbi*nbj*nbk blocks

    sA_reads_per_addr = nbj * nbk * nbi * Tj  # Wait, not quite
    # Correction: per (bi,bj,bk) block, sA(ii,kk) is read Tj times (one for each jj)
    # Total blocks = nbi * nbj * nbk
    # Total reads for sA(ii,kk) = Tj * nbi * nbj * nbk
    sA_total = Ti * Tk * Tj * nbi * nbj * nbk  # total across all sA addrs

    # sB(kk,jj) is read Ti times per block (one for each ii)
    sB_total = Tk * Tj * Ti * nbi * nbj * nbk

    # sC: read for add (not on first product which goes direct)
    # For each (bi,bj): sC(ii,jj) is accumulated nbk*Tk - 1 times (minus first direct mul)
    # Actually: nbk*Tk multiplications total, first one is direct (no sC read), rest are add+sC read
    # So sC reads = (nbk * Tk - 1) reads per (ii,jj) per (bi,bj)
    # Plus sC is read at end for output: Ti*Tj reads per (bi,bj)
    # Also reads of tmp: (nbk*Tk - 1) reads per (ii,jj) per (bi,bj)
    n_sC_reads = (nbk * Tk - 1) * Ti * Tj * nbi * nbj  # from add instructions
    n_sC_copy_reads = Ti * Tj * nbi * nbj  # copy to bulk C

    # tmp reads
    n_tmp_reads = (nbk * Tk - 1) * Ti * Tj * nbi * nbj  # from add instructions

    sA_cost = sum(addr_cost(sA_base + ii * Tk + kk)
                  for ii in range(Ti) for kk in range(Tk)) * nbj * nbk * nbi * Tj
    # Wait: sA(ii,kk) is read Tj times per block, total nbi*nbj*nbk blocks
    # Fix the formula: sA_total = nbi * nbj * nbk, each sA(ii,kk) read Tj times per block
    sA_cost_total = 0
    for ii in range(Ti):
        for kk in range(Tk):
            addr = sA_base + ii * Tk + kk
            sA_cost_total += addr_cost(addr) * Tj * nbi * nbj * nbk

    sB_cost_total = 0
    for kk in range(Tk):
        for jj in range(Tj):
            addr = sB_base + kk * Tj + jj
            sB_cost_total += addr_cost(addr) * Ti * nbi * nbj * nbk

    sC_cost_total = 0
    for ii in range(Ti):
        for jj in range(Tj):
            addr = sC_base + ii * Tj + jj
            # read (nbk*Tk - 1) times for accumulation + 1 for copy to bulk
            sC_cost_total += addr_cost(addr) * ((nbk * Tk - 1) + 1) * nbi * nbj

    tmp_cost_total = addr_cost(1) * (nbk * Tk - 1) * Ti * Tj * nbi * nbj

    return {
        'sA_base': sA_base, 'sB_base': sB_base, 'sC_base': sC_base,
        'scratch_size': Ti*Tk + Tk*Tj + Ti*Tj,
        'sA_cost': sA_cost_total, 'sB_cost': sB_cost_total,
        'sC_cost': sC_cost_total, 'tmp_cost': tmp_cost_total,
        'estimated_total': sA_cost_total + sB_cost_total + sC_cost_total + tmp_cost_total
    }


def sweep_nonsquare():
    """Try all valid non-square tile combinations and score them."""
    results = []
    tried = 0
    best_cost = RECORD
    best_config = None
    best_ir = None

    divisors_16 = [d for d in [1, 2, 4, 8, 16] if N % d == 0]

    print(f"Sweeping non-square tiles (Ti x Tj x Tk) where Ti, Tj, Tk divide {N}")
    print(f"Filtering: scratchpad size Ti*Tk + Tk*Tj + Ti*Tj <= 64")
    print(f"Record to beat: {RECORD:,}")
    print()

    for Ti in divisors_16:
        for Tj in divisors_16:
            for Tk in divisors_16:
                if Ti == Tj == Tk:
                    continue  # skip square (already tried)
                scratch_size = Ti * Tk + Tk * Tj + Ti * Tj
                if scratch_size > 64:
                    continue

                ir = generate_nonsquare(Ti, Tj, Tk)
                if ir is None:
                    continue
                tried += 1
                try:
                    cost = score_16x16(ir)
                    results.append((cost, Ti, Tj, Tk, ir))
                    delta = RECORD - cost
                    status = " *** BEATS RECORD!" if cost < RECORD else ""
                    print(f"  Ti={Ti} Tj={Tj} Tk={Tk}  scratch={scratch_size:>3}  cost={cost:>10,}  delta={delta:>+8,}{status}")
                    if cost < best_cost:
                        best_cost = cost
                        best_config = (Ti, Tj, Tk)
                        best_ir = ir
                except ValueError as e:
                    print(f"  Ti={Ti} Tj={Tj} Tk={Tk}  ERROR: {e}")

    print(f"\n{'='*60}")
    print(f"Tried {tried} non-square configurations")
    if best_config:
        Ti, Tj, Tk = best_config
        print(f"Best non-square: Ti={Ti} Tj={Tj} Tk={Tk}  cost={best_cost:,}  delta={RECORD-best_cost:+,}")
    else:
        print("No non-square config beat the record.")
    return results, best_cost, best_config, best_ir


def analyze_best_configs(results):
    """Detailed breakdown of the top configs."""
    sorted_results = sorted(results, key=lambda x: x[0])[:5]
    print("\nTop 5 configurations:")
    for cost, Ti, Tj, Tk, ir in sorted_results:
        print(f"\n  Ti={Ti} Tj={Tj} Tk={Tk}  cost={cost:,}")
        # Get breakdown
        from exp_explore import cost_breakdown
        tiers = cost_breakdown(ir)
        total_cost = sum(t["cost"] for t in tiers.values())
        for c in sorted(tiers):
            t = tiers[c]
            addrs = sorted(t["addrs"])
            pct = t["cost"] / total_cost * 100
            print(f"    cost={c:>2}  addrs={addrs[0]}-{addrs[-1]}  "
                  f"reads={t['reads']:>7,}  cost={t['cost']:>8,}  {pct:>5.1f}%")


if __name__ == "__main__":
    results, best_cost, best_config, best_ir = sweep_nonsquare()

    if results:
        analyze_best_configs(results)

    # Also try square T=4 as reference
    print("\n--- Reference: square T=4 ---")
    ir4 = generate_nonsquare(4, 4, 4)  # same as square T=4
    if ir4:
        cost4 = score_16x16(ir4)
        print(f"  Square T=4: cost={cost4:,}  delta={RECORD-cost4:+,}")

    if best_ir and best_cost < RECORD:
        Ti, Tj, Tk = best_config
        out_path = Path(__file__).parent / "ir" / f"new_record_{best_cost}.ir"
        out_path.parent.mkdir(exist_ok=True)
        out_path.write_text(best_ir + "\n")
        print(f"\nSaved new record IR to {out_path}")

    print(f"\nDone. Best non-square cost: {best_cost:,}  (record: {RECORD:,})")
