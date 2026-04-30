#!/usr/bin/env python3
"""
Exploration script for matmul optimization under the Dally cost model.

Run with:
    cd sutro-problems
    python3 matmul/exp_explore.py

Cost model recap:
  - Reading from address `addr` costs ceil(sqrt(addr))
  - Writes are FREE. Arithmetic (add/mul/sub) is FREE. Only reads cost.
  - Strategy: put values you read often at LOW addresses (scratchpad).
"""

import sys, math
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from matmul import score_16x16, score_4x4, generate_baseline_4x4, generate_baseline_16x16, generate_tiled_16x16, _parse

# ── Import opt1 from exp_layout_opt ──────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent))
from exp_layout_opt import generate_tiled_16x16_opt1


# ── Cost helpers ─────────────────────────────────────────────────────────────

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


def print_breakdown(label: str, ir: str, score: int):
    tiers = cost_breakdown(ir)
    total_cost = sum(t["cost"] for t in tiers.values())
    total_reads = sum(t["reads"] for t in tiers.values())
    print(f"\n{'='*60}")
    print(f"  {label}  —  score: {score:,}")
    print(f"{'='*60}")
    print(f"  {'Tier':>6}  {'Addr range':>16}  {'Reads':>7}  {'Cost':>8}  {'%':>6}")
    print(f"  {'-'*55}")
    for c in sorted(tiers):
        t = tiers[c]
        addrs = sorted(t["addrs"])
        rng = f"{addrs[0]}-{addrs[-1]}"
        pct = t["cost"] / total_cost * 100
        print(f"  cost={c:>2}  {rng:>16}  {t['reads']:>7,}  {t['cost']:>8,}  {pct:>5.1f}%")
    print(f"  {'TOTAL':>6}  {'':>16}  {total_reads:>7,}  {total_cost:>8,}")


# ── Baselines ─────────────────────────────────────────────────────────────────

def run_baselines():
    print("\n── Known baselines ──────────────────────────────────────────")
    results = [
        ("baseline_16x16 (naive)",    generate_baseline_16x16, score_16x16),
        ("tiled_16x16 (4×4 tiles)",   generate_tiled_16x16,    score_16x16),
        ("tiled_16x16_opt1 (tmp@1)",  generate_tiled_16x16_opt1, score_16x16),
    ]
    for label, gen, scorer in results:
        ir = gen()
        cost = scorer(ir)
        print(f"  {label:<38} {cost:>10,}")
    return results


# ── Ideas to try ──────────────────────────────────────────────────────────────
#
# Current best: tiled_16x16_opt1 = 110,743
#
# Cost breakdown shows:
#   - Scratchpad (addr 1-49) = 62% of total cost
#   - sC at addr 34-49 (cost 6-7)  = 22% alone  ← biggest target
#   - sA at addr 2-17 (cost 2-5)   = 12%
#   - sB at addr 18-33 (cost 5-6)  = 28%
#   - Bulk A/B reads (cost 8-24)   = 32%
#
# Opportunities:
#   1. Move sC to lower addresses (currently at 34-49, cost 6-7)
#      → Can't without bumping sA/sB up. Need a different layout.
#
#   2. Avoid redundant A-tile reloads.
#      Current loop: bi > bj > bk — A[bi,bk] is reloaded for every bj.
#      A only depends on (bi,bk), so it's loaded 4× unnecessarily.
#      → Reorder loop to bi > bk > bj, keep partial sC sums in bulk C.
#
#   3. Try different tile sizes (T=2, T=8).
#      Smaller T → cheaper scratchpad addresses, but more tile transitions.
#      Larger T → each bulk cell read once more, but addresses get expensive.
#
#   4. Eliminate tmp entirely.
#      Write mul result directly into sC (dest=sC) on first iteration,
#      then use sC as src on subsequent adds (in-place accumulate).


def generate_tiled_T(T: int) -> str:
    """Tiled 16x16 with configurable tile size T. Same layout logic as opt1."""
    n = 16
    assert n % T == 0, f"T={T} must divide n={n}"
    nb = n // T

    tmp = 1
    sA_base = 2
    sB_base = sA_base + T * T
    sC_base = sB_base + T * T

    A_base = sC_base + T * T
    B_base = A_base + n * n
    C_base = B_base + n * n

    sA = lambda ii, kk: sA_base + ii * T + kk
    sB = lambda kk, jj: sB_base + kk * T + jj
    sC = lambda ii, jj: sC_base + ii * T + jj
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
                for ii in range(T):
                    for jj in range(T):
                        for kk in range(T):
                            lines.append(f"mul {tmp},{sA(ii,kk)},{sB(kk,jj)}")
                            if bk == 0 and kk == 0:
                                lines.append(f"copy {sC(ii,jj)},{tmp}")
                            else:
                                lines.append(f"add {sC(ii,jj)},{tmp}")
            for ii in range(T):
                for jj in range(T):
                    lines.append(f"copy {C_at(bi*T+ii, bj*T+jj)},{sC(ii,jj)}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


def experiment_tile_sizes():
    """Try different tile sizes and compare costs."""
    print("\n── Experiment: tile size sweep ──────────────────────────────")
    print(f"  {'Tile size':>10}  {'sC addr range':>16}  {'sC cost':>8}  {'Total cost':>12}")
    print(f"  {'-'*55}")
    for T in [1, 2, 4, 8, 16]:
        n = 16
        if n % T != 0:
            continue
        nb = n // T
        sA_size = T * T
        sB_size = T * T
        sC_base = 2 + sA_size + sB_size
        sC_top  = sC_base + T * T - 1
        sC_cost_avg = sum(addr_cost(sC_base + i) for i in range(T*T)) / (T*T)

        ir = generate_tiled_T(T)
        cost = score_16x16(ir)
        print(f"  T={T:>2} ({nb}×{nb} blocks)  "
              f"sC@{sC_base}-{sC_top:>3}           "
              f"avg={sC_cost_avg:.1f}   {cost:>10,}")


# ── YOUR EXPERIMENTS GO HERE ─────────────────────────────────────────────────

def my_experiment() -> str:
    """
    Implement your optimized IR generator here.
    Return an IR string that passes score_16x16().

    Ideas:
    - Try eliminating the tmp register: use sC as direct mul destination
      on the first kk iteration (avoids 1 read per first-product).
    - Try reordering the inner ii/jj/kk loops.
    - Try a different address layout (e.g., interleave sA and sB).
    """
    # Placeholder: just return opt1 so the script always runs
    return generate_tiled_16x16_opt1()


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    run_baselines()
    experiment_tile_sizes()

    # Uncomment to see full cost breakdown for a specific IR:
    # ir = generate_tiled_16x16_opt1()
    # print_breakdown("opt1", ir, score_16x16(ir))

    # Run your experiment:
    ir = my_experiment()
    cost = score_16x16(ir)
    best = 110_743
    delta = best - cost
    print(f"\n── my_experiment ────────────────────────────────────────────")
    print(f"  cost:   {cost:,}")
    print(f"  best:   {best:,}")
    print(f"  delta:  {delta:+,}  ({'improvement' if delta > 0 else 'worse'})")
    if cost < best:
        print(f"\n  *** NEW RECORD! Save IR and update README. ***")
        # Uncomment to save:
        # Path("matmul/ir/my_new_record.ir").write_text(ir + "\n")
