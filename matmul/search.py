#!/usr/bin/env python3
"""
Automated search for cheaper 16×16 matmul IR under the Dally cost model.

Explores a parametrized space of tiled algorithms:
  - tile size T
  - outer loop order (bi, bj, bk permutations)
  - inner loop order (ii, jj, kk permutations)
  - scratchpad address permutation (sA/sB/sC slot order)
  - whether to skip tmp on the first product (mul direct into sC)

Run:
    cd sutro-problems
    python3 matmul/search.py
    python3 matmul/search.py --mode hill      # hill-climb from best known
    python3 matmul/search.py --mode exhaustive # try every combination (~1 728)
"""

import sys, math, itertools, random, argparse, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from matmul import score_16x16, _parse

BASELINE = 110_743   # tiled_16x16_opt1 — current record to beat
N = 16


# ── Address cost ──────────────────────────────────────────────────────────────

def addr_cost(addr: int) -> int:
    return math.isqrt(addr - 1) + 1


# ── Parametrized IR generator ─────────────────────────────────────────────────

def generate(
    T: int,
    outer_order: tuple,   # permutation of (0,1,2) → which of (bi,bj,bk) is outermost
    inner_order: tuple,   # permutation of (0,1,2) → which of (ii,jj,kk) is outermost
    slot_order: tuple,    # permutation of (0,1,2) → which slot gets (sA, sB, sC)
    direct_first: bool,   # mul directly into sC on first product (skips copy of tmp)
) -> str | None:
    """Generate an IR string for the given parameters, or None if invalid."""
    if N % T != 0:
        return None
    nb = N // T

    # Scratchpad layout:
    #   slot 0 = T² addrs starting at 2
    #   slot 1 = T² addrs starting at 2 + T²
    #   slot 2 = T² addrs starting at 2 + 2*T²
    # tmp is always at addr 1.
    tmp = 1
    slot_bases = [2 + slot_order[s] * T * T for s in range(3)]
    sA_base, sB_base, sC_base = slot_bases

    sA = lambda ii, kk: sA_base + ii * T + kk
    sB = lambda kk, jj: sB_base + kk * T + jj
    sC = lambda ii, jj: sC_base + ii * T + jj

    A_base = 2 + 3 * T * T        # bulk starts right after scratchpad
    B_base = A_base + N * N
    C_base = B_base + N * N

    A_at = lambda i, j: A_base + i * N + j
    B_at = lambda i, j: B_base + i * N + j
    C_at = lambda i, j: C_base + i * N + j

    inputs  = ([A_at(i, j) for i in range(N) for j in range(N)] +
               [B_at(i, j) for i in range(N) for j in range(N)])
    outputs = [C_at(i, j) for i in range(N) for j in range(N)]

    # Build loop ranges.  outer_order[k]=0 → bi, 1 → bj, 2 → bk
    outer_iters = [range(nb)] * 3
    inner_iters = [range(T)]  * 3

    lines = [",".join(map(str, inputs))]

    # Outer block loops
    for b0 in outer_iters[0]:
        for b1 in outer_iters[1]:
            for b2 in outer_iters[2]:
                bvals = [b0, b1, b2]
                bi = bvals[outer_order.index(0)]
                bj = bvals[outer_order.index(1)]
                bk = bvals[outer_order.index(2)]

                # Load A-tile
                for ii in range(T):
                    for kk in range(T):
                        lines.append(f"copy {sA(ii,kk)},{A_at(bi*T+ii, bk*T+kk)}")
                # Load B-tile
                for kk in range(T):
                    for jj in range(T):
                        lines.append(f"copy {sB(kk,jj)},{B_at(bk*T+kk, bj*T+jj)}")

                # Inner contraction
                first_set: set[tuple] = set()   # (ii,jj) pairs that have been initialized

                for i0 in inner_iters[0]:
                    for i1 in inner_iters[1]:
                        for i2 in inner_iters[2]:
                            ivals = [i0, i1, i2]
                            ii = ivals[inner_order.index(0)]
                            jj = ivals[inner_order.index(1)]
                            kk = ivals[inner_order.index(2)]

                            is_first = (bk == 0) and ((ii, jj) not in first_set)
                            if is_first:
                                first_set.add((ii, jj))

                            if is_first and direct_first:
                                # Write product directly into sC — no tmp read
                                lines.append(f"mul {sC(ii,jj)},{sA(ii,kk)},{sB(kk,jj)}")
                            else:
                                lines.append(f"mul {tmp},{sA(ii,kk)},{sB(kk,jj)}")
                                if is_first:
                                    lines.append(f"copy {sC(ii,jj)},{tmp}")
                                else:
                                    lines.append(f"add {sC(ii,jj)},{tmp}")

            # Write sC tile to bulk C (after all bk for this (bi,bj))
            for ii in range(T):
                for jj in range(T):
                    lines.append(f"copy {C_at(bi*T+ii, bj*T+jj)},{sC(ii,jj)}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


# ── Search strategies ─────────────────────────────────────────────────────────

PERMS3 = list(itertools.permutations(range(3)))

def random_params(T=None):
    return dict(
        T            = T or random.choice([1, 2, 4, 8]),
        outer_order  = random.choice(PERMS3),
        inner_order  = random.choice(PERMS3),
        slot_order   = random.choice(PERMS3),
        direct_first = random.choice([True, False]),
    )


def score_params(p: dict) -> int | None:
    ir = generate(**p)
    if ir is None:
        return None
    try:
        return score_16x16(ir)
    except ValueError:
        return None


def exhaustive_search():
    """Try all combinations in the parameter space (~1728 combos)."""
    print(f"Exhaustive search over parameter space...")
    best_cost = BASELINE
    best_params = None
    total = 0
    start = time.time()

    for T in [1, 2, 4, 8]:
        for outer in PERMS3:
            for inner in PERMS3:
                for slot in PERMS3:
                    for direct in [True, False]:
                        p = dict(T=T, outer_order=outer, inner_order=inner,
                                 slot_order=slot, direct_first=direct)
                        cost = score_params(p)
                        total += 1
                        if cost is not None and cost < best_cost:
                            best_cost = cost
                            best_params = p
                            print(f"  [{total:>5}] NEW BEST {cost:,}  params={p}")

    elapsed = time.time() - start
    print(f"\nDone. {total} combinations in {elapsed:.1f}s.")
    return best_cost, best_params


def hill_climb(n_iter: int = 5000, n_restarts: int = 5):
    """Hill climbing with random restarts."""
    print(f"Hill climbing: {n_iter} iters × {n_restarts} restarts")
    global_best_cost = BASELINE
    global_best_params = None
    start = time.time()

    for restart in range(n_restarts):
        # Start from a random point
        current = random_params()
        current_cost = score_params(current) or 10**9

        local_best_cost = current_cost
        local_best_params = current

        for i in range(n_iter):
            # Mutate one parameter
            candidate = dict(current)
            key = random.choice(list(candidate.keys()))
            if key == "T":
                candidate["T"] = random.choice([1, 2, 4, 8])
            elif key in ("outer_order", "inner_order", "slot_order"):
                candidate[key] = random.choice(PERMS3)
            else:
                candidate["direct_first"] = not candidate["direct_first"]

            cost = score_params(candidate)
            if cost is not None and cost < current_cost:
                current = candidate
                current_cost = cost
                if cost < local_best_cost:
                    local_best_cost = cost
                    local_best_params = candidate

        print(f"  restart {restart+1}: local best = {local_best_cost:,}")
        if local_best_cost < global_best_cost:
            global_best_cost = local_best_cost
            global_best_params = local_best_params
            print(f"    *** NEW GLOBAL BEST: {global_best_cost:,}  {global_best_params}")

    elapsed = time.time() - start
    print(f"\nDone. {elapsed:.1f}s elapsed.")
    return global_best_cost, global_best_params


def random_search(n_iter: int = 10_000):
    """Pure random sampling — good for surveying the landscape."""
    print(f"Random search: {n_iter} iterations")
    best_cost = BASELINE
    best_params = None
    start = time.time()

    for i in range(n_iter):
        p = random_params()
        cost = score_params(p)
        if cost is not None and cost < best_cost:
            best_cost = cost
            best_params = p
            print(f"  [{i:>6}] NEW BEST {cost:,}  params={p}")

        if (i + 1) % 1000 == 0:
            elapsed = time.time() - start
            print(f"  [{i+1:>6}] best so far: {best_cost:,}  ({elapsed:.1f}s)")

    elapsed = time.time() - start
    print(f"\nDone. {n_iter} samples in {elapsed:.1f}s.")
    return best_cost, best_params


# ── Main ──────────────────────────────────────────────────────────────────────

def save_best(cost: int, params: dict):
    if params is None:
        return
    ir = generate(**params)
    out = Path(__file__).parent / "ir" / f"search_best_{cost}.ir"
    out.write_text(ir + "\n")
    print(f"\nSaved: {out}")
    print(f"Best params: {params}")
    print(f"Best cost:   {cost:,}  (vs record {BASELINE:,}, delta {BASELINE-cost:+,})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["exhaustive", "hill", "random"],
                        default="exhaustive")
    parser.add_argument("--iters", type=int, default=5000)
    parser.add_argument("--restarts", type=int, default=5)
    args = parser.parse_args()

    print(f"Current record: {BASELINE:,}")
    print(f"Search mode:    {args.mode}\n")

    if args.mode == "exhaustive":
        cost, params = exhaustive_search()
    elif args.mode == "hill":
        cost, params = hill_climb(n_iter=args.iters, n_restarts=args.restarts)
    else:
        cost, params = random_search(n_iter=args.iters)

    if params and cost < BASELINE:
        save_best(cost, params)
    elif cost >= BASELINE:
        print(f"\nNo improvement found over {BASELINE:,}.")
