You are an optimization agent working on the Sutro Group matmul problem.

## The problem

Find the cheapest IR program that computes C = A @ B for 16×16 matrices
under the simplified Dally cost model:

  cost(read from address addr) = ceil(sqrt(addr))
  Writes are FREE. Arithmetic (add/mul/sub) is FREE. Only READS cost.

The key strategy: put values you read many times at LOW addresses.

## Files

Working directory: /Users/azhiboedova/PetProjects/SutroGroup/sutro-problems

- matmul/matmul.py     — scorer (score_16x16, score_4x4, _cost, _simulate)
- matmul/search.py     — existing parametrized search (already run to completion)
- matmul/exp_explore.py — cost breakdown utilities
- matmul/ir/           — saved IR files

## IR format

Line 1: comma-separated input addresses
Middle: one instruction per line (add/sub/mul/copy)
Last:   comma-separated output addresses

Instructions:
  mul dest, src1, src2  — mem[dest] = mem[src1] * mem[src2]; reads src1+src2
  add dest, src1, src2  — mem[dest] = mem[src1] + mem[src2]; reads src1+src2
  add dest, src         — in-place: mem[dest] += mem[src]; reads dest+src
  copy dest, src        — mem[dest] = mem[src]; reads src (1 read)
  exit line             — each output address pays one read at its cost

Inputs: A flattened row-major (256 values), then B row-major (256 values).
Outputs: C flattened row-major (256 values).
All addresses must be positive integers (>= 1).

## Current record to beat

110,487 — achieved by:
  T=4 tiled matmul (4 tiles of 4×4 per axis)
  Scratchpad layout (tmp@1, sA@2-17, sB@18-33, sC@34-49, bulk@50+)
  Inner loop: mul directly into sC on first product (no copy-through-tmp)

To score your IR: `python3 -c "import sys; sys.path.insert(0,'matmul'); from matmul import score_16x16; print(score_16x16(open('matmul/ir/candidate.ir').read()))"`

## What has already been tried (do NOT repeat)

1. Exhaustive sweep of T in {1,2,4,8} — T=4 wins
2. All 6 permutations of outer loop order (bi,bj,bk) — no difference (symmetric)
3. All 6 permutations of inner loop order (ii,jj,kk) — no difference (symmetric)
4. All 6 permutations of scratchpad slot assignment (sA/sB/sC) — no difference (symmetric)
5. direct_first=True/False — True saves exactly 256, already in the record

The symmetry analysis shows: within the T=4 tiling structure with fixed
bulk layout, all loop/slot permutations are equivalent because sA, sB, and
sC each have exactly 4,096 reads. Improvements must come from changing the
ALGORITHM, not just rearranging the parameters.

## Directions to explore (pick one, go deep, then try another)

**A. Strassen at the tile level**
   Strassen multiplies 2×2 matrices in 7 scalar multiplications instead of 8.
   Since mul is FREE in our model but READS are expensive, Strassen trades
   multiplications (free) for additions (also free) — the win comes from
   fewer reads of sA/sB per output element. Implement Strassen for the 2×2
   sub-blocks within each 4×4 tile computation.

**B. Two-level hierarchical tiling**
   Currently: 16×16 split into 4×4-cell tiles.
   Idea: add a second tiling level inside each 4×4 operation.
   Use a 2×2 inner scratchpad (4 cells at addr 2-5, cost 2) inside each
   4×4 tile step. More copies but each computation reads from cheaper addresses.

**C. Non-square / asymmetric tiles**
   Try tiles like 4×1 (one full column of A at a time), 1×4, or 2×8.
   Non-square tiles can change the read-count distribution across scratchpad cells.
   Implement a general generate(Ti, Tj, Tk) generator where Ti!=Tj!=Tk.

**D. Reduce bulk redundancy**
   In the current bi>bj>bk loop, A[bi,bk] is copied to sA once per bj
   iteration (4x more than necessary). Each A cell is read 4x from bulk.
   Fix: cache A tiles across bj iterations by storing multiple A tiles
   simultaneously in scratchpad. The extra scratchpad cells cost more per
   read but eliminate 3/4 of the bulk A reads (which cost 8-24 each).
   Calculate the break-even point and implement it.

**E. Reduce sC reads with a different accumulation strategy**
   Currently sC is read on every add (3,840 times at cost 6-7).
   Alternative: accumulate in a dedicated register at addr 1 and only
   write to sC at the end of each (ii,jj) computation.
   Or: unroll the kk loop to reduce the number of sC reads.

## Your loop

Iterate until you find a cost strictly below 110,487. For each idea:
1. Implement the generator in a new exp_*.py file
2. Score it with score_16x16()
3. Print the cost and the delta vs 110,487
4. If it is worse, analyze WHY (use cost_breakdown from exp_explore.py)
   and adjust before moving to the next idea
5. If it beats the record, save the IR to matmul/ir/ and stop

Keep a running log of every attempt and cost in your output so the human
can follow your reasoning.
