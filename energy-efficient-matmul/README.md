# energy-efficient-matmul

DeepMind's AlphaTensor https://github.com/google-deepmind/alphatensor discover a better 4x4 matrix multiplication algorithm in terms of arithmetic operations. 

What is best algorithm when we care about *energy* rather than *FLOP count*?

To estimate energy, use simplified version of Bill Dally's proposed *Parallel Explicit Communication
Model* [cybertronai/simplified-dally-model](https://github.com/cybertronai/simplified-dally-model).


## API

```python
from matmul import (
    score_1x1, score_4x4, score_16x16,
    generate_baseline_4x4, generate_baseline_16x16,
    generate_tiled_16x16,
)

# Verify your IR computes A @ B correctly and return its read-cost.
cost = matmul.score_1x1("1,2;mul 3,1,2;3")    # 5  (1+2 + 2)
cost = matmul.score_4x4(my_ir_text)
cost = matmul.score_16x16(my_ir_text)
```

The scorer parses the IR, simulates it against deterministic test
matrices `A[i,j] = i·n + j + 1`, `B[i,j] = j·n + i + 1`, verifies the
exit values match `A @ B`, and returns the total read cost. It raises
`ValueError` on a parse error, an uninitialized read, or a wrong
output.

The input convention is fixed: line 1 lists `2n²` addresses — the `n²`
elements of A row-major, then the `n²` elements of B row-major. The
exit line lists the `n²` addresses of C in row-major order.

### Baselines

```python
ir = generate_baseline_4x4()      # naive triple loop, 4×4
ir = generate_baseline_16x16()    # naive triple loop, 16×16
ir = generate_tiled_16x16()       # 4×4 scratchpad-cached tiles
```

## Current best results

| matrices | algorithm                           | cost     |
|----------|-------------------------------------|---------:|
| 4×4      | `generate_baseline_4x4` (naive)     |   1,316  |
| 16×16    | `generate_baseline_16x16` (naive)   | 340,704  |
| 16×16    | `generate_tiled_16x16` (4×4 tiles)  | 133,783  |

**Tiled wins by 2.55×** despite issuing 32 % more instructions —
the extra `mov` traffic loading 4×4 A/B tiles into addresses 1..32 is
more than paid back by the inner `mul`/`add` reads now hitting
distance-1..6 cells instead of distance-15..23 cells. This is the
energy-vs-FLOPs trade-off the problem is meant to expose.

Open: can you do better than 133,783 on 16×16? Strassen and the
AlphaTensor/AlphaEvolve schedules trade muls for extra adds. Under
this cost model the relative price of an add (2 reads at whatever
distances the schedule chose) versus a mul (also 2 reads) is no
longer 0:1 — every read counts. Submit a `generate_<your_method>()`
returning IR that beats the tiled baseline.

## Files

| File | Purpose |
|------|---------|
| `matmul.py` | Scorer, parser, simulator, and the three baseline generators (`generate_baseline_4x4`, `generate_baseline_16x16`, `generate_tiled_16x16`). |
| `README.md` | This file. |

## Notes

- The cost model has **no cache**. Every read pays its full Manhattan
  distance every time. Tiling helps because `mov` lets you copy a far
  cell into a near cell once and re-read the near cell many times.
- The IR is straight-line. Loops and branches at IR-emission time are
  fine (the baseline generators write out the unrolled instruction
  stream); the simulator just executes the resulting tape.
- Constants are not first-class. The baselines avoid `0` by using the
  first product of an inner sum to initialize the accumulator (no
  `add tmp, tmp, ZERO`). If your method needs literal constants, the
  natural extension is a `const dest, #literal` pseudo-op (zero cost,
  matches LLVM's immediate-operand convention).
