# energy-efficient-matmul

What is the cheapest algorithm to multiply two matrices when *energy*
is the cost function rather than *FLOP count*?

DeepMind's AlphaTensor / AlphaEvolve advertise their matrix-multiplication
algorithms by counting arithmetic operations. On modern hardware
arithmetic is essentially free; the cost is dominated by data movement.
This problem rescores the same algorithms under a model that prices
data movement explicitly and asks for the cheapest correct schedule.

## Cost model

A simplified version of Bill Dally's *Parallel Explicit Communication
Model* — see [cybertronai/simplified-dally-model](https://github.com/cybertronai/simplified-dally-model).

- Processor at the origin; memory is a 2D upper half-plane, linearly
  indexed. Cell `addr` sits at Manhattan distance `⌈√addr⌉` from the
  core.
- **Reads are priced** (one read = `⌈√addr⌉`).
- **Writes are free.**
- **Arithmetic is free.**
- **At the start of a call** the caller specifies where every input
  byte lives; placement is free.
- **At the end of a call** the caller specifies the output addresses;
  each output pays one standard read.

## IR

Three-address SSA-flavored code, `;`-or-newline separated:

| Instruction              | Effect                              | Reads |
|--------------------------|-------------------------------------|------:|
| `add dest, src1, src2`   | `mem[dest] = mem[src1] + mem[src2]` |  2    |
| `sub dest, src1, src2`   | `mem[dest] = mem[src1] - mem[src2]` |  2    |
| `mul dest, src1, src2`   | `mem[dest] = mem[src1] * mem[src2]` |  2    |
| `mov dest, src`          | `mem[dest] = mem[src]`              |  1    |

The 4th op (`mov`) is what makes scratchpad-style tiling possible in
this model — there is no implicit cache, so re-using a value cheaply
requires explicitly copying it to a low address. Two-operand short
form for the binary ops: `add dest, src` is sugar for
`add dest, dest, src`.

A program is the input-placement line, a sequence of instructions, and
the output-read line:

```
1,2                  ← inputs at addrs 1 and 2
mul 3,1,2            ← mem[3] = mem[1] * mem[2];  reads ⌈√1⌉ + ⌈√2⌉ = 1+2
3                    ← exit: read addr 3;        cost ⌈√3⌉ = 2
```

Total cost = `1 + 2 + 2 = 5`.

## Worked example (`myfunc(a,b,c,d,e) = a*b + c*d + e`)

```
1,2,3,4,5            ← a@1, b@2, c@3, d@4, e@5
mul 1,1,2            ← t1 = a*b → addr 1
mul 2,3,4            ← t2 = c*d → addr 2
add 1,1,2            ← s  = t1 + t2 → addr 1
add 1,5              ← r  = s + e   → addr 1     (in-place form)
1                    ← exit: read addr 1
```

| step | reads          | cost  |
|-----:|----------------|------:|
| 1    | `a@1`, `b@2`   | 1+2   |
| 2    | `c@3`, `d@4`   | 2+2   |
| 3    | `t1@1`, `t2@2` | 1+2   |
| 4    | `s@1`, `e@5`   | 1+3   |
| exit | `r@1`          | 1     |

Total: `(1+2) + (2+2) + (1+2) + (1+3) + 1 = 15`.

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

| matrices | algorithm                           | reads | cost     |
|----------|-------------------------------------|------:|---------:|
| 4×4      | `generate_baseline_4x4` (naive)     |   112 |   1,316  |
| 16×16    | `generate_baseline_16x16` (naive)   | 7,936 | 340,704  |
| 16×16    | `generate_tiled_16x16` (4×4 tiles)  |10,496 | 133,783  |

**Tiled wins by 2.55×** despite issuing **32 % more instructions** —
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
