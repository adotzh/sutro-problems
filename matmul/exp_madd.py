#!/usr/bin/env python3
"""
Experiment: hypothetical madd instruction for tiled 16x16 matmul

Hypothesis: adding a fused multiply-add instruction (madd dest, src1, src2 →
dest += src1*src2, reads src1+src2+dest, no tmp) eliminates all 4096 tmp reads
and allows sA/sB/sC to shift back to addresses 1-48, reducing total cost
significantly vs opt1 (110,743).

madd semantics:
  madd dest, src1, src2  →  mem[dest] = mem[dest] + mem[src1] * mem[src2]
  cost: _cost(src1) + _cost(src2) + _cost(dest)   [dest is read for accumulation]
  note: for first iteration we still use plain mul (no prior value in dest)

Layout (no tmp slot needed):
  sA at 1..16   (was 2..17 in opt1)
  sB at 17..32  (was 18..33)
  sC at 33..48  (was 34..49)
  A bulk: 49..304   (was 50..305)
  B bulk: 305..560  (was 306..561)
  C bulk: 561..816  (was 562..817)

Usage:
    cd ~/dev/research/sutro/sutro-problems
    python3 matmul/exp_madd.py
"""

import sys, math
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from matmul import _cost, _check_addrs, _matmul_test, _BINARY


# ---------------------------------------------------------------------------
# Extended scorer with madd support
# ---------------------------------------------------------------------------

def _simulate_madd(ir: str, inputs):
    text = ir.replace(";", "\n")
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    input_addrs  = [int(x) for x in lines[0].split(",")]
    output_addrs = [int(x) for x in lines[-1].split(",")]
    _check_addrs(input_addrs,  "input line")
    _check_addrs(output_addrs, "output line")
    ops = []
    for ln in lines[1:-1]:
        head, _, rest = ln.partition(" ")
        operands = [int(x) for x in rest.split(",")]
        _check_addrs(operands, f"`{head}` operands")
        ops.append((head, operands))

    mem = {a: v for a, v in zip(input_addrs, inputs)}
    cost = 0
    for op, oprs in ops:
        if op == "copy":
            dest, src = oprs
            cost += _cost(src)
            mem[dest] = mem[src]
        elif op == "madd":
            # madd dest, src1, src2: dest = dest + src1 * src2
            dest, s1, s2 = oprs
            cost += _cost(dest) + _cost(s1) + _cost(s2)
            mem[dest] = mem[dest] + mem[s1] * mem[s2]
        elif op == "mul":
            dest, s1, s2 = oprs
            cost += _cost(s1) + _cost(s2)
            mem[dest] = mem[s1] * mem[s2]
        elif op in _BINARY:
            if len(oprs) == 3:
                dest, s1, s2 = oprs
            else:
                dest, s2 = oprs; s1 = dest
            cost += _cost(s1) + _cost(s2)
            mem[dest] = _BINARY[op](mem[s1], mem[s2])
        else:
            raise ValueError(f"unknown op: {op!r}")
    outputs = []
    for a in output_addrs:
        cost += _cost(a)
        outputs.append(mem[a])
    return outputs, cost


def score_madd_16x16(ir: str) -> int:
    inputs, expected = _matmul_test(16)
    actual, cost = _simulate_madd(ir, inputs)
    if actual != expected:
        raise ValueError(f"correctness failed:\n  got {actual}\n  exp {expected}")
    return cost


# ---------------------------------------------------------------------------
# IR generator using madd
# ---------------------------------------------------------------------------

def generate_tiled_16x16_madd() -> str:
    """Tiled 16x16 using madd — no tmp register.

    Layout:
      sA at 1..16   sB at 17..32   sC at 33..48
      A bulk: 49..304   B bulk: 305..560   C bulk: 561..816
    """
    n, T = 16, 4
    sA = lambda ii, kk: 1 + ii * T + kk
    sB = lambda kk, jj: 1 + T * T + kk * T + jj
    sC = lambda ii, jj: 1 + 2 * T * T + ii * T + jj

    A_base = 1 + 3 * T * T   # 49
    B_base = A_base + n * n
    C_base = B_base + n * n
    A_at = lambda i, j: A_base + i * n + j
    B_at = lambda i, j: B_base + i * n + j
    C_at = lambda i, j: C_base + i * n + j

    inputs  = ([A_at(i, j) for i in range(n) for j in range(n)] +
               [B_at(i, j) for i in range(n) for j in range(n)])
    outputs = [C_at(i, j) for i in range(n) for j in range(n)]

    lines = [",".join(map(str, inputs))]
    nb = n // T
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
                            if bk == 0 and kk == 0:
                                lines.append(f"mul {sC(ii,jj)},{sA(ii,kk)},{sB(kk,jj)}")
                            else:
                                lines.append(f"madd {sC(ii,jj)},{sA(ii,kk)},{sB(kk,jj)}")
            for ii in range(T):
                for jj in range(T):
                    lines.append(f"copy {C_at(bi*T+ii, bj*T+jj)},{sC(ii,jj)}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


if __name__ == "__main__":
    from matmul import score_16x16, generate_tiled_16x16
    from matmul.exp_layout_opt import generate_tiled_16x16_opt1

    baseline  = score_16x16(generate_tiled_16x16())
    opt1      = score_16x16(generate_tiled_16x16_opt1())
    madd_ir   = generate_tiled_16x16_madd()
    madd_cost = score_madd_16x16(madd_ir)

    print(f"Results:")
    print(f"  tiled_16x16 (baseline):  {baseline:>10,}")
    print(f"  opt1 (tmp@1):            {opt1:>10,}")
    print(f"  madd (no tmp):           {madd_cost:>10,}")
    print(f"  vs baseline improvement: {(baseline - madd_cost)/baseline*100:.1f}%")
    print(f"  vs opt1 improvement:     {(opt1 - madd_cost)/opt1*100:.1f}%")

    ir_path = Path(__file__).parent / "ir" / "tiled_16x16_madd.ir"
    ir_path.parent.mkdir(exist_ok=True)
    ir_path.write_text(madd_ir + "\n")
    print(f"\n  Saved: {ir_path}")
