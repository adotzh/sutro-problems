#!/usr/bin/env python3
"""
Direction A: Strassen with M-values at very low addresses.

Previous Strassen attempt scored 196,026 because M1..M7 were placed at addr 50+
(cost 8), and those values get read 2-3 times each in the output combination step.

Fix: place M1..M7 at addr 1-7 (cost 1-3). Use a different overall layout.

Strassen 2x2 formulas for C = A @ B where A=[[a,b],[c,d]], B=[[e,f],[g,h]]:
  M1 = (a+d)*(e+h)   -> used in C[0,0] and C[1,1]
  M2 = (c+d)*e       -> used in C[1,0] and C[1,1]
  M3 = a*(f-h)       -> used in C[0,1] and C[1,1]
  M4 = d*(g-e)       -> used in C[0,0] and C[1,0]
  M5 = (a+b)*h       -> used in C[0,0] and C[0,1]
  M6 = (c-a)*(e+f)   -> used in C[1,1]
  M7 = (b-d)*(g+h)   -> used in C[0,0]
  C[0,0] = M1+M4-M5+M7
  C[0,1] = M3+M5
  C[1,0] = M2+M4
  C[1,1] = M1-M2+M3+M6

Tile 16x16 into 2x2 blocks of size 8x8 (nb=8 per axis).
For each (bi, bj) pair of 2x2 blocks:
  Accumulate Strassen products across all 8 bk blocks.

Layout:
  M1..M7   at addr 1-7    (cost 1-3, read 2-3 times each)
  sums/tmp at addr 8-14   (cost 3, used once for add/sub)
  sC[0,0..1,1] at addr 15-18  (cost 4, accumulated across bk)
  sA[0,0..1,1] at addr 19-22  (cost 5, read once per Strassen block)
  sB[0,0..1,1] at addr 23-26  (cost 5, read once per Strassen block)
  bulk A at addr 27..282   (cost 6-17)
  bulk B at addr 283..538  (cost 17-24)
  bulk C at addr 539..794  (output only, 1 exit read each)
"""

import os, sys, math
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from matmul import score_16x16, _parse

BEST = 73_602
N = 16
nb = 8  # 16/2 = 8 blocks per axis


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


def generate_strassen_2x8_v1():
    """
    Strassen for 2x2 block decomposition of 16x16 matmul.
    Block size is 8 (nb=8 blocks of size 2).

    Layout v1:
      M[0..6]  = addr 1..7   (cost 1-3)
      tmp      = addr 8       (cost 3)
      sC[0..3] = addr 9..12  (cost 3-4, 4 cells for 2x2 C block)
      sA[0..3] = addr 13..16 (cost 4)
      sB[0..3] = addr 17..20 (cost 5)
      A bulk   = addr 21..276 (cost 5-17)
      B bulk   = addr 277..532 (cost 17-24)
      C bulk   = addr 533..788 (output)
    """
    M = list(range(1, 8))  # 1..7
    tmp = 8
    sC_base = 9
    sA_base = 13
    sB_base = 17
    scratch_end = 21

    A_base = scratch_end
    B_base = A_base + N * N
    C_base = B_base + N * N

    A_at = lambda i, j: A_base + i * N + j
    B_at = lambda i, j: B_base + i * N + j
    C_at = lambda i, j: C_base + i * N + j

    # 2x2 block offsets: row r, col c in {0,1}
    sA = lambda r, c: sA_base + r * 2 + c
    sB = lambda r, c: sB_base + r * 2 + c
    sC = lambda r, c: sC_base + r * 2 + c

    inputs = ([A_at(i, j) for i in range(N) for j in range(N)] +
              [B_at(i, j) for i in range(N) for j in range(N)])
    outputs = [C_at(i, j) for i in range(N) for j in range(N)]

    lines = [",".join(map(str, inputs))]

    for bi in range(nb):
        for bj in range(nb):
            # Accumulate over all bk blocks
            for bk in range(nb):
                first = (bk == 0)

                # Load A block [bi, bk] (2x2 subblock)
                # A[2*bi : 2*bi+2, 2*bk : 2*bk+2]
                for r in range(2):
                    for c in range(2):
                        lines.append(f"copy {sA(r,c)},{A_at(2*bi+r, 2*bk+c)}")

                # Load B block [bk, bj] (2x2 subblock)
                # B[2*bk : 2*bk+2, 2*bj : 2*bj+2]
                for r in range(2):
                    for c in range(2):
                        lines.append(f"copy {sB(r,c)},{B_at(2*bk+r, 2*bj+c)}")

                # Now compute Strassen M1..M7 using sA[r,c] and sB[r,c]
                # sA: a=sA[0,0], b=sA[0,1], c=sA[1,0], d=sA[1,1]
                # sB: e=sB[0,0], f=sB[0,1], g=sB[1,0], h=sB[1,1]
                a = sA(0, 0); b = sA(0, 1); c = sA(1, 0); d = sA(1, 1)
                e = sB(0, 0); f = sB(0, 1); g = sB(1, 0); h = sB(1, 1)

                # M1 = (a+d)*(e+h) -> used in sC[0,0] and sC[1,1]
                lines.append(f"add {M[0]},{a},{d}")      # M1 = a+d
                lines.append(f"add {tmp},{e},{h}")        # tmp = e+h
                lines.append(f"mul {M[0]},{M[0]},{tmp}")  # M1 = (a+d)*(e+h)

                # M2 = (c+d)*e -> used in sC[1,0] and sC[1,1]
                lines.append(f"add {M[1]},{c},{d}")       # M2 = c+d
                lines.append(f"mul {M[1]},{M[1]},{e}")    # M2 = (c+d)*e

                # M3 = a*(f-h) -> used in sC[0,1] and sC[1,1]
                lines.append(f"sub {M[2]},{f},{h}")       # M3 = f-h
                lines.append(f"mul {M[2]},{a},{M[2]}")    # M3 = a*(f-h)

                # M4 = d*(g-e) -> used in sC[0,0] and sC[1,0]
                lines.append(f"sub {M[3]},{g},{e}")       # M4 = g-e
                lines.append(f"mul {M[3]},{d},{M[3]}")    # M4 = d*(g-e)

                # M5 = (a+b)*h -> used in sC[0,0] and sC[0,1]
                lines.append(f"add {M[4]},{a},{b}")       # M5 = a+b
                lines.append(f"mul {M[4]},{M[4]},{h}")    # M5 = (a+b)*h

                # M6 = (c-a)*(e+f) -> used in sC[1,1]
                lines.append(f"sub {M[5]},{c},{a}")       # M6 = c-a
                lines.append(f"add {tmp},{e},{f}")         # tmp = e+f
                lines.append(f"mul {M[5]},{M[5]},{tmp}")   # M6 = (c-a)*(e+f)

                # M7 = (b-d)*(g+h) -> used in sC[0,0]
                lines.append(f"sub {M[6]},{b},{d}")        # M7 = b-d
                lines.append(f"add {tmp},{g},{h}")          # tmp = g+h
                lines.append(f"mul {M[6]},{M[6]},{tmp}")    # M7 = (b-d)*(g+h)

                # Now combine into sC[r,c]
                # sC[0,0] = M1+M4-M5+M7
                # sC[0,1] = M3+M5
                # sC[1,0] = M2+M4
                # sC[1,1] = M1-M2+M3+M6
                c00 = sC(0, 0); c01 = sC(0, 1); c10 = sC(1, 0); c11 = sC(1, 1)

                if first:
                    # Initialize sC
                    # sC[0,0] = M1+M4-M5+M7
                    lines.append(f"add {c00},{M[0]},{M[3]}")  # c00 = M1+M4
                    lines.append(f"sub {c00},{M[4]}")          # c00 -= M5
                    lines.append(f"add {c00},{M[6]}")          # c00 += M7

                    # sC[0,1] = M3+M5
                    lines.append(f"add {c01},{M[2]},{M[4]}")  # c01 = M3+M5

                    # sC[1,0] = M2+M4
                    lines.append(f"add {c10},{M[1]},{M[3]}")  # c10 = M2+M4

                    # sC[1,1] = M1-M2+M3+M6
                    lines.append(f"sub {c11},{M[0]},{M[1]}")  # c11 = M1-M2
                    lines.append(f"add {c11},{M[2]}")          # c11 += M3
                    lines.append(f"add {c11},{M[5]}")          # c11 += M6
                else:
                    # Accumulate into sC
                    # sC[0,0] += M1+M4-M5+M7
                    lines.append(f"add {tmp},{M[0]},{M[3]}")   # tmp = M1+M4
                    lines.append(f"sub {tmp},{M[4]}")           # tmp -= M5
                    lines.append(f"add {tmp},{M[6]}")           # tmp += M7
                    lines.append(f"add {c00},{tmp}")            # c00 += tmp

                    # sC[0,1] += M3+M5
                    lines.append(f"add {tmp},{M[2]},{M[4]}")   # tmp = M3+M5
                    lines.append(f"add {c01},{tmp}")            # c01 += tmp

                    # sC[1,0] += M2+M4
                    lines.append(f"add {tmp},{M[1]},{M[3]}")   # tmp = M2+M4
                    lines.append(f"add {c10},{tmp}")            # c10 += tmp

                    # sC[1,1] += M1-M2+M3+M6
                    lines.append(f"sub {tmp},{M[0]},{M[1]}")   # tmp = M1-M2
                    lines.append(f"add {tmp},{M[2]}")           # tmp += M3
                    lines.append(f"add {tmp},{M[5]}")           # tmp += M6
                    lines.append(f"add {c11},{tmp}")            # c11 += tmp

            # Write 2x2 C block to bulk C
            for r in range(2):
                for c in range(2):
                    lines.append(f"copy {C_at(2*bi+r, 2*bj+c)},{sC(r,c)}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


def generate_strassen_2x8_v2():
    """
    Version 2: Optimize layout for lower sA/sB costs.

    Layout v2:
      M[0..6]  = addr 1..7   (cost 1-3)
      sC[0..3] = addr 8..11  (cost 3-4)
      tmp      = addr 12      (cost 4)
      sA[0..3] = addr 13..16 (cost 4)
      sB[0..3] = addr 17..20 (cost 5)
      A bulk   = addr 21..276
      B bulk   = addr 277..532
      C bulk   = addr 533..788

    Key: sC moves to addr 8-11 (lower cost = 3-4 vs 3-4 previously).
    tmp at 12 instead of 8.
    """
    M = list(range(1, 8))  # 1..7
    sC_base = 8
    tmp = 12
    sA_base = 13
    sB_base = 17
    scratch_end = 21

    A_base = scratch_end
    B_base = A_base + N * N
    C_base = B_base + N * N

    A_at = lambda i, j: A_base + i * N + j
    B_at = lambda i, j: B_base + i * N + j
    C_at = lambda i, j: C_base + i * N + j

    sA = lambda r, c: sA_base + r * 2 + c
    sB = lambda r, c: sB_base + r * 2 + c
    sC = lambda r, c: sC_base + r * 2 + c

    inputs = ([A_at(i, j) for i in range(N) for j in range(N)] +
              [B_at(i, j) for i in range(N) for j in range(N)])
    outputs = [C_at(i, j) for i in range(N) for j in range(N)]

    lines = [",".join(map(str, inputs))]

    for bi in range(nb):
        for bj in range(nb):
            for bk in range(nb):
                first = (bk == 0)

                for r in range(2):
                    for c in range(2):
                        lines.append(f"copy {sA(r,c)},{A_at(2*bi+r, 2*bk+c)}")
                for r in range(2):
                    for c in range(2):
                        lines.append(f"copy {sB(r,c)},{B_at(2*bk+r, 2*bj+c)}")

                a = sA(0, 0); b = sA(0, 1); c = sA(1, 0); d = sA(1, 1)
                e = sB(0, 0); f = sB(0, 1); g = sB(1, 0); h = sB(1, 1)

                # M1 = (a+d)*(e+h)
                lines.append(f"add {M[0]},{a},{d}")
                lines.append(f"add {tmp},{e},{h}")
                lines.append(f"mul {M[0]},{M[0]},{tmp}")

                # M2 = (c+d)*e
                lines.append(f"add {M[1]},{c},{d}")
                lines.append(f"mul {M[1]},{M[1]},{e}")

                # M3 = a*(f-h)
                lines.append(f"sub {M[2]},{f},{h}")
                lines.append(f"mul {M[2]},{a},{M[2]}")

                # M4 = d*(g-e)
                lines.append(f"sub {M[3]},{g},{e}")
                lines.append(f"mul {M[3]},{d},{M[3]}")

                # M5 = (a+b)*h
                lines.append(f"add {M[4]},{a},{b}")
                lines.append(f"mul {M[4]},{M[4]},{h}")

                # M6 = (c-a)*(e+f)
                lines.append(f"sub {M[5]},{c},{a}")
                lines.append(f"add {tmp},{e},{f}")
                lines.append(f"mul {M[5]},{M[5]},{tmp}")

                # M7 = (b-d)*(g+h)
                lines.append(f"sub {M[6]},{b},{d}")
                lines.append(f"add {tmp},{g},{h}")
                lines.append(f"mul {M[6]},{M[6]},{tmp}")

                c00 = sC(0, 0); c01 = sC(0, 1); c10 = sC(1, 0); c11 = sC(1, 1)

                if first:
                    lines.append(f"add {c00},{M[0]},{M[3]}")
                    lines.append(f"sub {c00},{M[4]}")
                    lines.append(f"add {c00},{M[6]}")
                    lines.append(f"add {c01},{M[2]},{M[4]}")
                    lines.append(f"add {c10},{M[1]},{M[3]}")
                    lines.append(f"sub {c11},{M[0]},{M[1]}")
                    lines.append(f"add {c11},{M[2]}")
                    lines.append(f"add {c11},{M[5]}")
                else:
                    lines.append(f"add {tmp},{M[0]},{M[3]}")
                    lines.append(f"sub {tmp},{M[4]}")
                    lines.append(f"add {tmp},{M[6]}")
                    lines.append(f"add {c00},{tmp}")

                    lines.append(f"add {tmp},{M[2]},{M[4]}")
                    lines.append(f"add {c01},{tmp}")

                    lines.append(f"add {tmp},{M[1]},{M[3]}")
                    lines.append(f"add {c10},{tmp}")

                    lines.append(f"sub {tmp},{M[0]},{M[1]}")
                    lines.append(f"add {tmp},{M[2]}")
                    lines.append(f"add {tmp},{M[5]}")
                    lines.append(f"add {c11},{tmp}")

            for r in range(2):
                for c in range(2):
                    lines.append(f"copy {C_at(2*bi+r, 2*bj+c)},{sC(r,c)}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


def generate_strassen_2x8_v3():
    """
    Version 3: Even more aggressive - sA and sB before sC to reduce sA/sB costs.

    Layout v3:
      M[0..6]  = addr 1..7   (cost 1-3)
      tmp      = addr 8       (cost 3)
      sA[0..3] = addr 9..12  (cost 3-4)   <- cheaper than v1!
      sB[0..3] = addr 13..16 (cost 4)
      sC[0..3] = addr 17..20 (cost 5)
      A bulk   = addr 21..276
      B bulk   = addr 277..532
      C bulk   = addr 533..788

    sA is cheaper (cost 3-4 instead of 4-5), but sC is more expensive (cost 5).
    sC is read less often (128 times per tile) vs sA (512 times per tile).
    Let's try it.
    """
    M = list(range(1, 8))  # 1..7
    tmp = 8
    sA_base = 9
    sB_base = 13
    sC_base = 17
    scratch_end = 21

    A_base = scratch_end
    B_base = A_base + N * N
    C_base = B_base + N * N

    A_at = lambda i, j: A_base + i * N + j
    B_at = lambda i, j: B_base + i * N + j
    C_at = lambda i, j: C_base + i * N + j

    sA = lambda r, c: sA_base + r * 2 + c
    sB = lambda r, c: sB_base + r * 2 + c
    sC = lambda r, c: sC_base + r * 2 + c

    inputs = ([A_at(i, j) for i in range(N) for j in range(N)] +
              [B_at(i, j) for i in range(N) for j in range(N)])
    outputs = [C_at(i, j) for i in range(N) for j in range(N)]

    lines = [",".join(map(str, inputs))]

    for bi in range(nb):
        for bj in range(nb):
            for bk in range(nb):
                first = (bk == 0)

                for r in range(2):
                    for c in range(2):
                        lines.append(f"copy {sA(r,c)},{A_at(2*bi+r, 2*bk+c)}")
                for r in range(2):
                    for c in range(2):
                        lines.append(f"copy {sB(r,c)},{B_at(2*bk+r, 2*bj+c)}")

                a = sA(0, 0); b = sA(0, 1); c_a = sA(1, 0); d = sA(1, 1)
                e = sB(0, 0); f = sB(0, 1); g = sB(1, 0); h = sB(1, 1)

                # M1 = (a+d)*(e+h)
                lines.append(f"add {M[0]},{a},{d}")
                lines.append(f"add {tmp},{e},{h}")
                lines.append(f"mul {M[0]},{M[0]},{tmp}")

                # M2 = (c+d)*e
                lines.append(f"add {M[1]},{c_a},{d}")
                lines.append(f"mul {M[1]},{M[1]},{e}")

                # M3 = a*(f-h)
                lines.append(f"sub {M[2]},{f},{h}")
                lines.append(f"mul {M[2]},{a},{M[2]}")

                # M4 = d*(g-e)
                lines.append(f"sub {M[3]},{g},{e}")
                lines.append(f"mul {M[3]},{d},{M[3]}")

                # M5 = (a+b)*h
                lines.append(f"add {M[4]},{a},{b}")
                lines.append(f"mul {M[4]},{M[4]},{h}")

                # M6 = (c-a)*(e+f)
                lines.append(f"sub {M[5]},{c_a},{a}")
                lines.append(f"add {tmp},{e},{f}")
                lines.append(f"mul {M[5]},{M[5]},{tmp}")

                # M7 = (b-d)*(g+h)
                lines.append(f"sub {M[6]},{b},{d}")
                lines.append(f"add {tmp},{g},{h}")
                lines.append(f"mul {M[6]},{M[6]},{tmp}")

                c00 = sC(0, 0); c01 = sC(0, 1); c10 = sC(1, 0); c11 = sC(1, 1)

                if first:
                    lines.append(f"add {c00},{M[0]},{M[3]}")
                    lines.append(f"sub {c00},{M[4]}")
                    lines.append(f"add {c00},{M[6]}")
                    lines.append(f"add {c01},{M[2]},{M[4]}")
                    lines.append(f"add {c10},{M[1]},{M[3]}")
                    lines.append(f"sub {c11},{M[0]},{M[1]}")
                    lines.append(f"add {c11},{M[2]}")
                    lines.append(f"add {c11},{M[5]}")
                else:
                    lines.append(f"add {tmp},{M[0]},{M[3]}")
                    lines.append(f"sub {tmp},{M[4]}")
                    lines.append(f"add {tmp},{M[6]}")
                    lines.append(f"add {c00},{tmp}")

                    lines.append(f"add {tmp},{M[2]},{M[4]}")
                    lines.append(f"add {c01},{tmp}")

                    lines.append(f"add {tmp},{M[1]},{M[3]}")
                    lines.append(f"add {c10},{tmp}")

                    lines.append(f"sub {tmp},{M[0]},{M[1]}")
                    lines.append(f"add {tmp},{M[2]}")
                    lines.append(f"add {tmp},{M[5]}")
                    lines.append(f"add {c11},{tmp}")

            for r in range(2):
                for c in range(2):
                    lines.append(f"copy {C_at(2*bi+r, 2*bj+c)},{sC(r,c)}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


def print_breakdown(ir: str, label: str):
    tiers = cost_breakdown(ir)
    total = sum(t["cost"] for t in tiers.values())
    print(f"  Cost breakdown for {label}:")
    for c in sorted(tiers):
        t = tiers[c]
        addrs = sorted(t["addrs"])
        pct = t["cost"] / total * 100
        print(f"    cost={c:>2}  addrs={addrs[0]:>3}-{addrs[-1]:>3}  "
              f"reads={t['reads']:>7,}  cost={t['cost']:>8,}  {pct:>5.1f}%")


if __name__ == "__main__":
    # Check for STOP_SIGNAL
    stop_path = Path(__file__).parent / "STOP_SIGNAL"
    if stop_path.exists():
        print("STOP_SIGNAL received — halting.")
        stop_path.unlink()
        sys.exit(0)

    print(f"Direction A: Strassen with M-values at low addresses")
    print(f"Record to beat: {BEST:,}")
    print()

    variants = [
        ("v1: M@1-7, tmp@8, sC@9-12, sA@13-16, sB@17-20", generate_strassen_2x8_v1),
        ("v2: M@1-7, sC@8-11, tmp@12, sA@13-16, sB@17-20", generate_strassen_2x8_v2),
        ("v3: M@1-7, tmp@8, sA@9-12, sB@13-16, sC@17-20", generate_strassen_2x8_v3),
    ]

    best_score = None
    best_ir = None
    best_label = None

    for label, gen in variants:
        try:
            ir = gen()
            score = score_16x16(ir)
            delta = BEST - score
            marker = " *** BEATS RECORD!" if score < BEST else ""
            print(f"  {label}")
            print(f"    score={score:,}  delta={delta:+,}{marker}")
            print_breakdown(ir, label)
            print()
            if best_score is None or score < best_score:
                best_score = score
                best_ir = ir
                best_label = label
        except Exception as e:
            print(f"  {label}: ERROR: {e}")
            print()

    print(f"Best Strassen score: {best_score:,}  ({best_label})")
    print(f"Delta vs record: {BEST - best_score:+,}")

    if best_score is not None and best_score < BEST:
        ir_path = Path(__file__).parent / "ir" / f"new_record_{best_score}.ir"
        ir_path.write_text(best_ir + "\n")
        print(f"*** NEW RECORD! Saved to {ir_path}")
