#!/usr/bin/env python3
"""
Direction B: Column streaming (Ti=16, Tj=1, Tk=1).

For each k (0..15):
  1. Copy all 16 values of A[:,k] into sA_col at addr 2-17 (cost 2-5).
     This column is reused for ALL 16 values of j.
  2. For each j (0..15):
     a. Copy B[k,j] to addr 1 (cost 1).
     b. For each i (0..15): mul/add using sA_col[i] and B_cache

  Total: A column read ONCE per k (into sA_col), then reused for all 16 j values.
  B[k,j] is cached at addr 1 for all 16 i values.

Layout:
  B_cache  = addr 1   (cost 1, 4096 reads)
  sA_col   = addr 2-17 (cost 2-5, 16 cells, read N*N=256 times each)
  sC       = addr 18-33 (cost 5-6, 16 cells for one output column)
  A bulk   = addr 34-289 (cost 6-17)
  B bulk   = addr 290-545 (cost 17-24)
  C bulk   = addr 546-801 (output)

Each A[i,k] is read from sA_col[i] for N=16 values of j.
sA_col[i] at addr 2+i (cost 2-5): read 16 times per k = 16*16 = 256 total per cell.
B_cache (addr 1) read 16 times per (k,j) = 16*16*16 = 4096 total.
sC at addr 18+i (cost 5-6): accumulated over N=16 k values per output column.

Variant: Process ROWS of C instead of columns.
For each i (0..15):
  Load A[i,:] into sA_row at addr 2-17.
  For each j (0..15):
    sC[j] accumulated over k.
"""

import os, sys, math
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from matmul import score_16x16, _parse

BEST = 73_602
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


def generate_col_stream_v1():
    """
    Column streaming: for each k, load A[:,k] into sA_col, then for each j
    cache B[k,j] at addr 1 and accumulate sC column.

    Layout:
      B_cache  = 1       (cost 1)
      sA_col   = 2-17   (cost 2-5)
      sC       = 18-33  (cost 5-6, 16 output cells for one i)
      A bulk   = 34-289
      B bulk   = 290-545
      C bulk   = 546-801

    Loop: for k: load sA_col; for j: load B_cache; accumulate sC[i] for each i
    At end of all k: write sC column j to C.

    Wait, this accumulates C[:,j] for a specific j across all k.
    But we process j in the inner loop... let me restructure.

    Better: for each j, compute the full column C[:,j].
      for j in range(N):
        for k in range(N):
          Load B[k,j] into B_cache = addr 1
          Load sA_col = A[:,k] (16 cells) into addr 2-17  -- but this repeats for each j!

    That's expensive. Instead:

    Approach: outer loop is k, inner loop is j.
      sA_col[i] = A[i,k]: loaded once per k, reused for all j.
      B_cache = B[k,j]: loaded once per (k,j), reused for all i.
      sC[j*N+i]: needs to be accessible, but that's N*N cells...

    Actually, let's think about it as: for each output COLUMN j (0..15):
      - Accumulate sC[i] for i=0..15 over k=0..15
      - For each k: load B[k,j] -> B_cache; read A[i,k] for each i
      - At end: write sC[0..15] to C[0..15, j]

    This processes ONE OUTPUT COLUMN at a time.
    For each j:
      for k:
        copy B_cache, B[k,j]         <- 1 bulk B read (expensive)
        for i:
          if k==0: mul sC[i], A[i,k], B_cache   <- read A[i,k] from bulk + read B_cache
          else: mul tmp, A[i,k], B_cache; add sC[i], tmp
      for i:
        copy C[i,j], sC[i]

    A[i,k] is read N times (once per j). Total A reads: N*N*N/N = N^2 = 256.
    But each A[i,k] is read once per j => N reads per A cell. 256 A cells * 16j = 4096 bulk A reads.
    vs current: 1024 bulk A reads (read 4 times each, not 16).

    This is WORSE for A.

    THE KEY INSIGHT: we need sA_col to cache A[:,k] so each A[i,k] is read once from bulk
    (into sA_col) and then re-read from sA_col for each j.

    So structure:
      for k in range(N):
        # Load A[:,k] into sA_col (read each A[i,k] ONCE from bulk)
        for i in range(N):
          copy sA_col[i], A[i,k]
        # Now for each j, use B[k,j] and sA_col
        for j in range(N):
          copy B_cache, B[k,j]   # read B[k,j] once from bulk
          for i in range(N):
            if k==0: mul sC[i,j], sA_col[i], B_cache
            else: mul tmp, sA_col[i], B_cache; add sC[i,j], tmp

    Problem: sC has N*N = 256 cells! That's too expensive.

    Alternative: swap j outer, k inner (so sC only has N cells for one column):
      for j in range(N):
        for k in range(N):
          copy B_cache, B[k,j]     # bulk B read: N*N = 256 total
          for i in range(N):
            # Need to read A[i,k] -- but sA_col not loaded yet
            if k==0: mul sC[i], A[i,k], B_cache   # bulk A read
            else: mul tmp, A[i,k], B_cache; add sC[i], tmp  # bulk A read
        for i:
          copy C[i,j], sC[i]

    This has N^3 = 4096 bulk A reads (each A[i,k] read N times).

    Correct approach: k outer, j inner, sC has N*N cells (or process j in blocks).

    OR: process j in BLOCKS of Tj at a time, with an sC scratchpad of N*Tj.

    For Tj=1 (one column at a time) and keeping sA_col:
      for k:
        load sA_col from A[:,k]  (N bulk A reads)
        for j:
          load B_cache from B[k,j]  (1 bulk B read)
          for i:
            acc C[i,j] using sA_col[i] and B_cache
            -- but C[i,j] is N*N cells in sC or bulk

    With Tj=1 and k outer:
      sC needs N*N=256 cells to hold the whole C matrix (all j at once).
      OR we process only ONE column j at a time with k as inner loop:
        for j:
          for k:
            load sA_col from A[:,k]
            load B_cache from B[k,j]
            accumulate sC[0..N-1]

    Here sA_col is reloaded for each j! That's N^2 copies from bulk A.

    Compromise: sA_col stays loaded for multiple j values... only works if j is inner.

    ACTUAL PROPOSAL in the prompt:
    For each k: load sA_col, then for each j: load B[k,j] to addr 1 and compute.
    This requires storing partial C somewhere.

    If we keep sC for ALL j: N*N=256 sC cells -> too expensive.
    If we keep sC for Tj=1: only N cells, but we can't do k outer with j inner.

    ACTUALLY: Process it as k outer, j inner, but write C DIRECTLY to bulk (no sC).
      for k:
        load sA_col from A[:,k]  (N bulk reads)
        for j:
          copy B_cache, B[k,j]   (1 bulk read)
          for i:
            if k==0: mul C[i,j], sA_col[i], B_cache  -- write directly to bulk C
            else: mul tmp, sA_col[i], B_cache; add C[i,j], tmp  -- read+write bulk C

    Bulk C reads (for k>0): (N-1) * N * N = 15 * 256 = 3840 reads at bulk C cost.
    Bulk C addresses: C_base + i*N + j.

    This is what the prompt calls approach B!
    """
    # Layout: B_cache=1, sA_col=2-17, tmp=18, then bulk
    B_cache = 1
    sA_base = 2   # sA_col[0..15] at addr 2-17
    tmp = 18
    scratch_end = 19

    A_base = scratch_end
    B_base = A_base + N * N   # 275
    C_base = B_base + N * N   # 531

    A_at = lambda i, j: A_base + i * N + j
    B_at = lambda i, j: B_base + i * N + j
    C_at = lambda i, j: C_base + i * N + j

    sA = lambda i: sA_base + i

    inputs = ([A_at(i, j) for i in range(N) for j in range(N)] +
              [B_at(i, j) for i in range(N) for j in range(N)])
    outputs = [C_at(i, j) for i in range(N) for j in range(N)]

    lines = [",".join(map(str, inputs))]

    for k in range(N):
        # Load A[:,k] into sA_col
        for i in range(N):
            lines.append(f"copy {sA(i)},{A_at(i, k)}")

        for j in range(N):
            # Load B[k,j] into B_cache
            lines.append(f"copy {B_cache},{B_at(k, j)}")

            for i in range(N):
                if k == 0:
                    # Initialize C[i,j] directly
                    lines.append(f"mul {C_at(i,j)},{sA(i)},{B_cache}")
                else:
                    # Accumulate: tmp = sA[i] * B_cache; C[i,j] += tmp
                    lines.append(f"mul {tmp},{sA(i)},{B_cache}")
                    lines.append(f"add {C_at(i,j)},{tmp}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


def generate_col_stream_v2():
    """
    Version 2: Add a small sC scratchpad for one OUTPUT ROW at a time.
    Process output in blocks: for each i, accumulate C[i,0..N-1] in sC.

    But for k outer, j inner: sC would need N*N = 256 cells.

    Alternative: process Ti=1 row at a time.
    for k:
      load sA_col  (but only 1 cell: sA[0..0])
      nope, that doesn't help.

    Best approach for sC: keep j outer, k inner, sC has N cells.
    for j:
      for k:
        copy B_cache, B[k,j]
        for i:
          copy sA_cache, A[i,k]    <- sA_cache at addr 1 (cost 1)
          if k==0: mul sC[i], sA_cache, B_cache
          else: mul tmp, sA_cache, B_cache; add sC[i], tmp
      for i:
        copy C[i,j], sC[i]

    This is the STANDARD approach but with Tj=1.
    sA_cache at addr 1: read 16 times per (j,k,i) = N^3 / N * 1 = N^2 per k per j...
    Actually: for each (j,k,i): read sA_cache once = N*N*N = 4096 total reads.
    sA_cache at addr 1 (cost 1): 4096 reads = 4096 cost.

    B_cache: for each (j,k): read N times (once per i) = N*N*N = 4096 total reads.
    B_cache at addr 2 (cost 2): 4096 reads = 8192 cost.

    sC[i] = addr 3+i (cost 2-5, 16 cells):
    For each (j,k,i) with k>0: read sC[i] once = N*(N-1)*N = 16*15*16 = 3840 reads.
    Plus exit read: N*N = 256 reads.

    bulk A: for each (j,k,i): read A[i,k] once = N^3 = 4096 reads at addr 19-274.

    This is Tj=1 of the current algorithm, which should score worse than Tj=4.
    Let me check.
    """
    sA_cache = 1
    B_cache = 2
    sC_base = 3    # sC[0..15] at addr 3-18
    tmp = 19
    scratch_end = 20

    A_base = scratch_end
    B_base = A_base + N * N
    C_base = B_base + N * N

    A_at = lambda i, j: A_base + i * N + j
    B_at = lambda i, j: B_base + i * N + j
    C_at = lambda i, j: C_base + i * N + j

    sC = lambda i: sC_base + i

    inputs = ([A_at(i, j) for i in range(N) for j in range(N)] +
              [B_at(i, j) for i in range(N) for j in range(N)])
    outputs = [C_at(i, j) for i in range(N) for j in range(N)]

    lines = [",".join(map(str, inputs))]

    for j in range(N):
        for k in range(N):
            lines.append(f"copy {B_cache},{B_at(k, j)}")
            for i in range(N):
                lines.append(f"copy {sA_cache},{A_at(i, k)}")
                if k == 0:
                    lines.append(f"mul {sC(i)},{sA_cache},{B_cache}")
                else:
                    lines.append(f"mul {tmp},{sA_cache},{B_cache}")
                    lines.append(f"add {sC(i)},{tmp}")
        for i in range(N):
            lines.append(f"copy {C_at(i,j)},{sC(i)}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


def generate_col_stream_v3():
    """
    Version 3: k outer, j inner, with sC as FULL N*N scratchpad (too expensive but let's see).
    Actually let's do a COLUMN sC: sC[i] for a specific j, process one j at a time but k outer.
    This requires sC to persist across all k iterations for ALL j simultaneously.

    OR: Use the true "column streaming" approach:
    Preload A[:,k] into sA_col ONCE, then for j in range(N): compute N elements of output column.
    Store partial sums in bulk C (accept the cost of re-reading bulk C).

    Layout:
      B_cache  = 1       (cost 1)
      sA_col   = 2-17   (cost 2-5, 16 cells)
      tmp      = 18      (cost 5)
      A bulk   = 19-274  (cost 5-17)
      B bulk   = 275-530 (cost 17-24)
      C bulk   = 531-786 (output, also used as accumulator for k>0)
    """
    B_cache = 1
    sA_base = 2
    tmp = 18
    scratch_end = 19

    A_base = scratch_end   # 19
    B_base = A_base + N * N   # 275
    C_base = B_base + N * N   # 531

    A_at = lambda i, j: A_base + i * N + j
    B_at = lambda i, j: B_base + i * N + j
    C_at = lambda i, j: C_base + i * N + j

    sA = lambda i: sA_base + i

    inputs = ([A_at(i, j) for i in range(N) for j in range(N)] +
              [B_at(i, j) for i in range(N) for j in range(N)])
    outputs = [C_at(i, j) for i in range(N) for j in range(N)]

    lines = [",".join(map(str, inputs))]

    for k in range(N):
        # Load A[:,k] into sA_col (16 reads from bulk A)
        for i in range(N):
            lines.append(f"copy {sA(i)},{A_at(i, k)}")

        for j in range(N):
            lines.append(f"copy {B_cache},{B_at(k, j)}")

            for i in range(N):
                if k == 0:
                    lines.append(f"mul {C_at(i,j)},{sA(i)},{B_cache}")
                else:
                    lines.append(f"mul {tmp},{sA(i)},{B_cache}")
                    lines.append(f"add {C_at(i,j)},{tmp}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


def generate_col_stream_v4():
    """
    Version 4: k outer, j in blocks of Tj=4, with sC scratchpad of N*Tj.

    This is: Ti=N=16, Tj=4, Tk=1, with the A column cached.
    sC has 16*4=64 cells at addr 19-82 (cost 5-10).

    Layout:
      B_cache  = 1       (cost 1)  <- B[k, bj*Tj+jj], read N times per (k,j)
      sA_col   = 2-17   (cost 2-5, 16 cells) <- A[:,k], read Tj times per (k,i)
      sC       = 18-81  (cost 5-9, 64 cells = 16*4)
      tmp      = 82      (cost 10)
      A bulk   = 83-338
      B bulk   = 339-594
      C bulk   = 595-850

    Loop: for k: load sA_col; for bj: for jj: load B[k,bj*4+jj] into B_j[jj];
          for i: for jj: acc sC[i*4+jj]

    Hmm, actually with Tj=4 we need to load 4 B values and have sC with 16*4 cells.
    The sC at addr 18-81 costs 5-9. That's 64 cells with high cost.

    Let me instead try: Tj=1, but with sC[i] as a scratchpad column (N cells).
    For one j at a time, accumulate over k using sA_col to cache A[:,k].
    But then sA_col needs to be reloaded for each j.

    Actually the cleanest approach (given the data):
    j outer, k inner, sA_cache=1 (single cell), sB_col cached across i.

    for j:
      for k:
        copy B_cache, B[k,j]  <- one B value
        for i:
          copy sA_cache, A[i,k]  <- one A value at addr 1
          acc sC[i]
      writeback sC[i] to C[i,j]

    This is exp_rank2_v11's Ti=16, Tj=1 variant. Let me try it.
    """
    sA_cache = 1
    sC_base = 2     # sC[0..15] at addr 2-17 (cost 2-5)
    B_cache = 18    # B cache at addr 18 (cost 5)
    tmp = 19        # tmp at addr 19 (cost 5)
    scratch_end = 20

    A_base = scratch_end
    B_base = A_base + N * N
    C_base = B_base + N * N

    A_at = lambda i, j: A_base + i * N + j
    B_at = lambda i, j: B_base + i * N + j
    C_at = lambda i, j: C_base + i * N + j

    sC = lambda i: sC_base + i

    inputs = ([A_at(i, j) for i in range(N) for j in range(N)] +
              [B_at(i, j) for i in range(N) for j in range(N)])
    outputs = [C_at(i, j) for i in range(N) for j in range(N)]

    lines = [",".join(map(str, inputs))]

    for j in range(N):
        for k in range(N):
            lines.append(f"copy {B_cache},{B_at(k, j)}")
            for i in range(N):
                lines.append(f"copy {sA_cache},{A_at(i, k)}")
                if k == 0:
                    lines.append(f"mul {sC(i)},{sA_cache},{B_cache}")
                else:
                    lines.append(f"mul {tmp},{sA_cache},{B_cache}")
                    lines.append(f"add {sC(i)},{tmp}")
        for i in range(N):
            lines.append(f"copy {C_at(i,j)},{sC(i)}")

    lines.append(",".join(map(str, outputs)))
    return "\n".join(lines)


def generate_col_stream_v5():
    """
    Version 5: The TRUE column streaming as described in the prompt.

    For each k (0..15):
      1. Copy A[:,k] (16 values) into sA_col at addr 2-17 ONCE.
      2. For each j (0..15):
         a. Copy B[k,j] to addr 1 (B_cache).
         b. For each i: accumulate C result.

    The C result goes into sC_col (one full column, N=16 cells),
    but we need N columns of sC (N*N cells total) or we process one column at a time.

    If we process one column j at a time, we need to:
    1. Load sA_col[:,k] for all k (repeated for each j) OR
    2. Accumulate into full C matrix (need sC with N*N=256 cells)

    ACTUAL EFFICIENT VERSION:
    - Process output by ROWS of C: for each row i, accumulate C[i,0..N-1].
    - sA_cache at 1 (scalar for A[i,k]), B_col at 2-17 (B[:,j] column).
    - This is the transpose of what we tried above.

    Wait, re-reading the prompt more carefully:
    "For each k: 1. Copy A[:,k] into sA_col at addr 2-17.
     2. For each j: a. Copy B[k,j] to addr 1. b. For each i: mul sA_col[i], B_cache."

    This means sC has N*N=256 cells (one for each C[i,j]).
    But those 256 cells at addr 35+ would be expensive.

    OR: Process output COLUMN BY COLUMN, j outer, k inner:
    - sA_col at 2-17 is reloaded for each j (not shared across j).
    - Actually no - the prompt says sA_col is used for ALL j in step 2.

    I think the prompt envisions sC as 256 cells accumulated throughout.
    Let me try that: sC at addr 19-274 (after A_col), then A at 275+, B at 531+, C at 787+.
    That's a very expensive sC...

    Actually the most compact interpretation:
    Keep sC as N cells for one OUTPUT ROW at a time.
    Process: for each output row i:
      - For each k: load A[i,k] into sA_cache at addr 1.
      - For each j: load B[k,j], accumulate into sC[j] (N cells).
      - Writeback sC to C[i,*].

    But A[i,k] changes each k, so sA_col (the column) is needed to reuse it across j.

    Most promising variant: ROWS outer, then k/j inner with sA cached scalar:
    for i:
      for j:
        for k:
          acc C[i,j] using A[i,k] * B[k,j]
    With k innermost: sA_cache at 1 (A[i,k] = 1 cell), B_cache at 2 (B[k,j] = 1 cell), sC at 3-18 (N cells for one row).
    """
    sA_cache = 1
    B_cache = 2
    sC_base = 3     # sC[0..15] at addr 3-18 (cost 2-5) -- one output ROW
    tmp = 19
    scratch_end = 20

    A_base = scratch_end
    B_base = A_base + N * N
    C_base = B_base + N * N

    A_at = lambda i, j: A_base + i * N + j
    B_at = lambda i, j: B_base + i * N + j
    C_at = lambda i, j: C_base + i * N + j

    sC = lambda j: sC_base + j

    inputs = ([A_at(i, j) for i in range(N) for j in range(N)] +
              [B_at(i, j) for i in range(N) for j in range(N)])
    outputs = [C_at(i, j) for i in range(N) for j in range(N)]

    lines = [",".join(map(str, inputs))]

    for i in range(N):
        for j in range(N):
            for k in range(N):
                lines.append(f"copy {sA_cache},{A_at(i, k)}")
                lines.append(f"copy {B_cache},{B_at(k, j)}")
                if k == 0:
                    lines.append(f"mul {sC(j)},{sA_cache},{B_cache}")
                else:
                    lines.append(f"mul {tmp},{sA_cache},{B_cache}")
                    lines.append(f"add {sC(j)},{tmp}")
        for j in range(N):
            lines.append(f"copy {C_at(i,j)},{sC(j)}")

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

    print(f"Direction B: Column streaming variants")
    print(f"Record to beat: {BEST:,}")
    print()

    variants = [
        ("v1: k-outer, j-inner, sA_col@2-17, direct-to-C", generate_col_stream_v1),
        ("v2: j-outer, k-inner, sA_cache@1, B_cache@2, sC@3-18", generate_col_stream_v2),
        ("v3: k-outer, j-inner, sA_col@2-17, tmp@18, direct-to-C", generate_col_stream_v3),
        ("v4: j-outer, k-inner, sA_cache@1, sC@2-17, B_cache@18", generate_col_stream_v4),
        ("v5: i-outer, j-mid, k-inner, sA@1, B@2, sC@3-18", generate_col_stream_v5),
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
            import traceback
            traceback.print_exc()
            print()

    print(f"Best column streaming score: {best_score:,}  ({best_label})")
    print(f"Delta vs record: {BEST - best_score:+,}")

    if best_score is not None and best_score < BEST:
        ir_path = Path(__file__).parent / "ir" / f"new_record_{best_score}.ir"
        ir_path.write_text(best_ir + "\n")
        print(f"*** NEW RECORD! Saved to {ir_path}")
